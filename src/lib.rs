use std::collections::HashSet;
use std::io::{self, Read};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// Variant type as classified by bcftools view -v/-V.
///
/// Classification follows bcftools' allele-level rules per ALT allele:
/// - SNP: REF and ALT have identical length AND exactly one base position differs.
/// - MNP: REF and ALT have identical length AND more than one base position differs.
/// - Indel: REF and ALT differ in length (insertion or deletion).
/// - Other: any ALT containing '<' (symbolic allele) or '*' (missing allele).
///
/// A record matches a type if ANY of its ALT alleles matches that type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VcfType {
    Snp,
    Indel,
    Mnp,
    Other,
}

impl VcfType {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "snps" | "snp" => Some(Self::Snp),
            "indels" | "indel" => Some(Self::Indel),
            "mnps" | "mnp" => Some(Self::Mnp),
            "other" => Some(Self::Other),
            _ => None,
        }
    }
}

fn allele_type(ref_allele: &[u8], alt: &[u8]) -> VcfType {
    // Symbolic or missing allele.
    if alt.contains(&b'<') || alt.contains(&b'*') {
        return VcfType::Other;
    }
    let rlen = ref_allele.len();
    let alen = alt.len();
    if rlen != alen {
        return VcfType::Indel;
    }
    // Same length: count differing positions.
    let diffs = ref_allele
        .iter()
        .zip(alt.iter())
        .filter(|(r, a)| r != a)
        .count();
    if diffs <= 1 {
        VcfType::Snp
    } else {
        VcfType::Mnp
    }
}

/// True if the record's ALT field has at least one allele matching any type in `types`.
fn record_matches_any_type(ref_allele: &[u8], alt_field: &[u8], types: &HashSet<VcfType>) -> bool {
    // ALT is comma-separated; a record matches if ANY allele matches.
    for alt in alt_field.split(|&b| b == b',') {
        if types.contains(&allele_type(ref_allele, alt)) {
            return true;
        }
    }
    false
}

/// True if FILTER field passes the filter-set constraint.
///
/// bcftools `view -f` semantics: keep the record if its FILTER value is in the
/// requested list. Both "." and "PASS" are treated as the PASS state:
/// - requesting "PASS" in the list matches records with FILTER == "PASS" or ".".
/// - requesting "." in the list matches records with FILTER == "PASS" or ".".
fn filter_passes(filter_field: &[u8], allowed: &HashSet<Vec<u8>>) -> bool {
    // FILTER may be multi-valued ("PASS;LowQual" is non-standard but semicolons appear;
    // the standard delimiter for multiple applied filters is ";").
    // For simplicity: split on ";" like bcftools does.
    let pass_bytes: &[u8] = b"PASS";
    let dot_bytes: &[u8] = b".";

    // Normalise: if the field is "." or "PASS", treat both as the PASS state.
    let wants_pass = allowed.contains(pass_bytes) || allowed.contains(dot_bytes);

    // If the FILTER field is "." or "PASS", the record has no filter applied.
    if filter_field == dot_bytes || filter_field == pass_bytes {
        return wants_pass;
    }

    // Otherwise split on ";" and check each filter token.
    for token in filter_field.split(|&b| b == b';') {
        if allowed.contains(token) {
            return true;
        }
    }
    false
}

pub struct ViewConfig {
    /// If Some, keep only records whose type is in this set.
    pub keep_types: Option<HashSet<VcfType>>,
    /// If Some, drop records whose type is in this set.
    pub exclude_types: Option<HashSet<VcfType>>,
    /// If Some, keep only records whose FILTER value is in this set.
    pub apply_filters: Option<HashSet<Vec<u8>>>,
    /// If Some, keep only these sample columns (by name, in requested order).
    pub samples: Option<Vec<String>>,
    /// Emit header lines only (no data records).
    pub header_only: bool,
    /// Suppress all header lines (data records only).
    pub no_header: bool,
    /// Strip FORMAT and all sample columns (picard MakeSitesOnlyVcf equivalent).
    /// ##FORMAT header lines are also suppressed.
    pub sites_only: bool,
}

pub struct ViewStats {
    pub total: u64,
    pub kept: u64,
}

/// Run the view/subset pass over `input`, writing to `output`.
pub fn view_vcf(input: &Path, output: &mut dyn io::Write, cfg: &ViewConfig) -> Result<ViewStats> {
    let raw = std::fs::read(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let data = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut d = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut d)
            .map_err(RsomicsError::Io)?;
        d
    } else {
        raw
    };

    // Resolve sample column indices and rewrite #CHROM if needed.
    let mut sample_indices: Option<Vec<usize>> = None;
    let mut chrom_line_rewritten: Option<Vec<u8>> = None;

    let all_lines: Vec<&[u8]> = data
        .split(|&b| b == b'\n')
        .map(|l| match l.last() {
            Some(b'\r') => &l[..l.len() - 1],
            _ => l,
        })
        .filter(|l| !l.is_empty())
        .collect();

    let header_end = all_lines.iter().position(|l| !l.starts_with(b"#"));
    let header_count = header_end.unwrap_or(all_lines.len());

    // Resolve sample column indices from the #CHROM line.
    if let Some(sample_names) = &cfg.samples {
        // Find the #CHROM line (last header line starting with a single '#').
        let chrom_line_idx = (0..header_count)
            .rfind(|&i| all_lines[i].starts_with(b"#CHROM") || all_lines[i].starts_with(b"#chrom"));
        if let Some(idx) = chrom_line_idx {
            let chrom_line = all_lines[idx];
            let cols: Vec<&[u8]> = chrom_line.split(|&b| b == b'\t').collect();
            // Fixed VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT sample…]
            // Sample names start at column index 9 (0-based).
            let fixed = 9usize;
            if cols.len() > fixed {
                let vcf_samples: Vec<&[u8]> = cols[fixed..].to_vec();
                // Map each requested name → its 0-based index within the sample columns.
                let mut indices: Vec<usize> = Vec::with_capacity(sample_names.len());
                for name in sample_names {
                    let pos = vcf_samples
                        .iter()
                        .position(|&s| s == name.as_bytes())
                        .ok_or_else(|| {
                            RsomicsError::InvalidInput(format!(
                                "sample '{}' not found in VCF header",
                                name
                            ))
                        })?;
                    indices.push(pos);
                }
                // Rewrite #CHROM line to retain only the requested samples.
                let mut new_cols: Vec<&[u8]> = cols[..fixed].to_vec();
                for &i in &indices {
                    new_cols.push(vcf_samples[i]);
                }
                let mut rewritten = new_cols.join(&b'\t');
                rewritten.push(b'\n');
                chrom_line_rewritten = Some(rewritten);
                sample_indices = Some(indices);
            } else {
                // No sample columns — nothing to subset.
                sample_indices = Some(vec![]);
            }
        }
    }

    // Write header lines.
    if !cfg.no_header {
        for &line in all_lines[..header_count].iter() {
            // --sites-only: drop ##FORMAT lines and truncate #CHROM at column 8.
            if cfg.sites_only {
                if line.starts_with(b"##FORMAT") {
                    continue;
                }
                if line.starts_with(b"#CHROM") || line.starts_with(b"#chrom") {
                    // Keep only the 8 fixed columns (CHROM POS ID REF ALT QUAL FILTER INFO).
                    let cols: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
                    let fixed = cols.len().min(8);
                    let truncated = cols[..fixed].join(&b'\t');
                    output.write_all(&truncated).map_err(RsomicsError::Io)?;
                    output.write_all(b"\n").map_err(RsomicsError::Io)?;
                    continue;
                }
            }
            // Replace #CHROM line if we rewrote it for sample subsetting.
            if sample_indices.is_some()
                && (line.starts_with(b"#CHROM") || line.starts_with(b"#chrom"))
                && let Some(ref rw) = chrom_line_rewritten
            {
                output.write_all(rw).map_err(RsomicsError::Io)?;
                continue;
            }
            output.write_all(line).map_err(RsomicsError::Io)?;
            output.write_all(b"\n").map_err(RsomicsError::Io)?;
        }
    }

    if cfg.header_only {
        return Ok(ViewStats { total: 0, kept: 0 });
    }

    let data_lines = &all_lines[header_count..];

    let mut stats = ViewStats { total: 0, kept: 0 };

    for &line in data_lines {
        stats.total += 1;

        let mut cols = line.splitn(9, |&b| b == b'\t');
        let _chrom = cols.next().unwrap_or(b"");
        let _pos = cols.next().unwrap_or(b"");
        let _id = cols.next().unwrap_or(b"");
        let ref_allele = cols.next().unwrap_or(b"");
        let alt_field = cols.next().unwrap_or(b"");
        let _qual = cols.next().unwrap_or(b"");
        let filter_field = cols.next().unwrap_or(b".");
        // INFO is col 7; FORMAT+samples are the rest (col 8 in splitn(9)).
        let _info = cols.next().unwrap_or(b"");
        let _rest = cols.next(); // FORMAT + all sample columns joined (we split further below)

        // -v / -V type filter.
        if let Some(ref keep) = cfg.keep_types
            && !record_matches_any_type(ref_allele, alt_field, keep)
        {
            continue;
        }
        if let Some(ref excl) = cfg.exclude_types
            && record_matches_any_type(ref_allele, alt_field, excl)
        {
            continue;
        }

        // -f filter.
        if let Some(ref allowed) = cfg.apply_filters
            && !filter_passes(filter_field, allowed)
        {
            continue;
        }

        stats.kept += 1;

        // --sites-only: emit only the 8 fixed columns (no FORMAT, no samples).
        if cfg.sites_only {
            let all_cols: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
            let fixed = all_cols.len().min(8);
            let truncated = all_cols[..fixed].join(&b'\t');
            output.write_all(&truncated).map_err(RsomicsError::Io)?;
            output.write_all(b"\n").map_err(RsomicsError::Io)?;
            continue;
        }

        // -s sample subsetting: reconstruct the line with only kept columns.
        if let Some(ref indices) = sample_indices {
            // Re-split the full line on tab to access all columns.
            let all_cols: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
            // Fixed columns 0-8: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
            // Sample columns start at 9.
            let fixed_end = 9usize.min(all_cols.len());
            let mut out_cols: Vec<&[u8]> = all_cols[..fixed_end].to_vec();
            let sample_cols = &all_cols[fixed_end..];
            for &i in indices {
                out_cols.push(sample_cols.get(i).copied().unwrap_or(b"."));
            }
            let reconstructed = out_cols.join(&b'\t');
            output.write_all(&reconstructed).map_err(RsomicsError::Io)?;
            output.write_all(b"\n").map_err(RsomicsError::Io)?;
        } else {
            output.write_all(line).map_err(RsomicsError::Io)?;
            output.write_all(b"\n").map_err(RsomicsError::Io)?;
        }
    }

    Ok(stats)
}

/// Parse a comma-or-space separated type list string into a HashSet of VcfType.
pub fn parse_type_list(s: &str) -> Result<HashSet<VcfType>> {
    let mut set = HashSet::new();
    for token in s.split([',', ' ']) {
        let t = token.trim();
        if t.is_empty() {
            continue;
        }
        let vt = VcfType::parse(t).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "unknown variant type '{}': expected snps, indels, mnps, other",
                t
            ))
        })?;
        set.insert(vt);
    }
    Ok(set)
}

/// Parse a comma-separated FILTER list into a HashSet of byte-vecs.
pub fn parse_filter_list(s: &str) -> HashSet<Vec<u8>> {
    s.split(',')
        .map(|t| t.trim().as_bytes().to_vec())
        .filter(|v| !v.is_empty())
        .collect()
}

/// Read a samples file (one name per line, ignoring blank lines and '#' comments).
pub fn read_samples_file(path: &Path) -> Result<Vec<String>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    Ok(content
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_owned)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn type_classification() {
        assert_eq!(allele_type(b"A", b"T"), VcfType::Snp);
        assert_eq!(allele_type(b"AT", b"CT"), VcfType::Snp);
        assert_eq!(allele_type(b"AT", b"CA"), VcfType::Mnp);
        assert_eq!(allele_type(b"A", b"AT"), VcfType::Indel);
        assert_eq!(allele_type(b"AT", b"A"), VcfType::Indel);
        assert_eq!(allele_type(b"A", b"<DEL>"), VcfType::Other);
        assert_eq!(allele_type(b"ACT", b"GTT"), VcfType::Mnp);
    }

    #[test]
    fn filter_pass_semantics() {
        let allowed: HashSet<Vec<u8>> = ["PASS".as_bytes().to_vec()].into();
        assert!(filter_passes(b"PASS", &allowed));
        assert!(filter_passes(b".", &allowed));
        assert!(!filter_passes(b"LowQual", &allowed));
    }

    #[test]
    fn filter_dot_semantics() {
        let allowed: HashSet<Vec<u8>> = [".".as_bytes().to_vec()].into();
        assert!(filter_passes(b"PASS", &allowed));
        assert!(filter_passes(b".", &allowed));
        assert!(!filter_passes(b"LowQual", &allowed));
    }

    #[test]
    fn filter_named_filter() {
        let allowed: HashSet<Vec<u8>> = ["LowQual".as_bytes().to_vec()].into();
        assert!(filter_passes(b"LowQual", &allowed));
        assert!(!filter_passes(b"PASS", &allowed));
        assert!(!filter_passes(b".", &allowed));
    }

    #[test]
    fn type_zero_diffs_is_snp() {
        // REF == ALT (degenerate): one diff count = 0 → SNP category (bcftools agrees).
        assert_eq!(allele_type(b"A", b"A"), VcfType::Snp);
    }

    #[test]
    fn sites_only_strips_format_and_samples() {
        use std::io::Cursor;
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let vcf = b"\
##fileformat=VCFv4.1\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12879\n\
chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t0/0\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t1/1\t0/1\n\
";
        let mut tmp = NamedTempFile::new().unwrap();
        tmp.write_all(vcf).unwrap();

        let mut out = Cursor::new(Vec::new());
        let cfg = ViewConfig {
            keep_types: None,
            exclude_types: None,
            apply_filters: None,
            samples: None,
            header_only: false,
            no_header: false,
            sites_only: true,
        };
        let stats = view_vcf(tmp.path(), &mut out, &cfg).unwrap();
        assert_eq!(stats.kept, 2);

        let result = String::from_utf8(out.into_inner()).unwrap();
        // ##FORMAT line must be absent.
        assert!(!result.contains("##FORMAT"));
        // #CHROM line must have exactly 8 fields (no FORMAT or sample columns).
        let chrom_line = result.lines().find(|l| l.starts_with("#CHROM")).unwrap();
        assert_eq!(chrom_line.split('\t').count(), 8, "#CHROM: {chrom_line}");
        // Data records must have exactly 8 fields.
        for data_line in result.lines().filter(|l| !l.starts_with('#')) {
            assert_eq!(data_line.split('\t').count(), 8, "data: {data_line}");
        }
    }
}
