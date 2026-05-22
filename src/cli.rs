use std::io::BufWriter;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_vcf_view::{
    ViewConfig, parse_filter_list, parse_type_list, read_samples_file, view_vcf,
};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-vcf-view",
    version,
    about,
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Input VCF file (plain or .vcf.gz).
    #[arg(value_name = "INPUT")]
    pub input: PathBuf,

    /// Output file (default stdout).
    #[arg(short = 'o', long = "output", default_value = "-")]
    pub output: String,

    /// Keep only these variant types: snps,indels,mnps,other (comma-separated).
    /// Long-only because -v is occupied by CommonFlags --verbose.
    #[arg(long = "types", value_name = "LIST")]
    pub types: Option<String>,

    /// Exclude these variant types: snps,indels,mnps,other (comma-separated).
    /// Long-only because -V conflicts with clap's auto --version flag.
    #[arg(long = "exclude-types", value_name = "LIST")]
    pub exclude_types: Option<String>,

    /// Keep records with FILTER in this list (comma-separated; PASS and . are equivalent).
    #[arg(short = 'f', long = "apply-filters", value_name = "LIST")]
    pub apply_filters: Option<String>,

    /// Comma-separated list of sample names to keep (genotype columns).
    /// INFO fields are NOT recomputed (equivalent to bcftools view -s ... -I).
    #[arg(short = 's', long = "samples", value_name = "LIST")]
    pub samples: Option<String>,

    /// File containing sample names to keep, one per line.
    #[arg(short = 'S', long = "samples-file", value_name = "FILE")]
    pub samples_file: Option<PathBuf>,

    /// Print header lines only (no data records).
    #[arg(long = "header-only")]
    pub header_only: bool,

    /// Suppress all header lines (print data records only).
    #[arg(short = 'H', long = "no-header")]
    pub no_header: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        if self.header_only && self.no_header {
            return Err(RsomicsError::InvalidInput(
                "--header-only and --no-header are mutually exclusive".into(),
            ));
        }

        let keep_types = self.types.as_deref().map(parse_type_list).transpose()?;
        let exclude_types = self
            .exclude_types
            .as_deref()
            .map(parse_type_list)
            .transpose()?;
        let apply_filters = self.apply_filters.as_deref().map(parse_filter_list);

        // Resolve sample list: -s overrides -S when both are given (match bcftools precedence).
        let samples: Option<Vec<String>> = if let Some(ref list) = self.samples {
            Some(
                list.split(',')
                    .map(|s| s.trim().to_owned())
                    .filter(|s| !s.is_empty())
                    .collect(),
            )
        } else if let Some(ref path) = self.samples_file {
            Some(read_samples_file(path)?)
        } else {
            None
        };

        let cfg = ViewConfig {
            keep_types,
            exclude_types,
            apply_filters,
            samples,
            header_only: self.header_only,
            no_header: self.no_header,
        };

        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                std::fs::File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };

        let stats = view_vcf(&self.input, &mut out, &cfg)?;

        if !self.common.quiet {
            eprintln!("{}/{} records kept", stats.kept, stats.total);
        }

        Ok(())
    }
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.execute()
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Subset and filter VCF records — type, FILTER, sample, and header subsetting.",
    origin: Some(Origin {
        upstream: "bcftools view",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/gigascience/giab008"),
    }),
    usage_lines: &["[OPTIONS] <INPUT.vcf>"],
    sections: &[
        Section {
            title: "SCOPE PARTITION",
            flags: &[FlagSpec {
                short: None,
                long: "NOTE",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "rsomics-vcf-view handles type/FILTER/sample/header subsetting. \
                              Expression-based filtering (-i/-e include/exclude expressions) \
                              is handled by rsomics-vcf-filter.",
                why_default: None,
            }],
        },
        Section {
            title: "OPTIONS",
            flags: &[
                FlagSpec {
                    short: None,
                    long: "INPUT",
                    aliases: &[],
                    value: Some("<path>"),
                    type_hint: Some("Path"),
                    required: true,
                    default: None,
                    description: "Input VCF file (plain or gzip-compressed).",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('o'),
                    long: "output",
                    aliases: &[],
                    value: Some("<path>"),
                    type_hint: Some("Path"),
                    required: false,
                    default: Some("-"),
                    description: "Output VCF file (default: stdout).",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "types",
                    aliases: &[],
                    value: Some("<LIST>"),
                    type_hint: Some("String"),
                    required: false,
                    default: None,
                    description: "Keep only these variant types. Comma-separated from: snps, indels, mnps, other. \
                                  (bcftools -v; long-only here because -v is --verbose in CommonFlags)",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "exclude-types",
                    aliases: &[],
                    value: Some("<LIST>"),
                    type_hint: Some("String"),
                    required: false,
                    default: None,
                    description: "Exclude these variant types. Same values as --types. \
                                  (bcftools -V; long-only here because -V conflicts with clap --version)",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('f'),
                    long: "apply-filters",
                    aliases: &[],
                    value: Some("<LIST>"),
                    type_hint: Some("String"),
                    required: false,
                    default: None,
                    description: "Keep records whose FILTER column is in this comma-separated list. \
                                  PASS and . are treated as equivalent.",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('s'),
                    long: "samples",
                    aliases: &[],
                    value: Some("<LIST>"),
                    type_hint: Some("String"),
                    required: false,
                    default: None,
                    description: "Comma-separated list of sample names to retain. \
                                  INFO fields are not recomputed (equivalent to bcftools -I).",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('S'),
                    long: "samples-file",
                    aliases: &[],
                    value: Some("<FILE>"),
                    type_hint: Some("Path"),
                    required: false,
                    default: None,
                    description: "File of sample names to retain, one per line.",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "header-only",
                    aliases: &[],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Print header lines only; suppress all data records.",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('H'),
                    long: "no-header",
                    aliases: &[],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Suppress header lines; print data records only.",
                    why_default: None,
                },
            ],
        },
    ],
    examples: &[
        Example {
            description: "Keep only SNPs",
            command: "rsomics-vcf-view --types snps input.vcf",
        },
        Example {
            description: "Exclude indels",
            command: "rsomics-vcf-view --exclude-types indels input.vcf",
        },
        Example {
            description: "Keep only PASS records",
            command: "rsomics-vcf-view -f PASS input.vcf",
        },
        Example {
            description: "Subset to two samples (INFO not recomputed, like bcftools -I)",
            command: "rsomics-vcf-view -s NA12878,NA12879 input.vcf",
        },
        Example {
            description: "Header lines only",
            command: "rsomics-vcf-view --header-only input.vcf",
        },
        Example {
            description: "Data records only (no header)",
            command: "rsomics-vcf-view -H input.vcf",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
