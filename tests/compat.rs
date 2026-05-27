/// Compatibility tests: rsomics-vcf-view data records must be byte-identical to
/// bcftools view (version 1.22, htslib 1.22.1) for the corresponding mode.
///
/// Tested modes:
///   -v snps         keep SNPs only
///   -V indels       exclude indels
///   -f PASS         keep PASS records only
///   -s sample1,sample2 -I   subset to two samples without recomputing INFO
///   -H              suppress header (data records only)
///   --header-only   header lines only (no data)
///
/// Header lines differ (bcftools stamps its own ##bcftools_viewCommand and
/// ##FILTER=<ID=PASS> lines). Comparison is on data records only.
///
/// Quirks matched:
///   - SNP classification: same-length REF+ALT with exactly one differing base.
///   - MNP: same-length REF+ALT with two or more differing bases.
///   - "." in FILTER treated as PASS for -f PASS matching.
///   - -s sample subsetting without -I would cause bcftools to recompute AC/AN;
///     we do NOT recompute (behaviour matches `bcftools view -s ... -I`).
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-vcf-view"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/mixed.vcf")
}

fn bcftools_path() -> Option<String> {
    let candidates = [
        "bcftools",
        "/opt/homebrew/bin/bcftools",
        "/opt/homebrew/Caskroom/miniforge/base/envs/imotif-pipeline/bin/bcftools",
        "/usr/bin/bcftools",
        "/usr/local/bin/bcftools",
    ];
    for candidate in &candidates {
        let ok = Command::new(candidate)
            .arg("--version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if ok {
            return Some(candidate.to_string());
        }
    }
    None
}

fn bcftools_version(path: &str) -> String {
    let out = Command::new(path).arg("--version").output().unwrap().stdout;
    String::from_utf8_lossy(&out)
        .lines()
        .next()
        .unwrap_or("")
        .to_owned()
}

/// Extract data (non-header) records as lines.
fn records(vcf: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(vcf)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

#[test]
fn types_snps_match_bcftools() {
    let Some(bcftools) = bcftools_path() else {
        eprintln!("skipping: bcftools not found");
        return;
    };
    eprintln!("{}", bcftools_version(&bcftools));

    let vcf = fixture();

    let ours = ours().args(["--types", "snps"]).arg(&vcf).output().unwrap();
    assert!(
        ours.status.success(),
        "rsomics-vcf-view failed: {}",
        String::from_utf8_lossy(&ours.stderr)
    );

    let theirs = Command::new(&bcftools)
        .args(["view", "-v", "snps"])
        .arg(&vcf)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .unwrap();
    assert!(theirs.status.success(), "bcftools view -v snps failed");

    assert_eq!(
        records(&ours.stdout),
        records(&theirs.stdout),
        "-v snps: data records differ"
    );
}

#[test]
fn exclude_indels_match_bcftools() {
    let Some(bcftools) = bcftools_path() else {
        eprintln!("skipping: bcftools not found");
        return;
    };

    let vcf = fixture();

    let ours = ours()
        .args(["--exclude-types", "indels"])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(ours.status.success());

    let theirs = Command::new(&bcftools)
        .args(["view", "-V", "indels"])
        .arg(&vcf)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .unwrap();
    assert!(theirs.status.success());

    assert_eq!(
        records(&ours.stdout),
        records(&theirs.stdout),
        "-V indels: data records differ"
    );
}

#[test]
fn apply_filters_pass_match_bcftools() {
    let Some(bcftools) = bcftools_path() else {
        eprintln!("skipping: bcftools not found");
        return;
    };

    let vcf = fixture();

    let ours = ours().args(["-f", "PASS"]).arg(&vcf).output().unwrap();
    assert!(ours.status.success());

    let theirs = Command::new(&bcftools)
        .args(["view", "-f", "PASS"])
        .arg(&vcf)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .unwrap();
    assert!(theirs.status.success());

    assert_eq!(
        records(&ours.stdout),
        records(&theirs.stdout),
        "-f PASS: data records differ"
    );
}

#[test]
fn samples_subset_match_bcftools() {
    let Some(bcftools) = bcftools_path() else {
        eprintln!("skipping: bcftools not found");
        return;
    };

    let vcf = fixture();

    // Use -I with bcftools to skip INFO recomputation, matching our behaviour.
    let ours = ours()
        .args(["-s", "sample1,sample2"])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(
        ours.status.success(),
        "rsomics-vcf-view -s failed: {}",
        String::from_utf8_lossy(&ours.stderr)
    );

    let theirs = Command::new(&bcftools)
        .args(["view", "-s", "sample1,sample2", "-I"])
        .arg(&vcf)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .unwrap();
    assert!(theirs.status.success(), "bcftools view -s -I failed");

    assert_eq!(
        records(&ours.stdout),
        records(&theirs.stdout),
        "-s sample1,sample2: data records differ"
    );
}

#[test]
fn no_header_match_bcftools() {
    let Some(bcftools) = bcftools_path() else {
        eprintln!("skipping: bcftools not found");
        return;
    };

    let vcf = fixture();

    let ours = ours().arg("-H").arg(&vcf).output().unwrap();
    assert!(ours.status.success());

    let theirs = Command::new(&bcftools)
        .args(["view", "-H"])
        .arg(&vcf)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .unwrap();
    assert!(theirs.status.success());

    // -H: no header lines in either output; compare directly.
    let ours_lines: Vec<String> = String::from_utf8_lossy(&ours.stdout)
        .lines()
        .filter(|l| !l.is_empty())
        .map(str::to_owned)
        .collect();
    let their_lines: Vec<String> = String::from_utf8_lossy(&theirs.stdout)
        .lines()
        .filter(|l| !l.is_empty())
        .map(str::to_owned)
        .collect();

    assert_eq!(ours_lines, their_lines, "-H: output lines differ");
}

#[test]
fn header_only_has_no_data_records() {
    let vcf = fixture();
    let ours = ours().arg("--header-only").arg(&vcf).output().unwrap();
    assert!(ours.status.success());

    let data = records(&ours.stdout);
    assert!(
        data.is_empty(),
        "--header-only: unexpected data records: {data:?}"
    );

    // Verify header lines are present.
    let header_lines: Vec<String> = String::from_utf8_lossy(&ours.stdout)
        .lines()
        .filter(|l| l.starts_with('#'))
        .map(str::to_owned)
        .collect();
    assert!(
        !header_lines.is_empty(),
        "--header-only: no header lines in output"
    );
}

#[test]
fn sites_only_matches_bcftools() {
    let Some(bcftools) = bcftools_path() else {
        eprintln!("skipping: bcftools not found");
        return;
    };

    let vcf = fixture();

    let ours = ours().arg("--sites-only").arg(&vcf).output().unwrap();
    assert!(
        ours.status.success(),
        "rsomics-vcf-view --sites-only failed: {}",
        String::from_utf8_lossy(&ours.stderr)
    );

    // bcftools view -G drops genotype columns (FORMAT + all samples).
    let theirs = Command::new(&bcftools)
        .args(["view", "-G"])
        .arg(&vcf)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .unwrap();
    assert!(theirs.status.success(), "bcftools view -G failed");

    assert_eq!(
        records(&ours.stdout),
        records(&theirs.stdout),
        "--sites-only: data records differ from bcftools view -G"
    );
}
