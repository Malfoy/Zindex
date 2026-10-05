use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use anyhow::{ensure, Context, Result};
use clap::ValueEnum;
use deko::read::AnyDecoder;
use helicase::input::{FromMmap, FromSlice};
use helicase::{Config, FastxParser, HelicaseParser, ParserOptions};
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};

pub mod findere;

pub const MAX_HASHES: usize = 32;
pub const DEFAULT_MINIMIZER_SIZE: usize = 21;
pub const DEFAULT_MODIMIZER_SAMPLING: u64 = 16;
pub const DEFAULT_INDEX_ZSTD_LEVEL: i32 = 0;
pub const ZSTD_FRAME_MAGIC: [u8; 4] = [0x28, 0xb5, 0x2f, 0xfd];

pub const PARSER_CONFIG: Config = ParserOptions::default()
    .ignore_headers()
    .ignore_quality()
    .dna_string()
    .split_non_actg()
    .return_record(false)
    .config();

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum IndexMode {
    /// Index every canonical k-mer.
    Kmers,
    /// Index canonical minimizers selected from each k-mer window.
    Minimizers,
    /// Index hash-sampled canonical k-mers using a FracMinHash-style rate.
    Modimizers,
    /// Index all (k-z)-mers and reconstruct k-mer membership with findere.
    Findere,
}

impl IndexMode {
    pub fn as_u8(self) -> u8 {
        match self {
            Self::Kmers => 0,
            Self::Minimizers => 1,
            Self::Modimizers => 2,
            Self::Findere => 3,
        }
    }

    pub fn from_u8(value: u8) -> Result<Self> {
        match value {
            0 => Ok(Self::Kmers),
            1 => Ok(Self::Minimizers),
            2 => Ok(Self::Modimizers),
            3 => Ok(Self::Findere),
            _ => anyhow::bail!("unsupported index mode {}", value),
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Kmers => "kmers",
            Self::Minimizers => "minimizers",
            Self::Modimizers => "modimizers",
            Self::Findere => "findere",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FeatureConfig {
    pub mode: IndexMode,
    pub minimizer_size: usize,
    pub modimizer_sampling: u64,
    pub findere_z: usize,
}

impl FeatureConfig {
    pub fn new(mode: IndexMode, minimizer_size: usize, modimizer_sampling: u64) -> Self {
        Self {
            mode,
            minimizer_size,
            modimizer_sampling,
            findere_z: findere::DEFAULT_Z,
        }
    }

    pub fn legacy_kmers() -> Self {
        Self {
            mode: IndexMode::Kmers,
            minimizer_size: DEFAULT_MINIMIZER_SIZE,
            modimizer_sampling: DEFAULT_MODIMIZER_SAMPLING,
            findere_z: findere::DEFAULT_Z,
        }
    }

    pub fn with_findere_z(mut self, z: usize) -> Self {
        self.findere_z = z;
        self
    }
}

#[derive(Clone, Debug)]
pub struct ColorInput {
    pub name: String,
    pub path: PathBuf,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct FastxStats {
    pub bases: u64,
    pub chunks: u64,
}

pub fn read_fof(path: &Path) -> Result<Vec<ColorInput>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut inputs = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields = trimmed.split('\t').collect::<Vec<_>>();
        let (name, input_path) = if fields.len() >= 2 {
            (
                fields[0].trim().to_string(),
                PathBuf::from(fields[1].trim()),
            )
        } else {
            let parts = trimmed.split_whitespace().collect::<Vec<_>>();
            if parts.len() >= 2 {
                (parts[0].to_string(), PathBuf::from(parts[1]))
            } else {
                let input_path = PathBuf::from(trimmed);
                (color_name_from_path(&input_path), input_path)
            }
        };
        ensure!(
            !name.is_empty(),
            "{}:{} has an empty color name",
            path.display(),
            line_no + 1
        );
        inputs.push(ColorInput {
            name,
            path: input_path,
        });
    }
    Ok(inputs)
}

pub fn color_name_from_path(path: &Path) -> String {
    let mut name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("color")
        .to_string();
    for ext in [
        ".zstd", ".zst", ".gzip", ".gz", ".fasta", ".fastq", ".fna", ".fa", ".fq", ".fas",
    ] {
        if let Some(stripped) = name.strip_suffix(ext) {
            name = stripped.to_string();
        }
    }
    if name.is_empty() {
        "color".to_string()
    } else {
        name
    }
}

pub fn create_index_writer(
    path: &Path,
) -> Result<zstd::stream::write::Encoder<'static, BufWriter<File>>> {
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    zstd::stream::Encoder::new(BufWriter::new(file), DEFAULT_INDEX_ZSTD_LEVEL)
        .with_context(|| format!("creating zstd encoder for {}", path.display()))
}

pub fn finish_index_writer(
    writer: zstd::stream::write::Encoder<'static, BufWriter<File>>,
) -> Result<()> {
    let mut inner = writer.finish().context("finishing zstd index stream")?;
    inner.flush().context("flushing index file")?;
    Ok(())
}

pub fn open_index_reader(path: &Path) -> Result<Box<dyn Read>> {
    let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut magic = [0u8; 4];
    let bytes = file
        .read(&mut magic)
        .with_context(|| format!("reading {}", path.display()))?;
    file.seek(SeekFrom::Start(0))
        .with_context(|| format!("seeking {}", path.display()))?;
    if bytes == magic.len() && magic == ZSTD_FRAME_MAGIC {
        let decoder = zstd::stream::Decoder::new(BufReader::new(file))
            .with_context(|| format!("opening zstd-compressed index {}", path.display()))?;
        Ok(Box::new(decoder))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

pub fn parse_features_from_file<F>(
    path: &Path,
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
    on_feature: F,
) -> Result<FastxStats>
where
    F: FnMut(u64),
{
    if is_compressed_path(path) {
        let file =
            File::open(path).with_context(|| format!("unable to open {}", path.display()))?;
        let mut decoder = AnyDecoder::new(file);
        let mut data = Vec::new();
        decoder
            .read_to_end(&mut data)
            .with_context(|| format!("decompressing {}", path.display()))?;
        ensure!(!data.is_empty(), "{} is empty", path.display());
        let parser = FastxParser::<PARSER_CONFIG>::from_slice(&data)
            .with_context(|| format!("parsing {}", path.display()))?;
        return parse_features_from_parser(parser, k, seed, feature_config, on_feature);
    }

    let parser = FastxParser::<PARSER_CONFIG>::from_file_mmap(path)
        .with_context(|| format!("unable to mmap {}", path.display()))?;
    parse_features_from_parser(parser, k, seed, feature_config, on_feature)
}

pub fn parse_features_from_parser<F>(
    mut parser: FastxParser<'_, PARSER_CONFIG>,
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
    mut on_feature: F,
) -> Result<FastxStats>
where
    F: FnMut(u64),
{
    let mut stats = FastxStats::default();
    while parser.next().is_some() {
        let seq = parser.get_dna_string();
        stats.bases += seq.len() as u64;
        stats.chunks += 1;
        scan_index_features(seq, k, seed, feature_config, &mut on_feature);
    }
    Ok(stats)
}

pub fn is_compressed_path(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .map(|extension| {
            matches!(
                extension.to_ascii_lowercase().as_str(),
                "gz" | "gzip" | "zst" | "zstd"
            )
        })
        .unwrap_or(false)
}

#[inline(always)]
pub fn scan_index_features<F>(
    seq: &[u8],
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
    on_feature: &mut F,
) where
    F: FnMut(u64),
{
    match feature_config.mode {
        IndexMode::Kmers => scan_canonical_kmers(seq, k, on_feature),
        IndexMode::Minimizers => {
            scan_simd_minimizers(seq, k, feature_config.minimizer_size, seed, on_feature)
        }
        IndexMode::Modimizers => {
            scan_modimizers(seq, k, seed, feature_config.modimizer_sampling, on_feature)
        }
        IndexMode::Findere => {
            if let Ok(s) = findere::indexed_length(k, feature_config.findere_z) {
                for run in seq.split(|&b| encode_base(b).is_none()) {
                    if run.len() >= k {
                        scan_canonical_kmers(run, s, on_feature);
                    }
                }
            }
        }
    }
}

#[inline(always)]
pub fn scan_canonical_kmers<F>(seq: &[u8], k: usize, on_kmer: &mut F)
where
    F: FnMut(u64),
{
    if seq.len() < k {
        return;
    }
    let mask = (1u64 << (2 * k)) - 1;
    let rc_shift = 2 * (k - 1);
    let mut forward = 0u64;
    let mut reverse = 0u64;
    let mut seen = 0usize;

    for &base in seq {
        let Some(code) = encode_base(base) else {
            forward = 0;
            reverse = 0;
            seen = 0;
            continue;
        };
        let code = code as u64;
        forward = ((forward << 2) | code) & mask;
        reverse = (reverse >> 2) | ((code ^ 0b11) << rc_shift);
        seen += 1;
        if seen >= k {
            on_kmer(forward.min(reverse));
        }
    }
}

pub fn scan_simd_minimizers<F>(
    seq: &[u8],
    k: usize,
    minimizer_size: usize,
    seed: u64,
    on_kmer: &mut F,
) where
    F: FnMut(u64),
{
    if seq.len() < k || minimizer_size == 0 || minimizer_size > k {
        return;
    }
    let window = k - minimizer_size + 1;
    let packed = PackedSeqVec::from_ascii(seq);
    let minimizer_seed = (seed ^ (seed >> 32)) as u32;
    let hasher =
        <simd_minimizers::seq_hash::NtHasher>::new_with_seed(minimizer_size, minimizer_seed);
    let mut positions = Vec::new();
    let _ = simd_minimizers::canonical_minimizers(minimizer_size, window)
        .hasher(&hasher)
        .run(packed.as_slice(), &mut positions);

    for position in positions {
        let position = position as usize;
        let Some(end) = position.checked_add(minimizer_size) else {
            continue;
        };
        if end > seq.len() {
            continue;
        }
        if let Some(encoded) = encode_full_kmer(&seq[position..end]) {
            on_kmer(canonical_encoded(encoded, minimizer_size));
        }
    }
}

#[inline(always)]
pub fn scan_modimizers<F>(seq: &[u8], k: usize, seed: u64, sampling: u64, on_kmer: &mut F)
where
    F: FnMut(u64),
{
    debug_assert!(sampling > 0);
    scan_canonical_kmers(seq, k, &mut |kmer| {
        if mixsplit(kmer, seed) % sampling == 0 {
            on_kmer(kmer);
        }
    });
}

#[inline(always)]
pub fn encode_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

pub fn encode_full_kmer(seq: &[u8]) -> Option<u64> {
    let mut value = 0u64;
    for &base in seq {
        value = (value << 2) | encode_base(base)? as u64;
    }
    Some(value)
}

pub fn canonical_encoded(encoded: u64, k: usize) -> u64 {
    let mut value = encoded;
    let mut rc = 0u64;
    for _ in 0..k {
        let code = value & 0b11;
        rc = (rc << 2) | (code ^ 0b11);
        value >>= 2;
    }
    encoded.min(rc)
}

#[inline(always)]
pub fn mixsplit(key: u64, seed: u64) -> u64 {
    splitmix64(key.wrapping_add(seed))
}

#[inline(always)]
pub fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn collect_kmers(seq: &[u8], k: usize) -> Vec<u64> {
        let mut out = Vec::new();
        scan_canonical_kmers(seq, k, &mut |kmer| out.push(kmer));
        out
    }

    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => *base,
            })
            .collect()
    }

    fn temp_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "zorindex_{name}_{}_{}",
            std::process::id(),
            splitmix64(name.len() as u64)
        ))
    }

    macro_rules! base_case {
        ($name:ident, $base:expr, $expected:expr) => {
            #[test]
            fn $name() {
                assert_eq!(encode_base($base), $expected);
            }
        };
    }

    base_case!(base_upper_a, b'A', Some(0));
    base_case!(base_upper_c, b'C', Some(1));
    base_case!(base_upper_g, b'G', Some(2));
    base_case!(base_upper_t, b'T', Some(3));
    base_case!(base_lower_a, b'a', Some(0));
    base_case!(base_lower_c, b'c', Some(1));
    base_case!(base_lower_g, b'g', Some(2));
    base_case!(base_lower_t, b't', Some(3));
    base_case!(base_n_is_invalid, b'N', None);
    base_case!(base_dash_is_invalid, b'-', None);
    base_case!(base_newline_is_invalid, b'\n', None);
    base_case!(base_space_is_invalid, b' ', None);

    macro_rules! color_name_case {
        ($name:ident, $path:expr, $expected:expr) => {
            #[test]
            fn $name() {
                assert_eq!(color_name_from_path(Path::new($path)), $expected);
            }
        };
    }

    color_name_case!(color_name_fa, "sample.fa", "sample");
    color_name_case!(color_name_fasta, "sample.fasta", "sample");
    color_name_case!(color_name_fna, "sample.fna", "sample");
    color_name_case!(color_name_fastq, "sample.fastq", "sample");
    color_name_case!(color_name_fq, "sample.fq", "sample");
    color_name_case!(color_name_fas, "sample.fas", "sample");
    color_name_case!(color_name_gz, "sample.fa.gz", "sample");
    color_name_case!(color_name_zst, "sample.fna.zst", "sample");
    color_name_case!(color_name_zstd, "sample.fastq.zstd", "sample");
    color_name_case!(color_name_nested, "/data/run/sample.fq.gz", "sample");
    color_name_case!(color_name_no_extension, "sample", "sample");
    color_name_case!(color_name_empty, "", "color");

    macro_rules! compressed_path_case {
        ($name:ident, $path:expr, $expected:expr) => {
            #[test]
            fn $name() {
                assert_eq!(is_compressed_path(Path::new($path)), $expected);
            }
        };
    }

    compressed_path_case!(compressed_gz, "a.fa.gz", true);
    compressed_path_case!(compressed_gzip, "a.fa.gzip", true);
    compressed_path_case!(compressed_zst, "a.fa.zst", true);
    compressed_path_case!(compressed_zstd, "a.fa.zstd", true);
    compressed_path_case!(compressed_uppercase_gz, "a.fa.GZ", true);
    compressed_path_case!(compressed_plain_fa, "a.fa", false);
    compressed_path_case!(compressed_plain_fastq, "a.fastq", false);
    compressed_path_case!(compressed_no_extension, "a", false);

    macro_rules! full_kmer_case {
        ($name:ident, $seq:expr, $expected:expr) => {
            #[test]
            fn $name() {
                assert_eq!(encode_full_kmer($seq), $expected);
            }
        };
    }

    full_kmer_case!(encode_kmer_a, b"A", Some(0));
    full_kmer_case!(encode_kmer_c, b"C", Some(1));
    full_kmer_case!(encode_kmer_g, b"G", Some(2));
    full_kmer_case!(encode_kmer_t, b"T", Some(3));
    full_kmer_case!(encode_kmer_ac, b"AC", Some(1));
    full_kmer_case!(encode_kmer_gt, b"GT", Some(11));
    full_kmer_case!(encode_kmer_acgt, b"ACGT", Some(27));
    full_kmer_case!(encode_kmer_lowercase, b"acgt", Some(27));
    full_kmer_case!(encode_kmer_invalid_n, b"AN", None);
    full_kmer_case!(encode_kmer_invalid_gap, b"A-", None);

    macro_rules! canonical_count_case {
        ($name:ident, $seq:expr, $k:expr, $expected:expr) => {
            #[test]
            fn $name() {
                assert_eq!(collect_kmers($seq, $k).len(), $expected);
            }
        };
    }

    canonical_count_case!(kmers_count_short, b"ACG", 4, 0);
    canonical_count_case!(kmers_count_exact_1, b"ACGT", 4, 1);
    canonical_count_case!(kmers_count_exact_2, b"AAAA", 4, 1);
    canonical_count_case!(kmers_count_five_k3, b"ACGTA", 3, 3);
    canonical_count_case!(kmers_count_eight_k4, b"ACGTACGT", 4, 5);
    canonical_count_case!(kmers_count_ten_k5, b"ACGTACGTAA", 5, 6);
    canonical_count_case!(kmers_count_n_resets_left, b"NACGT", 4, 1);
    canonical_count_case!(kmers_count_n_resets_middle, b"ACNGTAC", 3, 2);
    canonical_count_case!(kmers_count_n_resets_many, b"ACGTNNACGT", 4, 2);
    canonical_count_case!(kmers_count_lowercase, b"acgtac", 4, 3);
    canonical_count_case!(kmers_count_mixed_case, b"AcGtAc", 4, 3);
    canonical_count_case!(kmers_count_k1, b"ACGT", 1, 4);
    canonical_count_case!(kmers_count_k2, b"ACGT", 2, 3);
    canonical_count_case!(kmers_count_k3, b"ACGT", 3, 2);
    canonical_count_case!(kmers_count_k4, b"ACGT", 4, 1);
    canonical_count_case!(kmers_count_homopolymer, b"AAAAAAAA", 5, 4);

    macro_rules! reverse_complement_case {
        ($name:ident, $seq:expr, $k:expr) => {
            #[test]
            fn $name() {
                let mut left = collect_kmers($seq, $k);
                let rc = reverse_complement($seq);
                let mut right = collect_kmers(&rc, $k);
                left.sort_unstable();
                right.sort_unstable();
                assert_eq!(left, right);
            }
        };
    }

    reverse_complement_case!(rc_equivalence_01, b"ACGTACGT", 3);
    reverse_complement_case!(rc_equivalence_02, b"AAAACCCCGGGGTTTT", 4);
    reverse_complement_case!(rc_equivalence_03, b"TGCATGCATGCA", 5);
    reverse_complement_case!(rc_equivalence_04, b"ACACACACAC", 2);
    reverse_complement_case!(rc_equivalence_05, b"GGGGAAAATTTT", 6);
    reverse_complement_case!(rc_equivalence_06, b"TTTTCCCCAAAA", 3);
    reverse_complement_case!(rc_equivalence_07, b"ACGTTGCAACGT", 7);
    reverse_complement_case!(rc_equivalence_08, b"ATATCGCGATAT", 4);
    reverse_complement_case!(rc_equivalence_09, b"GATTACAGATTACA", 5);
    reverse_complement_case!(rc_equivalence_10, b"CCGGAATTCCGG", 6);
    reverse_complement_case!(rc_equivalence_11, b"ACGTNACGT", 4);
    reverse_complement_case!(rc_equivalence_12, b"AAAANNNNTTTT", 4);

    macro_rules! canonical_value_case {
        ($name:ident, $seq:expr) => {
            #[test]
            fn $name() {
                let encoded = encode_full_kmer($seq).unwrap();
                let rc = reverse_complement($seq);
                let rc_encoded = encode_full_kmer(&rc).unwrap();
                assert_eq!(
                    canonical_encoded(encoded, $seq.len()),
                    encoded.min(rc_encoded)
                );
            }
        };
    }

    canonical_value_case!(canonical_value_01, b"A");
    canonical_value_case!(canonical_value_02, b"C");
    canonical_value_case!(canonical_value_03, b"G");
    canonical_value_case!(canonical_value_04, b"T");
    canonical_value_case!(canonical_value_05, b"AC");
    canonical_value_case!(canonical_value_06, b"GT");
    canonical_value_case!(canonical_value_07, b"ACG");
    canonical_value_case!(canonical_value_08, b"CGT");
    canonical_value_case!(canonical_value_09, b"ACGT");
    canonical_value_case!(canonical_value_10, b"GATTACA");
    canonical_value_case!(canonical_value_11, b"CCCCAAAA");
    canonical_value_case!(canonical_value_12, b"TTTTGGGG");

    macro_rules! mode_roundtrip_case {
        ($name:ident, $mode:expr, $byte:expr, $label:expr) => {
            #[test]
            fn $name() {
                assert_eq!($mode.as_u8(), $byte);
                assert_eq!(IndexMode::from_u8($byte).unwrap(), $mode);
                assert_eq!($mode.label(), $label);
            }
        };
    }

    mode_roundtrip_case!(mode_kmers_roundtrip, IndexMode::Kmers, 0, "kmers");
    mode_roundtrip_case!(mode_findere_roundtrip, IndexMode::Findere, 3, "findere");
    mode_roundtrip_case!(
        mode_minimizers_roundtrip,
        IndexMode::Minimizers,
        1,
        "minimizers"
    );
    mode_roundtrip_case!(
        mode_modimizers_roundtrip,
        IndexMode::Modimizers,
        2,
        "modimizers"
    );

    #[test]
    fn mode_rejects_unknown_byte() {
        assert!(IndexMode::from_u8(99).is_err());
    }

    macro_rules! modimizer_sampling_one_case {
        ($name:ident, $seq:expr, $k:expr, $seed:expr) => {
            #[test]
            fn $name() {
                let mut all = Vec::new();
                scan_canonical_kmers($seq, $k, &mut |kmer| all.push(kmer));
                let mut sampled = Vec::new();
                scan_index_features(
                    $seq,
                    $k,
                    $seed,
                    FeatureConfig::new(IndexMode::Modimizers, DEFAULT_MINIMIZER_SIZE, 1),
                    &mut |kmer| sampled.push(kmer),
                );
                assert_eq!(sampled, all);
            }
        };
    }

    modimizer_sampling_one_case!(modimizer_sampling_one_01, b"ACGTACGT", 3, 1);
    modimizer_sampling_one_case!(modimizer_sampling_one_02, b"AAAACCCC", 4, 2);
    modimizer_sampling_one_case!(modimizer_sampling_one_03, b"TTTTGGGG", 5, 3);
    modimizer_sampling_one_case!(modimizer_sampling_one_04, b"GATTACA", 3, 4);
    modimizer_sampling_one_case!(modimizer_sampling_one_05, b"ACGTNACGT", 4, 5);
    modimizer_sampling_one_case!(modimizer_sampling_one_06, b"CCCCCCCC", 2, 6);

    macro_rules! minimizer_subset_case {
        ($name:ident, $seq:expr, $k:expr, $m:expr, $seed:expr) => {
            #[test]
            fn $name() {
                let mut emitted = Vec::new();
                scan_index_features(
                    $seq,
                    $k,
                    $seed,
                    FeatureConfig::new(IndexMode::Minimizers, $m, DEFAULT_MODIMIZER_SAMPLING),
                    &mut |kmer| emitted.push(kmer),
                );
                let mut candidates = Vec::new();
                scan_canonical_kmers($seq, $m, &mut |kmer| candidates.push(kmer));
                for kmer in emitted {
                    assert!(candidates.contains(&kmer));
                }
            }
        };
    }

    minimizer_subset_case!(minimizer_subset_01, b"ACGTACGT", 5, 3, 11);
    minimizer_subset_case!(minimizer_subset_02, b"AAAACCCCGGGG", 7, 4, 12);
    minimizer_subset_case!(minimizer_subset_03, b"TTTTGGGGCCCC", 7, 3, 13);
    minimizer_subset_case!(minimizer_subset_04, b"GATTACAGATTACA", 9, 5, 14);
    minimizer_subset_case!(minimizer_subset_05, b"ACGTACGTACGT", 9, 4, 15);
    minimizer_subset_case!(minimizer_subset_06, b"CCCCAAAATTTT", 7, 2, 16);

    macro_rules! splitmix_stability_case {
        ($name:ident, $input:expr, $expected:expr) => {
            #[test]
            fn $name() {
                assert_eq!(splitmix64($input), $expected);
            }
        };
    }

    splitmix_stability_case!(splitmix_stable_00, 0, 0xe220_a839_7b1d_cdaF);
    splitmix_stability_case!(splitmix_stable_01, 1, 0x910a_2dec_8902_5cc1);
    splitmix_stability_case!(splitmix_stable_02, 2, 0x9758_35de_1c97_56ce);
    splitmix_stability_case!(splitmix_stable_03, 3, 0x1d0b_14e4_db01_8fed);
    splitmix_stability_case!(splitmix_stable_04, 4, 0x6e73_e372_e233_8aca);
    splitmix_stability_case!(splitmix_stable_05, 5, 0x6303_3b0c_a389_c35a);
    splitmix_stability_case!(splitmix_stable_06, 6, 0xbd64_a5d9_adef_e000);
    splitmix_stability_case!(splitmix_stable_07, 7, 0x63cbe1e459320dd7);
    splitmix_stability_case!(splitmix_stable_08, 8, 0x9e5651b0ef953636);
    splitmix_stability_case!(splitmix_stable_09, 9, 0xaeaf52febe706064);

    macro_rules! mixsplit_relation_case {
        ($name:ident, $key:expr, $seed:expr) => {
            #[test]
            fn $name() {
                let key: u64 = $key;
                let seed: u64 = $seed;
                assert_eq!(mixsplit(key, seed), splitmix64(key.wrapping_add(seed)));
            }
        };
    }

    mixsplit_relation_case!(mixsplit_relation_01, 0, 0);
    mixsplit_relation_case!(mixsplit_relation_02, 1, 2);
    mixsplit_relation_case!(mixsplit_relation_03, 17, 99);
    mixsplit_relation_case!(mixsplit_relation_04, u64::MAX, 1);
    mixsplit_relation_case!(mixsplit_relation_05, 0x1234_5678, 0x9abc_def0);
    mixsplit_relation_case!(mixsplit_relation_06, 0xfeed_face_cafe_beef, 69);
    mixsplit_relation_case!(mixsplit_relation_07, 42, u64::MAX);
    mixsplit_relation_case!(mixsplit_relation_08, 123_456_789, 987_654_321);

    #[test]
    fn index_writer_writes_zstd_stream() {
        let path = temp_path("zstd_stream");
        {
            let mut writer = create_index_writer(&path).unwrap();
            writer.write_all(b"index-payload").unwrap();
            finish_index_writer(writer).unwrap();
        }

        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.starts_with(&ZSTD_FRAME_MAGIC));

        let mut reader = open_index_reader(&path).unwrap();
        let mut decoded = Vec::new();
        reader.read_to_end(&mut decoded).unwrap();
        assert_eq!(decoded, b"index-payload");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn index_reader_accepts_raw_stream() {
        let path = temp_path("raw_stream");
        std::fs::write(&path, b"raw-index-payload").unwrap();

        let mut reader = open_index_reader(&path).unwrap();
        let mut decoded = Vec::new();
        reader.read_to_end(&mut decoded).unwrap();
        assert_eq!(decoded, b"raw-index-payload");
        let _ = std::fs::remove_file(path);
    }
}
