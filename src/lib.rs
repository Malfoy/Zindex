use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use anyhow::{ensure, Context, Result};
use clap::ValueEnum;
use deko::read::AnyDecoder;
use helicase::input::{FromMmap, FromReader};
use helicase::{Config, FastxParser, HelicaseParser, ParserOptions};
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};

pub mod anchors;
pub mod findere;
pub mod open_closed_minimizers;
pub mod packed_fingerprints;
pub mod power_merge;

pub use anchors::{RandstrobeConfig, SpacedSeed};

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

#[derive(
    Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "kebab-case")]
pub enum MinimizerScheme {
    /// Standard canonical SIMD hash minimizers (also used by legacy indexes).
    #[default]
    Simd,
    /// Canonical open-closed mod minimizers, with a syncmer-prioritized anchor.
    OpenClosedMod,
}

impl MinimizerScheme {
    pub fn as_u8(self) -> u8 {
        match self {
            Self::Simd => 0,
            Self::OpenClosedMod => 1,
        }
    }

    pub fn from_u8(value: u8) -> Result<Self> {
        match value {
            0 => Ok(Self::Simd),
            1 => Ok(Self::OpenClosedMod),
            _ => anyhow::bail!("unsupported minimizer scheme {}", value),
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Simd => "simd",
            Self::OpenClosedMod => "open-closed-mod",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FeatureConfig {
    pub mode: IndexMode,
    pub minimizer_size: usize,
    pub minimizer_scheme: MinimizerScheme,
    pub homopolymer_compression: bool,
    pub seed_pattern: Option<SpacedSeed>,
    pub randstrobemers: Option<RandstrobeConfig>,
    pub modimizer_sampling: u64,
    pub findere_z: usize,
}

impl FeatureConfig {
    pub fn new(mode: IndexMode, minimizer_size: usize, modimizer_sampling: u64) -> Self {
        Self {
            mode,
            minimizer_size,
            minimizer_scheme: MinimizerScheme::Simd,
            homopolymer_compression: false,
            seed_pattern: None,
            randstrobemers: None,
            modimizer_sampling,
            findere_z: findere::DEFAULT_Z,
        }
    }

    pub fn legacy_kmers() -> Self {
        Self {
            mode: IndexMode::Kmers,
            minimizer_size: DEFAULT_MINIMIZER_SIZE,
            minimizer_scheme: MinimizerScheme::Simd,
            homopolymer_compression: false,
            seed_pattern: None,
            randstrobemers: None,
            modimizer_sampling: DEFAULT_MODIMIZER_SAMPLING,
            findere_z: findere::DEFAULT_Z,
        }
    }

    pub fn with_findere_z(mut self, z: usize) -> Self {
        self.findere_z = z;
        self
    }

    pub fn with_minimizer_scheme(mut self, scheme: MinimizerScheme) -> Self {
        self.minimizer_scheme = scheme;
        self
    }

    pub fn with_homopolymer_compression(mut self, enabled: bool) -> Self {
        self.homopolymer_compression = enabled;
        self
    }

    pub fn with_seed_pattern(mut self, pattern: Option<SpacedSeed>) -> Self {
        self.seed_pattern = pattern;
        self
    }

    pub fn with_randstrobemers(mut self, config: Option<RandstrobeConfig>) -> Self {
        self.randstrobemers = config;
        self
    }

    pub fn validate_anchors(self, k: usize) -> Result<()> {
        if self.seed_pattern.is_some() || self.randstrobemers.is_some() {
            ensure!(
                matches!(self.mode, IndexMode::Kmers | IndexMode::Minimizers),
                "spaced seeds and randstrobemers require --index-mode kmers or minimizers"
            );
        }
        let length = if self.mode == IndexMode::Minimizers {
            self.minimizer_size
        } else {
            k
        };
        if let Some(pattern) = self.seed_pattern {
            ensure!(
                pattern.width() == length,
                "seed pattern width {} must match {} ({})",
                pattern.width(),
                if self.mode == IndexMode::Minimizers {
                    "minimizer size m"
                } else {
                    "k"
                },
                length
            );
        }
        if let Some(config) = self.randstrobemers {
            config.validate()?;
        }
        Ok(())
    }
}

/// Collapse runs of identical DNA bases (case-insensitively). Ambiguous bases
/// remain separators: compression never creates a feature across such a break.
pub fn homopolymer_compress(seq: &[u8]) -> Vec<u8> {
    let mut result = Vec::with_capacity(seq.len());
    let mut previous = None;
    for &base in seq {
        let code = encode_base(base);
        if code.is_none() || code != previous {
            result.push(code.map_or(base, |c| b"ACGT"[c as usize]));
        }
        previous = code;
    }
    result
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

// Helicase's streaming reader panics on Read errors. Decode here and retain the
// original error while presenting EOF to the parser, then return it as Result.
// This also checks the compression trailer instead of accepting a damaged tail.
struct CheckedSequenceReader<R: Read> {
    decoder: AnyDecoder<R>,
    error: Arc<Mutex<Option<std::io::Error>>>,
}
impl<R: Read> Read for CheckedSequenceReader<R> {
    fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
        if self.error.lock().unwrap().is_some() {
            return Ok(0);
        }
        loop {
            match self.decoder.read(output) {
                Ok(n) => return Ok(n),
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(error) => {
                    *self.error.lock().unwrap() = Some(error);
                    return Ok(0);
                }
            }
        }
    }
}

fn with_streaming_parser<R: Read + Send + 'static, T>(
    reader: R,
    visit: impl FnOnce(FastxParser<'static, PARSER_CONFIG>) -> Result<T>,
) -> Result<T> {
    let error = Arc::new(Mutex::new(None));
    let parser = FastxParser::from_reader(CheckedSequenceReader {
        decoder: AnyDecoder::new(reader),
        error: Arc::clone(&error),
    });
    let result = parser.map_err(anyhow::Error::from).and_then(visit);
    if let Some(error) = error.lock().unwrap().take() {
        return Err(error).context("reading/decompressing sequence stream");
    }
    result
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
        ".zstd", ".zst", ".gzip", ".gz", ".xz", ".bz2", ".bz", ".bzip2", ".fasta", ".fastq",
        ".fna", ".fa", ".fq", ".fas",
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
        return parse_features_from_reader(file, k, seed, feature_config, on_feature)
            .with_context(|| format!("parsing compressed {}", path.display()));
    }

    let parser = FastxParser::<PARSER_CONFIG>::from_file_mmap(path)
        .with_context(|| format!("unable to mmap {}", path.display()))?;
    parse_features_from_parser(parser, k, seed, feature_config, on_feature)
}

/// Stream compressed or plain FASTA/FASTQ without retaining the entire input.
/// Record/ambiguity boundaries are handled by the same parser as mmap input;
/// feature vectors and sorting remain the caller's in-memory responsibility.
pub fn parse_features_from_reader<R: Read + Send + 'static, F: FnMut(u64)>(
    reader: R,
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
    on_feature: F,
) -> Result<FastxStats> {
    with_streaming_parser(reader, |parser| {
        parse_features_from_parser(parser, k, seed, feature_config, on_feature)
    })
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

/// Visit complete FASTA/FASTQ records, including compressed input. Used when
/// query scoring must keep the two strobemer orientations separate per read.
pub fn parse_sequences_from_file(
    path: &Path,
    mut on_sequence: impl FnMut(&[u8]),
) -> Result<FastxStats> {
    if is_compressed_path(path) {
        return with_streaming_parser(File::open(path)?, |parser| {
            visit_sequence_parser(parser, &mut on_sequence)
        })
        .with_context(|| format!("parsing compressed {}", path.display()));
    }
    visit_sequence_parser(
        FastxParser::<PARSER_CONFIG>::from_file_mmap(path)?,
        &mut on_sequence,
    )
}

fn visit_sequence_parser(
    mut parser: FastxParser<'_, PARSER_CONFIG>,
    on_sequence: &mut impl FnMut(&[u8]),
) -> Result<FastxStats> {
    let mut stats = FastxStats::default();
    while parser.next().is_some() {
        let seq = parser.get_dna_string();
        stats.bases += seq.len() as u64;
        stats.chunks += 1;
        on_sequence(seq);
    }
    Ok(stats)
}

pub fn is_compressed_path(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .map(|extension| {
            matches!(
                extension.to_ascii_lowercase().as_str(),
                "gz" | "gzip" | "zst" | "zstd" | "xz" | "bz2" | "bz" | "bzip2"
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
    if feature_config.homopolymer_compression {
        let compressed = homopolymer_compress(seq);
        scan_features_uncompressed(&compressed, k, seed, feature_config, on_feature);
    } else {
        scan_features_uncompressed(seq, k, seed, feature_config, on_feature);
    }
}

/// Reverse-complement a query without changing ambiguous-base separators.
pub use anchors::{reverse_complement, RANDSTROBE_SCHEME};

fn scan_features_uncompressed<F: FnMut(u64)>(
    seq: &[u8],
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
    on_feature: &mut F,
) {
    if feature_config.seed_pattern.is_some() || feature_config.randstrobemers.is_some() {
        anchors::scan(seq, k, seed, feature_config, on_feature);
        return;
    }
    scan_plain_features_uncompressed(seq, k, seed, feature_config, on_feature);
}

pub(crate) fn scan_plain_features_uncompressed<F: FnMut(u64)>(
    seq: &[u8],
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
    on_feature: &mut F,
) {
    match feature_config.mode {
        IndexMode::Kmers => scan_canonical_kmers(seq, k, on_feature),
        IndexMode::Minimizers => match feature_config.minimizer_scheme {
            MinimizerScheme::Simd => {
                scan_simd_minimizers(seq, k, feature_config.minimizer_size, seed, on_feature)
            }
            MinimizerScheme::OpenClosedMod => open_closed_minimizers::scan(
                seq,
                k,
                feature_config.minimizer_size,
                seed,
                on_feature,
            ),
        },
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
    scan_simd_minimizers_positions(seq, k, minimizer_size, seed, &mut |_, key| on_kmer(key));
}

pub fn scan_simd_minimizers_positions<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: usize,
    minimizer_size: usize,
    seed: u64,
    emit: &mut F,
) {
    let mut offset = 0;
    for run in seq.split(|&b| encode_base(b).is_none()) {
        scan_simd_minimizers_run(run, k, minimizer_size, seed, &mut |position, key| {
            emit(offset + position, key)
        });
        offset += run.len() + 1;
    }
}

fn scan_simd_minimizers_run<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: usize,
    minimizer_size: usize,
    seed: u64,
    on_kmer: &mut F,
) {
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
            on_kmer(position, canonical_encoded(encoded, minimizer_size));
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

pub use ocm_minimizers::{canonical_encoded, encode_base, encode_full_kmer, mixsplit, splitmix64};

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

    compressed_path_case!(compressed_xz, "a.fa.xz", true);
    compressed_path_case!(compressed_bz2, "a.fa.bz2", true);
    compressed_path_case!(compressed_bz, "a.fa.bz", true);
    compressed_path_case!(compressed_bzip2, "a.fa.bzip2", true);
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

    #[test]
    fn streaming_features_match_slice_across_refills_records_ambiguities_and_modes() {
        use deko::{
            write::{AnyEncoder, Compression},
            Format,
        };
        use helicase::input::FromSlice;
        use std::io::Cursor;
        // Long wrapped records cross multiple decoder/parser buffer refills.
        // Include CRLF, short records, lowercase, ambiguity and reverse strand.
        let dna: Vec<u8> = (0..350_017)
            .map(|i| b"ACGT"[(splitmix64(i) & 3) as usize])
            .collect();
        let mut fasta = b">long\r\n".to_vec();
        for line in dna.chunks(79) {
            fasta.extend_from_slice(line);
            fasta.extend_from_slice(b"\r\n");
        }
        fasta.extend_from_slice(b">short\nACGTNacgtNNAA\n>ambiguous\n");
        fasta.extend_from_slice(&dna[..100_019]);
        fasta.extend_from_slice(b"NNN");
        fasta.extend_from_slice(&reverse_complement(&dna[..85_021]));
        fasta.extend_from_slice(b"\n");
        let streams: Vec<_> = [
            Format::Verbatim,
            Format::Gz,
            Format::Zstd,
            Format::Xz,
            Format::Bz,
        ]
        .into_iter()
        .map(|format| {
            let mut encoder = AnyEncoder::new(Vec::new(), format, Compression::Default).unwrap();
            encoder.write_all(&fasta).unwrap();
            encoder.finish().unwrap()
        })
        .collect();
        for mode in [
            IndexMode::Kmers,
            IndexMode::Minimizers,
            IndexMode::Modimizers,
            IndexMode::Findere,
        ] {
            for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
                let config = FeatureConfig::new(mode, 21, 7).with_minimizer_scheme(scheme);
                let mut expected = Vec::new();
                let before = parse_features_from_parser(
                    FastxParser::from_slice(&fasta).unwrap(),
                    31,
                    69,
                    config,
                    |key| expected.push(key),
                )
                .unwrap();
                for stream in &streams {
                    let mut actual = Vec::new();
                    let after = parse_features_from_reader(
                        Cursor::new(stream.clone()),
                        31,
                        69,
                        config,
                        |key| actual.push(key),
                    )
                    .unwrap();
                    assert_eq!(actual, expected, "{mode:?}/{scheme:?}");
                    assert_eq!((before.bases, before.chunks), (after.bases, after.chunks));
                }
            }
        }
    }

    #[test]
    fn streaming_emits_before_eof_and_returns_reader_errors() {
        use std::io::{self, Cursor};
        use std::sync::atomic::{AtomicUsize, Ordering};
        let consumed = Arc::new(AtomicUsize::new(0));
        struct Reader {
            bytes: Cursor<Vec<u8>>,
            consumed: Arc<AtomicUsize>,
        }
        impl Read for Reader {
            fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
                let n = self.bytes.read(out)?;
                self.consumed.fetch_add(n, Ordering::SeqCst);
                Ok(n)
            }
        }
        let fasta = b">record\nACGTTGCAACGTTGCANN\n".repeat(100_000);
        let total = fasta.len();
        let mut early = false;
        parse_features_from_reader(
            Reader {
                bytes: Cursor::new(fasta),
                consumed: Arc::clone(&consumed),
            },
            7,
            69,
            FeatureConfig::legacy_kmers(),
            |_| {
                early |= consumed.load(Ordering::SeqCst) < total;
            },
        )
        .unwrap();
        assert!(
            early,
            "must emit features without decoding the entire file first"
        );
        struct Broken;
        impl Read for Broken {
            fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
                Err(io::Error::other("injected decoder failure"))
            }
        }
        let error =
            parse_features_from_reader(Broken, 7, 69, FeatureConfig::legacy_kmers(), |_| {})
                .unwrap_err();
        assert!(format!("{error:#}").contains("injected decoder failure"));
    }

    #[test]
    fn streaming_fastq_and_corrupt_gzip_tail() {
        use deko::{
            write::{AnyEncoder, Compression},
            Format,
        };
        use helicase::input::FromSlice;
        use std::io::Cursor;
        let dna: Vec<u8> = (0..180_101)
            .map(|i| b"ACGT"[(splitmix64(i ^ 42) & 3) as usize])
            .collect();
        let mut fastq = b"@long\n".to_vec();
        fastq.extend_from_slice(&dna);
        fastq.extend_from_slice(b"\n+\n");
        fastq.extend(std::iter::repeat_n(b'@', dna.len()));
        fastq.extend_from_slice(b"\n@short\nACGTNNacgtACGT\n+\nIIIIIIIIIIIIII\n");
        let config = FeatureConfig::new(IndexMode::Minimizers, 31, 1)
            .with_minimizer_scheme(MinimizerScheme::OpenClosedMod);
        let mut expected = Vec::new();
        parse_features_from_parser(
            FastxParser::from_slice(&fastq).unwrap(),
            49,
            69,
            config,
            |key| expected.push(key),
        )
        .unwrap();
        let mut encoder = AnyEncoder::new(Vec::new(), Format::Gz, Compression::Default).unwrap();
        encoder.write_all(&fastq).unwrap();
        let mut compressed = encoder.finish().unwrap();
        let mut actual = Vec::new();
        parse_features_from_reader(Cursor::new(compressed.clone()), 49, 69, config, |key| {
            actual.push(key)
        })
        .unwrap();
        assert_eq!(actual, expected);
        let last = compressed.len() - 1;
        compressed[last] ^= 0x80;
        assert!(
            parse_features_from_reader(Cursor::new(compressed), 49, 69, config, |_| {}).is_err()
        );
    }
}
