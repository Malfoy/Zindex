use super::*;

/// Per-record query output. Both formats retain queries without reported hits.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum RecordFormat {
    /// One TSV row per query/dataset match, with dataset names and matched fractions.
    #[default]
    Tsv,
    /// One row per query, with comma-separated dataset_id:count hits.
    Compact,
}

impl RecordFormat {
    pub(super) fn header(self) -> &'static str {
        match self {
            Self::Tsv => "query_id\tdataset_id\tdataset_name\tmatches\tqueried_seeds\tratio",
            Self::Compact => "query_id\tvalid_kmers\thits",
        }
    }

    pub(super) fn write_record(
        self,
        out: &mut impl Write,
        name: &[u8],
        total: u64,
        hits: impl Iterator<Item = (usize, u64)>,
        colors: &[Color],
    ) -> Result<()> {
        match self {
            Self::Compact => {
                out.write_all(name)?;
                write!(out, "\t{total}\t")?;
                for (i, (origin, count)) in hits.enumerate() {
                    if i != 0 {
                        out.write_all(b",")?;
                    }
                    write!(out, "{origin}:{count}")?;
                }
                writeln!(out)?;
            }
            Self::Tsv => {
                let mut any = false;
                for (origin, count) in hits {
                    any = true;
                    write_field(out, name)?;
                    write!(out, "\t{origin}\t")?;
                    write_field(out, colors[origin].name.as_bytes())?;
                    let ratio = if total == 0 {
                        0.0
                    } else {
                        count as f64 / total as f64
                    };
                    writeln!(out, "\t{count}\t{total}\t{ratio:.8}")?;
                }
                if !any {
                    write_field(out, name)?;
                    writeln!(out, "\t.\t.\t0\t{total}\t0.00000000")?;
                }
            }
        }
        Ok(())
    }
}

// TSV fields use standard CSV quoting, so names containing quotes or whitespace
// remain readable with csv.DictReader(delimiter='\t') and equivalent tools.
fn write_field(out: &mut impl Write, value: &[u8]) -> Result<()> {
    if value
        .iter()
        .any(|b| matches!(b, b'\t' | b'\r' | b'\n' | b'"'))
    {
        out.write_all(b"\"")?;
        for part in value.split_inclusive(|&b| b == b'"') {
            out.write_all(part)?;
            if part.last() == Some(&b'"') {
                out.write_all(b"\"")?;
            }
        }
        out.write_all(b"\"")?;
    } else {
        out.write_all(value)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tsv_quotes_names_and_retains_unmatched_records() {
        let colors = vec![Color {
            name: "sample\twith \"quotes\"\nand newline".into(),
            path: "".into(),
            occurrences: 0,
            unique_kmers: 0,
        }];
        let mut output = Vec::new();
        RecordFormat::Tsv
            .write_record(&mut output, b"q", 10, std::iter::once((0, 8)), &colors)
            .unwrap();
        assert_eq!(
            output,
            b"q\t0\t\"sample\twith \"\"quotes\"\"\nand newline\"\t8\t10\t0.80000000\n"
        );
        output.clear();
        RecordFormat::Tsv
            .write_record(&mut output, b"empty", 0, std::iter::empty(), &colors)
            .unwrap();
        assert_eq!(output, b"empty\t.\t.\t0\t0\t0.00000000\n");
    }
}
