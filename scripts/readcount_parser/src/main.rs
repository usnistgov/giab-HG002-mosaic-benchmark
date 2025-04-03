use clap::Parser;
use csv::ReaderBuilder;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// Command-line arguments structure
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path to the file containing the list of bam-readcount output files
    #[arg(short, long)]
    bam_readcount_list: String,

    /// Path to the variant list file
    #[arg(short, long)]
    variants: String,

    /// Path to the output file
    #[arg(short, long)]
    output: String,
}

/// Reads the variant list from a file and stores it in a HashMap
fn read_variants(
    variant_file: &str,
) -> Result<HashMap<(String, u32), (String, String)>, Box<dyn Error>> {
    let file = File::open(variant_file)?;
    let mut rdr = ReaderBuilder::new()
        .has_headers(true) // Indicates that the file has headers
        .delimiter(b'\t')  // Specifies tab-delimited fields
        .from_reader(file);

    let mut variants = HashMap::new();

    for result in rdr.records() {
        let record = result?;
        let chr = record.get(0).ok_or("Missing chromosome field")?.to_string();
        let pos: u32 = record
            .get(1)
            .ok_or("Missing position field")?
            .parse()
            .map_err(|_| "Invalid position value")?;
        let ref_base = record.get(2).ok_or("Missing ref_base field")?.to_string();
        let alt_base = record.get(3).ok_or("Missing alt_base field")?.to_string();
        variants.insert((chr, pos), (ref_base, alt_base));
    }
    Ok(variants)
}

/// Extracts dna_source and platform from the file path
fn extract_metadata(file_path: &str) -> Option<(String, String)> {
    let path = Path::new(file_path);
    let components: Vec<_> = path.components().collect();
    if components.len() >= 3 {
        let dna_source = components[components.len() - 2]
            .as_os_str()
            .to_string_lossy()
            .to_string();
        let filename = path.file_name()?.to_string_lossy();
        let platform = filename.split("_HG002").next()?.to_string();
        Some((dna_source, platform))
    } else {
        None
    }
}

/// Processes a single bam-readcount output file and writes the filtered results to the output writer
fn process_bam_readcount_file(
    bam_readcount_file: &str,
    variants: &HashMap<(String, u32), (String, String)>,
    dna_source: &str,
    platform: &str,
    writer: &mut impl Write,
) -> Result<(), Box<dyn Error>> {
    let file = File::open(bam_readcount_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chr = fields[0].to_string();
        let pos: u32 = fields[1].trim().parse().map_err(|e| {
            eprintln!("Error parsing position '{}' in bam-readcount file: {}", fields[1], e);
            e
        })?;
        let ref_base = fields[2].to_string();

        if let Some((ref_base_var, alt_base)) = variants.get(&(chr.clone(), pos)) {
            if ref_base != *ref_base_var {
                eprintln!(
                    "Reference base mismatch at {}:{} (expected {}, found {})",
                    chr, pos, ref_base_var, ref_base
                );
                continue;
            }

            let mut base_counts = HashMap::new();
            let mut depth = 0;

            for base_info in &fields[4..] {
                let base_fields: Vec<&str> = base_info.split(':').collect();
                if base_fields.len() < 2 {
                    continue;
                }
                let base = base_fields[0].to_uppercase();
                let count: u32 = base_fields[1].trim().parse().map_err(|e| {
                    eprintln!("Error parsing base count '{}' for base '{}' at {}:{}: {}", base_fields[1], base_fields[0], chr, pos, e);
                    e
                })?;
                base_counts.insert(base.clone(), count);
                depth += count;
            }

            let ref_count = base_counts.get(&ref_base).unwrap_or(&0);
            let alt_count = base_counts.get(alt_base).unwrap_or(&0);

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                dna_source, platform, chr, pos, ref_base, alt_base, ref_count, alt_count, depth
            )?;
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let variants = read_variants(&args.variants)?;

    let bam_readcount_list_file = File::open(&args.bam_readcount_list)?;
    let bam_readcount_list_reader = BufReader::new(bam_readcount_list_file);
    let bam_readcount_files: Vec<String> = bam_readcount_list_reader
        .lines()
        .collect::<Result<_, _>>()?;

    let mut output_file = File::create(&args.output)?;
    writeln!(
        output_file,
        "dna_source\tplatform\tCHROM\tPOS\tREF\tALT\tref_count\talt_count\tdepth"
    )?;

    for bam_readcount_file in bam_readcount_files {
        if let Some((dna_source, platform)) = extract_metadata(&bam_readcount_file) {
            process_bam_readcount_file(
                &bam_readcount_file,
                &variants,
                &dna_source,
                &platform,
                &mut output_file,
            )?;
        } else {
            eprintln!(
                "Failed to extract metadata from file path: {}",
                bam_readcount_file
            );
        }
    }

    Ok(())
}
