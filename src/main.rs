use anyhow::anyhow;
use anyhow::{Context, Result};
use clap::Parser;
use csv::ReaderBuilder;
use rust_htslib::{
    bam,
    bam::record::Aux,
    bam::{Format, Read, Record, Writer},
};
use seine::salmon::QuantRecord;
use std::cmp;
use std::fs::File;
use std::path::Path;
use std::vec::Vec;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    // Path to the SAM/BAM file whose alignments will be annotated
    #[clap(short, long)]
    alignments: Option<String>,

    // Path to the salmon quantification file
    #[clap(short, long)]
    quant: String,

    // Output path of modified alignments
    #[clap(short, long)]
    output: Option<String>,

    // Number of threads to use
    #[clap(short, long, default_value = "2")]
    num_threads: usize,
}

fn read_quants<P: AsRef<Path>>(p: P) -> Result<Vec<QuantRecord>, csv::Error> {
    // TODO error handling for csv::Error and io::Error
    let file = File::open(p)?;
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_reader(file);

    let mut quant_vec = Vec::<QuantRecord>::new();

    for quant_record in rdr.deserialize() {
        let quant_record: QuantRecord = quant_record?;
        quant_vec.push(quant_record);
    }

    Ok(quant_vec)
}

/// Using the quantification values present in `quants`, this function will compute a
/// posterior assignment probability for each of the alignments in `alns`, it will populate
/// the "ZW" field of the corresponding SAM/BAM record with this posterior probability estimate
/// and it will write the resulting (annotated) alignment records using `writer`.
fn process_alignment_group(alns: &mut Vec<Record>, quants: &Vec<QuantRecord>, writer: &mut Writer) {
    // tot_tpm will hold the denominator that we will
    // divide by to properly normalize the posterior probabilities
    let mut tot_tpm = 0.0;
    // an iterator over our alignment records for this read
    let mut ait = alns.iter();

    // get the TPM for all records for this read
    while let Some(a) = ait.next() {
        tot_tpm += quants[a.tid() as usize].tpm;
        if a.is_paired() {
            ait.next();
        }
    }

    assert!(tot_tpm > 0.0);
    let norm: f64 = 1.0 / tot_tpm;
    // now iterate over the records again, this time
    // computing the posterior probability and writing it
    // to the ZW tag.
    let mut a_mit = alns.iter_mut();
    while let Some(a) = a_mit.next() {
        let pp = quants[a.tid() as usize].tpm * norm;
        if let Ok(mut value) = a.aux(b"ZW") {
            value = Aux::Float(pp as f32);
        } else {
            a.push_aux(b"ZW", Aux::Float(pp as f32));
        }
        writer.write(a);
    }
}

fn assign_posterior_probabilities(args: &Args) -> Result<()> {
    // we will require to have at least 2 reading threads
    // and at least 1 writing thread.
    let read_threads: usize = cmp::max(2usize, &args.num_threads / 2);
    let write_threads: usize = if args.num_threads > read_threads {
        args.num_threads - read_threads
    } else {
        1
    };
    let _tot_threads = read_threads + write_threads;

    // open the input alignment file and set the number of reading threads
    let mut bam: bam::Reader = match &args.alignments {
        None => bam::Reader::from_stdin()
            .with_context(|| format!("failed to open SAM/BAM from STDIN"))?,
        Some(aln_file) => bam::Reader::from_path(aln_file)
            .with_context(|| format!("failed to open SAM/BAM file {}", aln_file))?,
    };
    bam.set_threads(read_threads);

    // we'll need to read the header
    let header = bam::Header::from_template(bam.header());

    // if no output option is provided, then the alignments will be written to STDOUT
    // in *BAM* format.  Otherwise, they will be written in *SAM* or *BAM* based on the
    // extension of the provided file.
    let mut writer = match &args.output {
        None => {
            // write to stdout
            Writer::from_stdout(&header, Format::Bam).unwrap()
        }
        Some(outname) => {
            if outname.ends_with(".bam") {
                Writer::from_path(&outname, &header, Format::Bam).unwrap()
            } else if outname.ends_with(".sam") {
                Writer::from_path(&outname, &header, Format::Sam).unwrap()
            } else {
                panic!("The provided output path was {}, but it must either be stdout (nothing), or end with .bam or .sam", outname)
            }
        }
    };

    // set the number of threads we will use for multithreaded writing
    writer.set_threads(write_threads);

    // read the input from the salmon quant.sf file
    if let Ok(quant_vec) = read_quants(&args.quant) {
        let hhm = header.to_hashmap();
        let nquant = quant_vec.len();
        // get the map of reference info
        let ref_map = &hhm["SQ"];
        // make sure the number of references in the sam/bam file
        // is equal to the number of quantified sequences
        // NOTE: what about decoys?
        if ref_map.len() != nquant {
            return Err(anyhow!(
                "quant file had {} records, alignment file had {} refs",
                nquant,
                ref_map.len()
            ));
        }

        // we require that every target / record in the quant.sf file has
        // a corresponding entry in the input SAM/BAM header.  Further, we
        // currently require that these entries are in the same order.
        for (i, e) in quant_vec.iter().enumerate() {
            let r = &ref_map[i];
            assert_eq!(&e.name, &r["SN"]);
        }

        // Now, we'll process the actual records.
        let mut curr_qname = Vec::<u8>::new();
        // to minimize unnecessary allocations, we'll store the
        // alignments for reach read in a vector of records, and
        // the records will be re-used between reads; growing
        // as necessary to accommodate reads with more alignments.
        let mut rvec = Vec::<Record>::new();
        let mut rcache = Vec::new();
        let mut record = Record::new();
        let mut first = true;

        while let Some(result) = bam.read(&mut record) {
            match result {
                Ok(_) => {
                    // if the current read name is different than the previous
                    // read name, then we are starting to process the alignments
                    // for a new read.
                    if &curr_qname[..] != record.qname() {
                        let _prev_name = String::from_utf8(curr_qname).unwrap();
                        curr_qname = record.qname().to_vec();
                        // if this wasn't the first read then
                        // process the alignments for the previous
                        // read.
                        if !first {
                            process_alignment_group(&mut rvec, &quant_vec, &mut writer);
                            for rec in rvec.drain(..) {
                                rcache.push(rec);
                            }
                        }
                        // this is no longer the first record (that needs to be treated specially).
                        first = false;
                    }
                    rvec.push(record);
                    // get the next available record to re-use from the cache, or allocate
                    // a new one if the cache is empty.
                    if let Some(rec) = rcache.pop() {
                        record = rec;
                    } else {
                        record = Record::new();
                    }
                }
                Err(_) => panic!("BAM parsing failed..."),
            }
        }
        // since we don't hit the conditional at the top of the above loop, we need
        // one more call to deal with the alignments for the last read.
        if rvec.len() > 0 {
            process_alignment_group(&mut rvec, &quant_vec, &mut writer);
        }
    }
    Ok(())
}

fn main() {
    let args = Args::parse();
    assign_posterior_probabilities(&args);
}
