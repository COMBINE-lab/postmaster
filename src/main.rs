use anyhow::anyhow;
use anyhow::{Context, Result};
use clap::Parser;
use csv::ReaderBuilder;
use rust_htslib::{bam, bam::Format, bam::Read, bam::Record, bam::record::Aux, bam::Writer};
use seine::salmon::QuantRecord;
use std::error::Error;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::vec::Vec;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    // Path to the SAM/BAM file whose alignments will be annotated
    #[clap(short, long)]
    alignments: String,

    // Path to the salmon quantification file
    #[clap(short, long)]
    quant: String,

    // Output path of modified alignments
    #[clap(short, long)]
    output: String,
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

fn process_alignment_group(alns: &mut Vec<Record>, quants: &Vec<QuantRecord>, writer: &mut Writer) {

    let mut tot_tpm = 0.0;
    let mut ait = alns.iter();

    while let Some(a) = ait.next() {
        tot_tpm += quants[a.tid() as usize].tpm;
        if (a.is_paired()) { ait.next(); }
    }

    let mut a_mit = alns.iter_mut();
    while let Some(a) = a_mit.next() {
        let pp = quants[a.tid() as usize].tpm / tot_tpm;
        
        if let Ok(mut value) = a.aux(b"ZW") {
            value = Aux::Float(pp as f32);
        } else {
            a.push_aux(b"ZW", Aux::Float(pp as f32));
        }
        writer.write(a);
    }
}

fn assign_posterior_probabilities(args: &Args) -> Result<()> {
    let mut bam = bam::Reader::from_path(&args.alignments)
        .with_context(|| format!("failed to open SAM/BAM file {}", &args.alignments))?;

    let header = bam::Header::from_template(bam.header());
    let mut writer = Writer::from_path(&args.output, &header, Format::Sam).unwrap();

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

        for (i, e) in quant_vec.iter().enumerate() {
            let r = &ref_map[i];
            assert_eq!(&e.name, &r["SN"]);
        }

        let mut curr_qname = Vec::<u8>::new();
        let mut rvec = Vec::<Record>::new();
        let mut rcache = Vec::new();
        let mut record = Record::new();
        let mut avail_record_ptr = 0;
        let mut first = true;

        while let Some(result) = bam.read(&mut record) {
            match result {
                Ok(_) => {
                    // starting a new read
                    if (&curr_qname[..] != record.qname()) {
                        let prev_name = String::from_utf8(curr_qname).unwrap();
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

                        first = false;
                    }
                    rvec.push(record);
                    if let Some(rec) = rcache.pop() {
                        record = rec;
                    } else {
                        record = Record::new();
                    }
                }
                Err(_) => panic!("BAM parsing failed..."),
            }
        }
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
