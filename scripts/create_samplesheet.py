#!/usr/bin/env python3
import sys
import logging
import argparse
from pathlib import Path


# Supported genome assembly file extensions
SUPPORTED_EXTENSIONS = ("*.fa", "*.fasta", "*.fna")
EXTENSIONS_TO_STRIP = (".fa", ".fasta", ".fna", ".gz")

HAPLOTYPE_TRANSLATION = {
    "paternal": "1", "maternal": "2",
    "hap1": "1", "hap2": "2",
    "h1": "1", "h2": "2",
    "1": "1", "2": "2"
}

logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def find_bam_dir(sample_id, bams_path):
    sample_dir = bams_path / sample_id
    if not sample_dir.exists() or not sample_dir.is_dir():
        return ""
    
    bam_files = list(sample_dir.rglob("*.bam"))
    
    if not bam_files:
        return ""
    
    # prefer directories containing "pass" in their name
    pass_bams = [f for f in bam_files if "pass" in str(f.parent).lower()]
    if pass_bams:
        return str(pass_bams[0].parent)
        
    # fallback: return the parent directory of the first BAM file found
    return str(bam_files[0].parent)


def generate_samplesheet(assemblies_dir="data/assemblies", bams_dir="data/bams", output_csv="data/samplesheet.csv",
                         allow_unphased=False, exclude_assemblies=False, exclude_bams=False):
    """
    Scans the directory containing the assemblies and bams and automatically generates a pipeline compatible samplesheet.
    """
    if exclude_assemblies and exclude_bams:
        logging.warning("Both assemblies and bams are excluded. Nothing to process.")
        sys.exit(1)
    
    assemblies_path = Path(assemblies_dir)
    bams_path = Path(bams_dir)
    output_path = Path(output_csv)
    
    if not exclude_assemblies and not assemblies_path.exists():
        logging.warning(f"Assemblies folder '{assemblies_dir}' does not exist.")
    if not exclude_bams and not bams_path.exists():
        logging.warning(f"BAMs folder '{bams_dir}' does not exist.")
    
    # csv header        
    lines = ["sample,haplotype,fasta,bam_dir"]
    
    # map sample_id to bam_dir
    bam_map = {}
    if not exclude_bams and bams_path.exists():
        for p_dir in bams_path.iterdir():
            if p_dir.is_dir():
                sample_id = p_dir.name
                bam_dir = find_bam_dir(sample_id, bams_path)
                if bam_dir:
                    bam_map[sample_id] = bam_dir
    
    
    # collect fasta files and their haplotypes
    fasta_records = []
    if not exclude_assemblies and assemblies_path.exists():
        filepaths = []
        for ext in SUPPORTED_EXTENSIONS:
            filepaths.extend(assemblies_path.glob(ext))
        filepaths = sorted(list(set(filepaths)))
        
        for filepath in filepaths:
            current_path = filepath
            while current_path.suffix.lower() in EXTENSIONS_TO_STRIP:
                current_path = current_path.with_suffix("")
                
            parts = current_path.name.split(".")
            sample = parts[0]
            
            haplotype = None
            if len(parts) > 1:
                haplotype = HAPLOTYPE_TRANSLATION.get(parts[1].lower())
            
            if not haplotype:
                if allow_unphased:
                    haplotype = "0"
                else:
                    logging.warning(f"Skipping {filepath.name}: Unrecognized haplotype.")
                    continue
            
            fasta_records.append({
                "sample": sample,
                "haplotype": haplotype,
                "fasta": str(filepath)
            })

    processed_samples = set()
    count = 0

    # process fasta samples
    for record in fasta_records:
        sample = record["sample"]
        bam_dir = bam_map.get(sample, "")
        lines.append(f"{sample},{record['haplotype']},{record['fasta']},{bam_dir}")
        processed_samples.add(sample)
        count += 1

    # process bam samples
    for sample, bam_dir in bam_map.items():
        if sample not in processed_samples:
            lines.append(f"{sample},,,{bam_dir}")
            count += 1
        
    # write csv output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")
    logging.info(f"Samplesheet successfully created at: {output_path}")
    logging.info(f"Total {count} entries indexed.")


def main():
    parser = argparse.ArgumentParser(description="Automated samplesheet generation for FASTA and BAM files.")
    parser.add_argument("--assemblies-dir", default="data/assemblies", help="Folder with the assemblies.")
    parser.add_argument("--bams-dir", default="data/bams", help="Folder with the BAM files.")
    parser.add_argument("--out", default="data/samplesheet.csv", help="Output path of the CSV.")
    parser.add_argument("--allow-unphased", action="store_true", help="Treat files with unrecognized haplotypes as unphased (assigned to '0') instead of skipping them.")
    parser.add_argument("--exclude-bams", action="store_true", help="Ignore BAM directories when generating the sheet.")
    parser.add_argument("--exclude-assemblies", action="store_true", help="Ignore assembly FASTA files when generating the sheet.")
    args = parser.parse_args()
    
    generate_samplesheet(assemblies_dir=args.assemblies_dir, bams_dir=args.bams_dir, output_csv=args.out,
                         allow_unphased=args.allow_unphased, exclude_bams=args.exclude_bams, exclude_assemblies=args.exclude_assemblies)
    
    sys.exit(0)

if __name__ == "__main__":
    main()