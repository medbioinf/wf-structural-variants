#!/usr/bin/env python3
import sys
import logging
import argparse
from pathlib import Path


# Supported genome assembly file extensions
SUPPORTED_EXTENSIONS = ("*.fa", "*.fasta", "*.fna")

# List of extensions to strip from filenames to get the base name
EXTENSIONS_TO_STRIP = (".fa", ".fasta", ".fna", ".gz")

# Mapping of haplotype indicators to 1 or 2
HAPLOTYPE_TRANSLATION = {
    "paternal": "1", "maternal": "2",
    "hap1": "1", "hap2": "2",
    "h1": "1", "h2": "2",
    "1": "1", "2": "2"
}

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def generate_samplesheet(assemblies_dir="data/assemblies", output_csv="data/samplesheet.csv", allow_unphased=False):
    """
    Scans the directory containing the assemblies and automatically generates
    an nf-core-compatible samplesheet based on HPRC naming conventions.
    """
    assemblies_path = Path(assemblies_dir)
    output_path = Path(output_csv)
    
    if not assemblies_path.exists():
        raise FileNotFoundError(f"Folder '{assemblies_dir}' does not exist.")
        
    # Samplesheet header
    lines = ["sample,haplotype,fasta"]
    
    # Supported file extensions
    filepaths = []
    for ext in SUPPORTED_EXTENSIONS:
        filepaths.extend(assemblies_path.glob(ext))
    
    # Sort lines and remove potential duplicates
    filepaths = sorted(list(set(filepaths)))
    
    if not filepaths:
        raise FileNotFoundError(f"No FASTA files found in '{assemblies_dir}' {SUPPORTED_EXTENSIONS}.")
        
    count = 0
    for filepath in filepaths:
        filename = filepath.name
        
        # Remove extensions to get the base name
        current_path = filepath
        while current_path.suffix.lower() in EXTENSIONS_TO_STRIP:
            current_path = current_path.with_suffix("")     # remove extension
            
        name_without_ext = current_path.name
        parts = name_without_ext.split(".")
            
        sample = parts[0]
        
        if len(parts) > 1:
            hap_indicator = parts[1].lower()
            haplotype = HAPLOTYPE_TRANSLATION.get(hap_indicator)
        else:
            hap_indicator = "none"
            haplotype = None
        
        if not haplotype:
            if allow_unphased:
                haplotype = "0"
                logging.info(f"No haplotype recognized for {filename}. Treating as unphased ('0').")
            else:
                logging.warning(f"Warning: No haplotype recognized for {filename}.")
                logging.error(f" -> Expected indicators like: paternal, maternal, hap1, h2, 1, etc.")
                logging.error(f" -> To allow unphased assemblies (fallback to '0'), run with --allow-unphased")
                logging.error(f" -> Skipping file to prevent downstream errors.")
                continue
                
        # Add to samplesheet
        lines.append(f"{sample},{haplotype},{filepath}")
        count += 1
        
    # CSV output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")
    logging.info(f"Samplesheet successfully created at: {output_path}")
    logging.info(f"Total {count} assemblies indexed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automated samplesheet generation for Assembly files.")
    parser.add_argument("--dir", default="data/assemblies", help="Folder with the assemblies.")
    parser.add_argument("--out", default="data/samplesheet.csv", help="Output path of the CSV.")
    parser.add_argument("--allow-unphased", action="store_true", 
                        help="Treat files with unrecognized haplotypes as unphased (assigned to '0') instead of skipping them.")
    args = parser.parse_args()
    
    generate_samplesheet(assemblies_dir=args.dir, output_csv=args.out, allow_unphased=args.allow_unphased)
