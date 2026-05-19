#!/usr/bin/env python3
import logging
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
    format='%(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def generate_samplesheet(assemblies_dir="data/assemblies", output_csv="data/samplesheet.csv"):
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
        
    # Sort and remove potential duplicates
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
        
        if len(parts) < 2:
            logging.warning(f"Skipping file with incorrect format: {filename}")
            continue
            
        sample = parts[0]
        hap_indicator = parts[1].lower()
        
        haplotype = HAPLOTYPE_TRANSLATION.get(hap_indicator)
        
        if not haplotype:
            if hap_indicator.isdigit() and hap_indicator in ("1", "2"):
                haplotype = hap_indicator
            else:
                logging.warning(f"Warning: Haplotye for {filename} ambiguous ('{parts[1]}'). Setting default '1'.")
                haplotype = "1"
                
        # Add to samplesheet
        lines.append(f"{sample},{haplotype},{filepath}")
        count += 1
        
    # CSV output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")
    logging.info(f"Samplesheet successfully created at: {output_path}")
    logging.info(f"Total {count} phased Assemblies indexed.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Automated samplesheet generation for HPRC Pangenome.")
    parser.add_argument("--dir", default="data/assemblies", help="Folder with the assemblies.")
    parser.add_argument("--out", default="data/samplesheet.csv", help="Output path of the CSV.")
    args = parser.parse_args()
    
    generate_samplesheet(assemblies_dir=args.dir, output_csv=args.out)
