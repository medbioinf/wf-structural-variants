"""
Based on the Swave software package (https://github.com/skandavb/Swave).
Originally licensed under the GPL-3.0.

Publication:
Wang, S., Xu, T., Zhang, P. & Ye, K. Population-level structural variant 
characterization using pangenome graphs. Nat Genet (2026). 
https://doi.org/10.1038/s41588-026-02538-6

Modified and refactored for Nextflow integration.
Copyright (c) 2026 Jonah Kapski <Jonah.Kapski@edu.ruhr-uni-bochum.de>
"""

import sys
import os
import re
import pysam
import logging
import numpy as np
import pickle
import gzip

from src.generate_dotplots_projections_mod.structures import Dotplot
from src.utils.seq_utils import calculate_stride_size


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def process_sample_alleles_to_matrices(alt_fasta_path, ref_fasta_path, pangenome_gfa_fasta_path, options=None):
    """
    Reads the extracted sample allele fasta (ALT) and the original reference genome
    fasta (REF) and extracts the alternative and reference sequences with padding.
    Stores the dotplots for each snarl into a dictionary and saves it as a pickle file.

    Args:
        alt_fasta_path (str): Path to the extracted sample allele fasta (ALT)
        ref_fasta_path (str): Path to the original reference genome fasta (REF)
        pangenome_gfa_fasta_path (str): Path to the GFA fasta file containing node sequences
        options: Additional options for processing
    """
    ref_file = pysam.FastaFile(ref_fasta_path)
    gfa_fasta_file = pysam.FastaFile(pangenome_gfa_fasta_path)
    
    with open(alt_fasta_path, 'r') as alt_file:
        total_snarl_count = sum(1 for line in alt_file if line.startswith('>'))
    
    snarl_count = 0
    current_header = None
    snarl_dotplot_dict = {}
    
    with open(alt_fasta_path, 'r') as alt_file:
        for line in alt_file:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                current_header = line[1:]   # e.g., "HG002_hap1|>s3>s5|chr1:3690-3752"
            else:
                snarl_count += 1
                raw_alt_seq = "" if line == "-" else line.upper()
                
                # header meta information
                header_parts = current_header.split('|')
                sample_id = header_parts[0]
                snarl_id = header_parts[1]
                coords = header_parts[2]  # "chr1:3690-3752"
                is_reversed_mapping = header_parts[3].split(':')[1] == "true"
                
                chrom, pos_range = coords.split(':')
                snarl_ref_start, snarl_ref_end = map(int, pos_range.split('-'))
                
                nodes_with_orients = re.findall(r'([><][a-zA-Z0-9]+)', snarl_id)
                snarl_start_node_with_orient = nodes_with_orients[0]    # e.g., ">s3"
                snarl_end_node_with_orient = nodes_with_orients[-1]     # e.g., ">s5"
                
                chrom_len = ref_file.get_reference_length(chrom)
                
                if is_reversed_mapping:
                    snarl_start_node_id = snarl_start_node_with_orient[1:]   # remove orient character, e.g., "s3"
                    snarl_end_node_id = snarl_end_node_with_orient[1:]
                    
                    # retrieve node lengths from GFA fasta
                    start_node_len = gfa_fasta_file.get_reference_length(snarl_start_node_id)
                    end_node_len = gfa_fasta_file.get_reference_length(snarl_end_node_id)
                    
                    # adjust reference coordinates by node lengths as in Swave code and ensure coordinates are within chromosome bounds
                    snarl_ref_start = max(0, snarl_ref_start - start_node_len)
                    snarl_ref_end = min(chrom_len, snarl_ref_end + end_node_len)
                
                raw_ref_seq = ref_file.fetch(chrom, snarl_ref_start, snarl_ref_end).upper()
                
                extend_len = min(10000, 2 * max([len(raw_ref_seq), len(raw_alt_seq), (snarl_ref_end - snarl_ref_start)]))
                
                left_padding_seq = ref_file.fetch(chrom, max(snarl_ref_start - extend_len, 0), snarl_ref_start).replace("N", "").upper()
                right_padding_seq = ref_file.fetch(chrom, snarl_ref_end, min(chrom_len, snarl_ref_end + extend_len)).replace("N", "").upper()

                final_ref_seq = left_padding_seq + raw_ref_seq + right_padding_seq
                final_alt_seq = left_padding_seq + raw_alt_seq + right_padding_seq
                
                process_and_plot_snarl(
                    snarl_id=snarl_id,
                    chrom=chrom,
                    snarl_ref_start=snarl_ref_start,
                    snarl_ref_end=snarl_ref_end,
                    final_ref_seq=final_ref_seq,
                    final_alt_seq=final_alt_seq,
                    left_padding_len=len(left_padding_seq),
                    right_padding_len=len(right_padding_seq),
                    snarl_dotplot_dict=snarl_dotplot_dict,
                    is_reversed_mapping=is_reversed_mapping,
                    snarl_count=snarl_count,
                    total_snarl_count=total_snarl_count,
                    options=options
                )
    
    ref_file.close()
    gfa_fasta_file.close()
    
    if snarl_dotplot_dict:
        ouput_pickle_path = f"{options.pkl_out_prefix}_dotplots.pkl.gz"
        
        with open(ouput_pickle_path, 'wb') as f:
            with gzip.GzipFile(fileobj=f) as gz:
                pickle.dump(snarl_dotplot_dict, gz)
                
        logging.info(f"Successfully saved {len(snarl_dotplot_dict)} snarl dotplot bundles.")


def process_and_plot_snarl(snarl_id, chrom, snarl_ref_start, snarl_ref_end, final_ref_seq,
                           final_alt_seq, left_padding_len, right_padding_len, snarl_dotplot_dict,
                           is_reversed_mapping, snarl_count, total_snarl_count, options):
    
    if snarl_count % 1000 == 0:
        logging.info(f"Generating {snarl_count}th snarl. In total {total_snarl_count} snarls")
    
    if is_reversed_mapping:
        if (snarl_ref_end - snarl_ref_start) > options.max_sv_size:
            return
    
    if options.spec_path is not None and options.spec_path not in snarl_id:
        return
            
    # Note: since Swave in its latest version is only using minigraph, resolving snarls into sub-snarls
    # is not necessary as minigraph snarls do not contain shared nodes in the middle of the snarl path.
    
    dotplot_ref_start = snarl_ref_start - left_padding_len
    dotplot_ref_end = snarl_ref_end + right_padding_len
    
    dotplot_id = "{}+++{}+++{}+++{}+++{}+++{}".format(
        snarl_id, snarl_ref_start, snarl_ref_end, chrom,
        dotplot_ref_start, dotplot_ref_end
    )
    
    dotplot_stride_size = calculate_stride_size(final_ref_seq, final_alt_seq)    
    dotplot_output_prefix = os.path.join(options.img_out_prefix, dotplot_id)
    
    dotplot_objects_bundle = generate_dotplots(
        ref_seq=final_ref_seq,
        alt_seq=final_alt_seq,
        dotplot_stride_size=dotplot_stride_size,
        dotplot_output_prefix=dotplot_output_prefix,
        options=options
    )
    
    snarl_dotplot_dict[dotplot_id] = dotplot_objects_bundle


def generate_dotplots(ref_seq, alt_seq, dotplot_stride_size, dotplot_output_prefix, options):
    """
    Generates the dotplot objects for the given reference and alternative sequence and saves matrix PNG visualizations if specified.
    """
    prefix_clean = dotplot_output_prefix.replace(">", "+").replace("<", "-").replace("*", "none")
    
    if options.save_dotplot_images:
        os.makedirs(os.path.dirname(prefix_clean), exist_ok=True)
    
    x2x_ref2ref = Dotplot(ref_seq, ref_seq, options.kmer_size, out_prefix=f"{prefix_clean}.ref2ref", stride_size=dotplot_stride_size)
    
    x2y_ref2alt = Dotplot(ref_seq, alt_seq, options.kmer_size, out_prefix=f"{prefix_clean}.ref2alt", stride_size=dotplot_stride_size, 
                          given_x_kmer_index=x2x_ref2ref.get_seq_x_kmer_index())
    
    x2x_alt2alt = Dotplot(alt_seq, alt_seq, options.kmer_size, out_prefix=f"{prefix_clean}.alt2alt", stride_size=dotplot_stride_size)
    
    if options.save_dotplot_images:
        x2x_ref2ref.to_png(out_img=True)
        x2y_ref2alt.to_png(out_img=True)
        
        if x2y_ref2alt.matrix_rev.size > 0 and np.max(x2y_ref2alt.matrix_rev) > 0:
            orig_prefix = x2y_ref2alt.out_prefix
            x2y_ref2alt.out_prefix = f"{prefix_clean}.ref2alt_reverse"
            x2y_ref2alt.to_png(reverse=True, out_img=True)
            x2y_ref2alt.out_prefix = orig_prefix
        
        x2x_alt2alt.to_png(out_img=True)
    
    return {
        "x2x_ref2ref": x2x_ref2ref,
        "x2y_ref2alt": x2y_ref2alt,
        "x2x_alt2alt": x2x_alt2alt,
        "stride_size": dotplot_stride_size
    }
