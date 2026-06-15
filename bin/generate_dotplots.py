#!/usr/bin/env python3
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
import argparse
import numpy as np
from matplotlib import pyplot as plt

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from seq_utils import reverse_complement_seq, calculate_stride_size


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


class Dotplot:
    def __init__(self, seq_x, seq_y, kmer_size, out_prefix, against="auto", stride_size=None, given_x_kmer_index=None, given_y_kmer_index=None, skip_forward=False, skip_reverse=False):

        self.seq_x = seq_x.upper()
        self.seq_y = seq_y.upper()

        self.seq_x_len = len(seq_x)
        self.seq_y_len = len(seq_y)

        self.kmer_size = kmer_size
        
        self.skip_forward = skip_forward
        self.skip_reverse = skip_reverse

        if stride_size is None:
            self.stride_size = calculate_stride_size(self.seq_x, self.seq_y)
        else:
            self.stride_size = stride_size

        self.out_prefix = out_prefix

        self.matrix = np.zeros((int(self.seq_y_len / self.stride_size) + 1, int(self.seq_x_len / self.stride_size) + 1))
        self.matrix_rev = np.zeros((int(self.seq_y_len / self.stride_size) + 1, int(self.seq_x_len / self.stride_size) + 1))

        if against == "auto":
            if self.seq_x_len >= self.seq_y_len:
                self.create_matrix_against_x(given_x_kmer_index)
            else:
                self.create_matrix_against_y(given_y_kmer_index)

        elif against == "x":
            self.create_matrix_against_x(given_x_kmer_index)

        elif against == "y":
            self.create_matrix_against_y(given_y_kmer_index)

        else:
            logging.error("No such against axis: {}. Choose from [auto, x, y]".format(against))

        self.matrix = self.matrix[:-1, :-1]
        self.matrix_rev = self.matrix_rev[:-1, :-1]

    def create_matrix_against_x(self, given_x_kmer_index=None):
        if given_x_kmer_index is None:
            self.seq_x_kmer_index = KmerIndex(self.seq_x, kmer_size=self.kmer_size, stride_size=self.stride_size)
        else:
            self.seq_x_kmer_index = given_x_kmer_index

        pos_on_y = 0
        while pos_on_y < self.seq_y_len:
            kmer_str = self.seq_y[pos_on_y: pos_on_y + self.kmer_size]

            index_on_y = int(pos_on_y / self.stride_size)

            # for original kmer
            if not self.skip_forward:
                indexes_on_x = self.seq_x_kmer_index.find_all(kmer_str)
                if indexes_on_x is not None:
                    for index_on_x in indexes_on_x:
                        self.matrix[index_on_y, index_on_x] = 1

            # for reversed kmer
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)

                indexes_on_x = self.seq_x_kmer_index.find_all(kmer_str_reversed)

                if indexes_on_x is not None:
                    for index_on_x in indexes_on_x:
                        self.matrix[index_on_y, index_on_x] = 1
                        self.matrix_rev[index_on_y, index_on_x] = 1

            pos_on_y += self.stride_size

    def create_matrix_against_y(self, given_y_kmer_index=None):
        if given_y_kmer_index is None:
            self.seq_y_kmer_index = KmerIndex(self.seq_y, kmer_size=self.kmer_size, stride_size=self.stride_size)
        else:
            self.seq_y_kmer_index = given_y_kmer_index

        pos_on_x = 0
        while pos_on_x < self.seq_x_len:

            kmer_str = self.seq_x[pos_on_x: pos_on_x + self.kmer_size]

            index_on_x = int(pos_on_x / self.stride_size)

            # for original kmer
            if not self.skip_forward:
                indexes_on_y = self.seq_y_kmer_index.find_all(kmer_str)

                if indexes_on_y is not None:
                    for index_on_y in indexes_on_y:
                        self.matrix[index_on_y, index_on_x] = 1

            # for reversed kmer
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)

                indexes_on_y = self.seq_y_kmer_index.find_all(kmer_str_reversed)

                if indexes_on_y is not None:
                    for index_on_y in indexes_on_y:
                        self.matrix[index_on_y, index_on_x] = 1
                        self.matrix_rev[index_on_y, index_on_x] = 1

            pos_on_x += self.stride_size

    def get_seq_x_kmer_index(self):
        return self.seq_x_kmer_index

    def get_seq_y_kmer_index(self):
        return self.seq_y_kmer_index

    def to_png(self, reverse=False, out_img=False):
        self.dotplot_file = self.out_prefix + ".dotplot.png"
        target_matrix = self.matrix_rev if reverse else self.matrix
        
        if target_matrix.size == 0:
            logging.warning(f"Skipping PNG generation for {self.dotplot_file}: Matrix is empty.")
            return
        
        if np.max(target_matrix) == 0:
            matrix_resize_norm = 255 * np.ones(np.shape(target_matrix))
        else:
            matrix_resize_norm = 255 * abs(target_matrix - np.max(target_matrix)) / (np.max(target_matrix) - np.min(target_matrix))

        if out_img:
            plt.imsave(self.dotplot_file, matrix_resize_norm, cmap='gray')
            plt.close()


class KmerIndex:
    def __init__(self, seq, kmer_size, stride_size):
        self.seq = seq
        self.index_table = {}

        pos_on_seq = 0
        while pos_on_seq < len(seq):
            index_on_seq = int(pos_on_seq / stride_size)

            for i in range(stride_size):
                kmer_str = seq[pos_on_seq + i: pos_on_seq + i + kmer_size]

                if kmer_str not in self.index_table:
                    self.index_table[kmer_str] = []

                self.index_table[kmer_str].append(index_on_seq)

            pos_on_seq += stride_size

    def find_all(self, seq_str):
        """
        Returns all matrix indices where this k-mer exists.
        """
        return self.index_table.get(seq_str)
    

def process_sample_alleles_to_matrices(alt_fasta_path, ref_fasta_path, pangenome_gfa_fasta_path, options=None):
    """
    Reads the extracted sample allele fasta (ALT) and the original reference genome
    fasta (REF) and extracts the alternative and reference sequences with padding.
    Stores the dotplot matrices for each snarl into a dictionary and saves it as an npz file.

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
        output_npz_path = f"{options.npz_out_prefix}_dotplot_matrices.npz"
        np.savez(output_npz_path, **snarl_dotplot_dict)
        logging.info(f"Successfully saved {len(snarl_dotplot_dict)} matrices.")


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
    
    matrices_dict = generate_dotplot_matrices(
        ref_seq=final_ref_seq,
        alt_seq=final_alt_seq,
        dotplot_stride_size=dotplot_stride_size,
        dotplot_output_prefix=dotplot_output_prefix,
        options=options
    )
    
    for matrix_type, matrix_array in matrices_dict.items():
        matrix_key = f"{dotplot_id}+++{matrix_type}+++{dotplot_stride_size}"
        snarl_dotplot_dict[matrix_key] = matrix_array


def generate_dotplot_matrices(ref_seq, alt_seq, dotplot_stride_size, dotplot_output_prefix, options):
    """
    Generates the dotplot matrices for the given reference and alternative sequences and saves PNG visualizations if specified.
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
        "x2x_ref2ref_matrix": x2x_ref2ref.matrix,
        "x2y_ref2alt_matrix": x2y_ref2alt.matrix,
        "x2y_ref2alt_matrix_rev": x2y_ref2alt.matrix_rev,
        "x2x_alt2alt_matrix": x2x_alt2alt.matrix,
        "stride_size": dotplot_stride_size
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates 2D-dotplot matrices (REF2REF, REF2ALT, ALT2ALT)"
                                     "from pangenome graph alleles and exports them into a compressed .npz archive.")
    
    parser.add_argument('--alt_fasta', required=True, help="Path to the extracted alleles fasta (ALT).")
    parser.add_argument('--ref_fasta', required=True, help="Path to the fasta file of the original reference genome (REF).")
    parser.add_argument('--gfa_fasta', required=True, help="Path to the pangenome gfa fasta file.")
    parser.add_argument('--npz_out_prefix', required=True, help="Output prefix for the generated .npz file and dotplot images.")
    parser.add_argument('--img_out_prefix', required=True, help="Output prefix for the generated dotplot images.")
    
    # Optional parameters
    parser.add_argument('--kmer_size', type=int, default=30, help="K-mer size (default: 30).")
    parser.add_argument('--max_sv_size', type=int, default=100000, help="Maximum SV size to process (default: 100000).")
    parser.add_argument('--spec_path', help="Process only a specific path or snarl ID.")
    parser.add_argument('--save_dotplot_images', action='store_true', help="Saves PNG dotplot images for each snarl.")
    parser.add_argument('--skip_forward', action='store_true', help="Skip forward k-mer matches in dotplot matrix generation.")
    parser.add_argument('--skip_reverse', action='store_true', help="Skip reverse k-mer matches in dotplot matrix generation.")
    
    options = parser.parse_args()
    
    logging.info(f"Starting dotplot matrix generation for {options.alt_fasta}...")
    
    process_sample_alleles_to_matrices(
        alt_fasta_path=options.alt_fasta,
        ref_fasta_path=options.ref_fasta,
        pangenome_gfa_fasta_path=options.gfa_fasta,
        options=options
    )
    
    logging.info(f"Dotplot matrix generation for {options.alt_fasta} completed.")
    
    sys.exit(0)