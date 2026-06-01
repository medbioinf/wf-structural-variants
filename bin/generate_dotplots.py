#!/usr/bin/env python3
"""
Based on the Swave software package (https://github.com/skandavb/Swave).
Originally licensed under the GPL-3.0.

Publication:
Wang, S., Xu, T., Zhang, P. & Ye, K. Population-level structural variant 
characterization using pangenome graphs. Nat Genet (2026). 
https://doi.org/10.1038/s41588-026-02538-6

Modified, refactored and optimized for Nextflow integration.
Copyright (c) 2026 Jonah Kapski <Jonah.Kapski@edu.ruhr-uni-bochum.de>
"""

import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import logging
import argparse
import numpy as np
from matplotlib import pyplot as plt
from seq_utils import reverse_complement_seq, calculate_stride_size


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


class Dotplot:
    def __init__(self, seq_x, seq_y, kmer_size, out_prefix, against="auto", stride_size=1, given_x_index_table=None, given_y_index_table=None, skip_forward=False, skip_reverse=False):
        self.seq_x = seq_x.upper()
        self.seq_y = seq_y.upper()
        self.seq_x_len = len(seq_x)
        self.seq_y_len = len(seq_y)

        self.kmer_size = kmer_size
        self.skip_forward = skip_forward
        self.stride_size = stride_size
        self.out_prefix = out_prefix
        self.skip_reverse = skip_reverse

        dim_y = int(self.seq_y_len / self.stride_size)
        dim_x = int(self.seq_x_len / self.stride_size)
        self.matrix = np.zeros((dim_y + 1, dim_x + 1))
        self.matrix_rev = np.zeros((dim_y + 1, dim_x + 1))

        if against == "auto":
            if self.seq_x_len >= self.seq_y_len:
                self.create_matrix_against_x(given_x_index_table)
            else:
                self.create_matrix_against_y(given_y_index_table)
        elif against == "x":
            self.create_matrix_against_x(given_x_index_table)
        elif against == "y":
            self.create_matrix_against_y(given_y_index_table)
        else:
            logging.error(f"No such against axis: {against}. Choose from [auto, x, y]")

    def rotate_to_alt2ref(self):
        self.matrix = np.fliplr(np.rot90(self.matrix, k=-1))
        self.matrix_rev = np.fliplr(np.rot90(self.matrix_rev, k=-1))

        # Switch seq x and y
        tmp_seq_x = self.seq_x
        self.seq_x = self.seq_y
        self.seq_y = tmp_seq_x

        # Switch seq x and y
        tmp_seq_x_len = self.seq_x_len
        self.seq_x_len = self.seq_y_len
        self.seq_y_len = tmp_seq_x_len

    def create_matrix_against_x(self, given_x_index_table=None):
        if given_x_index_table is None:
            self.seq_x_index_table = KmerIndex(self.seq_x, kmer_size=self.kmer_size, stride_size=self.stride_size)
        else:
            self.seq_x_index_table = given_x_index_table

        pos_on_y = 0
        while pos_on_y < self.seq_y_len - self.kmer_size:
            kmer_str = self.seq_y[pos_on_y: pos_on_y + self.kmer_size]
            index_on_y = int(pos_on_y / self.stride_size)
            
            # Check if index_on_y is within the bounds of the matrix
            if index_on_y >= self.matrix.shape[0]:
                break

            # Forward matches
            if not self.skip_forward:
                indexes_on_x = self.seq_x_index_table.find_all(kmer_str)
                if indexes_on_x is not None:
                    for index_on_x in indexes_on_x:
                        if index_on_x < self.matrix.shape[1]:
                            self.matrix[index_on_y, index_on_x] = 1

            # Reverse matches
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)
                indexes_on_x = self.seq_x_index_table.find_all(kmer_str_reversed)
                if indexes_on_x is not None:
                    for index_on_x in indexes_on_x:
                        if index_on_x < self.matrix.shape[1]:
                            self.matrix[index_on_y, index_on_x] = 1
                            self.matrix_rev[index_on_y, index_on_x] = 1

            pos_on_y += self.stride_size

    def create_matrix_against_y(self, given_y_index_table=None):
        if given_y_index_table is None:
            self.seq_y_index_table = KmerIndex(self.seq_y, kmer_size=self.kmer_size, stride_size=self.stride_size)
        else:
            self.seq_y_index_table = given_y_index_table

        pos_on_x = 0
        while pos_on_x < self.seq_x_len - self.kmer_size:
            kmer_str = self.seq_x[pos_on_x: pos_on_x + self.kmer_size]
            index_on_x = int(pos_on_x / self.stride_size)
            
            # Check if index_on_x is within the bounds of the matrix
            if index_on_x >= self.matrix.shape[1]:
                break

            # Forward matches
            if not self.skip_forward:
                indexes_on_y = self.seq_y_index_table.find_all(kmer_str)
                if indexes_on_y is not None:
                    for index_on_y in indexes_on_y:
                        if index_on_y < self.matrix.shape[0]:
                            self.matrix[index_on_y, index_on_x] = 1

            # Reverse matches
            if not self.skip_reverse:
                kmer_str_reversed = reverse_complement_seq(kmer_str)
                indexes_on_y = self.seq_y_index_table.find_all(kmer_str_reversed)
                if indexes_on_y is not None:
                    for index_on_y in indexes_on_y:
                        if index_on_y < self.matrix.shape[0]:
                            self.matrix[index_on_y, index_on_x] = 1
                            self.matrix_rev[index_on_y, index_on_x] = 1

            pos_on_x += self.stride_size

    def get_seq_x_index_table(self):
        return self.seq_x_index_table

    def to_png(self, reverse=False, out_img=False):
        suffix = "rev" if reverse else "fwd"
        self.dotplot_file = f"{self.out_prefix}.dotplot_{suffix}.png"

        target_matrix = self.matrix_rev if reverse else self.matrix

        if np.max(target_matrix) == 0:
            matrix_resize_norm = 255 * np.ones(np.shape(target_matrix))
        else:
            matrix_resize_norm = 255 * abs(target_matrix - np.max(target_matrix)) / (np.max(target_matrix) - np.min(target_matrix))

        if out_img:
            plt.imsave(self.dotplot_file, matrix_resize_norm, cmap='gray')

        return matrix_resize_norm


class KmerIndex:
    """
    A class to index a sequence by its k-mers for efficient O(1) lookup.
    """
    def __init__(self, seq, kmer_size, stride_size):
        self.seq = seq.upper()
        self.index_table = {}

        pos_on_seq = 0
        seq_len = len(seq)
        
        while pos_on_seq < seq_len:
            index_on_seq = int(pos_on_seq / stride_size)

            for i in range(stride_size):
                start_pos = pos_on_seq + i
                end_pos = start_pos + kmer_size
                
                # Check if end_pos is within the bounds of the sequence
                if end_pos > seq_len:
                    break
                
                kmer_str = self.seq[start_pos:end_pos]

                if kmer_str not in self.index_table:
                    self.index_table[kmer_str] = []

                self.index_table[kmer_str].append(index_on_seq)

            pos_on_seq += stride_size

    def find_all(self, seq_str):
        return self.index_table.get(seq_str)


def generate_dotplots(ref_seq, alt_seq, dotplot_output_prefix, options):
    """
    for the ref seq and alt seq, generate dotplots and projections, which are further segmented into matrixes
    """
    dotplot_stride_size = calculate_stride_size(len(ref_seq), len(alt_seq), max_matrix_dim=options.max_dotplot_dim)
    
    ref2ref_dp = Dotplot(ref_seq, ref_seq, options.kmer_size, f"{dotplot_output_prefix}_ref2ref", against="x", stride_size=dotplot_stride_size)
    ref2alt_dp = Dotplot(ref_seq, alt_seq, options.kmer_size, f"{dotplot_output_prefix}_ref2alt", against="x",
                         stride_size=dotplot_stride_size, given_x_index_table=ref2ref_dp.get_seq_x_index_table())
    
    alt2alt_dp = Dotplot(alt_seq, alt_seq, options.kmer_size, f"{dotplot_output_prefix}_alt2alt", against="x", stride_size=dotplot_stride_size)
    alt2ref_dp = Dotplot(alt_seq, ref_seq, options.kmer_size, f"{dotplot_output_prefix}_alt2ref", against="x",
                         stride_size=dotplot_stride_size, given_x_index_table=alt2alt_dp.get_seq_x_index_table())
    
    np.savez_compressed(
        f"{dotplot_output_prefix}_matrices.npz",
        stride=dotplot_stride_size,
        r2r_fwd=ref2ref_dp.matrix, r2r_rev=ref2ref_dp.matrix_rev,
        r2a_fwd=ref2alt_dp.matrix, r2a_rev=ref2alt_dp.matrix_rev,
        a2a_fwd=alt2alt_dp.matrix, a2a_rev=alt2alt_dp.matrix_rev,
        a2r_fwd=alt2ref_dp.matrix, a2r_rev=alt2ref_dp.matrix_rev
    )
    
    logging.info(f"Saved dotplot matrices to {dotplot_output_prefix}_matrices.npz")
    
    if options.save_dotplot_images:
        ref2ref_dp.to_png(reverse=False, out_img=True)
        ref2ref_dp.to_png(reverse=True, out_img=True)
        
        ref2alt_dp.to_png(reverse=False, out_img=True)
        ref2alt_dp.to_png(reverse=True, out_img=True)
        
        alt2alt_dp.to_png(reverse=False, out_img=True)
        alt2alt_dp.to_png(reverse=True, out_img=True)
        
        alt2ref_dp.to_png(reverse=False, out_img=True)
        alt2ref_dp.to_png(reverse=True, out_img=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate dotplots and their matrices for a reference and alternative sequence.")

    parser.add_argument('--ref_seq', required=True, help='Reference sequence in FASTA format')
    parser.add_argument('--alt_seq', required=True, help='Alternative sequence in FASTA format')
    parser.add_argument('--out_prefix', required=True, help='Output prefix for generated files')
    
    parser.add_argument('--kmer_size', type=int, default=30, help='K-mer size for dotplot generation (default: 30)')
    parser.add_argument('--max_dotplot_dim', type=int, default=1500, help='Maximum dimension for dotplot matrices (default: 1500)')
    parser.add_argument('--save_dotplot_images', action='store_true', help='Whether to save the generated dotplot images as PNG files (default: False)')
    
    options = parser.parse_args()
    
    generate_dotplots(
        ref_seq=options.ref_seq,
        alt_seq=options.alt_seq,
        dotplot_output_prefix=options.out_prefix,
        options=options
    )
    
    sys.exit(0)