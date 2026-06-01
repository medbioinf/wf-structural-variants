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


def is_kmer_similar(kmer1, kmer2, mismatch_thresh=0):
    """
    Compares two k-mers and returns True if they are similar within a specified mismatch threshold.
    """
    if len(kmer1) != len(kmer2):
        return False
    mismatch_cnt = sum(b1 != b2 for b1, b2 in zip(kmer1, kmer2))
    return mismatch_cnt <= mismatch_thresh


def calculate_stride_size(seq_len_x, seq_len_y, max_matrix_dim=1500):
    """
    Calculates the optimal stride size to ensure that the resulting 2D matrix
    does not exceed a defined maximum pixel dimension (default: 1500x1500).
    
    (Replaces Swave's hard thresholds with a smooth calculation to keep matrix dimensions within the desired target range.)
    """
    max_seq_len = max(seq_len_x, seq_len_y)
    
    # If the max seq length fits within the max matrix dimension, we can use a stride of 1 (no downsampling)
    if max_seq_len <= max_matrix_dim:
        return 1
    
    # Otherwise: divide the sequence length by the desired matrix dimension and round up (e.g., 150.000 bp / 1500 max dim = stride of 100)
    stride_size = int(max_seq_len / max_matrix_dim)
    
    return max(1, stride_size)


def reverse_complement_seq(seq):
    """
    Returns the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))