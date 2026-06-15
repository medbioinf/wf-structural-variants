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

def reverse_complement_seq(seq):
    """
    Returns the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def calculate_stride_size(seq_x, seq_y):
    """
    Calculates the stride size for dotplot alignment to keep the final matrix dimension manageable.
    """
    max_seq_len = max([len(seq_x), len(seq_y)])

    if max_seq_len <= 1500:
        stride_size = 1
    elif max_seq_len <= 15000:
        stride_size = 10
    elif max_seq_len <= 75000:
        stride_size = 30
    elif max_seq_len <= 100000:
        stride_size = 50
    else:
        stride_size = 100

    return stride_size


def is_kmer_similar(kmer1, kmer2, mismatch_thresh=0):
    """
    Compares two k-mers and returns True if they are similar within a specified mismatch threshold.
    """
    if len(kmer1) != len(kmer2):
        return False
    mismatch_cnt = sum(b1 != b2 for b1, b2 in zip(kmer1, kmer2))
    return mismatch_cnt <= mismatch_thresh