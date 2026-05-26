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
import logging


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


class Node:
    def __init__(self, node_id, length):
        self.id = node_id
        self.length = length


class Snarl:
    def __init__(self, start_node_id, start_node_orient, end_node_id, end_node_orient, ref_chrom, ref_start, ref_end, reversed_mapping=False):
        self.snarl_id = f"{start_node_orient}{start_node_id}{end_node_orient}{end_node_id}"     # e.g. ">s1>s3"
        
        self.start_node_id = start_node_id
        self.start_node_orient = start_node_orient
        self.end_node_id = end_node_id
        self.end_node_orient = end_node_orient
        
        self.ref_chrom = ref_chrom
        self.ref_start = ref_start
        self.ref_end = ref_end
        
        self.reversed_mapping = reversed_mapping
        self.ref_asm_path = None
        
        self.path_asm_dict = {}


def parse_gfa_nodes(gfa_path):
    """
    Reads a GFA file and extracts the ID and length for each node (S).
    Returns a dictionary { node_id: Node object }.
    
    Minigraph formats the GFA in a way that the node lines (S) contain the following fields (example):
    S	s1  ACGT    LN:i:4  SN:Z:chr1   S0:i:0  SR:i:0
    """
    logging.info(f"Scanning GFA for node metadata: {gfa_path}")
    nodes_dict = {}
    
    with open(gfa_path, 'r') as gfa_file:
        for line in gfa_file:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node_id = parts[1]
                sequence = parts[2]
                
                # use length from sequence if available, otherwise parse from LN tag
                if sequence != "*":
                    length = len(sequence)
                else:
                    length = 0
                    for tag in parts[3:]:
                        if tag.startswith("LN:i:"):
                            length = int(tag.split(':')[-1])
                            break

                nodes_dict[node_id] = Node(node_id, length)
    
    logging.info(f"Found {len(nodes_dict)} nodes in GFA.")
    return nodes_dict
