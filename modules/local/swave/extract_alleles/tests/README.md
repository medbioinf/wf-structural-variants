sample.bed: first 20 records of the pangenomesv pipeline's minigraph --call
output for HG002_hap1 (chm13v2.0 chr1 test dataset).

minigraph_gfa_fasta.fa: pangenome node fasta (chm13v2.0 chr1, minigraph
output), trimmed to only the node IDs referenced in sample.bed.

sample.vcf / ref.vcf: first 20 records of the pangenomesv pipeline's
vg-deconstruct output for the same chm13v2.0 chr1 test dataset, built with
pggb. sample.vcf is the per-sample split (HG002, with GT column), ref.vcf
is the reference split (no sample columns, -G).

pggb_gfa_fasta.fa: pangenome node fasta (chm13v2.0 chr1, pggb output),
trimmed to only the node IDs referenced in sample.vcf / ref.vcf.

Used to test SWAVE_EXTRACT_ALLELES with already-verified minigraph and
vg-deconstruct (pggb/cactus) pipeline output.
