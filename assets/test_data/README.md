# Pipeline Test Dataset

This directory contains a minimal, test dataset used for testing. 

To keep the pipeline lightweight and fast, the dataset represents a **~5 Megabase region (approx. 2% of Chromosome 1)**.

## Dataset Specifications

* **Reference Genome:** `mini_chm13v2.0.fa`
  * **Source:** T2T-CHM13v2.0
  * **Region:** `chr1:1-5000000` (5 MB)
* **Sample:** `HG002`
* **Haplotype 1 (Maternal):** `mini_HG002_hap1.fa`
  * **Original Contig:** `HG002#2#h2tg000171l` (Region: `1-5000000`)
  * *Note: This specific assembly contig matches the beginning of chromosome 1 in the reference.*
* **Haplotype 2 (Paternal):** `mini_HG002_hap2.fa`
  * **Original Contig:** `HG002#1#h1tg000219l` (Region: `1-5000000`)
  * *Note: This specific assembly contig matches the beginning of chromosome 1 in the reference.*

## Generation Command

The dataset was sliced using `samtools` (the sequence headers were cleaned using `sed` to prevent naming errors):

```bash
# Example for Haplotype 1
samtools faidx data/assemblies/HG002.maternal.f1_assembly_v2.fa HG002#2#h2tg000171l:1-5000000 | sed 's/>HG002#2#h2tg000171l:1-5000000/>HG002#2#h2tg000171l/' > assets/test_data/assemblies/mini_HG002_hap1.fa