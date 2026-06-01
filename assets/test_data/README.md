# Pipeline Test Dataset

This directory contains a minimal, optimized human test dataset used for CI/CD testing on GitHub Actions.

To keep the pipeline lightweight and fast, the dataset represents a **~1.2 Megabase region** at the beginning of Chromosome 1.

## Dataset Specifications

* **Reference Genome:** `mini_chm13v2.0.fa`
  * **Source:** T2T-CHM13v2.0
  * **Region:** `chr1:1-1200000` (~1.2 MB / 15,000 lines)
* **Sample:** `pangenome_test`
* **Haplotype 1 (Maternal):** `mini_HG002_hap1.fa`
  * **Original Contig:** `HG002#2#h2tg000171l` (Full contig, Region approx.: `1-1200000`)
  * *Note: This specific assembly contig perfectly matches the beginning of chromosome 1 in the reference.*
* **Haplotype 2 (Paternal):** `mini_HG002_hap2.fa`
  * **Original Contig:** `HG002#1#h1tg000219l` (Full contig, Region approx.: `1-1200000`)
  * *Note: This specific assembly contig perfectly matches the beginning of chromosome 1 in the reference.*

## Generation Command

The dataset was sliced using `samtools` to ensure all files share the exact same biological length. Sequence headers were cleaned using `sed` to prevent naming errors during the graph construction:

```bash
# Example for slicing the reference genome
samtools faidx data/references/chm13v2.0.fa chr1:1-1200000 | sed 's/>chr1:1-1200000/>chr1/' > assets/test_data/reference/mini_chm13v2.0.fa

# Example for slicing Haplotype 1
samtools faidx data/assemblies/HG002.maternal.f1_assembly_v2.fa HG002#2#h2tg000171l:1-1200000 | sed 's/>HG002#2#h2tg000171l:1-1200000/>HG002#2#h2tg000171l/' > assets/test_data/assemblies/mini_HG002_hap1.fa