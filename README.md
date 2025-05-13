# QIIME2 Paired-End Microbiome Analysis Workflow (with ANCOM Differential Abundance)

## Overview
This workflow performs paired-end 16S rRNA microbiome analysis using QIIME 2. It includes importing FASTQ files, optional trimming, denoising with DADA2, diversity metrics (alpha/beta), taxonomy classification, and differential abundance analysis using ANCOM at both ASV and genus levels.

## 0. Environment & Folder Setup
```bash
conda activate qiime2-amplicon-XXXX
export TMPDIR="$PWD/tmpQIIME2"

mkdir -p \
  raw_fastq/FQ \
  qc \
  qiime_artifacts/{import,cutadapt,denoise,phylogeny,core_metrics,taxonomy,ancom,reference_db} \
  qiime_visuals/{demux,cutadapt,denoise,core_metrics,taxonomy,alpha_beta_tests,ancom} \
  sample_info
```

## 1. Import Paired-End Reads
```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path raw_fastq/FQ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path qiime_artifacts/import/rawData.qza

qiime demux summarize \
  --i-data qiime_artifacts/import/rawData.qza \
  --o-visualization qiime_visuals/demux/rawData.qzv
```

## 2. Optional Primer & Quality Trimming
```bash
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences qiime_artifacts/import/rawData.qza \
  --p-minimum-length 200 \
  --p-quality-cutoff-5end 20 \
  --p-quality-cutoff-3end 20 \
  --o-trimmed-sequences qiime_artifacts/cutadapt/qcData.qza \
  --verbose &> qiime_artifacts/cutadapt/qcData.log

qiime demux summarize \
  --i-data qiime_artifacts/cutadapt/qcData.qza \
  --o-visualization qiime_visuals/cutadapt/qcData.qzv
```

## 3. DADA2 Denoising
```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime_artifacts/import/rawData.qza \  # or qcData.qza
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-trunc-len-f 0 --p-trunc-len-r 0 \
  --p-n-threads 0 \
  --output-dir qiime_artifacts/denoise

qiime feature-table summarize \
  --i-table qiime_artifacts/denoise/table.qza \
  --o-visualization qiime_visuals/denoise/trim_table.qzv

qiime metadata tabulate \
  --m-input-file qiime_artifacts/denoise/denoising_stats.qza \
  --o-visualization qiime_visuals/denoise/trim_stats.qzv
```

## 4. Phylogenetic Tree
```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences qiime_artifacts/denoise/representative_sequences.qza \
  --output-dir qiime_artifacts/phylogeny
```

## 5. Core Diversity Metrics
```bash
depth=<CHOSEN_DEPTH>

qiime diversity core-metrics-phylogenetic \
  --i-table qiime_artifacts/denoise/table.qza \
  --i-phylogeny qiime_artifacts/phylogeny/rooted_tree.qza \
  --p-sampling-depth $depth \
  --m-metadata-file sample_info/metadata.tsv \
  --output-dir qiime_artifacts/core_metrics

qiime diversity alpha \
  --i-table qiime_artifacts/denoise/table.qza \
  --p-metric simpson \
  --o-alpha-diversity qiime_artifacts/core_metrics/simpson_vector.qza

for metric in shannon simpson evenness faith_pd; do
  qiime diversity alpha-group-significance \
    --i-alpha-diversity qiime_artifacts/core_metrics/${metric}_vector.qza \
    --m-metadata-file sample_info/metadata.tsv \
    --o-visualization qiime_visuals/core_metrics/${metric}-group.qzv
done

for dist in unweighted_unifrac weighted_unifrac bray_curtis jaccard; do
  qiime diversity beta-group-significance \
    --i-distance-matrix qiime_artifacts/core_metrics/${dist}_distance_matrix.qza \
    --m-metadata-file sample_info/metadata.tsv \
    --m-metadata-column <GROUP_COLUMN> \
    --p-pairwise \
    --o-visualization qiime_visuals/alpha_beta_tests/${dist}_GROUP.qzv
done
```

## 6. Taxonomic Classification
```bash
qiime feature-classifier classify-sklearn \
  --i-reads qiime_artifacts/denoise/representative_sequences.qza \
  --i-classifier qiime_artifacts/reference_db/silva-*.qza \
  --p-n-jobs 32 \
  --o-classification qiime_artifacts/taxonomy/taxonomy.qza

qiime taxa barplot \
  --i-table qiime_artifacts/denoise/table.qza \
  --i-taxonomy qiime_artifacts/taxonomy/taxonomy.qza \
  --m-metadata-file sample_info/metadata.tsv \
  --o-visualization qiime_visuals/taxonomy/taxa_barplot.qzv
```

## 7. Differential Abundance with ANCOM
```bash
# Filter rare features
qiime feature-table filter-features \
  --i-table qiime_artifacts/denoise/table.qza \
  --p-min-frequency 100 --p-min-samples 5 \
  --o-filtered-table qiime_artifacts/ancom/table_abund.qza

# Add pseudocount
qiime composition add-pseudocount \
  --i-table qiime_artifacts/ancom/table_abund.qza \
  --o-composition-table qiime_artifacts/ancom/table_abund_comp.qza

# ANCOM (ASV)
qiime composition ancom \
  --i-table qiime_artifacts/ancom/table_abund_comp.qza \
  --m-metadata-file sample_info/metadata.tsv \
  --m-metadata-column <GROUP_COLUMN> \
  --o-visualization qiime_visuals/ancom/ancom_ASV_GROUP.qzv

# Collapse to genus
qiime taxa collapse \
  --i-table qiime_artifacts/ancom/table_abund.qza \
  --i-taxonomy qiime_artifacts/taxonomy/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table qiime_artifacts/ancom/genus_table.qza

qiime composition add-pseudocount \
  --i-table qiime_artifacts/ancom/genus_table.qza \
  --o-composition-table qiime_artifacts/ancom/genus_table_comp.qza

# ANCOM (Genus)
qiime composition ancom \
  --i-table qiime_artifacts/ancom/genus_table_comp.qza \
  --m-metadata-file sample_info/metadata.tsv \
  --m-metadata-column <GROUP_COLUMN> \
  --o-visualization qiime_visuals/ancom/ancom_Genus_GROUP.qzv
```

## 8. Review & Export
Open `.qzv` files in [QIIME 2 View](https://view.qiime2.org) for interactive plots and summaries.
