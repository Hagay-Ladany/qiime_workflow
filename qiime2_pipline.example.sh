#!/bin/bash
# QIIME2 Analysis Pipeline for <bla bla bla>
# Author: Hagay Ladany | Date: <DATE>
# Project: <PROJ_NAME>
# seqID: X samples <seqID>;<projectID>

# Working in an interactive Docker env (example)
 docker run \
 -it -v "$(pwd)":/data \
 -v "$HOME/docker_bash_history.txt":/root/.bash_history \
 -e HISTFILE=/root/.bash_history \
 quay.io/qiime2/amplicon:2025.4 bash

## -----------------------------------------------------------------------------
## Initialize Environment
## -----------------------------------------------------------------------------
export TMPDIR="./tmpQIIME2"
mkdir -p ./tmpQIIME2

THREADS=0 # Use all cores
CLASSIFIER="reference_db/silva-138-99-nb-classifier.qza"

# Directory Structure
mkdir -p qiime_artifacts/{00_raw,01_denoise,02_phylogeny,03_taxonomy,04_diversity,05_differential_abundance,06_longitudinal}
mkdir -p qiime_visuals/{00_raw,01_denoise,03_taxonomy,04_diversity,05_differential_abundance,06_longitudinal}

# Path Shortcuts
ART="qiime_artifacts"
VIS="qiime_visuals"
META="sample_info/metadata.tsv" # Main metadata
META_FACTORIAL="sample_info/metadata_factor.tsv" # Metadata with Time as T1-T4, Location, CW_type

# V4 primers
F_primer="GTGCCAGCMGCCGCGGTAA"
R_primer="GGACTACHVGGGTWTCTAAT"

## -----------------------------------------------------------------------------
## 0. Download Classifier
## -----------------------------------------------------------------------------
echo "=== 0. DOWNLOAD CLASSIFIER ==="
 mkdir -p reference_db
 wget -O "${CLASSIFIER}" "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
 mv silva-138-99-nb-classifier.qza reference_db/ # Adjust if CLASSIFIER path is already correct
## -----------------------------------------------------------------------------
## 1. Data Import & Quality Control
## -----------------------------------------------------------------------------
echo "=== 1. DATA IMPORT & QC ==="
# 1a. Import raw data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "raw_fastq/FQ" \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path "${ART}/00_raw/demux_paired_end.qza"

# 1b. Demux summary
qiime demux summarize \
  --i-data "${ART}/00_raw/demux_paired_end.qza" \
  --o-visualization "${VIS}/00_raw/demux_summary.qzv"
echo "---"

## -----------------------------------------------------------------------------
## 2. Primer Removal
## -----------------------------------------------------------------------------
echo "=== 2. PRIMER TRIMMING ==="
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "${ART}/00_raw/demux_paired_end.qza" \
  --p-front-f ${F_primer} \
  --p-front-r ${R_primer} \
  --p-discard-untrimmed \
  --p-error-rate 0.1 \
  --p-overlap 3 \
  --o-trimmed-sequences "${ART}/01_denoise/trimmed_seqs.qza" \
  --p-cores ${THREADS}
  
qiime demux summarize \
  --i-data "${ART}/01_denoise/trimmed_seqs.qza" \
  --o-visualization "${VIS}/01_denoise/trimmed_summary.qzv"
echo "---"

## -----------------------------------------------------------------------------
## 3. Denoising with DADA2 (Paired-End)
## -----------------------------------------------------------------------------
echo "=== 3. DENOISING WITH DADA2 (Paired-End) ==="
# due to 2x150 sequencing and minimal overlap and high sequencing quality I will not trim sequences.
TRUNC_LEN_F=0
TRUNC_LEN_R=0
OVERLAP=8

TABLE_PE_OUT="${ART}/01_denoise/paired_end/feature_table.F${TRUNC_LEN_F}R${TRUNC_LEN_R}overlap${OVERLAP}.qza"
REPSEQS_PE_OUT="${ART}/01_denoise/paired_end/rep_seqs.F${TRUNC_LEN_F}R${TRUNC_LEN_R}overlap${OVERLAP}.qza"
STATS_PE_OUT="${ART}/01_denoise/paired_end/denoising_stats.F${TRUNC_LEN_F}R${TRUNC_LEN_R}overlap${OVERLAP}.qza"

mkdir -p "${ART}/01_denoise/paired_end/" "${VIS}/01_denoise/paired_end/"
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${ART}/01_denoise/trimmed_seqs.qza" \
  --p-trunc-len-f ${TRUNC_LEN_F} \
  --p-trunc-len-r ${TRUNC_LEN_R} \
  --p-min-overlap ${OVERLAP} \
  --p-n-threads ${THREADS} \
  --o-table "${TABLE_PE_OUT}" \
  --o-representative-sequences "${REPSEQS_PE_OUT}" \
  --o-denoising-stats "${STATS_PE_OUT}" \
  --verbose

qiime feature-table summarize \
  --i-table "${TABLE_PE_OUT}" \
  --o-visualization "${VIS}/01_denoise/paired_end/feature_table_summary.qzv" \
  --m-sample-metadata-file "${META}"

qiime feature-table tabulate-seqs \
  --i-data "${REPSEQS_PE_OUT}" \
  --o-visualization "${VIS}/01_denoise/paired_end/rep_seqs_summary.qzv"

qiime metadata tabulate \
  --m-input-file "${STATS_PE_OUT}" \
  --o-visualization "${VIS}/01_denoise/paired_end/denoising_stats_summary.qzv"
echo "---"

## -----------------------------------------------------------------------------
## 4. Phylogenetic Analysis (Paired-End)
## -----------------------------------------------------------------------------
echo "=== 4. PHYLOGENETIC TREE (Paired-End) ==="
ROOTED_TREE_PE_OUT="${ART}/02_phylogeny/paired_end/rooted_tree.qza"
mkdir -p "${ART}/02_phylogeny/paired_end/"

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${REPSEQS_PE_OUT}" \
  --o-alignment "${ART}/02_phylogeny/paired_end/aligned_seqs.qza" \
  --o-masked-alignment "${ART}/02_phylogeny/paired_end/masked_alignment.qza" \
  --o-tree "${ART}/02_phylogeny/paired_end/unrooted_tree.qza" \
  --o-rooted-tree "${ROOTED_TREE_PE_OUT}" \
  --p-n-threads ${THREADS}
echo "---"

## -----------------------------------------------------------------------------
## 5. Taxonomic Classification & Bar Plots (Paired-End)
## -----------------------------------------------------------------------------
echo "=== 5. TAXONOMY & BAR PLOTS (Paired-End) ==="
TAXONOMY_PE_OUT="${ART}/03_taxonomy/paired_end/taxonomy.qza"
mkdir -p "${ART}/03_taxonomy/paired_end/" "${VIS}/03_taxonomy/paired_end/"

qiime feature-classifier classify-sklearn \
  --i-reads "${REPSEQS_PE_OUT}" \
  --i-classifier "${CLASSIFIER}" \
  --p-n-jobs ${THREADS} \
  --o-classification "${TAXONOMY_PE_OUT}"
  
qiime metadata tabulate \
  --m-input-file "${TAXONOMY_PE_OUT}" \
  --o-visualization "${VIS}/03_taxonomy/paired_end/taxonomy_summary.qzv"

qiime taxa barplot \
  --i-table "${TABLE_PE_OUT}" \
  --i-taxonomy "${TAXONOMY_PE_OUT}" \
  --m-metadata-file "${META}" \
  --o-visualization "${VIS}/03_taxonomy/paired_end/taxa_barplot.qzv"
echo "---"

## -----------------------------------------------------------------------------
## 6. Core Diversity Analysis (Paired-End)
## -----------------------------------------------------------------------------
echo "=== 6. CORE DIVERSITY ANALYSIS (Paired-End) ==="
SAMPLING_DEPTH=70000

CORE_METRICS_PE_DIR="${ART}/04_diversity/paired_end/core_metrics"
mkdir -p "${VIS}/04_diversity/paired_end/" "${CORE_METRICS_PE_DIR}"

qiime diversity alpha-rarefaction \
  --i-table "${TABLE_PE_OUT}" \
  --i-phylogeny "${ROOTED_TREE_PE_OUT}" \
  --p-max-depth ${SAMPLING_DEPTH} \
  --p-metrics observed_features \
  --p-metrics shannon \
  --p-metrics faith_pd \
  --p-metrics simpson \
  --m-metadata-file "${META}" \
  --o-visualization "${VIS}/04_diversity/paired_end/alpha_rarefaction.qzv"

qiime diversity core-metrics-phylogenetic \
  --i-table "${TABLE_PE_OUT}" \
  --i-phylogeny "${ROOTED_TREE_PE_OUT}" \
  --p-sampling-depth ${SAMPLING_DEPTH} \
  --m-metadata-file "${META}" \
  --p-n-jobs-or-threads ${THREADS} \
  --output-dir "${CORE_METRICS_PE_DIR}"

for METRIC_VECTOR_QZA in "${CORE_METRICS_PE_DIR}"/*_vector.qza; do
  if [ -f "${METRIC_VECTOR_QZA}" ]; then
    METRIC_NAME=$(basename "${METRIC_VECTOR_QZA}" _vector.qza)
    qiime diversity alpha-group-significance \
        --i-alpha-diversity "${METRIC_VECTOR_QZA}" \
        --m-metadata-file "${META}" \
        --o-visualization "${VIS}/04_diversity/paired_end/${METRIC_NAME}_alpha_significance.qzv"
  fi
done

echo "--- Overall Beta Diversity PERMANOVA (Paired-End) ---"
UNWEIGHTED_UNIFRAC_DM_PAIRED="${CORE_METRICS_PE_DIR}/unweighted_unifrac_distance_matrix.qza"
WEIGHTED_UNIFRAC_DM_PAIRED="${CORE_METRICS_PE_DIR}/weighted_unifrac_distance_matrix.qza"

qiime diversity beta-group-significance \
  --i-distance-matrix "${UNWEIGHTED_UNIFRAC_DM_PAIRED}" \
  --m-metadata-file "${META}" \
  --m-metadata-column "CW_type" \
  --p-method permanova \
  --p-permutations 9999 \
  --o-visualization "${VIS}/04_diversity/paired_end/beta_permanova_overall_CW_type_unweighted.qzv"

for COLUMN in CW_type Location; do
  qiime diversity beta-group-significance \
    --i-distance-matrix "${WEIGHTED_UNIFRAC_DM_PAIRED}" \
    --m-metadata-file "${META}" \
    --m-metadata-column "${COLUMN}" \
    --p-method permanova \
    --p-permutations 9999 \
    --o-visualization "${VIS}/04_diversity/paired_end/beta_permanova_overall_${COLUMN}_weighted.qzv"
done

qiime diversity beta-group-significance \
  --i-distance-matrix "${WEIGHTED_UNIFRAC_DM_PAIRED}" \
  --m-metadata-file "${META_FACTORIAL}" \
  --m-metadata-column "Time" \
  --p-method permanova \
  --p-permutations 9999 \
  --o-visualization "${VIS}/04_diversity/paired_end/beta_permanova_overall_Time_weighted.qzv"
echo "---"
## -----------------------------------------------------------------------------
## 7. Stratified Beta Diversity Visualizations & Group Significance (Paired-End Unweighted UniFrac)
## -----------------------------------------------------------------------------
echo "=== 7. STRATIFIED BETA DIVERSITY & GROUP SIGNIFICANCE (Paired-End Unweighted UniFrac) ==="

VIS_REQUEST_1="${VIS}/04_diversity/paired_end/emperor_plots_Request1_PooledByTreatment"
ART_REQUEST_1="${ART}/04_diversity/paired_end/stratified_Req1_PooledByTreatment_artifacts"
mkdir -p "${VIS_REQUEST_1}" "${ART_REQUEST_1}"

VIS_REQUEST_2="${VIS}/04_diversity/paired_end/stratified_Req2_TreatmentByTime"
ART_REQUEST_2="${ART}/04_diversity/paired_end/stratified_Req2_TreatmentByTime_artifacts"
mkdir -p "${VIS_REQUEST_2}" "${ART_REQUEST_2}"

VIS_REQUEST_3="${VIS}/04_diversity/paired_end/stratified_Req3_PooledByLocation"
ART_REQUEST_3="${ART}/04_diversity/paired_end/stratified_Req3_PooledByLocation_artifacts"
mkdir -p "${VIS_REQUEST_3}" "${ART_REQUEST_3}"

TREATMENT_TYPES=("Control" "Treatment")
TIME_POINTS=("T1" "T2" "T3" "T4")
LOCATIONS=("A" "B" "C" "D" "E" "F" "G" "H" "I" "J") 

echo "--- Generating Custom Emperor Plots and PERMANOVA for Specific Groups ---"

# --- Request 1: PCoA Plots & Stats Pooled by Treatment ---
echo ""
echo "## Processing Request 1: Samples pooled by Treatment ##"
for treatment in "${TREATMENT_TYPES[@]}"; do
    echo "--- Group (Request 1): ${treatment} (All Time Points & Locations) ---"
    GROUP_NAME_TAG_REQ1="Req1_${treatment}_AllTime_AllLocations"
    FILTERED_DM_REQ1_OUT="${ART_REQUEST_1}/unweighted_unifrac_dm_${GROUP_NAME_TAG_REQ1}.qza"
    PCOA_REQ1_OUT="${ART_REQUEST_1}/unweighted_unifrac_pcoa_${GROUP_NAME_TAG_REQ1}.qza"
    EMPEROR_PLOT_REQ1_OUT="${VIS_REQUEST_1}/unweighted_unifrac_emperor_${GROUP_NAME_TAG_REQ1}.qzv"

    WHERE_CLAUSE_REQ1="[CW_type]='${treatment}'"
    echo "Applying filter: ${WHERE_CLAUSE_REQ1}"

    qiime diversity filter-distance-matrix \
        --i-distance-matrix "${UNWEIGHTED_UNIFRAC_DM_PAIRED}" \
        --m-metadata-file "${META_FACTORIAL}" \
        --p-where "${WHERE_CLAUSE_REQ1}" \
        --o-filtered-distance-matrix "${FILTERED_DM_REQ1_OUT}"
    
    qiime diversity pcoa \
        --i-distance-matrix "${FILTERED_DM_REQ1_OUT}" \
        --o-pcoa "${PCOA_REQ1_OUT}"
    
    qiime emperor plot \
        --i-pcoa "${PCOA_REQ1_OUT}" \
        --m-metadata-file "${META_FACTORIAL}" \
        --o-visualization "${EMPEROR_PLOT_REQ1_OUT}"
    
    echo "--- PERMANOVA for Time effect in ${treatment} group (Unweighted UniFrac) ---"
    qiime diversity beta-group-significance \
      --i-distance-matrix "${FILTERED_DM_REQ1_OUT}" \
      --m-metadata-file "${META_FACTORIAL}" \
      --m-metadata-column "Time" \
      --o-visualization "${VIS_REQUEST_1}/permanova_${GROUP_NAME_TAG_REQ1}_vs_Time.qzv" \
      --p-method permanova \
      --p-pairwise \
      --p-permutations 9999

    echo "--- PERMANOVA for Location effect in ${treatment} group (Unweighted UniFrac) ---"
    qiime diversity beta-group-significance \
      --i-distance-matrix "${FILTERED_DM_REQ1_OUT}" \
      --m-metadata-file "${META_FACTORIAL}" \
      --m-metadata-column "Location" \
      --o-visualization "${VIS_REQUEST_1}/permanova_${GROUP_NAME_TAG_REQ1}_vs_Location.qzv" \
      --p-method permanova \
      --p-pairwise \
      --p-permutations 9999
    echo "-----------------------------------------------------"
done
echo "## Finished Request 1 ##"


# --- Request 2: PCoA Plots Stratified by Treatment AND Time Point ---
echo ""
echo "## Processing Request 2: PCoA for each Treatment-Time stratum ##"
for treatment in "${TREATMENT_TYPES[@]}"; do
    for time_point in "${TIME_POINTS[@]}"; do
        echo "--- Group (Request 2): ${treatment} at ${time_point} ---"
        GROUP_NAME_TAG_REQ2="Req2_${treatment}_${time_point}"
        FILTERED_DM_REQ2_OUT="${ART_REQUEST_2}/unweighted_unifrac_dm_${GROUP_NAME_TAG_REQ2}.qza"
        PCOA_REQ2_OUT="${ART_REQUEST_2}/unweighted_unifrac_pcoa_${GROUP_NAME_TAG_REQ2}.qza"
        EMPEROR_PLOT_REQ2_OUT="${VIS_REQUEST_2}/unweighted_unifrac_emperor_${GROUP_NAME_TAG_REQ2}.qzv"

        # Note: This clause filters to a single Treatment and Time, useful for visualizing that specific group's spread
        WHERE_CLAUSE_REQ2="[CW_type]='${treatment}' AND [Time]='${time_point}'"
        echo "Applying filter: ${WHERE_CLAUSE_REQ2}"

        qiime diversity filter-distance-matrix \
            --i-distance-matrix "${UNWEIGHTED_UNIFRAC_DM_PAIRED}" \
            --m-metadata-file "${META_FACTORIAL}" \
            --p-where "${WHERE_CLAUSE_REQ2}" \
            --o-filtered-distance-matrix "${FILTERED_DM_REQ2_OUT}"
        
        # Only proceed if the filtered matrix is not empty
        if qiime tools validate "${FILTERED_DM_REQ2_OUT}" &> /dev/null; then
            qiime diversity pcoa \
                --i-distance-matrix "${FILTERED_DM_REQ2_OUT}" \
                --o-pcoa "${PCOA_REQ2_OUT}"
            
            qiime emperor plot \
                --i-pcoa "${PCOA_REQ2_OUT}" \
                --m-metadata-file "${META_FACTORIAL}" \
                --o-visualization "${EMPEROR_PLOT_REQ2_OUT}"
        else
            echo "WARNING: No samples found for filter '${WHERE_CLAUSE_REQ2}'. Skipping PCoA and Emperor plot."
        fi
        echo "-----------------------------------------------------"
    done
done
echo "## Finished Request 2 ##"


# --- Request 2b: Comparison Stratified ONLY by Time (Unweighted & Weighted UniFrac) ---
echo ""
echo "## Processing Request 2b: Control vs. Treatment comparison within each Time Point ##"

# Create dedicated directories for this new analysis
VIS_REQUEST_2B="${VIS}/04_diversity/paired_end/stratified_Req2b_TimeOnly"
ART_REQUEST_2B="${ART}/04_diversity/paired_end/stratified_Req2b_TimeOnly"
mkdir -p "${VIS_REQUEST_2B}" "${ART_REQUEST_2B}"

# Define the diversity metrics to loop through
BETA_METRICS=("unweighted_unifrac" "weighted_unifrac")

for metric in "${BETA_METRICS[@]}"; do
    echo "=== Processing Metric: ${metric} ==="
    
    # Select the correct input distance matrix based on the metric
    if [ "$metric" == "unweighted_unifrac" ]; then
        INPUT_DM="${UNWEIGHTED_UNIFRAC_DM_PAIRED}"
    else
        INPUT_DM="${WEIGHTED_UNIFRAC_DM_PAIRED}"
    fi

    for time_point in "${TIME_POINTS[@]}"; do
        echo "--- Group (Request 2b): All samples at ${time_point} using ${metric} ---"
        GROUP_NAME_TAG_REQ2B="Req2b_${time_point}_${metric}"
        FILTERED_DM_REQ2B_OUT="${ART_REQUEST_2B}/dm_${GROUP_NAME_TAG_REQ2B}.qza"
        PCOA_REQ2B_OUT="${ART_REQUEST_2B}/pcoa_${GROUP_NAME_TAG_REQ2B}.qza"
        EMPEROR_PLOT_REQ2B_OUT="${VIS_REQUEST_2B}/emperor_${GROUP_NAME_TAG_REQ2B}_by_CW_type.qzv"

        # Filter for all samples belonging to the current time point
        WHERE_CLAUSE_REQ2B="[Time]='${time_point}'"
        echo "Applying filter: ${WHERE_CLAUSE_REQ2B} for ${metric}"

        qiime diversity filter-distance-matrix \
            --i-distance-matrix "${INPUT_DM}" \
            --m-metadata-file "${META_FACTORIAL}" \
            --p-where "${WHERE_CLAUSE_REQ2B}" \
            --o-filtered-distance-matrix "${FILTERED_DM_REQ2B_OUT}"
        
        echo "--- Generating PCoA plot for ${time_point} (colored by CW_type) ---"
        qiime diversity pcoa \
            --i-distance-matrix "${FILTERED_DM_REQ2B_OUT}" \
            --o-pcoa "${PCOA_REQ2B_OUT}"
        
        qiime emperor plot \
            --i-pcoa "${PCOA_REQ2B_OUT}" \
            --m-metadata-file "${META_FACTORIAL}" \
            --o-visualization "${EMPEROR_PLOT_REQ2B_OUT}"
            
        echo "--- PERMANOVA for CW_type effect within ${time_point} (${metric}) ---"
        qiime diversity beta-group-significance \
          --i-distance-matrix "${FILTERED_DM_REQ2B_OUT}" \
          --m-metadata-file "${META_FACTORIAL}" \
          --m-metadata-column "CW_type" \
          --o-visualization "${VIS_REQUEST_2B}/permanova_${GROUP_NAME_TAG_REQ2B}_vs_CW_type.qzv" \
          --p-method permanova \
          --p-permutations 9999

        echo "-----------------------------------------------------"
    done
done
echo "## Finished Request 2b ##"






# --- Request 3: PCoA Plots Stratified by Location (pooling Control & Treatment) ---
echo ""
echo "## Processing Request 3: Samples pooled by Location (Control & Treatment combined) ##"
for location_val in "${LOCATIONS[@]}"; do
    location_val_fn=$(echo "${location_val}" | sed 's/[^a-zA-Z0-9_-]/_/g')
    
    echo "--- Group (Request 3): Location '${location_val}' (Control & Treatment pooled) ---"
    GROUP_NAME_TAG_REQ3="Req3_Location_${location_val_fn}_CandTpooled"
    FILTERED_DM_REQ3_OUT="${ART_REQUEST_3}/unweighted_unifrac_dm_${GROUP_NAME_TAG_REQ3}.qza"
    PCOA_REQ3_OUT="${ART_REQUEST_3}/unweighted_unifrac_pcoa_${GROUP_NAME_TAG_REQ3}.qza"
    EMPEROR_PLOT_REQ3_OUT="${VIS_REQUEST_3}/unweighted_unifrac_emperor_${GROUP_NAME_TAG_REQ3}.qzv"

    WHERE_CLAUSE_REQ3="[Location]='${location_val}'"
    echo "Applying filter: ${WHERE_CLAUSE_REQ3}"

    qiime diversity filter-distance-matrix \
        --i-distance-matrix "${UNWEIGHTED_UNIFRAC_DM_PAIRED}" \
        --m-metadata-file "${META_FACTORIAL}" \
        --p-where "${WHERE_CLAUSE_REQ3}" \
        --o-filtered-distance-matrix "${FILTERED_DM_REQ3_OUT}"
        
    qiime diversity pcoa \
        --i-distance-matrix "${FILTERED_DM_REQ3_OUT}" \
        --o-pcoa "${PCOA_REQ3_OUT}"

    qiime emperor plot \
        --i-pcoa "${PCOA_REQ3_OUT}" \
        --m-metadata-file "${META_FACTORIAL}" \
        --o-visualization "${EMPEROR_PLOT_REQ3_OUT}"
    echo "-----------------------------------------------------"
done
echo "## Finished Request 3 ##"
echo "---"

## -----------------------------------------------------------------------------
## 8. Advanced Analyses - Specific Research Questions
## -----------------------------------------------------------------------------
echo "=== 8. ADVANCED ANALYSES - SPECIFIC RESEARCH QUESTIONS ==="

SUB_VIS_PAIRED="${VIS}/04_diversity/paired_end/advanced_beta_stratified_analysis"
SUB_ART_PAIRED="${ART}/04_diversity/paired_end/advanced_beta_stratified_analysis_artifacts"
LONG_VIS_PAIRED="${VIS}/06_longitudinal/paired_end/"
LONG_ART_PAIRED="${ART}/06_longitudinal/paired_end/"

mkdir -p "${SUB_VIS_PAIRED}" "${SUB_ART_PAIRED}"
mkdir -p "${LONG_VIS_PAIRED}" "${LONG_ART_PAIRED}"

TREATMENT_META_FACTOR_ADV="sample_info/treatment_metadata_factor_advanced.tsv" # New name to avoid conflict
echo "--- Preparing treatment-only metadata with categorical Time (for Step 8) ---"
awk -F'\t' 'BEGIN{OFS="\t"} NR==1 || $2=="Treatment" {print}' "${META_FACTORIAL}" > "${TREATMENT_META_FACTOR_ADV}"
echo "Created ${TREATMENT_META_FACTOR_ADV}"

TREATMENT_META_CONTINUOUS="sample_info/treatment_metadata_continuous_time.tsv"
echo "--- Preparing treatment-only metadata with continuous Time (from main META file) ---"
head -n 1 "${META}" > "${TREATMENT_META_CONTINUOUS}"
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $2=="Treatment" {print}' "${META}" >> "${TREATMENT_META_CONTINUOUS}"
echo "Created ${TREATMENT_META_CONTINUOUS}"

TREATMENT_META_MIN_TWO_GROUPS="sample_info/treatment_metadata_factor_with_min_two_groups.tsv"
echo "--- Preparing treatment-only metadata with dummy AnalysisGroup for pairwise-distances ---"
awk -F'\t' 'BEGIN{OFS="\t"}
            NR==1 {print $0, "AnalysisGroup"; next} 
            {
                current_location = $4; 
                group_val = "OtherDummyGroup"; 
                if (current_location ~ /^[A-E]$/) { group_val = "DummyGroup1"; }
                else if (current_location ~ /^[F-J]$/) { group_val = "DummyGroup2"; }
                print $0, group_val;
            }' "${TREATMENT_META_FACTOR_ADV}" > "${TREATMENT_META_MIN_TWO_GROUPS}" # Use the newly created one
echo "Created ${TREATMENT_META_MIN_TWO_GROUPS}"
echo ""

echo "--- Q.a: Beta Diversity - Effect of Time within Treatment & Control Groups (Paired-End, Weighted UniFrac) ---"
FILTERED_DIST_PAIRED_CW_TREATMENT="${SUB_ART_PAIRED}/dist_paired_weighted_unifrac_CW_Treatment.qza"
qiime diversity filter-distance-matrix \
  --i-distance-matrix "${WEIGHTED_UNIFRAC_DM_PAIRED}" \
  --m-metadata-file "${META}" \
  --p-where "[CW_type]='Treatment'" \
  --o-filtered-distance-matrix "${FILTERED_DIST_PAIRED_CW_TREATMENT}"

qiime diversity beta-group-significance \
  --i-distance-matrix "${FILTERED_DIST_PAIRED_CW_TREATMENT}" \
  --m-metadata-file "${META_FACTORIAL}" \
  --m-metadata-column "Time" \
  --p-method permanova \
  --p-pairwise \
  --p-permutations 9999 \
  --o-visualization "${SUB_VIS_PAIRED}/beta_permanova_Time_within_Treatment_weighted.qzv"

FILTERED_DIST_PAIRED_CW_CONTROL="${SUB_ART_PAIRED}/dist_paired_weighted_unifrac_CW_Control.qza"
qiime diversity filter-distance-matrix \
  --i-distance-matrix "${WEIGHTED_UNIFRAC_DM_PAIRED}" \
  --m-metadata-file "${META}" \
  --p-where "[CW_type]='Control'" \
  --o-filtered-distance-matrix "${FILTERED_DIST_PAIRED_CW_CONTROL}"

qiime diversity beta-group-significance \
  --i-distance-matrix "${FILTERED_DIST_PAIRED_CW_CONTROL}" \
  --m-metadata-file "${META_FACTORIAL}" \
  --m-metadata-column "Time" \
  --p-method permanova \
  --p-pairwise \
  --p-permutations 9999 \
  --o-visualization "${SUB_VIS_PAIRED}/beta_permanova_Time_within_Control_weighted.qzv"

echo "--- Q.b: ANCOM-BC - Microbes changing over Categorical Time in Treatment (Paired-End) ---"
FEATURE_TABLE_TREATMENT_ONLY="${SUB_ART_PAIRED}/feature_table_TreatmentOnly.qza"
qiime feature-table filter-samples \
  --i-table "${TABLE_PE_OUT}" \
  --m-metadata-file "${TREATMENT_META_FACTOR_ADV}" \
  --o-filtered-table "${FEATURE_TABLE_TREATMENT_ONLY}"
  
qiime composition ancombc \
  --i-table "${FEATURE_TABLE_TREATMENT_ONLY}" \
  --m-metadata-file "${TREATMENT_META_FACTOR_ADV}" \
  --p-formula "Time" \
  --o-differentials "${SUB_ART_PAIRED}/ancombc_Treatment_CategoricalTime.qza"
  
qiime composition tabulate \
  --i-data "${SUB_ART_PAIRED}/ancombc_Treatment_CategoricalTime.qza" \
  --o-visualization "${SUB_VIS_PAIRED}/ancombc_Treatment_CategoricalTime_tabular.qzv"

echo "--- Q.c: Longitudinal Beta Diversity - Plant Resilience/Proneness (Paired-End, Weighted UniFrac) ---"
qiime longitudinal pairwise-distances \
  --i-distance-matrix "${FILTERED_DIST_PAIRED_CW_TREATMENT}" \
  --m-metadata-file "${TREATMENT_META_MIN_TWO_GROUPS}" \
  --p-state-column Time \
  --p-individual-id-column Location \
  --p-state-1 T1 \
  --p-state-2 T4 \
  --p-group-column AnalysisGroup \
  --o-visualization "${LONG_VIS_PAIRED}/pairwise_distances_Treatment_T1vsT4_by_Location_weighted.qzv"

FEATURE_TABLE_TREATMENT_ONLY_RELATIVE="${SUB_ART_PAIRED}/feature_table_TreatmentOnly_relative.qza"
qiime feature-table relative-frequency \
  --i-table "${FEATURE_TABLE_TREATMENT_ONLY}" \
  --o-relative-frequency-table "${FEATURE_TABLE_TREATMENT_ONLY_RELATIVE}"

qiime longitudinal volatility \
  --i-table "${FEATURE_TABLE_TREATMENT_ONLY_RELATIVE}" \
  --m-metadata-file "${TREATMENT_META_FACTOR_ADV}" \
  --p-state-column Time \
  --p-individual-id-column Location \
  --p-default-metric 'weighted_unifrac' \
  --p-default-group-column Location \
  --o-visualization "${LONG_VIS_PAIRED}/volatility_Treatment_Time_by_Location.qzv"

# Model weighted-UniFrac change across ALL time points (T1-T4)
# distance from each sample to its own T1 (baseline)
qiime longitudinal first-distances \
  --i-distance-matrix  "qiime_artifacts/04_diversity/paired_end/beta_stratified_analysis_artifacts/dist_paired_CW_Treatment.qza" \
  --m-metadata-file    "${TREATMENT_META_CONTINUOUS}" \
  --p-state-column     Time \
  --p-individual-id-column Location \
  --p-baseline         1 \
  --o-first-distances  "${LONG_ART_PAIRED}/wuf_to_T1.qza"

# LME on weighted-UniFrac distance-to-baseline (all four times)
qiime longitudinal linear-mixed-effects \
  --m-metadata-file  "${LONG_ART_PAIRED}/wuf_to_T1.qza" \
  --m-metadata-file  "${TREATMENT_META_CONTINUOUS}" \
  --p-metric         Distance \
  --p-state-column   Time \
  --p-individual-id-column Location \
  --o-visualization  "${LONG_VIS_PAIRED}/lme_wuf_all_times.qzv"


echo "--- General ANCOM (All Samples, Paired-End) ---"
TABLE_PAIR_ANCOM_FILTERED="${ART}/05_differential_abundance/paired_end/table_ancom_filtered_full.qza"
COMPOSITION_TABLE_PAIR_FULL="${ART}/05_differential_abundance/paired_end/composition_table_full.qza"
mkdir -p "${ART}/05_differential_abundance/paired_end/" "${VIS}/05_differential_abundance/paired_end/"

qiime feature-table filter-features \
  --i-table "${TABLE_PE_OUT}" \
  --p-min-frequency 100 \
  --p-min-samples 2 \
  --o-filtered-table "${TABLE_PAIR_ANCOM_FILTERED}"
  
qiime composition add-pseudocount \
  --i-table "${TABLE_PAIR_ANCOM_FILTERED}" \
  --o-composition-table "${COMPOSITION_TABLE_PAIR_FULL}"

qiime composition ancom \
  --i-table "${COMPOSITION_TABLE_PAIR_FULL}" \
  --m-metadata-file "${META}" \
  --m-metadata-column CW_type \
  --o-visualization "${VIS}/05_differential_abundance/paired_end/ancom_overall_CW_type.qzv"

qiime composition ancom \
  --i-table "${COMPOSITION_TABLE_PAIR_FULL}" \
  --m-metadata-file "${META}" \
  --m-metadata-column Location \
  --o-visualization "${VIS}/05_differential_abundance/paired_end/ancom_overall_Location.qzv"

qiime composition ancom \
  --i-table "${COMPOSITION_TABLE_PAIR_FULL}" \
  --m-metadata-file "${META_FACTORIAL}" \
  --m-metadata-column Time \
  --o-visualization "${VIS}/05_differential_abundance/paired_end/ancom_overall_Time_categorical.qzv"
echo "---"

echo "=== FULL ANALYSIS SCRIPT LOGIC COMPLETED ==="
