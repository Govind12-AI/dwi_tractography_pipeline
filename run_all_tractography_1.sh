#!/bin/bash


# ============================================================================
# NATIVE SPACE TRACTOGRAPHY & CONNECTOME PIPELINE - PRODUCTION VERSION
# 
# Purpose: Generate anatomically-constrained tractography and structural 
#          connectivity matrices for 35 HCP-style subjects in native DWI space
#
# Key Features:
# - Multi-tissue intensity normalization for accurate FOD amplitudes
# - Processes data entirely in native DWI space (no MNI space tractography)
# - Uses 5 million streamlines per subject for dense connectivity mapping
# - SIFT2 filtering for biologically accurate connection strengths
# - Parallel processing of 3 subjects simultaneously
# - Labeled connectivity matrices with Morel atlas regions
# - Centralized connectivity matrix storage for easy access
# 
# Pipeline Overview:
# 1. Brain extraction and 5TT generation from T1
# 2. Multi-tissue intensity normalization (mtnormalise)
# 3. Registration: T1->DWI and MNI->T1 transformations
# 4. Warp atlas and seed regions to native DWI space
# 5. Anatomically-constrained tractography (ACT) with 5M streamlines
# 6. SIFT2 filtering for biological accuracy
# 7. Generate labeled weighted connectivity matrices
# 8. Copy matrices to central repository
# ============================================================================

# Enable strict error handling for robust pipeline execution
set -euxo pipefail

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================

# Base paths - adjust these to match your system
BASE_DIR="/home/kaushik/HCP_DWI_data"
DATA_DIR="${BASE_DIR}/data"                    # Raw HCP data location
PREPROCESSED_DIR="${BASE_DIR}/outputs"          # DWI preprocessing outputs
TRACTOGRAPHY_BASE="${BASE_DIR}/Tractography_analysis"  # Main analysis directory

# Atlas and seed files (in MNI space - will be warped to native space)
ATLAS_FILE="${BASE_DIR}/morel_unified_atlas/1mm/Morel_Unified_Atlas_1mm.nii.gz"
SEED_FILE="${BASE_DIR}/seed_mask_pain.nii"
ATLAS_LABELS="${BASE_DIR}/morel_unified_atlas/morel_labels_simple.txt"

# Subject list - all 35 subjects from mgh_1001 to mgh_1035
SUBJECTS=()
for i in $(seq -f "%04g" 1001 1035); do
    SUBJECTS+=("mgh_${i}")
done

# Processing parameters - optimized for 24-core system with 62GB RAM
PARALLEL_JOBS=3      # Number of subjects to process simultaneously
NTHREADS=7           # Threads per subject (3*7=21 cores used, 3 for system)
STREAMLINES=3000000  # 3 million streamlines per subject

# ============================================================================
# SETUP MAIN ANALYSIS DIRECTORY STRUCTURE
# ============================================================================

echo "============================================================================"
echo "NATIVE SPACE CONNECTIVITY ANALYSIS - PRODUCTION RUN"
echo "Analysis directory: ${TRACTOGRAPHY_BASE}"
echo "Number of subjects: ${#SUBJECTS[@]}"
echo "Streamlines per subject: ${STREAMLINES}"
echo "Parallel jobs: ${PARALLEL_JOBS}"
echo "Threads per subject: ${NTHREADS}"
echo "System: 24 cores, 62GB RAM"
echo "============================================================================"

# Create organized directory structure
mkdir -p "${TRACTOGRAPHY_BASE}"/{scripts,logs,reference,connectomes,subjects}

# Define key directories for easy reference
SCRIPTS_DIR="${TRACTOGRAPHY_BASE}/scripts"
LOGS_DIR="${TRACTOGRAPHY_BASE}/logs"
REF_DIR="${TRACTOGRAPHY_BASE}/reference"
CONNECTOMES_DIR="${TRACTOGRAPHY_BASE}/connectomes"
SUBJECTS_DIR="${TRACTOGRAPHY_BASE}/subjects"


# ============================================================================
# CHECK FOR ALREADY COMPLETED SUBJECTS (RESUME CAPABILITY)
# ============================================================================

echo "============================================================================"
echo "CHECKING FOR PREVIOUSLY COMPLETED SUBJECTS"
echo "============================================================================"

# Arrays to track subject status
COMPLETED_SUBJECTS=()
REMAINING_SUBJECTS=()

# Check which subjects are already done
for subject in "${SUBJECTS[@]}"; do
    # Check if connectivity matrix exists (main indicator of completion)
    if [[ -f "${CONNECTOMES_DIR}/${subject}_connectivity_matrix.csv" ]] && \
       [[ -f "${CONNECTOMES_DIR}/${subject}_connectivity_matrix_labeled.csv" ]]; then
        COMPLETED_SUBJECTS+=("${subject}")
        echo "✓ ${subject} - Already completed, skipping"
    else
        REMAINING_SUBJECTS+=("${subject}")
        echo "○ ${subject} - Will process"
    fi
done

echo ""
echo "Summary:"
echo "  - Completed: ${#COMPLETED_SUBJECTS[@]}/${#SUBJECTS[@]} subjects"
echo "  - Remaining: ${#REMAINING_SUBJECTS[@]} subjects to process"

# Exit if all subjects are done
if [[ ${#REMAINING_SUBJECTS[@]} -eq 0 ]]; then
    echo "✓ All subjects have been successfully processed!"
    exit 0
fi

# ============================================================================
# PREPARE REFERENCE FILES
# ============================================================================

echo "Setting up reference files for analysis..."

# Copy MNI template to reference directory for future visualization/overlay purposes
FSL_MNI_PATH="${HOME}/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"
MNI_TEMPLATE_TARGET="${REF_DIR}/MNI152_T1_1mm_brain.nii.gz"

if [[ -f "${FSL_MNI_PATH}" ]]; then
    echo "✓ Found 1mm MNI template. Copying to reference directory."
    cp -f "${FSL_MNI_PATH}" "${MNI_TEMPLATE_TARGET}"
else
    echo "ERROR: 1mm MNI template not found at ${FSL_MNI_PATH}"
    echo "Please verify FSL installation and path."
    exit 1
fi

# Copy atlas, seed, and label files to reference directory
cp -f "${ATLAS_FILE}" "${REF_DIR}/Morel_Atlas.nii.gz"
cp -f "${SEED_FILE}" "${REF_DIR}/seed_mask.nii"
cp -f "${ATLAS_LABELS}" "${REF_DIR}/morel_labels.txt"

# Create node lookup file in MRtrix format (without header for tck2connectome)
tail -n +2 "${ATLAS_LABELS}" | awk '{print $1 "\t" $2}' > "${REF_DIR}/morel_nodes.txt"

echo "✓ Reference files prepared in ${REF_DIR}"

# ============================================================================
# CREATE THE PER-SUBJECT PROCESSING SCRIPT
# ============================================================================

echo "Creating subject processing script..."

cat > "${SCRIPTS_DIR}/process_subject_native.sh" << 'SCRIPT_EOF'
#!/bin/bash
set -euxo pipefail

# ============================================================================
# PER-SUBJECT NATIVE SPACE PROCESSING SCRIPT
# This script is called in parallel for each subject
# ============================================================================

# Parse command-line arguments
SUBJECT=$1
DATA_DIR=$2
PREPROCESSED_DIR=$3
SUBJECTS_DIR=$4
CONNECTOMES_DIR=$5
NTHREADS=$6
STREAMLINES=$7
MNI_TEMPLATE=$8
ATLAS_FILE=$9
SEED_FILE=${10}
LOGS_DIR=${11}
NODE_LABELS=${12}

# ============================================================================
# DEFINE SUBJECT-SPECIFIC PATHS
# ============================================================================

# Input directories
DWI_PREPROC_DIR="${PREPROCESSED_DIR}/${SUBJECT}/dwi"
T1_DIR="${DATA_DIR}/${SUBJECT}/anat/T1"

# Output directories for this subject
SUB_DIR="${SUBJECTS_DIR}/${SUBJECT}"
REG_DIR="${SUB_DIR}/registration"
TRACT_DIR="${SUB_DIR}/tractography"
CONN_DIR="${SUB_DIR}/connectome"

# Create subject directory structure
mkdir -p "${REG_DIR}" "${TRACT_DIR}" "${CONN_DIR}"

# Setup logging for this subject
LOG_FILE="${LOGS_DIR}/${SUBJECT}_processing.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================================================"
echo "[${SUBJECT}] Starting NATIVE SPACE processing at $(date)"
echo "[${SUBJECT}] Streamlines target: ${STREAMLINES}"
echo "[${SUBJECT}] Threads allocated: ${NTHREADS}"
echo "============================================================================"

# Safety check: Skip if already completed 
if [[ -f "${CONNECTOMES_DIR}/${SUBJECT}_connectivity_matrix.csv" ]] && \
   [[ -f "${CONNECTOMES_DIR}/${SUBJECT}_connectivity_matrix_labeled.csv" ]]; then
    echo "[${SUBJECT}] ✓ Already completed! Skipping to prevent duplicate processing."
    exit 0
fi

# ============================================================================
# STEP 1: PREPARE NATIVE SPACE IMAGES
# ============================================================================

echo "[${SUBJECT}] Step 1/7: Preparing native space images..."
echo "[${SUBJECT}] - Extracting mean b0 image from DWI data"
echo "[${SUBJECT}] - Brain extraction for b0 and T1 images"
echo "[${SUBJECT}] - Generating 5-tissue-type (5TT) image for ACT"
cd "${REG_DIR}"

# Extract mean b0 image from preprocessed DWI data
DWI_FINAL_PREPROC="${DWI_PREPROC_DIR}/dwi_unbiased.mif"
dwiextract "${DWI_FINAL_PREPROC}" - -bzero | mrmath - mean b0_mean.nii.gz -axis 3 -force

# Brain extraction on b0 image (conservative settings for DWI)
bet b0_mean.nii.gz b0_mean_brain.nii.gz -f 0.2 -g 0 -m

# Brain extraction on T1 image (slightly more aggressive for better 5TT)
bet "${T1_DIR}/T1.nii.gz" T1_brain.nii.gz -f 0.3 -g 0 -B

# Generate 5TT image from T1 for anatomically-constrained tractography
# This creates tissue maps for: GM, subcortical GM, WM, CSF, and pathological tissue
5ttgen fsl T1_brain.nii.gz 5tt_native.mif -premasked -force -nthreads "${NTHREADS}"

echo "[${SUBJECT}] ✓ Native space images prepared"

# ============================================================================
# STEP 2: MULTI-TISSUE INTENSITY NORMALIZATION
# ============================================================================

echo "[${SUBJECT}] Step 2/7: Performing multi-tissue intensity normalization..."
echo "[${SUBJECT}] - Normalizing FOD amplitudes across tissues (WM, GM, CSF)"
echo "[${SUBJECT}] - This ensures comparability across subjects"
cd "${DWI_PREPROC_DIR}"

# Check if normalized FODs already exist to avoid redundant processing
if [[ ! -f "wmfod_norm.mif" ]] || [[ ! -f "gmfod_norm.mif" ]] || [[ ! -f "csffod_norm.mif" ]]; then
    # Perform multi-tissue normalization
    # This balances the FOD amplitudes based on the tissue properties
    mtnormalise wmfod.mif wmfod_norm.mif \
                gmfod.mif gmfod_norm.mif \
                csffod.mif csffod_norm.mif \
                -mask "${DWI_PREPROC_DIR}/mask.mif" \
                -force -nthreads "${NTHREADS}"
    echo "[${SUBJECT}] ✓ FOD normalization complete"
else
    echo "[${SUBJECT}] ✓ Normalized FODs already exist, skipping normalization"
fi

# ============================================================================
# STEP 3: CALCULATE REGISTRATION TRANSFORMATIONS
# ============================================================================

echo "[${SUBJECT}] Step 3/7: Computing registration transformations..."
echo "[${SUBJECT}] - T1 to DWI registration (rigid transformation)"
echo "[${SUBJECT}] - MNI to T1 registration (non-linear SyN)"
cd "${REG_DIR}"

# Register T1 to DWI space (rigid registration for within-subject alignment)
# Output: T1_to_DWI_0GenericAffine.mat
antsRegistrationSyNQuick.sh -d 3 \
    -f b0_mean_brain.nii.gz \
    -m T1_brain.nii.gz \
    -o T1_to_DWI_ \
    -t r

# Register MNI template to native T1 space (non-linear for atlas mapping)
# Outputs: MNI_to_T1_0GenericAffine.mat, MNI_to_T1_1Warp.nii.gz
antsRegistrationSyN.sh -d 3 \
    -f T1_brain.nii.gz \
    -m "${MNI_TEMPLATE}" \
    -o MNI_to_T1_ \
    -t s \
    -n "${NTHREADS}"

echo "[${SUBJECT}] ✓ Registration transformations computed"


# ============================================================================
# STEP 4: WARP ANATOMICAL STRUCTURES TO DWI SPACE AND COMBINE SEED WITH ATLAS
# ============================================================================

echo "[${SUBJECT}] Step 4/7: Warping anatomical structures to native DWI space..."

# --- 4a. Transform 5TT image (handle each tissue type separately) ---
echo "[${SUBJECT}] - Transforming 5TT tissue maps to DWI space"

# Split 5TT into individual tissue volumes for transformation
mrconvert 5tt_native.mif -coord 3 0 - | mrconvert - tissue0.nii.gz -force  # GM
mrconvert 5tt_native.mif -coord 3 1 - | mrconvert - tissue1.nii.gz -force  # Subcortical GM
mrconvert 5tt_native.mif -coord 3 2 - | mrconvert - tissue2.nii.gz -force  # WM
mrconvert 5tt_native.mif -coord 3 3 - | mrconvert - tissue3.nii.gz -force  # CSF
mrconvert 5tt_native.mif -coord 3 4 - | mrconvert - tissue4.nii.gz -force  # Pathological

# Transform each tissue volume to DWI space
for i in 0 1 2 3 4; do
    antsApplyTransforms -d 3 \
        -i tissue${i}.nii.gz \
        -r b0_mean_brain.nii.gz \
        -o tissue${i}_coreg.nii.gz \
        -t T1_to_DWI_0GenericAffine.mat \
        -n Linear
done

# Recombine transformed tissue maps into single 5TT image
mrcat tissue0_coreg.nii.gz tissue1_coreg.nii.gz tissue2_coreg.nii.gz \
      tissue3_coreg.nii.gz tissue4_coreg.nii.gz 5tt_coreg.mif -axis 3 -force

# Clean up temporary tissue files
rm -f tissue*.nii.gz

# --- 4b. Transform MNI atlas to DWI space (chain transformations) ---
echo "[${SUBJECT}] - Warping MNI atlas to DWI space"
# Chain: MNI -> T1 -> DWI using both transformations
antsApplyTransforms -d 3 \
    -i "${ATLAS_FILE}" \
    -r b0_mean.nii.gz \
    -o atlas_coreg.nii.gz \
    -t T1_to_DWI_0GenericAffine.mat \
    -t MNI_to_T1_1Warp.nii.gz \
    -t MNI_to_T1_0GenericAffine.mat \
    -n NearestNeighbor  # Preserve atlas labels

# --- 4c. Transform seed mask to DWI space ---
echo "[${SUBJECT}] - Warping seed mask to DWI space"
antsApplyTransforms -d 3 \
    -i "${SEED_FILE}" \
    -r b0_mean.nii.gz \
    -o seed_coreg.nii.gz \
    -t T1_to_DWI_0GenericAffine.mat \
    -t MNI_to_T1_1Warp.nii.gz \
    -t MNI_to_T1_0GenericAffine.mat \
    -n NearestNeighbor  # Preserve binary mask

# --- 4d. Combine seed with atlas as region 77 ---
echo "[${SUBJECT}] - Combining seed region as node 77 in atlas"
# First, zero out any atlas values where seed mask is present (to avoid overlap)
# Then add seed as region 77
# This ensures the seed region is tracked as a node in connectivity analysis
mrcalc atlas_coreg.nii.gz \
       seed_coreg.nii.gz 0 -gt \
       -mult \
       seed_coreg.nii.gz 77 -mult \
       -add \
       atlas_with_seed.nii.gz \
       -force

# Verify the combined atlas
echo "[${SUBJECT}] - Verifying combined atlas has 77 regions:"
mrstats atlas_with_seed.nii.gz -mask atlas_with_seed.nii.gz -output max

echo "[${SUBJECT}] ✓ All anatomical structures warped to DWI space and seed added as region 77"

# ============================================================================
# STEP 5: ANATOMICALLY-CONSTRAINED TRACTOGRAPHY (ACT)
# ============================================================================

echo "[${SUBJECT}] Step 5/7: Performing anatomically-constrained tractography..."
echo "[${SUBJECT}] - Generating ${STREAMLINES} streamlines"
echo "[${SUBJECT}] - Using seed mask and ACT constraints"
echo "[${SUBJECT}] - Using normalized WM FOD for improved accuracy"
cd "${TRACT_DIR}"

# Generate tractography with anatomical constraints using normalized FOD
# Key parameters:
#   -act: Use 5TT image for anatomical constraints
#   -seed_image: Start streamlines from seed region
#   -mask: Restrict to brain mask
#   -select: Target number of valid streamlines
#   -backtrack: Allow backtracking for difficult regions
#   -cutoff: FOD amplitude cutoff threshold
tckgen "${DWI_PREPROC_DIR}/wmfod_norm.mif" tracks_3M.tck \
    -act "${REG_DIR}/5tt_coreg.mif" \
    -seed_image "${REG_DIR}/seed_coreg.nii.gz" \
    -mask "${REG_DIR}/b0_mean_brain_mask.nii.gz" \
    -select ${STREAMLINES} \
    -backtrack \
    -cutoff 0.06 \
    -force -nthreads "${NTHREADS}"

echo "[${SUBJECT}] ✓ Tractography complete: ${STREAMLINES} streamlines generated"

# ============================================================================
# STEP 6: SIFT2 FILTERING FOR BIOLOGICAL ACCURACY
# ============================================================================

echo "[${SUBJECT}] Step 6/7: Applying SIFT2 filtering for biological accuracy..."
echo "[${SUBJECT}] - Computing streamline weights based on normalized FOD"
echo "[${SUBJECT}] - This preserves all streamlines but assigns weights"

# SIFT2: Spherical-deconvolution Informed Filtering of Tractograms
# Assigns weights to streamlines to match fiber density from normalized FOD
# Output is a text file with weights, not a new tractogram
tcksift2 tracks_3M.tck "${DWI_PREPROC_DIR}/wmfod_norm.mif" sift2_weights.txt \
    -act "${REG_DIR}/5tt_coreg.mif" \
    -force -nthreads "${NTHREADS}"

echo "[${SUBJECT}] ✓ SIFT2 weights computed"

# ============================================================================
# STEP 7: GENERATE LABELED CONNECTIVITY MATRIX WITH SEED INCLUDED
# ============================================================================

echo "[${SUBJECT}] Step 7/7: Building labeled structural connectivity matrix..."
cd "${CONN_DIR}"

# IMPORTANT: Now using atlas_with_seed.nii.gz which includes seed as region 77
# Generate weighted connectivity matrix using SIFT2 weights with node labels
# This will create a 77x77 matrix where row/column 77 is the seed region
tck2connectome "${TRACT_DIR}/tracks_3M.tck" \
    "${REG_DIR}/atlas_with_seed.nii.gz" \
    connectivity_matrix.csv \
    -tck_weights_in "${TRACT_DIR}/sift2_weights.txt" \
    -symmetric \
    -zero_diagonal \
    -scale_invnodevol \
    -out_assignment assignments.csv \
    -force -nthreads "${NTHREADS}"

# For labeled version, we need to update the node names to include seed
# First create a temporary label file with seed added
echo "77 Seed_Pain" >> "${NODE_LABELS}_with_seed.txt"

# Generate labeled version with node names for easier interpretation
tck2connectome "${TRACT_DIR}/tracks_3M.tck" \
    "${REG_DIR}/atlas_with_seed.nii.gz" \
    connectivity_matrix_labeled.csv \
    -tck_weights_in "${TRACT_DIR}/sift2_weights.txt" \
    -symmetric \
    -zero_diagonal \
    -scale_invnodevol \
    -out_assignment assignments_labeled.csv \
    -config WriteCsvLabels true \
    -force -nthreads "${NTHREADS}"

echo "[${SUBJECT}] ✓ Connectivity matrices generated (77x77 including seed as region 77)"

# ============================================================================
# COPY FILES TO CENTRAL REPOSITORY
# ============================================================================

echo "[${SUBJECT}] Copying connectivity matrices to central repository..."

# Copy with descriptive filenames to central connectomes folder
cp connectivity_matrix.csv "${CONNECTOMES_DIR}/${SUBJECT}_connectivity_matrix.csv"
cp connectivity_matrix_labeled.csv "${CONNECTOMES_DIR}/${SUBJECT}_connectivity_matrix_labeled.csv"
cp assignments.csv "${CONNECTOMES_DIR}/${SUBJECT}_assignments.csv"

echo "[${SUBJECT}] ✓ Files copied to central repository"

# ============================================================================
# PROCESSING COMPLETE
# ============================================================================

echo "============================================================================"
echo "[${SUBJECT}] ✓ NATIVE SPACE PROCESSING COMPLETE at $(date)"
echo "[${SUBJECT}] Output locations:"
echo "[${SUBJECT}]   - Registration: ${REG_DIR}"
echo "[${SUBJECT}]   - Tractography: ${TRACT_DIR}"
echo "[${SUBJECT}]   - Connectivity: ${CONN_DIR}"
echo "[${SUBJECT}]   - Central copy: ${CONNECTOMES_DIR}/${SUBJECT}_connectivity_matrix.csv"
echo "============================================================================"
SCRIPT_EOF

# Make the processing script executable
chmod +x "${SCRIPTS_DIR}/process_subject_native.sh"

echo "✓ Subject processing script created"

# ============================================================================
# LAUNCH PARALLEL PROCESSING FOR ALL SUBJECTS
# ============================================================================

echo "============================================================================"
echo "LAUNCHING PARALLEL PROCESSING"
echo "Processing ${#SUBJECTS[@]} subjects with ${PARALLEL_JOBS} parallel jobs"
echo "Using ${NTHREADS} threads per subject (21/24 cores utilized)"
echo "============================================================================"

# Create main processing log
MAIN_LOG="${LOGS_DIR}/main_processing_$(date +%Y%m%d_%H%M%S).log"

# Run subjects in parallel using GNU parallel
# This will process PARALLEL_JOBS subjects simultaneously
printf "%s\n" "${REMAINING_SUBJECTS[@]}" | \
    parallel -j ${PARALLEL_JOBS} --eta \
    --joblog "${LOGS_DIR}/parallel_run.log" \
    "${SCRIPTS_DIR}/process_subject_native.sh" {} \
        "${DATA_DIR}" \
        "${PREPROCESSED_DIR}" \
        "${SUBJECTS_DIR}" \
        "${CONNECTOMES_DIR}" \
        "${NTHREADS}" \
        "${STREAMLINES}" \
        "${MNI_TEMPLATE_TARGET}" \
        "${ATLAS_FILE}" \
        "${SEED_FILE}" \
        "${LOGS_DIR}" \
        "${REF_DIR}/morel_nodes.txt" \
    2>&1 | tee "${MAIN_LOG}"

# ============================================================================
# POST-PROCESSING: VERIFY COMPLETENESS
# ============================================================================

echo "============================================================================"
echo "VERIFYING PROCESSING COMPLETION"
echo "============================================================================"

# Count successful completions
completed_count=0
failed_subjects=()

for subject in "${SUBJECTS[@]}"; do
    if [[ -f "${CONNECTOMES_DIR}/${subject}_connectivity_matrix.csv" ]]; then
        ((completed_count++))
    else
        failed_subjects+=("${subject}")
    fi
done

echo "✓ Successfully processed: ${completed_count}/${#SUBJECTS[@]} subjects"

if [[ ${#failed_subjects[@]} -gt 0 ]]; then
    echo "⚠ Failed subjects: ${failed_subjects[*]}"
    echo "Check individual logs in: ${LOGS_DIR}"
else
    echo "✓ All subjects processed successfully!"
fi

# ============================================================================
# GENERATE SUMMARY REPORT
# ============================================================================

echo "============================================================================"
echo "GENERATING PROCESSING SUMMARY"
echo "============================================================================"

SUMMARY_FILE="${LOGS_DIR}/processing_summary.txt"
{
    echo "TRACTOGRAPHY PROCESSING SUMMARY"
    echo "================================"
    echo "Date: $(date)"
    echo "Total subjects: ${#SUBJECTS[@]}"
    echo "Completed: ${completed_count}"
    echo "Failed: ${#failed_subjects[@]}"
    echo ""
    echo "Configuration:"
    echo "  - Streamlines per subject: ${STREAMLINES}"
    echo "  - Parallel jobs: ${PARALLEL_JOBS}"
    echo "  - Threads per subject: ${NTHREADS}"
    echo "  - FOD normalization: Applied (mtnormalise)"
    echo "  - Atlas: Morel Unified Atlas (76 regions)"
    echo "  - Connectivity weighting: SIFT2"
    echo ""
    if [[ ${#failed_subjects[@]} -gt 0 ]]; then
        echo "Failed subjects:"
        printf '  - %s\n' "${failed_subjects[@]}"
    fi
} > "${SUMMARY_FILE}"

cat "${SUMMARY_FILE}"

# ============================================================================
# FINAL SUMMARY
# ============================================================================

echo "============================================================================"
echo "✓ TRACTOGRAPHY PIPELINE COMPLETE!"
echo "============================================================================"
echo "Analysis directory: ${TRACTOGRAPHY_BASE}"
echo ""
echo "Key outputs:"
echo "  - Connectivity matrices: ${CONNECTOMES_DIR}/"
echo "    • Unlabeled: *_connectivity_matrix.csv"
echo "    • Labeled: *_connectivity_matrix_labeled.csv"
echo "  - Individual subject data: ${SUBJECTS_DIR}/"
echo "  - Processing logs: ${LOGS_DIR}/"
echo "  - Reference files: ${REF_DIR}/"
echo ""
echo "Processing details:"
echo "  - FOD normalization: Applied"
echo "  - Atlas labels: Included (76 Morel regions)"
echo "  - System utilization: 21/24 cores (87.5%)"
echo ""
echo "Next steps:"
echo "  1. Quality check connectivity matrices"
echo "  2. Run group-level statistical analyses"
echo "  3. Visualize connectivity patterns"
echo "  4. Generate TDI maps if needed (separate script)"
echo "============================================================================"
