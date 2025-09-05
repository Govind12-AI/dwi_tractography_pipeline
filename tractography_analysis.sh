#!/bin/bash

# ============================================================================
# NATIVE SPACE TRACTOGRAPHY & CONNECTOME PIPELINE (Final Corrected Version)
#
# FIXES: 
# 1. transformconvert syntax corrected to use proper itk_import format
# 2. File path corrected in tck2connectome command
# ============================================================================

# Stricter error handling for a more robust pipeline
set -euxo pipefail

# --- CONFIGURATION ---
BASE_DIR="/home/kaushik/HCP_DWI_data"
DATA_DIR="${BASE_DIR}/data"
PREPROCESSED_DIR="${BASE_DIR}/outputs"
TRACTOGRAPHY_BASE="${BASE_DIR}/tractography"
ATLAS_FILE="${BASE_DIR}/morel_unified_atlas/1mm/Morel_Unified_Atlas_1mm.nii.gz"
SEED_FILE="${BASE_DIR}/seed_mask_pain.nii"
TEST_SUBJECTS=("mgh_1001" "mgh_1002" "mgh_1003")
PARALLEL_JOBS=3
NTHREADS=20
# --- END CONFIGURATION ---

# --- SETUP ANALYSIS DIRECTORY AND REFERENCE FILES ---
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
ANALYSIS_NAME="native_space_analysis_final_${TIMESTAMP}"
ANALYSIS_DIR="${TRACTOGRAPHY_BASE}/${ANALYSIS_NAME}"

echo "============================================================================"
echo "NATIVE SPACE CONNECTIVITY ANALYSIS"
echo "Analysis directory: ${ANALYSIS_DIR}"
echo "============================================================================"

mkdir -p "${ANALYSIS_DIR}"/{00_scripts,01_logs,02_reference}
REF_DIR="${ANALYSIS_DIR}/02_reference"
FSL_MNI_PATH="${HOME}/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"
MNI_TEMPLATE_TARGET="${REF_DIR}/MNI152_T1_1mm_brain.nii.gz"

if [[ -f "${FSL_MNI_PATH}" ]]; then
    echo "Found 1mm MNI template. Copying to reference directory."
    cp "${FSL_MNI_PATH}" "${MNI_TEMPLATE_TARGET}"
else
    echo "FATAL: 1mm MNI template not found at ${FSL_MNI_PATH}. Please check the path."
    exit 1
fi
echo "? Reference files prepared."

# --- CREATE THE PER-SUBJECT PROCESSING SCRIPT ---
cat > "${ANALYSIS_DIR}/00_scripts/process_subject_native.sh" << 'SCRIPT_EOF'
#!/bin/bash
set -euxo pipefail

# --- ARGUMENT PARSING ---
SUBJECT=$1
DATA_DIR=$2
PREPROCESSED_DIR=$3
ANALYSIS_DIR=$4
NTHREADS=$5
MNI_TEMPLATE=$6
ATLAS_FILE=$7
SEED_FILE=$8
# --- END ARGUMENT PARSING ---

# --- PATH DEFINITIONS ---
DWI_PREPROC_DIR="${PREPROCESSED_DIR}/${SUBJECT}/dwi"
T1_DIR="${DATA_DIR}/${SUBJECT}/anat/T1"
SUB_DIR="${ANALYSIS_DIR}/${SUBJECT}"
REG_DIR="${SUB_DIR}/01_registration"
TRACT_DIR="${SUB_DIR}/02_tractography"
CONN_DIR="${SUB_DIR}/03_connectome"
LOG_FILE="${ANALYSIS_DIR}/01_logs/${SUBJECT}_processing.log"

mkdir -p "${REG_DIR}" "${TRACT_DIR}" "${CONN_DIR}"

exec > >(tee -a "${LOG_FILE}") 2>&1
echo "============================================================================"
echo "[${SUBJECT}] Starting NATIVE SPACE processing at $(date)"
echo "============================================================================"

# --- STEP 1: PREPARE NATIVE SPACE IMAGES ---
echo "[${SUBJECT}] Step 1: Preparing native images (b0, T1 brain, 5TT)..."
cd "${REG_DIR}"

DWI_FINAL_PREPROC="${DWI_PREPROC_DIR}/dwi_unbiased.mif"
dwiextract "${DWI_FINAL_PREPROC}" - -bzero | mrmath - mean b0_mean.nii.gz -axis 3 -force
bet b0_mean.nii.gz b0_mean_brain.nii.gz -f 0.2 -g 0 -m
bet "${T1_DIR}/T1.nii.gz" T1_brain.nii.gz -f 0.3 -g 0 -B

5ttgen fsl T1_brain.nii.gz 5tt_native.mif -premasked -force -nthreads "${NTHREADS}"

# --- STEP 2: CALCULATE TRANSFORMATIONS WITH ANTS ---
echo "[${SUBJECT}] Step 2: Calculating transformations (T1->DWI and MNI->T1)..."
antsRegistrationSyNQuick.sh -d 3 -f b0_mean_brain.nii.gz -m T1_brain.nii.gz -o T1_to_DWI_ -t r
antsRegistrationSyN.sh -d 3 -f T1_brain.nii.gz -m "${MNI_TEMPLATE}" -o MNI_to_T1_ -t s -n "${NTHREADS}"

# --- STEP 3: WARP ANATOMICALS INTO NATIVE DWI SPACE ---
echo "[${SUBJECT}] Step 3: Warping 5TT, Atlas, and Seed into native DWI space..."

# a) Warp the 5TT image - need to handle each tissue volume separately
echo "[${SUBJECT}] ... handling 5TT image (transforming each tissue type)..."
# Split the 5TT into individual tissue volumes
mrconvert 5tt_native.mif -coord 3 0 - | mrconvert - tissue0.nii.gz -force
mrconvert 5tt_native.mif -coord 3 1 - | mrconvert - tissue1.nii.gz -force
mrconvert 5tt_native.mif -coord 3 2 - | mrconvert - tissue2.nii.gz -force
mrconvert 5tt_native.mif -coord 3 3 - | mrconvert - tissue3.nii.gz -force
mrconvert 5tt_native.mif -coord 3 4 - | mrconvert - tissue4.nii.gz -force

# Transform each tissue volume separately
for i in 0 1 2 3 4; do
    antsApplyTransforms -d 3 -i tissue${i}.nii.gz -r b0_mean_brain.nii.gz \
        -o tissue${i}_coreg.nii.gz \
        -t T1_to_DWI_0GenericAffine.mat \
        -n Linear
done

# Recombine into a single 5TT image
mrcat tissue0_coreg.nii.gz tissue1_coreg.nii.gz tissue2_coreg.nii.gz \
      tissue3_coreg.nii.gz tissue4_coreg.nii.gz 5tt_coreg.mif -axis 3 -force

# Clean up temporary files
rm -f tissue*.nii.gz

# b) Warp the MNI Atlas into DWI space by chaining the MNI->T1 and T1->DWI transforms
echo "[${SUBJECT}] ... warping MNI atlas to DWI space..."
antsApplyTransforms -d 3 -i "${ATLAS_FILE}" -r b0_mean.nii.gz \
    -o atlas_coreg.nii.gz \
    -t T1_to_DWI_0GenericAffine.mat \
    -t MNI_to_T1_1Warp.nii.gz \
    -t MNI_to_T1_0GenericAffine.mat \
    -n NearestNeighbor

# c) Warp the MNI Seed into DWI space using the same chained transform
echo "[${SUBJECT}] ... warping MNI seed to DWI space..."
antsApplyTransforms -d 3 -i "${SEED_FILE}" -r b0_mean.nii.gz \
    -o seed_coreg.nii.gz \
    -t T1_to_DWI_0GenericAffine.mat \
    -t MNI_to_T1_1Warp.nii.gz \
    -t MNI_to_T1_0GenericAffine.mat \
    -n NearestNeighbor

# --- STEP 4: NATIVE SPACE TRACTOGRAPHY ---
echo "[${SUBJECT}] Step 4: Performing Anatomically-Constrained Tractography..."
cd "${TRACT_DIR}"

tckgen "${DWI_PREPROC_DIR}/wmfod.mif" tracks_1M.tck \
    -act "${REG_DIR}/5tt_coreg.mif" \
    -seed_image "${REG_DIR}/seed_coreg.nii.gz" \
    -mask "${REG_DIR}/b0_mean_brain_mask.nii.gz" \
    -select 1000000 \
    -backtrack \
    -cutoff 0.06 \
    -force -nthreads "${NTHREADS}"

# --- STEP 5: FILTER TRACTOGRAM WITH SIFT2 ---
echo "[${SUBJECT}] Step 5: Filtering tractogram with SIFT2 for biological accuracy..."
# SIFT2 assigns weights to streamlines (doesn't reduce their number)
# Output is a weights file (.txt), not a new tractogram
tcksift2 tracks_1M.tck "${DWI_PREPROC_DIR}/wmfod.mif" sift2_weights.txt \
    -act "${REG_DIR}/5tt_coreg.mif" \
    -force -nthreads "${NTHREADS}"

# --- STEP 6: GENERATE THE CONNECTIVITY MATRIX ---
echo "[${SUBJECT}] Step 6: Building the final connectivity matrix..."
cd "${CONN_DIR}"

# Use the original tractogram with SIFT2 weights for connectivity
tck2connectome "${TRACT_DIR}/tracks_1M.tck" "${REG_DIR}/atlas_coreg.nii.gz" connectivity_matrix.csv \
    -tck_weights_in "${TRACT_DIR}/sift2_weights.txt" \
    -symmetric -zero_diagonal \
    -scale_invnodevol \
    -out_assignment assignments.csv \
    -force -nthreads "${NTHREADS}"

echo "[${SUBJECT}] ? NATIVE SPACE PROCESSING COMPLETE at $(date)"
SCRIPT_EOF

chmod +x "${ANALYSIS_DIR}/00_scripts/process_subject_native.sh"

# ============================================================================
# RUN PROCESSING IN PARALLEL
# ============================================================================
echo "--> Launching parallel processing for ${#TEST_SUBJECTS[@]} subjects..."
cd "${ANALYSIS_DIR}"

printf "%s\n" "${TEST_SUBJECTS[@]}" | \
    parallel -j ${PARALLEL_JOBS} --eta \
    --joblog "${ANALYSIS_DIR}/01_logs/parallel_run.log" \
    "${ANALYSIS_DIR}/00_scripts/process_subject_native.sh" {} \
        "${DATA_DIR}" \
        "${PREPROCESSED_DIR}" \
        "${ANALYSIS_DIR}" \
        "${NTHREADS}" \
        "${MNI_TEMPLATE_TARGET}" \
        "${ATLAS_FILE}" \
        "${SEED_FILE}"

echo "============================================================================"
echo "? PIPELINE COMPLETE!"
echo "Final connectivity matrices are located in: ${ANALYSIS_DIR}/{subject_id}/03_connectome/"
echo "Check logs in: ${ANALYSIS_DIR}/01_logs/"
echo "============================================================================"
