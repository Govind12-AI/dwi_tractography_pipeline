#!/bin/bash

# ==========================================================================================
#  GPU-Optimized and Parallelized dMRI Preprocessing and Modeling Pipeline
# ==========================================================================================
#
# DESCRIPTION:
# This script provides a comprehensive, start-to-finish pipeline for processing
# diffusion-weighted MRI (dMRI) data. It takes raw NIfTI files and produces
# high-quality Fibre Orientation Distributions (FODs) using the Multi-Shell,
# Multi-Tissue Constrained Spherical Deconvolution (MSMT-C  SD) algorithm.
#
# FEATURES:
#   - Robust Preprocessing: Includes denoising, Gibbs ringing removal, and FSL's
#     'eddy' correction for head motion and eddy current distortions.
#   - Advanced Modeling: Implements the recommended MSMT-CSD workflow, which
#     accurately models crossing fibers even in the presence of gray matter and CSF.
#   - GPU Acceleration: Automatically detects and uses 'eddy_cuda' if available,
#     significantly speeding up the most time-consuming preprocessing step.
#   - Parallel Processing: Utilizes GNU Parallel to process multiple subjects
#     simultaneously, maximizing hardware efficiency.
#   - Customizable: Designed to be easily adapted for new projects and datasets
#     by modifying a few key variables.
#
# DEPENDENCIES:
# This script requires the following software packages to be installed and accessible
# in the system's PATH or correctly sourced in the environment setup below.
#   1. MRtrix3 (https://www.mrtrix.org/): Core dMRI processing tools.
#   2. FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki): Essential for 'eddy' and other utilities.
#   3. ANTs (https://antsx.github.io/): Required for the 'dwibiascorrect ants' step.
#   4. FreeSurfer (https://surfer.nmr.mgh.harvard.edu/): Sourced for compatibility, though not
#      actively used in this specific script.
#   5. Conda / Miniconda (https://docs.conda.io/): For managing the Python environment for MRtrix3.
#   6. GNU Parallel (https://www.gnu.org/software/parallel/): For parallel execution.
#
# ---
#
# CUSTOMIZATION FOR YOUR PROJECT:
# To use this script for your data, please review and edit all sections marked with:
#   ---> EDIT THIS: <---
#
# CUSTOMIZATION FOR SINGLE-SHELL DATA:
# This pipeline is optimized for multi-shell data (required for MSMT-CSD).
# If you only have single-shell data, you will need to switch to a single-shell
# CSD workflow. The main changes would be:
#   1. In `dwi2response`, use an algorithm like `dwi2response tournier`.
#   2. In `dwi2fod`, switch from `msmt_csd` to `csd`.
#   This would require modifying the `process_subject` function.
#
# ==========================================================================================

# --- 1. ENVIRONMENT SETUP ---
# This section ensures that all necessary software (Conda, FSL, FreeSurfer)
# is available in the script's environment.
echo "--- Loading Processing Environment ---"
# ---> EDIT THIS: Update the path if your Conda installation is in a different location.
CONDA_BASE_DIR=$(conda info --base)
source "${CONDA_BASE_DIR}/etc/profile.d/conda.sh"
# ---> EDIT THIS: Activate the specific Conda environment that has MRtrix3 installed.
conda activate base
# ---> EDIT THIS: Set the FSLDIR variable to your FSL installation directory.
export FSLDIR=~/HCP_DWI_data/software/fsl
source "${FSLDIR}/etc/fslconf/fsl.sh"
# ---> EDIT THIS: Set the FREESURFER_HOME variable to your FreeSurfer installation directory.
export FREESURFER_HOME=~/HCP_DWI_data/software/freesurfer
source "${FREESURFER_HOME}/SetUpFreeSurfer.sh"
echo "--- Environment Loaded Successfully. Starting Preprocessing... ---"

# --- 2. GLOBAL CONFIGURATION ---
# These variables define the core paths and settings for the entire project.
# Modifying these is the primary way to customize the script for your data.

# ---> EDIT THIS: Set this to the main directory for your project.
export PROJECT_DIR=~/HCP_DWI_data
# ---> EDIT THIS: Set this to the base directory containing your raw subject data
# The script expects a structure like: .../data/SUBJECT_ID/diff/raw/...
export RAW_DATA_BASE_DIR="$PROJECT_DIR/data"
# ---> EDIT THIS: Set this to the directory where all processed outputs will be saved.
export PREPROC_OUTPUT_DIR="$PROJECT_DIR/outputs"
# ---> EDIT THIS: Set the maximum number of subjects to process in parallel.
# For a system with a powerful GPU and ample RAM (e.g., >32GB), you can often
# increase this value. For a 24GB VRAM GPU, 4-6 jobs is a safe starting point.
# For a high-end system with 48GB+ VRAM and 128GB+ RAM, you could try 8 or more.
# Monitor your 'nvidia-smi' and 'htop' output to find the optimal number.
export MAX_PARALLEL_JOBS=4


# --- 3. SUBJECT PROCESSING FUNCTION ---
# This function contains the entire step-by-step pipeline for a single subject.
# It is designed to be executed in parallel for each subject.
process_subject() {
    # --- A. Initialization ---
    local subj_id=$1
    local nthreads=$(nproc) # Use all available CPU cores for multi-threaded steps.
    echo "[$(date)] --- Starting Subject: ${subj_id} ---"
    # Create a dedicated output directory for this subject's DWI processing.
    local subj_preproc_dir="${PREPROC_OUTPUT_DIR}/${subj_id}/dwi"
    mkdir -p "${subj_preproc_dir}"
    cd "${subj_preproc_dir}" || return 1
    # Define paths to the raw data files for this subject.
    local raw_dwi_nii="${RAW_DATA_BASE_DIR}/${subj_id}/diff/raw/mri/diff.nii.gz"
    local raw_bvec="${RAW_DATA_BASE_DIR}/${subj_id}/diff/raw/bvecs_fsl.txt"
    local raw_bval="${RAW_DATA_BASE_DIR}/${subj_id}/diff/raw/bvals.txt"

    # Input data check.
    if [ ! -f "$raw_dwi_nii" ]; then echo "ERROR: Raw NIfTI not found for ${subj_id}."; return 1; fi

    # --- B. Core Preprocessing and Modeling Pipeline ---
    # The '&& \' chain ensures that the script will stop processing this subject
    # if any single command fails, preventing errors from propagating.

    # Step 1: Convert raw data from NIfTI to MRtrix format (.mif), incorporating gradient information.
    mrconvert "$raw_dwi_nii" dwi_initial.mif -fslgrad "${raw_bvec}" "${raw_bval}" -force -nthreads "$nthreads" && \

    # Step 2: Denoise the DWI data. This removes thermal noise and improves model fits.
    dwidenoise dwi_initial.mif dwi_den.mif -force -nthreads "$nthreads" && \

    # Step 3: Remove Gibbs Ringing artifacts. These are common artifacts near tissue boundaries.
    mrdegibbs dwi_den.mif dwi_den_unr.mif -axes 0,1 -force -nthreads "$nthreads" && \

    # Step 4: Eddy current, head motion, and EPI distortion correction using FSL's 'eddy'.
    # This is the most computationally intensive step. It will automatically use 'eddy_cuda' if found.
    dwifslpreproc dwi_den_unr.mif dwi_preproc.mif -rpe_none -pe_dir AP -eddy_options " --slm=linear --data_is_shelled --very_verbose" -force && \

    # Step 5: B1 field inhomogeneity correction using ANTs. This corrects for spatial intensity variations.
    dwibiascorrect ants dwi_preproc.mif dwi_unbiased.mif -bias bias.mif -force && \

    # Step 6: Create a brain mask from the final preprocessed DWI data.
    dwi2mask dwi_unbiased.mif mask.mif -force -nthreads "$nthreads" && \

    # Step 7: Estimate tissue response functions for White Matter, Gray Matter, and CSF.
    # The 'dhollander' algorithm is specifically designed for multi-shell data.
    dwi2response dhollander dwi_unbiased.mif wm_response.txt gm_response.txt csf_response.txt -force -nthreads "$nthreads" && \

    # Step 8: Perform Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution (MSMT-CSD).
    # This is the core modeling step that estimates the orientation of fiber bundles in each voxel.
    dwi2fod msmt_csd dwi_unbiased.mif wm_response.txt wmfod.mif gm_response.txt gmfod.mif csf_response.txt csffod.mif -mask mask.mif -force -nthreads "$nthreads"
    
   
    # --- C. Final Status Check ---
    # Check the exit code of the last command to confirm success.
    if [ $? -eq 0 ]; then
        echo "[$(date)] --- Subject ${subj_id} Finished Successfully ---"
    else
        echo "[$(date)] --- Subject ${subj_id} FAILED at one of the steps. See log for details. ---"
        return 1
    fi
}
# Export the function and key variables so they are accessible to GNU Parallel.
export -f process_subject
export RAW_DATA_BASE_DIR PREPROC_OUTPUT_DIR

# --- 4. MAIN EXECUTION LOGIC ---
# This is the master controller that runs the pipeline.
mkdir -p "$PREPROC_OUTPUT_DIR"
echo "Starting parallel batch preprocessing with up to ${MAX_PARALLEL_JOBS} jobs."

# ---> EDIT THIS: Define the list of subject IDs to be processed.
# The `seq` command is useful for numerical IDs. For alphanumeric IDs, you could
# use: SUBJECT_LIST=("sub-01" "sub-02" "sub-05")
# The current script assumes your subject folders are named e.g., 'mgh_1001', 'mgh_1002'.
SUBJECT_LIST=$(seq 1001 1035)

# Use GNU Parallel to robustly schedule and execute the processing for each subject.
# -j: Sets the number of parallel jobs.
# --joblog: Creates a log file tracking the start/end time and success/failure of each job.
# The 'bash -c' part ensures that each parallel job starts with a clean shell and
# properly sources the environment before calling our function.
# Each subject's standard output and error are redirected to a dedicated log file.
parallel -j $MAX_PARALLEL_JOBS --joblog "${PREPROC_OUTPUT_DIR}/parallel_joblog.txt" \
    'bash -c "source ~/.bashrc; process_subject mgh_{}" > "${PREROC_OUTPUT_DIR}/mgh_{}_main.log" 2>&1' ::: $SUBJECT_LIST

echo "All subjects have been dispatched for processing. See joblog for status."
