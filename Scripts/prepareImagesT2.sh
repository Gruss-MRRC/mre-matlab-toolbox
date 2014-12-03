#!/usr/bin/env sh
# Prepare images for MR elastography reconstruction
#
# Dependencies:
# convertAndRescale script
# getTagValueFromNIFTIHeader.script
# FSL
# Matlab
# ImageJ + Miji plugin
#        + Phase_Tools plugin
#
# Author: Alex Krolick <amk283@cornell.edu>

# -----------------------------------------------------------------------------
# USAGE
# -----------------------------------------------------------------------------

if [ $# -le 0 ] # If run without arguments, print some help
then
  echo "Prepare images for MR elastography reconstruction"
  echo "Usage:"
  echo "prepareImagesT2 <MRE> <outputDir>"
  echo "  <MRE> is a 3D-sensitized MRE file in DICOM format"
  echo "  <outputDir> output directory (optional, defaults to NIFTI)"
  exit 2
fi

# -----------------------------------------------------------------------------
# INPUTS
# -----------------------------------------------------------------------------

mreFile=$1

if [ $# == 1 ]; then
  outputDir=NIFTI
else 
  outputdir=$2
fi

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------

scripts="/gmrrc/ALLSHARE/EAS/SOFTWARE/SCRIPTS"

# -----------------------------------------------------------------------------
# STEP 1: Convert DICOMs to NIFTI
# -----------------------------------------------------------------------------

echo "Converting DICOM to NIFTI..."

# Create output dir if it does not exist
[ -e ${outputDir} ] || mkdir ${outputDir}

mreName=$(basename "$mreFile") # Strip path
mreName="${mreName%.*}" # Strip file extension

convertAndRescale $outputDir/$mreName $mreFile

# -----------------------------------------------------------------------------
# STEP 2: Extract phase images
# -----------------------------------------------------------------------------

echo "Extracting phase images..."

cd $outputDir

# Split into magnitude and phase images
# First half of the volumes are mag; second are phase
# Get number of volumes, NVOL

NVOL=`$scripts/getTagValueFromNIFTIHeader.script $mreName dim4`

fslroi ${mreName}.nii ${mreName}_PHASE.nii $((NVOL/2)) $((NVOL/2))
fslroi ${mreName}.nii ${mreName}_MAG.nii   0           $((NVOL/2))

# Correct for bulk motion; register subsequent volumes to first volume
# motioncorrect=$scripts/dtiEddyCorrection.script
# $motioncorrect ${mreName}_PHASE.nii ${mreName}_PHASE_.nii
# $motioncorrect ${mreName}_MAG.nii ${mreName}_MAG_.nii

# -----------------------------------------------------------------------------
# STEP 3: White/gray matter separation
# -----------------------------------------------------------------------------

echo "Separating white & gray matter..."

mreBrain=${mreName}_MAG_BRAIN.nii.gz

# Run brain extraction
bet ${mreName}_MAG.nii.gz $mreBrain -R -f 0.25

# Partial volume extraction on skull-stripped magnitude
# This is used to compare white and gray matter stiffness
# -t flag sets mode to T2 or T1
fast -t 2 -g $mreBrain
      
# -----------------------------------------------------------------------------
# STEP 4: Phase unwrapping
# -----------------------------------------------------------------------------

# Unzip everything so that FIJI can open it
gunzip *.nii.gz

# Load Matlab and run an unwrapper in ImageJ
matlab=/apps1/matlab2013a/bin/matlab
i1="'${mreName}_MAG_BRAIN.nii'"
i2="'${mreName}_MAG.nii'"
i3="'${mreName}_PHASE.nii'"
i4="'${mreName}'"
$matlab -nodesktop -nosplash -r "unwrapAndNormalize($i1,$i2,$i3,$i4);exit"

gzip *.nii

# -----------------------------------------------------------------------------
# STEP 5: Housekeeping
# -----------------------------------------------------------------------------

# for some inscrutable reason, the files created by Matlab 
# always end up one directory level above where Matlab was launched
cd .. 
[ -e UNWRAP ] || mkdir UNWRAP
[ -e RECON ] || mkdir RECON

mv ${mreName}_unwrap* UNWRAP/

# -----------------------------------------------------------------------------
echo "Done"




