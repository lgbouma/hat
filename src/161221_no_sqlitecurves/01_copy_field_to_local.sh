# This is a script to copy a HAT field to local computer for processing.
# Usage:
# > bash 01_copy_field_to_local.sh

photometry_dir="/H/BIGPROJ/hatuser/2007_hatnet_phot/"
field_dir="G154"
bls_subdir="/BLS"

copy_dir=$photometry_dir$field_dir$bls_subdir

if [ ! -d "../data/HATpipe/$field_dir$bls_subdir" ]; then
  mkdir ../data/HATpipe/$field_dir
  mkdir ../data/HATpipe/$field_dir$bls_subdir 
fi
scp -r lbouma@phn12:$copy_dir ../data/HATpipe/$field_dir/.
