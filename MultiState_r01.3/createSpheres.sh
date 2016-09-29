#!/bin/bash
#
# input X Y Z radius

function convert {
    x=$1; y=$2; z=$3
}

MNIC="$1 $2 $3"
radius=$4
VOI=$5
# convert MNI to voxel coordinates
VoxC=$(echo $MNIC | std2imgcoord  -img $FSLDIR/data/standard/MNI152_T1_2mm_brain -std $FSLDIR/data/standard/MNI152_T1_2mm_brain -vox)
convert $VoxC
#create point mask
fslmaths $FSLDIR/data/standard/MNI152_T1_2mm -mul 0 -add 1 -roi $x 1 $y 1 $z 1 0 1 VOIs/VOIfiles/point -odt float
#create sphere from point
fslmaths VOIs/VOIfiles/point -kernel sphere $radius -fmean VOIs/$VOI -odt float
fslmaths VOIs/$VOI -bin VOIs/$VOI
