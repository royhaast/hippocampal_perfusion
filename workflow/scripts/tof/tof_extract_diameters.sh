#!/bin/bash

source activate neuro
export PATH=/usr/lib/afni/bin:$PATH
export PATH=/opt/braincharter:$PATH
export AFNI_NIFTI_TYPE_WARN=NO

###############################################################################
#                                   Input
###############################################################################
image=$1
ext=$2

scriptpath=/opt/braincharter
printf $scriptpath

###############################################################################
#                               Execution
###############################################################################

pushd "$(dirname "${image}")" 
image="$(basename "${image}" .nii.gz)"  

printf "\n+-+- +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 0. Binarize segmentations.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_bin.${ext} ]; then
    echo "Binarize segmentations"
    fsl5.0-fslmaths ${image}.${ext} -thr 255 -bin ${image}_bin.${ext}
fi
image=${image}_bin

printf "\n+-+- +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 1. Diameters extraction.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_diameters.${ext} ]; then
	echo "Diameter Extraction"
	#halfsize=`echo "scale=3; ${smalldim} / 2.0" | bc`
	halfsize=0.15
	3dresample -overwrite -dxyz ${halfsize} ${halfsize} ${halfsize} -rmode Cu -prefix ${image}_HALF.${ext} -inset ${image}.${ext}
	${scriptpath}/ExtractDiameter.py ${image}_HALF.${ext} ${image}_diameters.${ext}
	3dresample -overwrite -master ${image}.${ext} -rmode Cu -prefix ${image}_diameters.${ext} -inset ${image}_diameters.${ext}
	3dcalc -overwrite -a ${image}_diameters.${ext} -b ${image}.${ext} -expr "step(a)*step(b)*a" -prefix ${image}_diameters.${ext}
	rm -rf *HALF*
else
	printf "Diameters file already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 2. Centerlines extraction.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_skel.${ext} ]; then
	echo "Centerlines Extraction"
	${scriptpath}/ExtractCenterline.py ${image}.${ext} ${image}_skel.${ext}
	echo "Diameters in centerlines only"
	3dcalc -overwrite  -a ${image}_skel.${ext} -b ${image}_diameters.${ext}  -expr "step(a)*b" -prefix ${image}_centerdia.${ext} -datum float
else
	printf "Centerline file already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 3. Diameter extraction (cleaner).\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_centerdia_clean.${ext} ]; then
	echo "Centerlines clean"
	3dresample -overwrite -master ${image}.${ext} -rmode Cu -prefix ${image}_diameters.${ext} -inset ${image}_diameters.${ext}
	3dBlurInMask -prefix ${image}_diameters_clean.${ext} -mask ${image}.${ext} -FWHM 4 ${image}_diameters.${ext}
	3dcalc -overwrite  -a ${image}_skel.${ext} -b ${image}_diameters.${ext}  -expr "step(a)*b" -prefix ${image}_centerdia.${ext} -datum float
	3dcalc -overwrite  -a ${image}_skel.${ext} -b ${image}_diameters_clean.${ext}  -expr "step(a)*b" -prefix ${image}_centerdia_clean.${ext} -datum float
else
	printf "Centerline file already exists for this subject.\n"
fi

for file in `ls *_bin_*.nii.gz` ; do 
	fsl5.0-fslcpgeom ${image}.${ext} $file 
done

fsl5.0-fslmaths ${image}_diameters_clean.${ext} -dilM ${image}_diameters_clean_dil.${ext} 

printf "Pipeline process completed.\n\n"
unset AFNI_NIFTI_TYPE_WARN
popd