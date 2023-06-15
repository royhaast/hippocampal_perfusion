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
imgType=$3
scale_min=$4
scale_max=$5

scriptpath=/opt/braincharter
echo $scriptpath

###############################################################################
#                               Parse argument
###############################################################################

export AFNI_NIFTI_TYPE_WARN=NO

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 0. Parse the arguments.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"

if [ "${1}" = "-h" ]; then
	printf "extract_vessels.sh filename ext type(TOF or SWI or OTHER) [-m --mask mask_filename] [-p --phase phase_filename] [-s --std denoised_STD_value] [-f --factor upsampleFactor] [-c --centerline]"
	exit
fi

if [ -z "$image" ]; then
	printf "Input image name is required"
	exit
fi

if [ -z "$ext" ]; then
	printf "Extension is required"
	exit
fi

if [ -z "$imgType" ]; then
	printf "Image type {TOF|SWI} is required"
	exit
else
	if [ "${imgType}" = "TOF" ]; then
		printf "Time-of-flight is selected. \n"
	elif [ "${imgType}" = "SWI" ]; then
		printf "SWI is selected. \n"
	elif [ "${imgType}" = "OTHER" ]; then
		printf "SWI is selected. \n"
	else
		printf "ImageType Error, must be {TOF|SWI}"
		exit
	fi
fi

############################### Optional #####################################

function checkArg {
	if [ -z "$1" ]; then
		printf "Must provide a value after key.\n"
		exit
	fi
	case $1 in
		-m|--mask)
			printf "Value must not be a key.\n"
			exit
			;;
		-p|--phase)
			printf "Value must not be a key.\n"
			exit
			;;
		-f|--factor)
			printf "Value must not be a key.\n"
			exit
			;;
		-c|--centerline)
			printf "Value must not be a key.\n"
			exit
			;;
		-d|--diameters)
			printf "Value must not be a key.\n"
			exit
			;;
		-s|--scales)
			printf "Value must not be a key.\n"
			exit
			;;
		*)
			;;
	esac
}

if(( $(echo "$# > 3" | bc -l) )); then

	# Pass the 3 required arguments
	shift; shift; shift;

	while [[ $# -gt 0 ]]
	do
		key="$1"
		case $key in
			-m|--mask)
				mask="$2"
				checkArg ${image}_mask
				shift # past argument
				;;
			-p|--phase)
				phase="$2"
				checkArg ${phase}
				shift # past argument
				;;
			-f|--factor)
				upsampleFactor="$2"
				checkArg ${upsampleFactor}
				shift # past argument
				;;
			-c|--centerline)
				getCenterline=true
				shift # past argument
				;;
			-d|--diameters)
				getDiameters=true
				shift # past argument
				;;
			-s|--scales)
				scale_min=$2
				scale_max=$3
				# past argument
				;;
			*)
				#echo "unknow args!"
				#exit
				;;
		esac
		shift # past argument or value
	done


else
	printf "No optional args. \n"
fi

if [ -z "$upsampleFactor" ]; then
	printf "No upsample factor has been passed as argument.\nData won't be affect except a simple resample to isotropic for diameter extraction.\n"
else
	if(( $(echo "${upsampleFactor} < 2" | bc -l) )); then
		printf "Error: Your upsample factor must be greater than 1.\n"
		exit
	fi
fi

echo INPUT FILE          = "${image}"
echo EXTENSION           = "${ext}"
echo IMAGE TYPE          = "${imgType}"
echo FILE MASK           = "${image}_mask"
echo FILE PHASE          = "${phase}"
echo STD                 = "${stdDenoised}"
echo UPSAMPLED FACTOR    = "${upsampleFactor}"
echo GENERATE DIAMETERS = "${getDiameters}"
echo GENERATE CENTERLINE = "${getCenterline}"

###############################################################################
#                               Execution
###############################################################################

pushd "$(dirname "${image}")" 
image="$(basename "${image}" .nii.gz)"  

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 0. Reslicing to make isotropic. \n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"

min() {
	printf "%s\n" "${@:2}" | sort "$1" | head -n1
}

xdim=`fsl5.0-fslval $image.$ext pixdim1`
ydim=`fsl5.0-fslval $image.$ext pixdim2`
zdim=`fsl5.0-fslval $image.$ext pixdim3`

if [ "$xdim" = "$ydim" ] && [ "$ydim" = "$zdim" ] ; then
	printf "Data already isotropic, no reslicing necessary\n"
else
	mv $image.$ext ${image}_orig.$ext
	minsize=`min -g $xdim $ydim $zdim`
	echo fsl5.0-flirt -in ${image}_orig.$ext -ref ${image}_orig.$ext -applyisoxfm $minsize -out $image.$ext
	fsl5.0-flirt -in ${image}_orig.$ext -ref ${image}_orig.$ext -applyisoxfm $minsize -out $image.$ext
	printf "Data resliced to $minsize isotropic voxel size\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 2. SkullStripping. \n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_mask.${ext} ]; then
	fsl5.0-fslmaths ${image}.${ext} -div ${image}.${ext} ${image}_mask.${ext}
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 4. Vessel Extracting with VED.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"

if [ "${imgType}" = "SWI" ]; then
	if [ ! -f ${image}_Contrasted.${ext} ]; then
		${scriptpath}/InvertContrast.py ${image}_std${stdDenoised}_denoised.${ext} ${image}_mask.${ext} ${image}_Contrasted.${ext}
	else
		printf "Inverted contrast file already exists for this subject.\n"
	fi
else
	cp -rf ${image}.${ext} ${image}_Contrasted.${ext}
	printf "No needs to invert the constrast for TOF.\n"
fi

# Checking for smallest pix dim.
smalldim=99999.0

xdim=`3dinfo -adi ${image}_Contrasted.${ext}`
ydim=`3dinfo -adj ${image}_Contrasted.${ext}`
zdim=`3dinfo -adk ${image}_Contrasted.${ext}`

smalldim=`min -g $xdim $ydim $zdim`

# Hack to change smallest diameter value (added to original)
# smalldim=0.3

echo "Smallest dim : ${smalldim} \n"

if [ ! -z "$scale_min" ] && [ ! -z "$scale_max" ] ; then
	echo "Scales inputted."
	small_scale=$scale_min
	large_scale=$scale_max
	large_scale_clarity=$scale_max
else
	echo "Scales not inputted.. will be based on smalldims = ${smalldim}"
	a=0.3
	b=`echo "scale=3; ${smalldim} * 1.05" | bc`
	#small_scale=$(( a > b ? a : b ))
	small_scale=$b
	#small_scale=`echo "scale=3; ${small_scale} * 0.8" | bc`
	large_scale=`echo "scale=3; ${small_scale} * 2.3" | bc`
	large_scale_clarity=`echo "scale=3; ${small_scale} * 1.8" | bc`
fi

echo "small scales = ${small_scale}, large scales = ${large_scale}"

if [ ! -f ${image}_Ved.${ext} ]; then
	if [ "${imgType}" = "TOF" ]; then
		3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_upsampled.${ext} -inset ${image}_Contrasted.${ext}
		#${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${small_scale} -O -M ${large_scale} -t 1 -n 20 -s 2 -w 90 -I --out_folder "./${image}_iterations"
		${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${small_scale} -O -M 6 -t 1 -n 20 -s 2 -w 90 -I --out_folder "./${image}_iterations"
		#${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${smalldim} -M 6 -t 18 -n 10 -s 5 -w 25 #--generate_scale -D 'scales'
	elif [ "${imgType}" = "SWI" ]; then
		3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_upsampled.${ext} -inset ${image}_Contrasted.${ext}
		${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${small_scale} -O -M ${large_scale} -t 1 -n 20 -s 2 -w 90 -I --out_folder "./${image}_iterations"
		#${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${smalldim} -M 6 -t 18 -n 10 -s 5 -w 25
	elif [ "${imgType}" = "OTHER" ]; then
		3dresample -overwrite -dxyz ${smalldim} ${smalldim} ${smalldim} -rmode Cu -prefix ${image}_upsampled.${ext} -inset ${image}_Contrasted.${ext}
		${scriptpath}/ComputeVED.py ${image}_upsampled.${ext} ${image}_Ved.${ext} -m ${small_scale} -O -M ${large_scale_clarity} -t 1 -n 15 -s 2 -w 90 -I --out_folder "./${image}_iterations"
	fi
	rm -rf ./${image}_iterations
else
	printf "Vessel enhancing diffusion file already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 5. VED post-processing \n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_newVed_unscaled.${ext} ]; then
	#post-processing of Ved

	if [ ! -d ${image}_Ved_scales_nofilters ] ; then	mkdir ${image}_Ved_scales_nofilters ; else rm -rf ${image}_Ved_scales_nofilters/*	; fi
	mv -f Scale_NOW* ${image}_Ved_scales_nofilters/
	3dTcat -prefix ${image}_Ved_scales_nofilters/${image}_Ved_scales.${ext}  ${image}_Ved_scales_nofilters/Scale_NOW*
	#rm -rf Scale_NOW*
	mv -f ${image}_Ved_scales_nofilters/${image}_Ved_scales.${ext} ./
	3dTstat -overwrite -max -prefix ${image}_Ved.${ext} ${image}_Ved_scales.${ext}

	if [ ! -d ${image}_Ved_scales ] ; then mkdir ${image}_Ved_scales ; else	rm -rf ${image}_Ved_scales/* ; fi
	mv -f Scale_rescaled_* ${image}_Ved_scales/
	3dTcat -prefix ${image}_Ved_scales/${image}_newVed_scales.${ext} ${image}_Ved_scales/Scale_rescaled_*
	#rm -rf Scale_NOW*
	mv -f ${image}_Ved_scales/${image}_newVed_scales.${ext} ./
	3dTstat -overwrite -max -prefix ${image}_newVed.${ext} ${image}_newVed_scales.${ext}

	if [ ! -d ${image}_Ved_rescales ] ; then mkdir ${image}_Ved_rescales ; else	rm -rf ${image}_Ved_rescales/* ; fi
	mv -f Scale_processed_* ${image}_Ved_rescales/
	3dTcat -prefix ${image}_Ved_rescales/${image}_newVed_notscaled_scales.${ext} ${image}_Ved_rescales/Scale_processed_*
	#rm -rf Scale_NOW*
	mv -f ${image}_Ved_rescales/${image}_newVed_notscaled_scales.${ext} ./
	3dTstat -overwrite -max -prefix ${image}_newVed_unscaled.${ext} ${image}_newVed_notscaled_scales.${ext}

	if [ "${imgType}" = "SWI" ]; then
		3dresample -overwrite -master ${image}_newVed_scales.${ext} -rmode Cu -prefix ${image}_mask_up.${ext} -inset ${image}_mask.${ext}
		3dmask_tool -overwrite -input ${image}_mask.${ext} -prefix ${image}_newmask.${ext} -dilate_input -2
	elif [ "${imgType}" = "ToF" ]; then
		3dresample -overwrite -master ${image}_newVed_scales.${ext} -rmode Cu -prefix ${image}_mask_up.${ext} -inset ${image}_mask.${ext}
		3dmask_tool -overwrite -input ${image}_mask.${ext} -prefix ${image}_newmask.${ext} -dilate_input -1
	else
		3dresample -master ${image}_newVed_scales.${ext} -rmode Cu -prefix ${image}_newmask.${ext} -inset ${image}_mask.${ext}
	fi

	3dcalc -overwrite -a ${image}_newVed.${ext} -b ${image}_newmask.${ext} -expr "step(b) * a " -prefix ${image}_newVed_sqrt.${ext}
	3dcalc -overwrite -a ${image}_newVed_unscaled.${ext} -b ${image}_newmask.${ext} -expr "step(b) * ( (1-astep(a, 1000))*sqrt(a) + astep(a, 10)*10 ) " -prefix ${image}_newVed_unscaled_sqrt.${ext}
	3dcalc -overwrite -a ${image}_Ved.${ext} -b ${image}_newmask.${ext} -expr "step(b) * ( (1-astep(a, 1000))*sqrt(a) + astep(a, 1000)*1000 ) " -prefix ${image}_Ved_sqrt.${ext}

	3dcalc -overwrite -a ${image}_Ved_scales.${ext} -b ${image}_newmask.${ext} -expr "step(b) * ( (1-astep(a, 1000))*sqrt(a) + astep(a, 1000)*1000 ) " -prefix ${image}_Ved_scales_sqrt.${ext}
	3dcalc -overwrite -a ${image}_newVed_notscaled_scales.${ext} -b ${image}_newmask.${ext} -expr "step(b) * ( (1-astep(a, 1000))*sqrt(a) + astep(a, 1000)*1000 ) " -prefix ${image}_newVed_notscaled_scales_sqrt.${ext}
	3dcalc -overwrite -a ${image}_newVed_scales.${ext} -b ${image}_newmask.${ext} -expr "step(b) * a " -prefix ${image}_newVed_scales_sqrt.${ext}

else
	printf "Processed VED files already exists for this subject.\n"
fi

printf "\n+-+- +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 6. Revise scales (post-processing) \n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_Ved_corrected.${ext} ] || [ -f ${image}_Ved_scales_sqrt.${ext} ]; then

	3dTstat -overwrite -max -prefix ${image}_Ved_corrected.${ext} ${image}_Ved_scales_sqrt.${ext}[0..10]
	3dTstat -overwrite -max -prefix ${image}_newVed_corrected.${ext} ${image}_newVed_scales_sqrt.${ext}[0..10]
	3dTstat -overwrite -max -prefix ${image}_newVed_unscaled_corrected.${ext} ${image}_newVed_notscaled_scales_sqrt.${ext}[0..10]

	3dcalc -overwrite -a ${image}_Ved_corrected.${ext} -b ${image}_newmask.${ext} -expr "step(b)*astep(a,1.0)" -prefix ${image}_Ved_Thr.${ext} -datum short
	3dcalc -overwrite -a ${image}_newVed_corrected.${ext} -b ${image}_newmask.${ext} -expr "step(b)*astep(a,40.0)" -prefix ${image}_newVed_Thr.${ext} -datum short
	3dcalc -overwrite -a ${image}_newVed_unscaled_corrected.${ext} -b ${image}_newmask.${ext} -expr "step(b)*astep(a,1.0)" -prefix ${image}_newVed_unscaled_Thr.${ext} -datum short

	3dcalc -overwrite -a ${image}_Ved_corrected.${ext} -b ${image}_newmask.${ext} -expr "step(b)*astep(a,1.0)*log(a)*step(log(a))" -prefix ${image}_Ved_log.${ext}
	3dcalc -overwrite -a ${image}_newVed_corrected.${ext} -b ${image}_newmask.${ext} -expr "step(b)*astep(a,40.0)*log(a)*step(log(a))" -prefix ${image}_newVed_log.${ext}
	3dcalc -overwrite -a ${image}_newVed_unscaled_corrected.${ext} -b ${image}_newmask.${ext} -expr "step(b)*astep(a,1.0)*log(a)*step(log(a))" -prefix ${image}_newVed_unscaled_log.${ext}

	3dcalc -overwrite -a ${image}_Ved_corrected.${ext} -b ${image}_newmask.${ext} -expr "step(b)*astep(a,1.0)*log(log(a))*step(log(a))" -prefix ${image}_Ved_loglog.${ext}
	3dcalc -overwrite -a ${image}_newVed_corrected.${ext} -b ${image}_newmask.${ext} -expr "step(b)*astep(a,40.0)*log(log(a))*step(log(a))" -prefix ${image}_newVed_loglog.${ext}
	3dcalc -overwrite -a ${image}_newVed_unscaled_corrected.${ext} -b ${image}_newmask.${ext} -expr "step(b)*astep(a,1.0)*log(log(a))*step(log(a))" -prefix ${image}_newVed_unscaled_loglog.${ext}

	3dmask_tool -overwrite -dilate_inputs 1 -1 -input ${image}_Ved_Thr.${ext} -prefix ${image}_Ved_Thr_opened.${ext}
	3dmask_tool -overwrite -dilate_inputs 1 -1 -input ${image}_newVed_Thr\.${ext} -prefix ${image}_newVed_Thr_opened.${ext}
	3dmask_tool -overwrite -dilate_inputs 1 -1 -input ${image}_newVed_unscaled_Thr.${ext} -prefix ${image}_newVed_unscaled_Thr_opened.${ext}

	3dmerge -overwrite -dxyz=1 -isovalue -1clust 1.01 60 -prefix ${image}_Ved_Thr_clean.${ext} ${image}_Ved_Thr_opened.${ext}
	3dmerge -overwrite -dxyz=1 -isovalue -1clust 1.01 60 -prefix ${image}_newVed_Thr_clean.${ext} ${image}_newVed_Thr_opened.${ext}
	3dmerge -overwrite -dxyz=1 -isovalue -1clust 1.01 60 -prefix ${image}_newVed_unscaled_Thr_clean.${ext} ${image}_newVed_unscaled_Thr_opened.${ext}

else
	printf "Post-processed files already for this subject.\n"
fi

printf "\n+-+- +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 7. Diameters extraction.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_diameters.${ext} ]; then
	echo "Diameter Extraction"
	halfsize=`echo "scale=3; ${smalldim} / 2.0" | bc`
	3dresample -overwrite -dxyz ${halfsize} ${halfsize} ${halfsize} -rmode Cu -prefix ${image}_newVed_unscaled_Thr_HALF.${ext} -inset ${image}_newVed_unscaled_Thr_clean.${ext}
	${scriptpath}/ExtractDiameter.py ${image}_newVed_unscaled_Thr_HALF.${ext} ${image}_diameters.${ext}
	3dresample -overwrite -master ${image}_newVed_unscaled_Thr_clean.${ext} -rmode Cu -prefix ${image}_diameters.${ext} -inset ${image}_diameters.${ext}
	3dcalc -overwrite -a ${image}_diameters.${ext} -b ${image}_newVed_unscaled_Thr_clean.${ext} -expr "step(a)*step(b)*a" -prefix ${image}_diameters.${ext}
	rm -rf *HALF*
else
	printf "Diameters file already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 8. Centerlines extraction.\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_skel.${ext} ]; then
	echo "Centerlines Extraction"
	${scriptpath}/ExtractCenterline.py ${image}_newVed_unscaled_Thr_clean.${ext} ${image}_skel.${ext}
	echo "Diameters in centerlines only"
	3dcalc -overwrite  -a ${image}_skel.${ext} -b ${image}_diameters.${ext}  -expr "step(a)*b" -prefix ${image}_centerdia.${ext} -datum float
else
	printf "Centerline file already exists for this subject.\n"
fi

printf "\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
printf "Step 9. Diameter extraction (cleaner).\n"
printf "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n"
if [ ! -f ${image}_centerdia_clean.${ext} ]; then
	echo "Centerlines clean"
	3dresample -overwrite -master ${image}_newVed_unscaled_Thr_clean.${ext} -rmode Cu -prefix ${image}_diameters.${ext} -inset ${image}_diameters.${ext}
	3dBlurInMask -prefix ${image}_diameters_clean.${ext} -mask ${image}_newVed_unscaled_Thr_clean.${ext} -FWHM 4 ${image}_diameters.${ext}
	3dcalc -overwrite  -a ${image}_skel.${ext} -b ${image}_diameters.${ext}  -expr "step(a)*b" -prefix ${image}_centerdia.${ext} -datum float
	3dcalc -overwrite  -a ${image}_skel.${ext} -b ${image}_diameters_clean.${ext}  -expr "step(a)*b" -prefix ${image}_centerdia_clean.${ext} -datum float
else
	printf "Centerline file already exists for this subject.\n"
fi

printf "Pipeline process completed.\n\n"
unset AFNI_NIFTI_TYPE_WARN
popd