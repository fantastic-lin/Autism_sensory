#!/bin/bash
operation_dir="/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git" ## Change to your operation pathway
subject_dir="FS_git" ## Change to your subject directory

export SUBJECTS_DIR="$operation_dir/$subject_dir"

## Assign values from the volume of AAL to the surface vertex of standard brain template
mri_vol2surf --mov $operation_dir/AAL116_1mm.nii --mni152reg --hemi rh --out_type mgh --o rh.aal3.mgh
mri_vol2surf --mov $operation_dir/AAL116_1mm.nii --mni152reg --hemi lh --out_type mgh --o lh.aal3.mgh

## AAL index: odd number ==> left hemisphere; even number ==> right hemisphere
olf_ids_l=(5 9 15 21 25 29 37 39 41) ## 5:Frontal_Sup_Orb_L 9:Frontal_Mid_Orb_L 15:Frontal_Inf_Orb_L 21:Olfactory_L 25:Frontal_Med_Orb_L 29:Insula_L 37:Hippocampus_L 39:ParaHippocampal_L 41:Amygdala_L

olf_ids_R=(6 10 16 22 26 30 38 40 42) ## 6:Frontal_Sup_Orb_R 10:Frontal_Mid_Orb_R 16:Frontal_Inf_Orb_R 22:Olfactory_R 26:Frontal_Med_Orb_R 30:Insula_R 38:Hippocampus_R 40:ParaHippocampal_R 42:Amygdala_R

## Note: the left brain index of AAL actually corresponds to the right brain of Freesurfer. This can be confirmed by two examples: 
## (1) using the index of the right brain of AAL to run mris_vol2label (--h rh) will result in an error;
## (2) using the index of the left brain of AAL to run mris_vol2label (--h rh) will not result in error, and using freeview can correctly sgement the corresponding right brain region.

## Traverse the left and right hemisphere
hemis=(lh rh)
for he in "${hemis[@]}";do
	if [ $he == lh ];then
		olf_index=("${olf_ids_R[@]}")
	else
		olf_index=("${olf_ids_l[@]}")
	fi
	echo "此次执行的是${he}, 其脑区index为${olf_index[@]}"
## Create lh/rh.num{}.label in freesurfer
	for id in "${olf_index[@]}"; do
		mri_vol2label --i ${he}.aal3.mgh \
			--id "$id" \
			--l "${he}.num${id}.label" \
			--surf fsaverage ${he}
	done

## Merge lh/rh.num{}.label into a lh/rh.aal_olfactory.annot file 
	mris_label2annot --s fsaverage --h ${he} \
		 --ctab $operation_dir/${he}_aal_olfactory_RGB_note_1.txt \
		 --a aal_olfactory \
		 --l $SUBJECTS_DIR/fsaverage/label/${he}.num${olf_index[0]}.label \
		 --l $SUBJECTS_DIR/fsaverage/label/${he}.num${olf_index[1]}.label \
		 --l $SUBJECTS_DIR/fsaverage/label/${he}.num${olf_index[2]}.label \
		 --l $SUBJECTS_DIR/fsaverage/label/${he}.num${olf_index[3]}.label \
		 --l $SUBJECTS_DIR/fsaverage/label/${he}.num${olf_index[4]}.label \
		 --l $SUBJECTS_DIR/fsaverage/label/${he}.num${olf_index[5]}.label \
                 --l $SUBJECTS_DIR/fsaverage/label/${he}.num${olf_index[6]}.label \
                 --l $SUBJECTS_DIR/fsaverage/label/${he}.num${olf_index[7]}.label \
                 --l $SUBJECTS_DIR/fsaverage/label/${he}.num${olf_index[8]}.label
done

