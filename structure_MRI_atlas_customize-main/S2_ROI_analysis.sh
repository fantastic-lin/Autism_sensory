#!/bin/bash
#PBS -N ROI_analysis
#PBS -l nodes=1:ppn=6
#PBS -j oe
#PBS -q batch

operation_dir="/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git" ## Change to your operation pathway
subject_dir="FS_git" ## Change to your subject directory
export SUBJECTS_DIR="${operation_dir}/${subject_dir}"

json_file="${operation_dir}/FS_successful_download.json"
subject_list=($(jq -r '.[]' "$json_file"))
echo "subject_list:${subject_list[@]}"

## Generate structual measurement file for each subject
declare -a error_subjects_file=() ## Define an empty array for storing error files
declare -a fail_extract_subj=() ## Define an empty array for storing fail genearted subjects

counter_success_extract=0
hemis=(lh rh)
for subj in "${subject_list[@]}";do
        he_number=0
        for he in "${hemis[@]}";do
		rm $SUBJECTS_DIR/${subj}/stats/${subj}_${he}.txt
                if [ -f "$SUBJECTS_DIR/${subj}/stats/${subj}_${he}.txt" ]; then
			((he_number++))
                else
			echo "************ first step: mri_surf2surf *******************"
			if mri_surf2surf --hemi ${he} --srcsubject fsaverage --sval-annot $SUBJECTS_DIR/fsaverage/label/${he}.aal_olfactory.annot --trgsubject ${subj} --trgsurfval $SUBJECTS_DIR/${subj}/label/${he}.aal_olf_${subj}.annot; then
				echo "***************** second step: mris_anatomical_stats ******************"
				if mris_anatomical_stats -a $SUBJECTS_DIR/${subj}/label/${he}.aal_olf_${subj}.annot -f $SUBJECTS_DIR/${subj}/stats/${subj}_${he}.txt -b ${subj} ${he}; then
					echo "Anatomical stats generation successful for ${subj} ${he}"
					((he_number++))
				else
					echo "!!!!Error in mris_anatomical_stats for ${subj} ${he}"
					error_subjects_file+=("${subj}_${he}.txt")
				fi
			fi
		fi
	done
	if [ "$he_number" -eq 2 ];then
		((counter_success_extract++))
		echo "${subj} generated successfully"
	else
		fail_extract_subj+=("${subj}")
	fi
done

echo "error files: ${error_subjects_file[@]}"

echo "fail genearted subjects: ${fail_extract_subj[@]}"

echo "successful subjects: ${counter_success_extract}
