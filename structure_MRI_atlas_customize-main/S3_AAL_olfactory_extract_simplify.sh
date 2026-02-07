#!bin/bash
operation_dir="/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git" ## Change to your subject directory
subject_dir="FS_git" ## Change to your subject directory

export SUBJECTS_DIR="${operation_dir}/${subject_dir}"

json_file="${operation_dir}/FS_successful_download.json"
subject_list=($(jq -r '.[]' "$json_file"))
echo "subject_list:${subject_list[@]}"

missing_files=FS_git_without_*h_subj.csv

hemis=(lh rh)
for subj in "${subject_list[@]}";do
        echo $subj
	for he in "${hemis[@]}"; do
		if [ -f "$SUBJECTS_DIR/${subj}/stats/${subj}_${he}_1.txt" ]; then
			echo "$SUBJECTS_DIR/${subj}/stats/${subj}_${he}_1.txt already exists"
			rm "$SUBJECTS_DIR/${subj}/stats/${subj}_${he}_1.txt"
		else   
			echo "$SUBJECTS_DIR/${subj}/stats/${subj}_${he}_1.txt does not exists"
			if [ -f "$SUBJECTS_DIR/${subj}/stats/${subj}_${he}.txt" ]; then
				sed -n '61,66p' $SUBJECTS_DIR/${subj}/stats/${subj}_${he}.txt | awk -F ' ' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > $SUBJECTS_DIR/${subj}/stats/${subj}_${he}_1.txt
			else
				echo "${subj}_${he}.txt" >> "$missing_files"
				echo "File ${subj}_${he}.txt does not exist. Skipping..."
			fi
		fi
	done
done
