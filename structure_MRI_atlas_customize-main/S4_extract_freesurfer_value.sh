#!/bin/bash
operation_dir="/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git" ## Change to your operation pathway
subject_dir="FS_git" ## Change to your subject directory

export SUBJECTS_DIR="$operation_dir/$subject_dir"

## Surface Area
output_file="$operation_dir/FS_git_surface_area.csv"

for subj in $(find $SUBJECTS_DIR/* -maxdepth 0 -type d -not -name fsaverage -printf "%f\n"); do
	echo $subj
	lh_sum=$(sed -n '54,87p' $SUBJECTS_DIR/${subj}/stats/lh.aparc.stats | awk -F' ' '{sum+=$3} END {print sum}')
	echo "Sum of lh_surfer_area: $lh_sum"
	rh_sum=$(sed -n '54,87p' $SUBJECTS_DIR/${subj}/stats/rh.aparc.stats | awk -F' ' '{sum+=$3} END {print sum}')
	echo "Sum of rh_surfer_area: $rh_sum"
	total_sum=$(($lh_sum+$rh_sum))
	echo "Total area: $total_sum"
	echo "$subj,$total_sum" >> "$output_file"
done


## Volume
volume_output_file="$operation_dir/FS_git_gray_matter_volume.csv"
for subj in $(find $SUBJECTS_DIR/* -maxdepth 0 -type d -not -name fsaverage -printf "%f\n"); do
        echo $subj
        TotalGrayVol=$(grep 'Measure TotalGray, TotalGrayVol, Total gray matter volume' $SUBJECTS_DIR/${subj}/stats/aseg.stats | awk -F',' '{print $4}')
        echo "TotalGrayVol: $TotalGrayVol"
	Left_Hippocampus_vol=$(grep 'Left-Hippocampus' $SUBJECTS_DIR/${subj}/stats/aseg.stats | awk -F' ' '{print $4}')
	Right_Hippocampus_vol=$(grep 'Right-Hippocampus' $SUBJECTS_DIR/${subj}/stats/aseg.stats | awk -F' ' '{print $4}')
	Left_Amygdala_vol=$(grep 'Left-Amygdala' $SUBJECTS_DIR/${subj}/stats/aseg.stats | awk -F' ' '{print $4}')
        Right_Amygdala_vol=$(grep 'Right-Amygdala' $SUBJECTS_DIR/${subj}/stats/aseg.stats | awk -F' ' '{print $4}')
	echo "$subj,$TotalGrayVol,$Left_Hippocampus_vol,$Right_Hippocampus_vol,$Left_Amygdala_vol,$Right_Amygdala_vol" >> "$volume_output_file"
done
