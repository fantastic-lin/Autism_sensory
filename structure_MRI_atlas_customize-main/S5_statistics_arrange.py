#!/usr/bin/python
#coding=utf-8
import pandas as pd
import os
import json
import csv
aspect_columns_surface_area = ["Frontal_Sup_Orb_L_5","Frontal_Sup_Orb_R_6","Frontal_Mid_Orb_L_9","Frontal_Mid_Orb_R_10","Frontal_Inf_Orb_L_15","Frontal_Inf_Orb_R_16","Olfactory_L_21","Olfactory_R_22","Frontal_Mid_Orb_L_25","Frontal_Mid_Orb_R_26","Insula_L_29","Insula_R_30"]
aspect_columns_gray_matter_volume = ["Frontal_Sup_Orb_L_5","Frontal_Sup_Orb_R_6","Frontal_Mid_Orb_L_9","Frontal_Mid_Orb_R_10","Frontal_Inf_Orb_L_15","Frontal_Inf_Orb_R_16","Olfactory_L_21","Olfactory_R_22","Frontal_Mid_Orb_L_25","Frontal_Mid_Orb_R_26","Insula_L_29","Insula_R_30",'TotalGrayVol','Left_Hippocampus_vol','Right_Hippocampus_vol','Left_Amygdala_vol','Right_Amygdala_vol']

SUBJECTS_DIR='/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/FS_git'

Subdirectories=[name for name in os.listdir(SUBJECTS_DIR) if os.path.isdir(os.path.join(SUBJECTS_DIR, name))]

remove_element=['fsaverage']
Subdirectories_1=[x for x in Subdirectories if x not in remove_element]

## Abtain the total_surface_area value of each subject
def regression_dict(filename):
	regression_aspect=pd.read_csv(filename, index_col=False, header=None)
	regression_aspect.columns = ['subj_ID','aspect']
	regression_aspect_dict=regression_aspect.set_index('subj_ID')['aspect'].to_dict()
	return regression_aspect_dict

## Abtain the total_gray_matter_volume value of each subject
def regression_dict_vol(filename):
	regression_aspect=pd.read_csv(filename, index_col=False, header=None)
	regression_aspect.columns = ['subj_ID','TotalGrayVol','Left_Hippocampus_vol','Right_Hippocampus_vol','Left_Amygdala_vol','Right_Amygdala_vol']
	regression_aspect_dict=regression_aspect.set_index('subj_ID').to_dict(orient='index')
	return regression_aspect_dict

## Generate dataframe for each subject's measurement values
def subject_df_generate(file_path,character_col_index,subject_id):
	subject_dict = {}
	df = pd.read_csv(file_path,sep=" ",header=None).iloc[:,[0,character_col_index]]	
	df.columns=["key","value"]
	result_dict = df.set_index("key").to_dict()["value"]
	subject_dict[subject_id]=result_dict
	subject_df=pd.DataFrame(subject_dict)
	return subject_df

## Aspect
#Surface area=2
#Gray matter volume=3

## Integrate all subjects'data into one dataframe
subject_without_htxt=[] ## Recode subjects without "_lh/rh_1.txt"
def subj_statistic_df_generated(aspect,aspect_columns,aspect_dict):
	subject_dict={}
	aspect_df=pd.DataFrame(columns=aspect_columns)
	for i in Subdirectories_1:
		file_name_lh=i+"_lh_1.txt"
		file_name_rh=i+"_rh_1.txt"
		file_path_lh=os.path.join(SUBJECTS_DIR,i,"stats",file_name_lh)
		file_path_rh=os.path.join(SUBJECTS_DIR,i,"stats",file_name_rh)
		if os.path.exists(file_path_lh) and os.path.exists(file_path_rh):
			subject_df_lh=subject_df_generate(file_path_lh,aspect,i)
			subject_df_rh=subject_df_generate(file_path_rh,aspect,i)
			subject_df_merge=pd.concat([subject_df_lh, subject_df_rh])
			subject_dict[subject_df_merge.columns[0]]={} 
			for index,row in subject_df_merge.iterrows():
				subject_dict[subject_df_merge.columns[0]][index]=row[subject_df_merge.columns[0]] ##  Dictionary structure: {key:{key1:value1,key2:value2}}
				if aspect==2:
					subject_dict[subject_df_merge.columns[0]]['Total_surface_area']=aspect_dict[subject_df_merge.columns[0]]
				else:
					subject_dict[subject_df_merge.columns[0]].update(aspect_dict[subject_df_merge.columns[0]])
		else:
			subject_without_htxt.append(i)
	## Fill data to dataframe
	for key, values in subject_dict.items():
		aspect_df.loc[key]={key: float(value) for key, value in values.items()}
	aspect_df.reset_index(inplace=True)
	aspect_df.rename(columns={'index':'Subject_ID'},inplace=True)
	return aspect_df


## (1) Surface area
total_surface_area_dict=regression_dict('/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/FS_git_surface_area.csv')

surface_area_df=subj_statistic_df_generated(2,aspect_columns_surface_area,total_surface_area_dict)
surface_area_df.to_csv('FS_git_surface_area_df.csv', index=False)

## (2) Gray matter volume
gray_matter_volume_dict=regression_dict_vol('/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/FS_git_gray_matter_volume.csv')

gray_matter_volume_df=subj_statistic_df_generated(3,aspect_columns_gray_matter_volume,gray_matter_volume_dict)
gray_matter_volume_df.to_csv('FS_git_gray_matter_volume_df.csv',index=False)

