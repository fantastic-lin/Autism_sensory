#!/usr/bin/python
#coding=utf-8
import os
import requests
import pandas as pd
import csv
import json

## (1)Count the number of files in each subject's directory
def countFiles(root_path):
        assert os.path.exists(root_path)
        total_files = 0
        item_list = os.listdir(root_path)
        if len(item_list) == 0:
                return 0
        for item in item_list:
                next_path = os.path.join(root_path, item)
                if os.path.isfile(next_path):
                        total_files += 1
                else:
                        total_files += countFiles(next_path)
        return total_files

## (2)Download files
def download_file(url, filename):
    ## Detect the size of local files and decide where to continue downloading 
    try:
        file_size = os.path.getsize(filename)
    except FileNotFoundError:
        file_size = 0
    headers = {'Range': f'bytes={file_size}-'} if file_size else None
    with requests.get(url, headers=headers, stream=True) as r:
        with open(filename, 'ab') as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

## (3)Retrieve all file names for a specific directory and obtain the absolute path
def get_all_files_in_directory(directory):
    all_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if root == directory:
                all_files.append(file)
            else:
                subdir = os.path.relpath(root, directory)
                file_path_sub_dir = os.path.join(subdir, file)
                all_files.append(file_path_sub_dir)
    return all_files

## (4)Output files that are not in the target list
def find_extra_elements(subject_list, target_list):
    set1 = set(subject_list)
    set2 = set(target_list)
    if set2.issubset(set1):
        return set()
    else:
        return set2.difference(set1)

## (5)Build a new directory
def ensure_directory_exists(directory):
    ## Check if the file exist
    if not os.path.exists(directory):
        ## If the file does not exist, create it
        os.makedirs(directory)
        print(f"Directory '{directory}' created.")
    else:
        print(f"Directory '{directory}' already exists.")



## File storage path
prefix_path='/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/FS_git' ## Convert to your file storage path
## File download path
s3_prefix='https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/freesurfer/5.1'
## Acquire subject_ID
Phenotypic_file=pd.read_csv('/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/Phenotypic_V1_0b_preprocessed1.csv',header=0,index_col=None) ## Pay attention to changing the file storage path

## Create a dictionary for the subfiles to be downloaded
with open("/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/each_sub_necessary_files.txt", "r") as f:
        file_name = f.read()
        file_name1=file_name.split("]")
        parsed_data = {}
for i in file_name1:
        key_value_pair = i.split(':')
        if len(key_value_pair) >= 2:
                key = key_value_pair[0].strip()
                value = key_value_pair[1].strip()
                value1 = value.replace('\n','').replace(' ','')
                value2 = value1.strip("[]")
                value3 = value2.split(",")
                parsed_data[key] = value3
        else:
                continue

file_name_paths = [f"{key}/{value}" for key, values in parsed_data.items() for value in values]

## Retrieve the existing subject folder in FS_git
downloaded_file=next(os.walk('/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/FS_git')) [1] ## Convert to your file storage path

## Loop of download files
error_load_dict=dict()
successed_download=0
file_identifier_number=0
file_identifier_list=[]
for index, file_identifier in Phenotypic_file['FILE_ID'].items():
	## For testing purposes, only the first five subjects will be downloaded. To analyze all data, please delete the following two lines
	if index >= 5:
        	break
	if file_identifier != "no_filename":
		file_identifier_number+=1
		file_ID_path= os.path.join(prefix_path,file_identifier)
		subject_file_list=get_all_files_in_directory(file_ID_path)
		unload_file=find_extra_elements(subject_file_list,file_name_paths)
		if unload_file:
			print(f'the necessary file that the subject did not dowload: {unload_file}')
			for file in list(unload_file):
				file_store_path=os.path.join(file_ID_path,file)
				file_store_path_dir=os.path.dirname(file_store_path)
				ensure_directory_exists(file_store_path_dir)
				file_download_address='/'.join([s3_prefix,file_identifier,file])
				try:
					download_file(file_download_address,file_store_path)
					print(f"file {file_download_address} have been download to local path {file_store_path}")
				except Exception as exc: 
					print(f"an error occurred while downloading the file: {exc}")
					error_load_dict[file_store_path]=str(exc)
		else:
			print("this subject have downloaded all neccessory file")
			successed_download+=1
		file_identifier_list.append(file_identifier)

print(file_identifier_list)

## Save the fully downloaded subject_ID into a JSON file
with open('FS_successful_download.json', 'w') as f:
    json.dump(file_identifier_list, f)
