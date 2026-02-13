#!/usr/bin/env python3
# coding=utf-8

import os
import json
import requests
import pandas as pd
from typing import List, Set, Dict, Tuple

# ============================================================
# (1) Count the number of files in each subject's directory
# ============================================================
def count_files(root_path: str) -> int:
    if not os.path.exists(root_path):
        return 0
    total_files = 0
    for item in os.listdir(root_path):
        next_path = os.path.join(root_path, item)
        if os.path.isfile(next_path):
            total_files += 1
        else:
            total_files += count_files(next_path)
    return total_files


# ============================================================
# (2) Download files (supports resume safely)
# ============================================================
def download_file(url: str, filename: str, chunk_size: int = 8192, timeout: int = 60) -> Tuple[bool, str]:
    """
    Returns (ok, message). Handles Range resume safely.
    - If local file exists and server supports Range (206), append.
    - If server does NOT support Range (200) but local file exists, overwrite to avoid corruption.
    """
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    local_size = os.path.getsize(filename) if os.path.exists(filename) else 0
    headers = {"Range": f"bytes={local_size}-"} if local_size > 0 else {}

    try:
        with requests.get(url, headers=headers, stream=True, timeout=timeout) as r:
            # If server ignores Range and returns full content (200),
            # appending would corrupt the file. Overwrite instead.
            if r.status_code == 200 and local_size > 0:
                local_size = 0

            r.raise_for_status()

            mode = "ab" if (local_size > 0 and r.status_code == 206) else "wb"
            with open(filename, mode) as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)

        return True, "ok"

    except Exception as exc:
        return False, str(exc)


# ============================================================
# (3) Retrieve all file names for a specific directory (relative paths)
# ============================================================
def get_all_files_in_directory(directory: str) -> List[str]:
    if not os.path.exists(directory):
        return []
    all_files: List[str] = []
    for root, _, files in os.walk(directory):
        for file in files:
            if root == directory:
                all_files.append(file)
            else:
                subdir = os.path.relpath(root, directory)
                all_files.append(os.path.join(subdir, file))
    return all_files


# ============================================================
# (4) Output files that are missing from subject_list compared to target_list
# ============================================================
def find_missing_elements(subject_list: List[str], target_list: List[str]) -> Set[str]:
    set_subject = set(subject_list)
    set_target = set(target_list)
    return set_target.difference(set_subject)


# ============================================================
# (5) Build a new directory
# ============================================================
def ensure_directory_exists(directory: str) -> None:
    os.makedirs(directory, exist_ok=True)


# ============================================================
# Parse each_sub_necessary_files.txt (your original format)
# ============================================================
def parse_required_files(txt_path: str) -> Dict[str, List[str]]:
    """
    Expects a format roughly like:
      key: [a,b,c]
      key2: [d,e]
    Your original code split by ']' and ':', so we replicate a safer version.
    """
    parsed: Dict[str, List[str]] = {}
    with open(txt_path, "r", encoding="utf-8") as f:
        raw = f.read()

    # Split by ']' blocks; each block may contain "key: [.."
    blocks = raw.split("]")
    for blk in blocks:
        blk = blk.strip()
        if not blk or ":" not in blk:
            continue
        key, rest = blk.split(":", 1)
        key = key.strip()
        rest = rest.strip()

        # remove leading '[' if present and cleanup
        rest = rest.replace("\n", "").replace(" ", "")
        rest = rest.lstrip("[")

        if not key:
            continue
        values = [v for v in rest.split(",") if v]
        if values:
            parsed[key] = values

    return parsed


def build_required_paths(parsed_data: Dict[str, List[str]]) -> List[str]:
    """
    Build relative paths like f"{key}/{value}" for each value under each key.
    """
    out: List[str] = []
    for key, values in parsed_data.items():
        for v in values:
            out.append(os.path.join(key, v))
    return out


# ============================================================
# Main
# ============================================================
def main():
    # ---------- Config ----------
    prefix_path = "/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/FS_git"
    s3_prefix = "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/freesurfer/5.1"

    phenotypic_csv = "/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/Phenotypic_V1_0b_preprocessed1.csv"
    required_txt = "/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/each_sub_necessary_files.txt"

    # For testing only: set to None to download all subjects
    max_subjects = 5

    # Output
    out_success_json = os.path.join(prefix_path, "FS_successful_download.json")
    out_error_json = os.path.join(prefix_path, "FS_download_errors.json")

    # ---------- Load data ----------
    ensure_directory_exists(prefix_path)

    Phenotypic_file = pd.read_csv(phenotypic_csv, header=0, index_col=None)

    parsed_data = parse_required_files(required_txt)
    file_name_paths = build_required_paths(parsed_data)

    # ---------- Download loop ----------
    error_log: Dict[str, str] = {}
    successful_ids: List[str] = []
    attempted_ids: List[str] = []

    file_identifier_number = 0

    for index, file_identifier in Phenotypic_file["FILE_ID"].items():
        if max_subjects is not None and index >= max_subjects:
            break

        if file_identifier == "no_filename":
            continue

        file_identifier_number += 1
        attempted_ids.append(file_identifier)

        file_ID_path = os.path.join(prefix_path, file_identifier)
        ensure_directory_exists(file_ID_path)

        subject_file_list = get_all_files_in_directory(file_ID_path)
        missing_files = find_missing_elements(subject_file_list, file_name_paths)

        if missing_files:
            print(f"[{file_identifier}] Missing {len(missing_files)} required files. Downloading...")
            for rel_path in sorted(missing_files):
                file_store_path = os.path.join(file_ID_path, rel_path)
                ensure_directory_exists(os.path.dirname(file_store_path))

                file_download_address = "/".join([s3_prefix, file_identifier, rel_path.replace("\\", "/")])

                ok, msg = download_file(file_download_address, file_store_path)
                if ok:
                    print(f"  OK: {file_download_address} -> {file_store_path}")
                else:
                    print(f"  ERR: {file_download_address} | {msg}")
                    error_log[file_store_path] = msg

            # Re-check after downloads
            subject_file_list2 = get_all_files_in_directory(file_ID_path)
            missing_after = find_missing_elements(subject_file_list2, file_name_paths)
            if not missing_after and not any(p.startswith(file_ID_path + os.sep) for p in error_log.keys()):
                successful_ids.append(file_identifier)
                print(f"[{file_identifier}] All required files present ✅")
            else:
                print(f"[{file_identifier}] Still missing {len(missing_after)} files or had errors ❌")

        else:
            # Already complete
            successful_ids.append(file_identifier)
            print(f"[{file_identifier}] Already has all required files ✅")

    print("\nAttempted subject IDs:")
    print(attempted_ids)

    # ---------- Save outputs ----------
    with open(out_success_json, "w", encoding="utf-8") as f:
        json.dump(successful_ids, f, indent=2)

    with open(out_error_json, "w", encoding="utf-8") as f:
        json.dump(error_log, f, indent=2)

    print("\nSaved:")
    print(f"  Successful IDs: {out_success_json} ({len(successful_ids)})")
    print(f"  Errors:         {out_error_json} ({len(error_log)})")
    print(f"  Total attempted: {len(attempted_ids)}")


if __name__ == "__main__":
    main()
