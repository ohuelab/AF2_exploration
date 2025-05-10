#!/bin/bash

# PATH
MGL_PATH=/home/0/uw02030/uchikawa/apps/mgltools_x86_64Linux2_1.5.7

# INPUT
tgt_protein=$1
optimize_id=$2
log_ind_dir=$3

# log_ind_dirに*_trim_align_prep.pdbqtが存在しないなら
if ! ls "$log_ind_dir"/*_trim_align_prep.pdbqt 1> /dev/null 2>&1; then
    echo "[prep_rcp] don't exist trim_align_prep.pdbqt, starting the process!"

    python scripts/trim_models_plddt.py $tgt_protein $optimize_id $log_ind_dir
    echo "[prep_rcp] finish triming"
    sleep 1s

    # align
    trim_rcp_path=`find $log_ind_dir -type f -name *_trim.pdb`
    if [ -z "$trim_rcp_path" ]; then
        echo "Error: No matching file found in $log_ind_dir" >&2
        exit 1
    fi
    python scripts/align.py $tgt_protein $trim_rcp_path
    echo "[prep_rcp] finish align"
    sleep 1s

    # prep receptor
    input_rcp_path=`find $log_ind_dir -type f -name *_trim_align.pdb`
    prep_rcp_path=${input_rcp_path/.pdb/_prep.pdbqt}

    if [ -z "$input_rcp_path" ]; then
        echo "Error: No matching file found in $log_ind_dir" >&2
        exit 1
    fi

    echo "[prep_rcp] prepare_receptor4"
    LD_LIBRARY_PATH=${MGL_PATH}/lib ${MGL_PATH}/bin/pythonsh \
    ${MGL_PATH}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
    -r $input_rcp_path \
    -o $prep_rcp_path
else
    echo "[prep_rcp] exist trim_align_prep.pdbqt, finish the process!"
fi
