#!/bin/bash

# INPUT
tgt_protein=$1
log_ind_dir=$2
split_num=1
lig_dir=$4
output_dir=$5
gpu_num=4
docking_option="balance"

prep_rcp_path=`find $log_ind_dir -type f -name *_prep.pdbqt`

lig_num=`find $lig_dir -type f | wc -l`
total_integer=$lig_num

start_list=()
end_list=()
divisor=$((split_num * $gpu_num))
quotient=$((total_integer / divisor))
remainder=$((total_integer % divisor))
for ((i = 1; i <= divisor; i++)); do
    if ((i <= remainder)); then
        start=$(( (i-1) * (quotient+1) + 1 ))
        end=$(( i * (quotient+1) ))
    else
        start=$(( (i-1) * quotient + remainder + 1 ))
        end=$(( i * quotient + remainder ))
    fi
    start_list+=($start)
    end_list+=($end)
done

# make outputdir
mkdir -p $output_dir

batch=$(find $lig_dir -type f -name "*_prep.pdbqt")

# docking
unidock \
    --receptor $prep_rcp_path \
    --gpu_batch $batch \
    --config scripts/docking/config/config_${tgt_protein}.txt \
    --search_mode balance \
    --scoring vina \
    --dir $output_dir \
    --seed 1
