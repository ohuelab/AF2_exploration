#!/bin/sh

# jaxのseed固定
export TF_CUDNN_DETERMINISTIC=1

tgt_protein=$1
max_seq=$2
max_extra_seq=$3
num_recycle=$4
random_seed=$5
num_models=$6
model_order=$7
msa_mutation_status_str=$8
output_dir=$9
optimize_id=${10}

INPUTFILE="dat/tgt_protein/$tgt_protein/input.fasta"

mkdir -p $output_dir
# msaのcopy
cp -r dat/af2_msa_env/$tgt_protein/* $output_dir
tgt_ind_msa_dir=$(echo $output_dir/*_env )
# msa mutation
python scripts/mutate_msa.py $tgt_protein $optimize_id $tgt_ind_msa_dir "$msa_mutation_status_str"

#source ~/.bashrc
chmod +x ${INPUTFILE}

# colabfold_batch \
#     --num-recycle ${num_recycle} \
#     --random-seed ${random_seed} \
#     --max-seq ${max_seq} \
#     --max-extra-seq ${max_extra_seq} \
#     --num-models ${num_models} \
#     --model-order ${model_order} \
#     ${INPUTFILE} \
#     ${output_dir}

