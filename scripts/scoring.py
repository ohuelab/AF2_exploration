import sys, os
import subprocess
import re
from glob import glob
from datetime import datetime
import numpy as np
import time
import pandas as pd

from utils import submit_job, submit_job_once

from config_loader import ConfigLoader

if len(sys.argv) > 1:
    TGT_PROTEIN = sys.argv[1]
    OPTIMIZE_ID = sys.argv[2]
    LOG_ROOT_DIR = f"dat/optimize_log/{TGT_PROTEIN}/{OPTIMIZE_ID}"
    LOGFILE_PATH = f"{LOG_ROOT_DIR}/progress_log.json"
    POP_LOGFILE_PATH = f"{LOG_ROOT_DIR}/population_log.csv"
    config = ConfigLoader(LOGFILE_PATH)
else:
    print("CAN NOT CONFIG LOAD!")


def exploration():
    log_dir = f"{LOG_ROOT_DIR}/loop_{config.now_loop}"

    i = int(sys.argv[3])
    print(f"[exploration] >> individual_{i}")
    tgt_dir = f"{log_dir}/population/individual_{i}"

    lig_tgt = f"exploration_{config.active_size}_{config.split_id}"

    # prep
    if not any(file.endswith("_prep.pdbqt") for file in os.listdir(tgt_dir)):
        command = f"scripts/prep_rcp.sh {config.tgt_protein} {config.optimize_id} {tgt_dir}"
        submit_job_once(command)

    # docking
    lig_tgt_dir = f"dat/dude_cv/{config.tgt_protein}/{lig_tgt}/train"
    docking_res_dir = f"{tgt_dir}/{lig_tgt}"
    if not os.path.exists(f"{docking_res_dir}/vina.csv"):
        command = f"scripts/unidock_all.sh {config.tgt_protein} {tgt_dir} 1 {lig_tgt_dir} {docking_res_dir}"
        submit_job_once(command)

        # rocauc
        command = f"python scripts/all_docking_rocauc.py {config.tgt_protein} {tgt_dir} {lig_tgt_dir} {docking_res_dir}"
        submit_job_once(command)

    if i == config.population_size-1:
        # wait
        docking_finished_num = 0
        while docking_finished_num != config.population_size:
            docking_finished_num = 0
            for i in range(config.population_size):
                log_ind_dir = f"{log_dir}/population/individual_{i}/exploration_{config.active_size}_{config.split_id}"
                docking_finished_num += len(glob(f"{log_ind_dir}/rocauc.txt"))
            time.sleep(5)

        rocauc_scores = get_overall_rocauc(log_dir, f"exploration_{config.active_size}_{config.split_id}")
        top_score_index = rocauc_scores.index(max(rocauc_scores))
        tgt_dir = f"{log_dir}/population/individual_{top_score_index}"
        lig_tgt = f"exploration_{config.active_size}_{config.split_id}"
        docking_res_dir = f"{tgt_dir}/{lig_tgt}_test"

        if not os.path.exists(f"{docking_res_dir}/vina.csv"):
            command = f"python scripts/top_score_rocauc.py {config.tgt_protein} {config.optimize_id}"
            submit_job_once(command)


        command = f"python scripts/record_rocauc.py {config.tgt_protein} {config.optimize_id}"
        submit_job_once(command)


    return "[scoring] rocauc_unidock is running..."

def get_overall_rocauc(log_dir, tgt_method):
    rocauc_scores = []
    for i in range(config.population_size):
        log_ind_dir = f"{log_dir}/population/individual_{i}/{tgt_method}"
        try:
            with open(f"{log_ind_dir}/rocauc.txt") as f:
                vina_score = float(f.read())
        except:
            vina_score = 0

        rocauc_scores.append(vina_score)
    return rocauc_scores

def unidock_all():
    log_dir = f"{LOG_ROOT_DIR}/loop_{config.now_loop}"

    i = int(sys.argv[3])
    print(f"[unidock_all_nonparallel] >> individual_{i}")
    tgt_dir = f"{log_dir}/population/individual_{i}"

    lig_tgt = f"all"

    # prep
    if not any(file.endswith("_prep.pdbqt") for file in os.listdir(tgt_dir)):
        command = f"scripts/prep_rcp.sh {config.tgt_protein} {config.optimize_id} {tgt_dir}"
        submit_job_once(command)

    # docking
    lig_tgt_dir = f"dat/dude_cv/{config.tgt_protein}/{lig_tgt}"
    docking_res_dir = f"{tgt_dir}/unidock_{lig_tgt}"
    if not os.path.exists(f"{docking_res_dir}/vina.csv"):
        command = f"scripts/unidock_all.sh {config.tgt_protein} {tgt_dir} 1 {lig_tgt_dir} {docking_res_dir}"
        submit_job_once(command)


        # rocauc
        command = f"python scripts/all_docking_rocauc.py {config.tgt_protein} {tgt_dir} {lig_tgt_dir} {docking_res_dir}"
        submit_job_once(command)


    if i == config.population_size-1:
        docking_finished_num = 0
        while docking_finished_num != config.population_size:
            docking_finished_num = 0
            for i in range(config.population_size):
                log_ind_dir = f"{log_dir}/population/individual_{i}/unidock_all"
                docking_finished_num += len(glob(f"{log_ind_dir}/rocauc.txt"))
            time.sleep(5)

        command = f"./scripts/record_rocauc.sh {config.tgt_protein} {config.optimize_id}"
        submit_job_once(command)


    return "[scoring] rocauc_unidock is running..."


def default_case():
    return "invalid option"




def switch_case(option):
    switch_dict = {
        "unidock_all": unidock_all,
        "exploration": exploration,
    }
    selected_case = switch_dict.get(option, default_case)
    return selected_case()

def scoring_core():
    now = datetime.now()
    datetime_string = now.strftime("%Y-%m-%d_%H%M")

    print(f"[scoring] optimize_id : {config.optimize_id}")
    print(f"[scoring] loop : {config.now_loop}")
    print(f"[scoring] scoring_method : {config.scoring_method}")
    print(f"[scoring] date : {datetime_string}")

    result = switch_case(config.scoring_method)
    print(result)

if __name__ == "__main__":
    scoring_core()
