import os
import sys
import json
import subprocess
from datetime import datetime
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

def main():

    now = datetime.now()
    datetime_string = now.strftime("%Y-%m-%d_%H%M")

    print(f"[record_rocauc] optimize_id : {config.optimize_id}")
    print(f"[record_rocauc] loop : {config.now_loop}")
    print(f"[record_rocauc] date : {datetime_string}")

    log_dir = f"{LOG_ROOT_DIR}/loop_{config.now_loop}"

    # top scoreの個体のindexを取得
    rocauc_scores = get_overall_rocauc(log_dir, f"exploration_{config.active_size}_{config.split_id}")
    top_score_index = rocauc_scores.index(max(rocauc_scores))

    tgt_dir = f"{log_dir}/population/individual_{top_score_index}"
    print("tgt_dir:", tgt_dir)

    # prep
    if not any(file.endswith("_prep.pdbqt") for file in os.listdir(tgt_dir)):
        command = f"scripts/prep_rcp.sh {config.tgt_protein} {config.optimize_id} {tgt_dir}"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        print(command)

    # docking
    lig_tgt = f"exploration_{config.active_size}_{config.split_id}"
    lig_tgt_dir = f"dat/dude_cv/{config.tgt_protein}/{lig_tgt}/test"
    docking_res_dir = f"{tgt_dir}/{lig_tgt}_test"

    if not os.path.exists(f"{docking_res_dir}/vina.csv"):
        command = f"scripts/unidock_all.sh {config.tgt_protein} {tgt_dir} {config.t3_dock_array_num} {lig_tgt_dir} {docking_res_dir}"
        print(command)
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

        command = f"python scripts/all_docking_rocauc.py {config.tgt_protein} {tgt_dir} {lig_tgt_dir} {docking_res_dir}"
        print(command)
        result = subprocess.run(command, shell=True, capture_output=True, text=True)


if __name__ == "__main__":
    main()


