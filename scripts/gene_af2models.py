import sys, os
from datetime import datetime
import pandas as pd
import time
import shutil
from glob import glob

from utils import submit_job_once
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

def generate_af2model():
    log_dir = f"{LOG_ROOT_DIR}/loop_{config.now_loop}"
    population_df = pd.read_csv(f"{log_dir}/population.csv")
    pop_log_df = pd.read_csv(POP_LOGFILE_PATH)

    now = datetime.now()
    datetime_string = now.strftime("%Y-%m-%d_%H%M")

    print(f"[gene_af2models] optimize_id : {config.optimize_id}")
    print(f"[gene_af2models] loop : {config.now_loop}")
    print(f"[gene_af2models] date : {datetime_string}")


    # generate
    hold_jid = ""
    for i, individual in population_df.iterrows():
        print(f"[gene_af2models] individual_{i}: {individual}")

        matching_mask = (pop_log_df[config.all_columns] == individual[config.all_columns]).all(axis=1)
        matching_row = pop_log_df[matching_mask].iloc[0]
        # このパラメータがまだ評価されていない新規パラメータなとき
        if (matching_row["status"] == "pending"):
            print(f"[gene_af2models] individual_{i} is pending!  preparing a job!")

            # params
            max_seq = individual["max_seq"]
            max_extra_seq = individual["max_extra_seq"]
            num_recycle = individual["num_recycle"]
            model_order = individual["model_order"]
            random_seed = 1
            num_models = 1
            msa_mutation_status_str = individual["msa_mutation"]

            # job
            if len(glob(f"{log_dir}/population/individual_{i}/*alphafold2*.pdb")) == 0:
                command = f"./scripts/run_localcolab.sh " \
                        f"{config.tgt_protein} {max_seq} {max_extra_seq} {num_recycle} " \
                        f'{random_seed} {num_models} {model_order} "{msa_mutation_status_str}" {log_dir}/population/individual_{i} {config.optimize_id}'
                submit_job_once(command)

            command = f"python scripts/scoring.py {config.tgt_protein} {config.optimize_id} {i}"
            submit_job_once(command)

        else:
            print(f"[gene_af2models] individual_{i} is finished!  reusing past results!")

            source_dir = matching_row["path"]
            destination_dir = f"{log_dir}/population/individual_{i}"
            if not os.path.exists(destination_dir):
                print(f"[gene_af2models] copy : {source_dir} -> {destination_dir}")
                shutil.copytree(source_dir, destination_dir, dirs_exist_ok=True)
            else:
                print(f"[gene_af2models] exist (do not copy) : {source_dir}")

            # tmp
            # params
            max_seq = individual["max_seq"]
            max_extra_seq = individual["max_extra_seq"]
            num_recycle = individual["num_recycle"]
            model_order = individual["model_order"]
            random_seed = 1
            num_models = 1
            msa_mutation_status_str = individual["msa_mutation"]

            # job
            command = f"./scripts/run_localcolab.sh " \
                        f"{config.tgt_protein} {max_seq} {max_extra_seq} {num_recycle} " \
                        f'{random_seed} {num_models} {model_order} "{msa_mutation_status_str}" {log_dir}/population/individual_{i} {config.optimize_id}'
            submit_job_once(command)


            command = f"python scripts/scoring.py {config.tgt_protein} {config.optimize_id} {i}"
            submit_job_once(command)

        time.sleep(0.5)




if __name__ == "__main__":
    generate_af2model()