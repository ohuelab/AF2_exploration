import pandas as pd
import numpy as np
import sys
import json
from datetime import datetime
from scipy.spatial.distance import pdist, squareform
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

def scale_rocauc(rocauc):
    return rocauc

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
    log_dir = f"{LOG_ROOT_DIR}/loop_{config.now_loop}"

    now = datetime.now()
    datetime_string = now.strftime("%Y-%m-%d_%H%M")

    print(f"[record_rocauc] optimize_id : {config.optimize_id}")
    print(f"[record_rocauc] loop : {config.now_loop}")
    print(f"[record_rocauc] date : {datetime_string}")

    if "exploration" == config.scoring_method:
        rocauc_scores = get_overall_rocauc(log_dir, f"exploration_{config.active_size}_{config.split_id}")
        scaled_rocaucs = [scale_rocauc(x) for x in rocauc_scores]

        # average rocauc
        ave_score = np.mean(rocauc_scores)

        # rocaucをpopulation.csvに追記
        pop_df = pd.read_csv(f"{log_dir}/population.csv")
        pop_df["score"] = [scale_rocauc(x) for x in rocauc_scores]
        pop_df.to_csv(f"{log_dir}/population.csv", index=False)

        # population_log.csvに保存
        pop_log_df = pd.read_csv(POP_LOGFILE_PATH)
        for i, individual in pop_df.iterrows():
            matching_mask = (pop_log_df[config.all_columns] == individual[config.all_columns]).all(axis=1)
            matching_idx = matching_mask.idxmax()
            pop_log_df.loc[matching_idx, "score"] = scaled_rocaucs[i]
            pop_log_df.loc[matching_idx, config.scoring_method] = rocauc_scores[i]
            pop_log_df.loc[matching_idx, "status"] = "finished"
            pop_log_df.loc[matching_idx, "path"] = f"{log_dir}/population/individual_{i}"

        pop_log_df.to_csv(POP_LOGFILE_PATH, index=False)
        pop_log_df.to_csv(f"{log_dir}/population_log.csv", index=False)

        # 今世代のMSA mutationの平均ハミング距離を計算し収束度を見る（ランダムで0.5）
        msa_df = pop_df["msa_mutation"]
        convert_to_list = lambda x: list(map(int, x.split(',')))
        bit_arrays = msa_df.apply(convert_to_list).tolist()
        bit_matrix = np.array(bit_arrays)
        distances = pdist(bit_matrix, metric='hamming')
        distance_matrix = squareform(distances)
        average_hamming_distance = np.mean(distance_matrix[np.triu_indices(len(msa_df), k=1)])

        # top rocaucの取得
        top_score_index = rocauc_scores.index(max(rocauc_scores))
        rocauc_txt = f"{log_dir}/population/individual_{top_score_index}/exploration_{config.active_size}_{config.split_id}_test/rocauc.txt"
        with open(rocauc_txt, "r") as f:
            top_rocauc = float(f.read())

        # progress.jsonを更新
        max_ser = pop_df.loc[pop_df["score"].idxmax()]
        max_param_dict = max_ser.to_dict()
        max_param_dict[f"{config.scoring_method}_ave_score"] = ave_score
        max_param_dict[f"{config.scoring_method}_score"] = max(rocauc_scores)
        max_param_dict["run_date"] = datetime_string
        max_param_dict["ind"] = pop_df["score"].idxmax()
        max_param_dict["average_hamming_distance"] = average_hamming_distance
        max_param_dict["top_rocauc"] = top_rocauc

        with open(LOGFILE_PATH, "r") as f:
            params = json.load(f)
        params[config.now_loop] = max_param_dict
        params["now_loop"] += 1
        with open(LOGFILE_PATH, "w") as f:
            json.dump(params, f, indent=4)

    else:
        if "unidock_all" in config.scoring_method:
            docking_res_dirname = "unidock_all"
        rocauc_scores = get_overall_rocauc(log_dir, docking_res_dirname)
        scaled_rocaucs = [scale_rocauc(x) for x in rocauc_scores]

        # average rocauc
        ave_score = np.mean(rocauc_scores)

        # rocaucをpopulation.csvに追記
        pop_df = pd.read_csv(f"{log_dir}/population.csv")
        pop_df["score"] = [scale_rocauc(x) for x in rocauc_scores]
        pop_df.to_csv(f"{log_dir}/population.csv", index=False)

        # population_log.csvに保存
        pop_log_df = pd.read_csv(POP_LOGFILE_PATH)
        for i, individual in pop_df.iterrows():
            matching_mask = (pop_log_df[config.all_columns] == individual[config.all_columns]).all(axis=1)
            matching_idx = matching_mask.idxmax()
            pop_log_df.loc[matching_idx, "score"] = scaled_rocaucs[i]
            pop_log_df.loc[matching_idx, docking_res_dirname] = rocauc_scores[i]
            pop_log_df.loc[matching_idx, "status"] = "finished"
            pop_log_df.loc[matching_idx, "path"] = f"{log_dir}/population/individual_{i}"
        pop_log_df.to_csv(POP_LOGFILE_PATH, index=False)

        # 今世代のMSA mutationの平均ハミング距離を計算し収束度を見る（ランダムで0.5）
        msa_df = pop_df["msa_mutation"]
        convert_to_list = lambda x: list(map(int, x.split(',')))
        bit_arrays = msa_df.apply(convert_to_list).tolist()
        bit_matrix = np.array(bit_arrays)
        distances = pdist(bit_matrix, metric='hamming')
        distance_matrix = squareform(distances)
        average_hamming_distance = np.mean(distance_matrix[np.triu_indices(len(msa_df), k=1)])

        # top rocaucの取得
        top_score_index = rocauc_scores.index(max(rocauc_scores))
        rocauc_txt = f"{log_dir}/population/individual_{top_score_index}/unidock_all/rocauc.txt"
        with open(rocauc_txt, "r") as f:
            top_rocauc = float(f.read())

        # progress.jsonを更新
        max_ser = pop_df.loc[pop_df["score"].idxmax()]
        max_param_dict = max_ser.to_dict()
        max_param_dict[f"{docking_res_dirname}_ave_score"] = ave_score
        max_param_dict[f"{docking_res_dirname}_score"] = max(rocauc_scores)
        max_param_dict["run_date"] = datetime_string
        max_param_dict["ind"] = pop_df["score"].idxmax()
        max_param_dict["average_hamming_distance"] = average_hamming_distance
        max_param_dict["top_rocauc"] = top_rocauc

        with open(LOGFILE_PATH, "r") as f:
            params = json.load(f)
        params[config.now_loop] = max_param_dict
        params["now_loop"] += 1
        with open(LOGFILE_PATH, "w") as f:
            json.dump(params, f, indent=4)


    if config.generations > config.now_loop:
        command = f"python scripts/gene_next_population.py {config.tgt_protein} {config.optimize_id}"
        submit_job_once(command)

    else:
        now = datetime.now()
        datetime_string = now.strftime("%Y-%m-%d_%H%M")
        print(datetime_string)

if __name__ == "__main__":
    main()


