import os, sys
from datetime import datetime
import yaml
import json

from utils import submit_job, submit_job_once
from config_loader import ConfigLoader

CONFIG_PATH = "config.yaml"

def main():
    # 初回実行
    if len(sys.argv) == 1:
        print(f"[optimize_core] initial optimization!")

        # optimize_idの作成
        now = datetime.now()
        datetime_string = now.strftime("%Y-%m-%d_%H%M%S")
        optimize_id = datetime_string

        # configのひな型作成
        now_loop = 1
        params = {"id": datetime_string,
                "now_loop": now_loop}

        with open(CONFIG_PATH, 'r') as f:
            yaml_config = yaml.safe_load(f)
        params.update(yaml_config)

        # LOG_ROOT_DIR, LOGFILE_PATHの作成
        TGT_PROTEIN = params["TGT_PROTEIN"]
        LOG_ROOT_DIR = f"dat/optimize_log/{TGT_PROTEIN}/{optimize_id}"
        LOGFILE_PATH = f"{LOG_ROOT_DIR}/progress_log.json"

        # log(progress_log.json)の作成
        print(f"[optimize_core] make {LOGFILE_PATH}...")
        os.makedirs(LOG_ROOT_DIR, exist_ok=True)
        with open(LOGFILE_PATH, "w") as f:
            json.dump(params, f, indent=4)

        # configのload
        config = ConfigLoader(LOGFILE_PATH)

        # 初期個体の生成
        print(f"[optimize_core] generate init population...")
        command = f"python scripts/gene_next_population.py {config.tgt_protein} {config.optimize_id}"

        submit_job_once(command)

        if config.alphafold_version == 2:
            print(f"[optimize_core] generate af2models...")
            command = f"python scripts/gene_af2models.py {config.tgt_protein} {config.optimize_id}"
            submit_job_once(command)

    else:
        print(f"[optimize_core] optimization is continuing! ")
        # 実行
        print(f"[optimize_core] generate next population & af2_models...")
        TGT_PROTEIN = sys.argv[1]
        OPTIMIZE_ID = sys.argv[2]

        LOG_ROOT_DIR = f"dat/optimize_log/{TGT_PROTEIN}/{OPTIMIZE_ID}"
        LOGFILE_PATH = f"{LOG_ROOT_DIR}/progress_log.json"

        # configのload
        config = ConfigLoader(LOGFILE_PATH)

        command = f"python scripts/gene_af2models.py {config.tgt_protein} {config.optimize_id}"

        submit_job_once(command)


if __name__ == "__main__":
    main()