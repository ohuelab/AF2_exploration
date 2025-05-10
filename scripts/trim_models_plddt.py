import sys
import glob
import os
import json
import re
from Bio.PDB import *
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

# alphafold出力ディレクトリが対象
# ディレクトリ内のpdbを全て対応するscore_jsonでtrim

SCORE_THRESH = 50

class NotLowPlddt(Select):
    def __init__(self, plddt_scores):
        self.plddt_scores = plddt_scores
    def accept_residue(self, residue):
        if self.plddt_scores[residue.id[1] - 1] >= SCORE_THRESH:
            return True
        else:
            return False

def main():
    tgt_dir = sys.argv[3]

    pathlist = glob.glob(f"{tgt_dir}/*.pdb")

    for path in pathlist:
        if config.alphafold_version == 2:
            match = re.search(r"rank_(\d+)_.*?model_(\d+)", path)
            rank, model = match.groups() if match else (None, None)
            scores_pathlist = glob.glob(f"{os.path.dirname(path)}/*scores*.json")
            matching_score_path = [path for path in scores_pathlist if re.search(rf"rank_{rank}_.*?model_{model}", path)][0]
            print(path)
            print(matching_score_path)

            try:
                with open(matching_score_path) as f:
                    d = json.load(f)
            except Exception as e:
                print(f'[ERROR] {type(e)}:{str(e)} : {matching_score_path}')
                exit(1)

            plddt_scores = d["plddt"]

            if sum(plddt_scores)/len(plddt_scores) < SCORE_THRESH:
                exit(1)

            # pdb load
            parser = PDBParser()
            struc = parser.get_structure("tgt_protein", path)

            # trimming by plddt
            io=PDBIO()
            io.set_structure(struc)
            out_path = f"{os.path.splitext(path)[0]}_trim.pdb"
            io.save(out_path, select=NotLowPlddt(plddt_scores))

        elif config.alphafold_version == 3:
            matching_score_path = f"{os.path.dirname(path)}/alphafold3/alphafold3_confidences.json"
            print(path)
            print(matching_score_path)

            try:
                with open(matching_score_path) as f:
                    d = json.load(f)
            except Exception as e:
                print(f'[ERROR] {type(e)}:{str(e)} : {matching_score_path}')
                exit(1)

            plddt_scores = d["atom_plddts"]

            # pdb load
            parser = PDBParser()
            struc = parser.get_structure("tgt_protein", path)

            # trimming by plddt
            io=PDBIO()
            io.set_structure(struc)
            out_path = f"{os.path.splitext(path)[0]}_trim.pdb"
            io.save(out_path)


if __name__ == '__main__':
    main()