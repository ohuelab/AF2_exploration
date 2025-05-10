from sklearn.metrics import roc_auc_score
from glob import glob
from tqdm import tqdm
import os
import pandas as pd
import sys
import zipfile

def compress_files_with_extension(directory, extension, zip_filename):
    with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(directory):
            for file in files:
                if file.endswith(extension):
                    file_path = os.path.join(root, file)
                    zipf.write(file_path, os.path.relpath(file_path, directory))

def count_non_hydrogen_atoms(pdbqt_file):
    non_hydrogen_count = 0
    with open(pdbqt_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_type = line[12:16].strip()
                if not atom_type.startswith("H"):
                    non_hydrogen_count += 1
    return non_hydrogen_count

def main():
    tgt_protein = sys.argv[1]
    log_ind_dir = sys.argv[2]
    lig_tgt_dir = sys.argv[3]
    docking_res_dir = sys.argv[4]

    inpaths=glob(f'{lig_tgt_dir}/*.pdbqt')
    paths=glob(f"{docking_res_dir}/*.pdbqt")

    # すでに処理済みならexit
    if os.path.exists(f"{docking_res_dir}/vina.csv"):
        exit(1)

    print(len(paths), len(inpaths))
    if len(paths) != len(inpaths):
        with open(f'{docking_res_dir}/rocauc.txt',mode='w') as f:
            f.write("0")
        assert len(paths) == len(inpaths)

    scores={}
    non_h_atoms = []
    for path in tqdm(paths):
        with open(path) as f:
            l=f.readline()
            l=f.readline()
        score=float(l.split()[3])
        ligname=os.path.basename(path).replace('_out.pdbqt','')
        # docking_score
        scores[ligname]=score
        # ligand efficiency
        lig_pdbqt = os.path.join(lig_tgt_dir, os.path.basename(path).replace('_out.pdbqt','.pdbqt'))
        non_h_atoms.append(count_non_hydrogen_atoms(lig_pdbqt))


    scores_df = pd.DataFrame(scores.values(), index=scores.keys(), columns=['vina score'])
    scores_df['label'] = scores_df.index.str.contains('actives')
    scores_df['non_h_atoms'] = non_h_atoms
    scores_df['le'] = scores_df['vina score'] / scores_df['non_h_atoms']

    scores_df.to_csv(f'{docking_res_dir}/vina.csv')

    vina_score = roc_auc_score(scores_df['label'], -scores_df['vina score'])
    vina_score_le = roc_auc_score(scores_df['label'], -scores_df['le'])
    print(f"vina_score: {vina_score}")

    with open(f'{docking_res_dir}/rocauc.txt',mode='w') as f:
        f.write(str(vina_score))

    # ファイルを圧縮
    output_zip_filename = f"{docking_res_dir}/compressed.zip"
    compress_files_with_extension(docking_res_dir, ".pdbqt", output_zip_filename)
    print(f".pdbqt compressed to {output_zip_filename}")

    # .pdbqtを削除
    for path in paths:
        try:
            os.remove(path)
        except Exception as e:
            print(f"{path} の削除中にエラーが発生しました: {e}")

if __name__ == "__main__":
    main()


