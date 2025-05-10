import sys, os
from glob import glob

from bindingsites import PROTEIN_BINDING_SITES

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

# MSAを変異させる関数mutate_msaを含む

# mutate_msa( a3m_lines, { 15: "A", 155: "A" } )
def mutate_msa(a3m_lines, pos_res) -> str:

    for target_res in pos_res.values():
        assert len(target_res) == 1

    output = []

    # Iterate over alignment lines
    for line in a3m_lines.split("\n"):
        if line.startswith(">"):
            output.append(line)
        elif len(line) > 1:
            line = list(line)
            for pos, res in pos_res.items():
                if line[pos-1] in "ACDEFGHIKLMNPQRSTVWY":
                    line[pos-1] = res
            output.append("".join(line))
        else:
            output.append(line)

    return "\n".join(output)

def main():
    args = sys.argv
    tgt_dir = args[3]
    msa_mutation_status_str = args[4]

    if msa_mutation_status_str == "none":
        print("[mutate_msa] MSA mutation was not executed.")
        exit(1)
    
    msa_mutation_status = eval(msa_mutation_status_str)

    # 変異する残基のdictを作成する
    resi = PROTEIN_BINDING_SITES.get(config.tgt_protein)[config.msa_mutation_target]

    # 残基番号と変異するしないを対応付ける
    mapping = dict(zip(resi, msa_mutation_status))

    # pos_res: 変異対象の残基辞書 {残基番号: "A"}
    pos_res = {}
    for key, value in mapping.items():
        if value == 0:
            continue  # valueが0の場合は何もせずに次の要素へ
        elif value == 1:
            pos_res[key] = "A"  # valueが1の場合は"A"に変更
    print("pos_res:", pos_res)

    # a3mファイルのパスリストを取得
    tgt_a3m_pathlist = glob(f"{tgt_dir}/**/*.a3m", recursive=True)

    for tgt_a3m_path in tgt_a3m_pathlist:
        # a3mを読み込む
        with open (tgt_a3m_path, "r") as f:
            a3m_lines = f.read()

        # mutate MSA
        mutated_msa = mutate_msa(a3m_lines, pos_res)

        # output
        with open (tgt_a3m_path, "w") as f:
            f.write(mutated_msa)

if __name__ == "__main__":
    main()
