import os, sys
import numpy as np
import pandas as pd

from bindingsites import PROTEIN_BINDING_SITES
from utils import submit_job
from genetic_algo import GeneticAlgorithm

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

def gene_af2_optimization_population():
    ga = GeneticAlgorithm(LOGFILE_PATH)
    next_population = ga.genetic_algorithm(opt_target="AF2_PARAMS")
    return next_population

def gene_msa_optimization_population():
    ga = GeneticAlgorithm(LOGFILE_PATH)
    next_population = ga.genetic_algorithm(opt_target="MSA_MUTATION")
    return next_population

def gene_af2_random_population():
    # af2params
    af2_population = []
    individuals_t = []
    for param_range in config.af2_params_space:
        rand_vals = np.arange(param_range[0], param_range[1]+1, step=param_range[2])
        individuals_t.append(np.random.choice(rand_vals, config.population_size))
    af2_population = np.transpose(individuals_t)

    return pd.DataFrame(af2_population, columns=config.af2_columns)

def gene_msa_random_population():
    # msa mutation
    msa_population = []
    resi = PROTEIN_BINDING_SITES.get(config.tgt_protein)[config.msa_mutation_target]
    # 0: そのまま, 1: 変異
    for i in range(config.population_size):
        msa_mutation_binary = np.random.choice([0, 1], size=len(resi), p=[1-config.msa_mutation_p, config.msa_mutation_p])
        msa_population.append(','.join(map(str, msa_mutation_binary)))

    return pd.DataFrame(msa_population, columns=config.msa_columns)

def gene_init_msa_random_population():
    rand_msa_pop_df = gene_msa_random_population()

    for col in config.msa_columns:
        ind0_msa_gene = rand_msa_pop_df.loc[0, col]
        if col == "msa_mutation":
            rand_msa_pop_df.loc[0, col] = ind0_msa_gene.replace("1", "0")

    return rand_msa_pop_df


def gene_af2_default_population():
    af2_population = [config.default_af2_params for _ in range(config.population_size)]
    return pd.DataFrame(af2_population, columns=config.af2_columns)

def gene_msa_default_population():
    if config.default_msa_mutation == "none":
        resi = PROTEIN_BINDING_SITES.get(config.tgt_protein)[config.msa_mutation_target]
        default_msa_gene = ",".join(["0"] * len(resi))
        msa_population = [default_msa_gene for _ in range(config.population_size)]
    else:
        msa_population = [config.default_msa_mutation for _ in range(config.population_size)]
    return pd.DataFrame(msa_population, columns=config.msa_columns)

def export_next_population(next_population):
    log_dir = f"{LOG_ROOT_DIR}/loop_{config.now_loop}"
    os.makedirs(log_dir, exist_ok=True)
    next_population.to_csv(f"{log_dir}/population.csv", index=False)

    # population_log.csvに追記
    if os.path.exists(POP_LOGFILE_PATH):
        prev_pop_log_df = pd.read_csv(POP_LOGFILE_PATH)
        next_population["status"] = "pending"
        next_population["path"] = "none"
        pop_log_df = pd.concat([prev_pop_log_df, next_population], axis=0, ignore_index=True)
        pop_log_df = pop_log_df.drop_duplicates(subset=config.all_columns, keep="first")
        pop_log_df.to_csv(POP_LOGFILE_PATH, index=False)
        print(pop_log_df)
    else:
        # ファイルが存在しない場合、新しくファイルを作成
        pop_log_df = next_population
        pop_log_df["status"] = "pending"
        pop_log_df["path"] = "none"
        pop_log_df.to_csv(POP_LOGFILE_PATH, index=False)

# 初期実行用の関数
def gene_init_population():
    af2_gene_functions = {
        "optimize": gene_af2_random_population,
        "random": gene_af2_random_population,
        "default": gene_af2_default_population
    }

    msa_gene_functions = {
        "optimize": gene_init_msa_random_population,
        "random": gene_msa_random_population,
        "default": gene_msa_default_population
    }

    population_elems = []

    af2_population = get_population(config.optimization_option_af2, af2_gene_functions)
    population_elems.append(af2_population)

    msa_population = get_population(config.optimization_option_msa, msa_gene_functions)
    population_elems.append(msa_population)

    next_population = pd.concat(population_elems, axis=1)

    print(next_population)

    export_next_population(next_population)

def get_population(optimization_option, gene_functions):
    return gene_functions.get(optimization_option, gene_functions["default"])()

def gene_next_population():
    af2_gene_functions = {
        "optimize": gene_af2_optimization_population,
        "random": gene_af2_random_population,
        "default": gene_af2_default_population
    }

    msa_gene_functions = {
        "optimize": gene_msa_optimization_population,
        "random": gene_msa_random_population,
        "default": gene_msa_default_population
    }

    population_elems = []

    af2_population = get_population(config.optimization_option_af2, af2_gene_functions)
    population_elems.append(af2_population)

    msa_population = get_population(config.optimization_option_msa, msa_gene_functions)
    population_elems.append(msa_population)

    next_population = pd.concat(population_elems, axis=1)

    print("next_population")
    print(next_population)

    export_next_population(next_population)

    command = f"python scripts/gene_af2models.py {config.tgt_protein} {config.optimize_id}"
    submit_job(command)


if __name__ == "__main__":
    if config.now_loop == 1:
        gene_init_population()
    else:
        gene_next_population()