import numpy as np
import pandas as pd
from datetime import datetime
from config_loader import ConfigLoader

class GeneticAlgorithm:
    def __init__(self, LOGFILE_PATH):
        self.config = ConfigLoader(LOGFILE_PATH)

        self.tgt_protein = self.config.tgt_protein
        self.optimize_id = self.config.optimize_id
        self.now_loop = self.config.now_loop
        self.population_size = self.config.population_size
        self.af2_mutation_rate = self.config.af2_mutation_rate
        self.msa_mutation_rate = self.config.msa_mutation_rate
        self.af2_crossover_p = self.config.af2_crossover_p
        self.msa_crossover_p = self.config.msa_crossover_p
        self.af2_params_space = self.config.af2_params_space
        self.af2_columns = self.config.af2_columns
        self.msa_columns = self.config.msa_columns
        self.LOG_ROOT_DIR = f"dat/optimize_log/{self.tgt_protein}/{self.optimize_id}"
        self.POP_LOGFILE_PATH = f"{self.LOG_ROOT_DIR}/population_log.csv"

    def id_to_seed(self):
        """
        Optimize ID から乱数シードを生成
        """
        date_part, time_part = self.optimize_id.split('_')
        year, month, day = [int(part) for part in date_part.split('-')]
        hour, minute = int(time_part[:2]), int(time_part[2:])
        seed = year * 100000000 + month * 1000000 + day * 10000 + hour * 100 + minute
        return seed % (2 ** 32)

    def roulette_selection(self, pop_df):
        """
        ルーレット選択
        """
        if self.now_loop == 2:
            pop_df = pop_df.sort_values(by='score', ascending=False).head(5)
            pop_df = pop_df.reset_index(drop=True)
            population = [pop_df.drop(columns=['score']).iloc[[i]].reset_index(drop=True) for i in range(len(pop_df))]
            return population
        else:
            population = [pop_df.drop(columns=['score']).iloc[[i]].reset_index(drop=True) for i in range(len(pop_df))]
            fitness_values = pop_df["score"].to_list()

            scaled_fitness = [fv**2 for fv in fitness_values]
            total_scaled_fitness = sum(scaled_fitness)

            # スケーリング後の確率を計算
            probabilities = [sf / total_scaled_fitness for sf in scaled_fitness]

            selected_population = []
            for _ in range(len(population)):
                selected_idx = np.random.choice(range(len(population)), size=1, p=probabilities)
                selected_population.append(population[selected_idx[0]])
            np.random.shuffle(selected_population)
            return selected_population

    def tournament_selection(self, pop_df):
        k = 5
        # 現在の集団（個体のリスト）
        population = [pop_df.drop(columns=['score']).iloc[[i]].reset_index(drop=True) for i in range(len(pop_df))]

        # 適応度（fitness）のリスト
        fitness_values = pop_df["score"].to_list()

        # トーナメント選択による新しい集団の生成
        selected_population = []
        for _ in range(len(population)):
            # トーナメントの実施: ランダムにk個体を選ぶ
            tournament_indices = np.random.choice(range(len(population)), size=k, replace=False)
            # トーナメント内で最も適応度が高い個体を選択
            best_idx = max(tournament_indices, key=lambda idx: fitness_values[idx])
            selected_population.append(population[best_idx])

        # 選択された個体をシャッフル
        np.random.shuffle(selected_population)
        return selected_population

    def crossover(self, parent1, parent2, target):
        """
        交叉操作
        """
        if target == "AF2_PARAMS":
            child1, child2 = self.crossover_af2params(parent1, parent2)
        elif target == "MSA_MUTATION":
            child1, child2 = self.crossover_msa(parent1, parent2)
        return child1, child2

    def crossover_af2params(self, parent1, parent2):
        parent1 = parent1[self.af2_columns].iloc[0].tolist()
        parent2 = parent2[self.af2_columns].iloc[0].tolist()
        if np.random.uniform() < self.af2_crossover_p:
            child1, child2 = self.uniform_crossover(parent1, parent2)
        else:
            child1, child2 = parent1, parent2
        child1 = pd.DataFrame([child1], columns=self.af2_columns)
        child2 = pd.DataFrame([child2], columns=self.af2_columns)
        return child1, child2

    def crossover_msa(self, parent1, parent2):
        parent1_str = parent1[self.msa_columns].iloc[0].tolist()[0]
        parent2_str = parent2[self.msa_columns].iloc[0].tolist()[0]
        parent1 = list(eval(parent1_str))
        parent2 = list(eval(parent2_str))
        if np.random.uniform() < self.msa_crossover_p:
            child1, child2 = self.two_point_crossover(parent1, parent2)
        else:
            child1, child2 = parent1, parent2
        child1 = pd.DataFrame([str(child1)], columns=self.msa_columns)
        child2 = pd.DataFrame([str(child2)], columns=self.msa_columns)
        return child1, child2

    def uniform_crossover(self, parent1, parent2):
        child1, child2 = [], []
        for gene1, gene2 in zip(parent1, parent2):
            if np.random.uniform() < 0.5:
                child1.append(gene1)
                child2.append(gene2)
            else:
                child1.append(gene2)
                child2.append(gene1)
        return child1, child2

    def two_point_crossover(self, parent1, parent2):
        length = len(parent1)
        crossover_points = np.sort(np.random.choice(range(1, length), 2, replace=False))
        child1 = parent1[:crossover_points[0]] + parent2[crossover_points[0]:crossover_points[1]] + parent1[crossover_points[1]:]
        child2 = parent2[:crossover_points[0]] + parent1[crossover_points[0]:crossover_points[1]] + parent2[crossover_points[1]:]
        return child1, child2

    def mutate(self, individual, target):
        # 各パラメータに対して、一定の確率で突然変異を実行
        if target == "AF2_PARAMS":
            mutated_individual = self.mutate_af2params(individual)
            mutated_individual = pd.DataFrame([mutated_individual], columns=self.af2_columns)
        elif target == "MSA_MUTATION":
            mutated_individual = self.mutate_msa(individual)
            mutated_individual = pd.DataFrame([mutated_individual], columns=self.msa_columns)
        return mutated_individual

    def mutate_af2params(self, individual):
        ind_af2 = individual[self.af2_columns].iloc[0].tolist()
        mutated_af2 = []
        for param, param_range in zip(ind_af2, self.af2_params_space):
            if np.random.uniform() < self.af2_mutation_rate:
                rand_vals = np.arange(param_range[0], param_range[1] + 1, step=param_range[2])
                mutated_af2.append(np.random.choice(rand_vals))
            else:
                mutated_af2.append(param)
        return mutated_af2

    def mutate_msa(self, individual):
        ind_msa_str = individual[self.msa_columns].iloc[0].tolist()[0]
        ind_msa = list(eval(ind_msa_str))
        mutated_msa = ind_msa[:]
        for i in range(len(mutated_msa)):
            if np.random.uniform() < self.msa_mutation_rate:
                mutated_msa[i] = 1 - mutated_msa[i]
        return ','.join(map(str, mutated_msa))

    def genetic_algorithm(self, opt_target):
        now = datetime.now()
        datetime_string = now.strftime("%Y-%m-%d_%H%M")
        print(f"[new_genes] optimize_id : {self.optimize_id}")
        print(f"[new_genes] loop : {self.now_loop}")
        print(f"[new_genes] date : {datetime_string}")

        prev_pop_path = f"{self.LOG_ROOT_DIR}/loop_{self.now_loop-1}/population.csv"
        prev_pop = pd.read_csv(prev_pop_path)

        selected_population = self.tournament_selection(prev_pop)

        # 最適化対象に応じてdfを抽出
        if opt_target == "AF2_PARAMS":
            selected_population = [df[self.af2_columns] for df in selected_population]
        elif opt_target == "MSA_MUTATION":
            selected_population = [df[self.msa_columns] for df in selected_population]
        print("selected_population : \n", selected_population)

        # 交叉を実行
        offspring = []
        if self.now_loop == 2:
            for i in range(0, self.population_size-2, 2):
                idx1, idx2 = np.random.choice(len(selected_population), size=2, replace=False)
                child1, child2 = self.crossover(selected_population[idx1], selected_population[idx2], target=opt_target)
                offspring.extend([child1, child2])
        else:
            for i in range(0, self.population_size-2, 2):
                child1, child2 = self.crossover(selected_population[i], selected_population[i + 1], target=opt_target)
                offspring.extend([child1, child2])
        print("offspring : \n", offspring)

        mutated_offspring = [self.mutate(ind, target=opt_target) for ind in offspring]
        print("mutated_offspring : \n", mutated_offspring)

        # elite個体を追加
        if opt_target == "AF2_PARAMS":
            selected_prev_pop = prev_pop[self.af2_columns+["score"]]
        elif opt_target == "MSA_MUTATION":
            selected_prev_pop = prev_pop[self.msa_columns+["score"]]
        prev_pop_sorted = selected_prev_pop.sort_values(by="score", ascending=False)
        prev_pop_sorted = prev_pop_sorted.drop_duplicates(subset=prev_pop_sorted.columns[:-1])
        prev_pop_sorted = prev_pop_sorted.head(2)
        elite_ind = prev_pop_sorted.iloc[:, :-1]
        for i in range(2):
            mutated_offspring.append(elite_ind.iloc[[i]].reset_index(drop=True))
            print(elite_ind.iloc[[i]].reset_index(drop=True))
        print("mutated_offspring : \n",mutated_offspring)

        new_pop_df = pd.concat(mutated_offspring, axis=0)
        # df化, 連番のindexをつける
        new_pop_df.reset_index(drop=True, inplace=True)
        print("next population : \n", new_pop_df)

        return new_pop_df