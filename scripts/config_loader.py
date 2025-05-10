import json


class ConfigLoader:
    def __init__(self, config_path):
        self._params = self._load_config(config_path)

    def _load_config(self, config_path):
        with open(config_path, "r") as f:
            return json.load(f)

    @property
    def optimize_id(self):
        return self._params["id"]

    @property
    def now_loop(self):
        return self._params["now_loop"]

    @property
    def tgt_protein(self):
        return self._params["TGT_PROTEIN"]

    @property
    def scoring_method(self):
        return self._params["SCORING_METHOD"]

    @property
    def alphafold_version(self):
        return self._params["ALPHAFOLD_VERSION"]

    @property
    def population_size(self):
        if self.now_loop == 1:
            return self._params["INIT_POPULATION_SIZE"]
        return self._params["POPULATION_SIZE"]

    @property
    def generations(self):
        return self._params["GENERATIONS"]

    @property
    def active_size(self):
        return self._params["ACTIVE_SIZE"]

    @property
    def split_id(self):
        return self._params["SPLIT_ID"]

    @property
    def parallel(self):
        return self._params["PARALLEL"]

    @property
    def af2_mutation_rate(self):
        return self._params["AF2_MUTATION_RATE"]

    @property
    def msa_mutation_rate(self):
        return self._params["MSA_MUTATION_RATE"]

    @property
    def af2_crossover_p(self):
        return self._params["AF2_CROSSOVER_P"]

    @property
    def msa_crossover_p(self):
        return self._params["MSA_CROSSOVER_P"]

    @property
    def optimization_option_af2(self):
        return self._params["OPTIMIZATION_OPTION_AF2"]

    @property
    def optimization_option_msa(self):
        return self._params["OPTIMIZATION_OPTION_MSA"]

    @property
    def msa_mutation_p(self):
        return self._params["MSA_MUTATION_P"]

    @property
    def msa_mutation_target(self):
        return self._params["MSA_MUTATION_TARGET"]

    @property
    def af2_params_space(self):
        return self._params["AF2_PARAMS_SPACE"]

    @property
    def all_columns(self):
        return self._params["ALL_COLUMNS"]

    @property
    def af2_columns(self):
        return self._params["AF2_COLUMNS"]

    @property
    def msa_columns(self):
        return self._params["MSA_COLUMNS"]

    @property
    def default_af2_params(self):
        return self._params["DEFAULT_AF2_PARAMS"]

    @property
    def default_msa_mutation(self):
        return self._params["DEFAULT_MSA_MUTATION"]
