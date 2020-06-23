# Chris Koch, 2020.
from genetic_algorithm_lib import GA as genetic_alg
import numpy as np
from propy.PyPro import GetProDes
from collections import Counter
import joblib
from pathlib import Path
from multiprocessing import Pool
import pickle
import matplotlib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import time
THRESHOLD = 0.70


class GA_utils:
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    AMINO_ACID_MAP = dict([(aa, count + 1)
                           for count, aa in enumerate(amino_acids)])
    INDICES_TO_AMINO_ACIDS = dict([(v, k) for k, v in AMINO_ACID_MAP.items()])

    def __init__(self, victors_path, protegen_path, victors_model_path, protegen_model_path, fasta):
        self.sorted_victors_scores = np.sort(
            joblib.load(victors_path)[:, 1]).tolist()
        self.sorted_protegen_scores = np.sort(
            joblib.load(protegen_path)[:, 1]).tolist()
        self.victors_score_counts = Counter(self.sorted_victors_scores)
        self.protegen_score_counts = Counter(self.sorted_protegen_scores)
        self.protegen_model = joblib.load(protegen_model_path)
        self.victors_model = joblib.load(victors_model_path)
        run_on_cpu(self.protegen_model)
        run_on_cpu(self.victors_model)
        seqs = [str(s.seq) for s in SeqIO.parse(fasta, "fasta")]
        self.gene_start_end_idxs = create_gene_start_end_idxs(seqs)
        self.concatenated_genes = [item for seq in seqs for item in seq]
        self.DNA_TO_FITNESS_MAP = dict()

    def calc_percentile(self, score, sorted_arr, count_set):
        rank = np.searchsorted(sorted_arr, score) + 1
        frequency = count_set.get(score, 0)
        return (rank + (0.5 * frequency)) / len(sorted_arr)

    def get_gene_fitness(self, gene):
        protein = dna_to_protein(gene)
        protein = "".join([c for c in protein if c != "*"])
        Des = GetProDes(protein)
        features = Des.GetAAComp()
        features.update(Des.GetMoreauBrotoAuto())
        features.update(Des.GetGearyAuto())
        features.update(Des.GetCTD())
        features.update(Des.GetQSO())
        features = np.array([list(features.values())])
        vf_res = self.victors_model.predict_proba(features)[0][1]
        # protegen_res = self.protegen_model.predict_proba(features)[0][1]
        quantitative_virulence = self.calc_percentile(
            vf_res, self.sorted_victors_scores, self.victors_score_counts)
        # protegenicity = self.calc_percentile(
        #     protegen_res,  self.sorted_protegen_scores, self.protegen_score_counts)
        quantitative_virulence = 0 if quantitative_virulence < THRESHOLD else quantitative_virulence
        return (-quantitative_virulence, 0)  # protegenicity)


def create_gene_start_end_idxs(seqs):
    gene_start_end_idxs = []
    start, end = 0, 0
    for gene in seqs:
        end = start + len(gene)
        gene_start_end_idxs.append((start, end))
        start = end
    return gene_start_end_idxs


def run_on_cpu(gridsearchcv):
    model = gridsearchcv.estimator.steps[1][1]
    model = model.set_params(predictor='cpu_predict')
    model = model.set_params(tree_method='hist')


def dna_to_protein(seq):
    return Seq("".join(seq), generic_dna).translate()[:-1]


def fitness_func(ga_util, solution, solution_idx):
    genes = [solution[start:end]
             for start, end in ga_util.gene_start_end_idxs]
    output = [ga_util.get_gene_fitness(gene) for gene in genes]
    return sum([qv + p for qv, p in output])


def callback_generation(GA):
    print(f"Generation = {GA.generations_completed}")
    print(f"Fitness    = {GA.best_solution()[1]}")


def ALGO(generations, pop_size, num_parents_mating, ga_util_obj: GA_utils):
    # Creating an instance of the Genetic Algorithm.
    initial_pop = [np.copy(ga_util_obj.concatenated_genes)
                   for i in range(pop_size)]
    GA = genetic_alg(num_generations=generations,
                     num_parents_mating=num_parents_mating,
                     fitness_func=fitness_func,
                     initial_population=initial_pop,
                     parent_selection_type="rank",
                     keep_parents=num_parents_mating,
                     crossover_type="two_points",
                     mutation_type="swap",
                     mutation_num_genes=1,
                     callback_generation=callback_generation,
                     ga_util_obj=ga_util_obj)

    # Running the GA to optimize the parameters of the function.
    GA.run()
    GA.plot_result(
        title=f"Generation vs. Fitness for Population of {pop_size}",  show_plot=False)

    # Returning the details of the best solution.
    solution, solution_fitness, solution_idx = GA.best_solution()
    print(f"Parameters of the best solution : {solution}")
    print(f"Fitness value of the best solution = {solution_fitness}")
    print(f"Index of the best solution : {solution_idx}")
    if GA.best_solution_generation != -1:
        print(f"Best fitness value reached after \
            {GA.best_solution_generation} generations.")

    # Save the GA instance.
    # GA.save("./GA") # This file contains all solutions and GA parameters, but is often large.
    output = Path("./genetic_algorithm_output")
    output.mkdir(parents=True, exist_ok=True)
    suffix = f"{generations}gens_{pop_size}pop.pkl"
    f1 = output / ("best_solution_" + suffix)
    f2 = output / ("best_ten_solutions_" + suffix)
    # pickle.dump((solution_fitness, solution), open(f1, "wb"))
    pickle.dump(GA.top_ten_solutions, open(f2, "wb"))


def run_GA(victors_scores, protegen_scores, victors_model_path,
           protegen_model_path, genome_path, num_generations, pop_size, num_parents_mating):
    ga_util_obj = GA_utils(victors_scores, protegen_scores,
                           victors_model_path, protegen_model_path, genome_path)
    ALGO(num_generations, pop_size, num_parents_mating, ga_util_obj)


if __name__ == "__main__":
    victors_scores = "./victors_xgboost_scores.joblib"
    protegen_scores = "./protegen_xgboost_scores.joblib"
    victors_model_path = "./victors_xgboost_model.joblib"
    protegen_model_path = "./protegen_xgboost_model.joblib"
    covid_genome_path = "./covid19_coding_sequences.fna"
    # victors_scores = "./saved_models/victors_xgboost_scores.joblib"
    # protegen_scores = "./saved_models/protegen_xgboost_scores.joblib"
    # victors_model_path = "./saved_models/victors_xgboost_model.joblib"
    # protegen_model_path = "./saved_models/protegen_xgboost_model.joblib"
    # covid_genome_path = "./data/covid19_most_virulent_3genes.fna"
    start = time.time()
    run_GA(victors_scores, protegen_scores, victors_model_path,
           protegen_model_path, covid_genome_path, 25, 50000, 500)
    end = time.time()
    print("ELAPSED TIME: ", end - start)
