import numpy as np
from Bio import pairwise2
import pickle
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from genetic_algorithm import create_gene_start_end_idxs, GA_utils
import matplotlib.pyplot as plt


def filter_stop_codons(gene):
    codons = ["".join(gene[i:i+3]) for i in range(0, len(gene), 3)]
    temp = [c for c in codons[:-1] if c not in ["TAG", "TAA", "TGA"]]
    return "".join(temp + [codons[-1]])


def substitute_stop_codons(gene):
    codons = ["".join(gene[i:i+3]) for i in range(0, len(gene), 3)]
    temp = ["***" if c in ["TAG", "TAA", "TGA"] else c for c in codons[:-1]]
    res = "".join(temp + [codons[-1]])
    return res


def decode_solution_to_genes(gene_start_end_idxs, solution):
    genes = [solution[start:end] for start, end in gene_start_end_idxs]
    return [filter_stop_codons(gene) for gene in genes]


def decode_solution_to_unfiltered_genes(gene_start_end_idxs, solution):
    genes = [solution[start:end] for start, end in gene_start_end_idxs]
    return [substitute_stop_codons(gene) for gene in genes]


def analyze_mutations(original_genes, mutated_genes, fasta):
    ga_util_obj = create_GA_util_obj(fasta)
    og = [ga_util_obj.get_gene_fitness(gene) for gene in original_genes]
    mut = [ga_util_obj.get_gene_fitness(gene) for gene in mutated_genes]
    plot(og, mut)


def plot(original_fitness_vals, mutated_fitness_vals):
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.bar(list(range(len(original_fitness_vals))), original_fitness_vals)
    plt.show()


def write_solutions(base, fasta, top_ten_solutions):
    seqs = [seq for fid, seq in SimpleFastaParser(open(fasta))]
    ids = [fid for fid, seq in SimpleFastaParser(open(fasta))]
    gene_start_end_idxs = create_gene_start_end_idxs(seqs)
    count = 1
    for fitness, solution in top_ten_solutions:
        genes = decode_solution_to_genes(gene_start_end_idxs, solution)
        records = [SeqRecord(Seq(gene), fid, '', '')
                   for fid, gene in zip(ids, genes)]
        f = round(fitness, 5)
        path = base / \
            f"vaccine_candidates/solution_rank_{count}_fitness_{f}.fasta"
        SeqIO.write(records, open(path, "w"), "fasta")
        count += 1
    return top_ten_solutions


def write_report(top_ten, fasta):
    original_genes = [seq for fid, seq in SimpleFastaParser(open(fasta))]
    gene_start_end_idxs = create_gene_start_end_idxs(original_genes)
    avg_mutations, avg_stop_codons = 0, 0
    for _, solution in top_ten:
        mutated_genes = decode_solution_to_unfiltered_genes(
            gene_start_end_idxs, solution)
        num_mutations = []
        stop_codons = []
        for g1, g2 in zip(original_genes, mutated_genes):
            stop_codons.append(g2.count("*") / 3)
            mutations = sum([1 for c1, c2 in zip(g1, g2) if c1 != c2])
            num_mutations.append(mutations)
        avg_mutations += sum(num_mutations)
        avg_stop_codons += sum(stop_codons)
        analyze_mutations(original_genes, mutated_genes, fasta)
    avg_mutations /= 10
    avg_stop_codons /= 10
    print(avg_mutations, avg_stop_codons)


def analyze_stop_codons():
    pass


def create_GA_util_obj(genome_path):
    victors_scores = "./saved_models/victors_xgboost_scores.joblib"
    protegen_scores = "./saved_models/protegen_xgboost_scores.joblib"
    victors_model_path = "./saved_models/victors_xgboost_model.joblib"
    protegen_model_path = "./saved_models/protegen_xgboost_model.joblib"
    return GA_utils(victors_scores, protegen_scores,
                    victors_model_path, protegen_model_path, genome_path)


if __name__ == "__main__":
    base = Path("./genetic_algorithm_output")
    fasta = Path("./data/covid19_coding_sequences.fna")
    input_path = base / "best_ten_solutions_25gens_10000pop.pkl"
    top_ten_solutions = pickle.load(open(input_path, "rb"))
    # top_ten_solutions = write_solutions(base, fasta, input_path)
    write_report(top_ten_solutions, fasta)
    # output = Path("./genetic_algorithm_output")
    # output.mkdir(parents=True, exist_ok=True)
