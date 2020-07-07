import numpy as np
import pickle
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from genetic_algorithm import create_gene_start_end_idxs, GA_utils
import matplotlib.pyplot as plt
from matplotlib import pylab
import os


def filter_stop_codons(gene):
    codons = ["".join(gene[i:i+3]) for i in range(0, len(gene), 3)]
    temp = [c for c in codons[:-1] if c not in ["TAG", "TAA", "TGA"]]
    return "".join(temp + [codons[-1]])


def substitute_stop_codons(gene):
    codons = ["".join(gene[i:i+3]) for i in range(0, len(gene), 3)]
    temp = [("***" if c in {"TAG", "TAA", "TGA"} else c)
            for c in codons[:-1]]
    res = "".join(temp + [codons[-1]])
    return res


def decode_solution_to_genes(gene_start_end_idxs, solution):
    genes = [solution[start:end] for start, end in gene_start_end_idxs]
    return [filter_stop_codons(gene) for gene in genes]


def decode_solution_to_unfiltered_genes(gene_start_end_idxs, solution):
    genes = [solution[start:end] for start, end in gene_start_end_idxs]
    return [substitute_stop_codons(gene) for gene in genes]


def plot_quantitative_virulence(original_fitness_vals, mutated_fitness_vals):
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.35
    inds = np.arange(len(original_fitness_vals))
    ax.bar(inds, [-x for x, y in original_fitness_vals], width, color="blue")
    ax.bar(inds + width,
           [-x for x, y in mutated_fitness_vals], width, color="green")
    ax.legend(("Original Sequence", "Mutated Sequence"))
    fig.suptitle(
        "Quantitative Virulence Before and After Attenuation by Genetic Algorithm")
    pylab.ylabel("Quantitative Virulence")
    plt.xticks(range(len(original_fitness_vals)),
               range(1, 1 + len(original_fitness_vals)))
    plt.xlabel(
        "Coding DNA Sequences in Order of Appearance in Genome")
    plt.show()
    print("Quantitative virulence values after mutation: ", mutated_fitness_vals)


def plot_intial_quantitative_virulence(original_fitness_vals):
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.35
    inds = np.arange(len(original_fitness_vals))
    ax.bar(inds, [-x for x, y in original_fitness_vals], width, color="blue")
    fig.suptitle(
        "Initial Quantitative Virulence of Covid-19 Coding DNA Sequences")
    plt.xlabel(
        "Coding DNA Sequences in Order of Appearance in Genome")
    pylab.ylabel("Quantitative Virulence")
    plt.xticks(range(len(original_fitness_vals)),
               range(1, 1 + len(original_fitness_vals)))
    plt.show()
    print("Initial quantitative virulence values: ", original_fitness_vals)


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


def write_blast_report(base, dir1, f1):
    f1 = "./data/covid19_most_virulent_3genes.fna"
    f2 = Path(dir1) / os.listdir(dir1)[-1]
    cmd = f"blastn -query {f1} -subject {f2} -outfmt 7 -out {base}/blast.txt"
    os.system(cmd)


def generate_quantitative_virulence_graph1(fasta):
    original_genes = [seq for fid, seq in SimpleFastaParser(open(fasta))]
    ga_util_obj = create_GA_util_obj(fasta)
    og = [ga_util_obj.get_gene_fitness(gene) for gene in original_genes]
    plot_intial_quantitative_virulence(og)


def generate_quantitative_virulence_graph2(solution, fasta):
    original_genes = [seq for fid, seq in SimpleFastaParser(open(fasta))]
    gene_start_end_idxs = create_gene_start_end_idxs(original_genes)
    mutated_genes = decode_solution_to_unfiltered_genes(
        gene_start_end_idxs, solution)
    ga_util_obj = create_GA_util_obj(fasta)
    og = [ga_util_obj.get_gene_fitness(gene) for gene in original_genes]
    mut = [ga_util_obj.get_gene_fitness(
        gene.replace("*", "")) for gene in mutated_genes]
    plot_quantitative_virulence(og, mut)


def create_GA_util_obj(genome_path):
    victors_scores = "./saved_models/victors_xgboost_scores.joblib"
    victors_model_path = "./saved_models/victors_xgboost_model.joblib"
    protegen_scores = "./saved_models/protegen_xgboost_scores.joblib"
    protegen_model_path = "./saved_models/protegen_xgboost_model.joblib"
    return GA_utils(victors_scores, protegen_scores,
                    victors_model_path, protegen_model_path, genome_path)


if __name__ == "__main__":
    base = Path("./genetic_algorithm_output")
    most_virulent_fasta = Path("./data/covid19_most_virulent_3genes.fna")
    full_fasta = Path("./data/covid19_coding_sequences.fna")
    input_path = base / "best_ten_solutions_150gens_10000pop_1000mates.pkl"
    output_dir = base / "vaccine_candidates"
    top_ten_solutions = pickle.load(open(input_path, "rb"))
    best_sol = top_ten_solutions[0][1]
    write_solutions(base, most_virulent_fasta, top_ten_solutions)
    write_blast_report(base, output_dir, most_virulent_fasta)
    generate_quantitative_virulence_graph1(full_fasta)
    generate_quantitative_virulence_graph2(best_sol, most_virulent_fasta)
