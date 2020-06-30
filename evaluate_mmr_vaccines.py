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
from evaluate_ga_output import create_GA_util_obj


def fn(x, y):
    return -x  # if -x < 0.5 else 0


def plot_quantitative_virulence(original_fitness_vals, mutated_fitness_vals):
    original_fitness_vals = [fn(x, y) for x, y in original_fitness_vals]
    mutated_fitness_vals = [fn(x, y) for x, y in mutated_fitness_vals]
    print("QV: ", sum(original_fitness_vals), sum(mutated_fitness_vals))
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.35
    inds = np.arange(len(original_fitness_vals))
    ax.bar(inds, original_fitness_vals, width, color="blue")
    ax.bar(inds + width, mutated_fitness_vals, width, color="green")
    ax.legend(("Original Sequence", "Mutated Sequence"))
    fig.suptitle(
        "Quantitative Virulence Before and After Attenuation by Genetic Algorithm")
    pylab.ylabel("Quantitative Virulence")
    plt.xticks(range(len(original_fitness_vals)),
               range(1, 1 + len(original_fitness_vals)))
    plt.xlabel(
        "Coding DNA Sequences in Order of Appearance in Genome")
    plt.show()


def substitute_ambiguous_symbols(gene):
    invalid_chars = gene.replace("A", "").replace(
        "C", "").replace("G", "").replace("T", "")
    if invalid_chars != "":
        print(f"NON-STANDARD CHARACTERS FOUND IN DNA: {invalid_chars}")
    return gene.replace("X", "A").replace("N", "A").replace("Y", "C").replace("K", "G").replace("R", "A").replace("W", "T")


def generate_quantitative_virulence_graph1(fasta1, fasta2):
    og = [seq for fid, seq in SimpleFastaParser(open(fasta1))]
    mut1 = [seq for fid, seq in SimpleFastaParser(open(fasta2))]
    min_len = min(len(og), len(mut1))
    ga_util_obj = create_GA_util_obj(fasta1)
    og = [ga_util_obj.get_gene_fitness(gene)for gene in og][:min_len]
    mut1 = [ga_util_obj.get_gene_fitness(
        substitute_ambiguous_symbols(gene)) for gene in mut1][:min_len]
    # no_of_decreased_virulence_seqs = sum(
    #     [1 for x, y in zip(og, mut1) if fn(x[0], x[1]) > fn(y[0], y[1])])
    # no_of_increased_virulence_seqs = sum(
    #     [1 for x, y in zip(og, mut1) if fn(x[0], x[1]) < fn(y[0], y[1])])
    # print("Number of decreased: ", no_of_decreased_virulence_seqs)
    # print("Number of increased: ", no_of_increased_virulence_seqs)
    plot_quantitative_virulence(og, mut1)


def rubella():
    rubella_wild_type = Path("./data/rubella/wild_type.fna")
    matsuba = Path("./data/rubella/matsuba_vaccine.fna")
    matsuura = Path("./data/rubella/matsuura_vaccine.fna")
    takahashi = Path("./data/rubella/takahashi_vaccine.fna")
    TCRB19 = Path("./data/rubella/TCRB19_vaccine.fna")
    TO_336 = Path("./data/rubella/TO_336_vaccine.fna")

    # Graph rubella strains.
    generate_quantitative_virulence_graph1(rubella_wild_type, matsuba)
    generate_quantitative_virulence_graph1(rubella_wild_type, matsuura)
    generate_quantitative_virulence_graph1(rubella_wild_type, takahashi)
    generate_quantitative_virulence_graph1(rubella_wild_type, TCRB19)
    generate_quantitative_virulence_graph1(rubella_wild_type, TO_336)


def mumps():
    mumps_ODATE3 = Path("./data/mumps/mumps_ODATE3.fna")
    jeryl_lynn_minor = Path("./data/mumps/jeryl_lynn_minor.fna")
    jeryl_lynn_major = Path("./data/mumps/jeryl_lynn_major.fna")
    S79_minor_vaccine = Path("./data/mumps/S79_minor_vaccine.fna")
    S79_major_vaccine = Path("./data/mumps/S79_major_vaccine.fna")
    zagreb_vaccine = Path("./data/mumps/zagreb_vaccine.fna")

    # Graph mumps strains.
    generate_quantitative_virulence_graph1(mumps_ODATE3, jeryl_lynn_minor)
    generate_quantitative_virulence_graph1(mumps_ODATE3, jeryl_lynn_major)
    generate_quantitative_virulence_graph1(mumps_ODATE3, S79_minor_vaccine)
    generate_quantitative_virulence_graph1(mumps_ODATE3, S79_major_vaccine)
    generate_quantitative_virulence_graph1(mumps_ODATE3, zagreb_vaccine)


def measles():
    measles_wild_type = Path("./data/measles/edmonston_wild_type.fna")
    AIK_C_vaccine = Path("./data/measles/AIK_C_vaccine.fna")
    moraten = Path("./data/measles/moraten_vaccine.fna")
    rubeovax = Path("./data/measles/rubeovax_vaccine.fna")
    schwarz = Path("./data/measles/schwarz_vaccine.fna")
    zagreb = Path("./data/measles/zagreb_vaccine.fna")
    CAM_70 = Path("./data/measles/CAM_70.fna")
    leningrad_4 = Path("./data/measles/leningrad_4.fna")
    shanghai_191 = Path("./data/measles/shanghai_191.fna")
    changchun_70 = Path("./data/measles/changchun_70.fna")

    # Graph measles strains.
    generate_quantitative_virulence_graph1(measles_wild_type, AIK_C_vaccine)
    generate_quantitative_virulence_graph1(measles_wild_type, moraten)
    generate_quantitative_virulence_graph1(measles_wild_type, rubeovax)
    generate_quantitative_virulence_graph1(measles_wild_type, schwarz)
    generate_quantitative_virulence_graph1(measles_wild_type, zagreb)
    generate_quantitative_virulence_graph1(measles_wild_type, CAM_70)
    generate_quantitative_virulence_graph1(measles_wild_type, leningrad_4)
    generate_quantitative_virulence_graph1(measles_wild_type, shanghai_191)
    generate_quantitative_virulence_graph1(measles_wild_type, changchun_70)


if __name__ == "__main__":
    measles()
    mumps()
    rubella()
