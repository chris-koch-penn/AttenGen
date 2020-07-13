# Chris Koch, 2020.
import pickle
from pathlib import Path
import random
from Bio import SeqIO
import csv


def save_victors(uniprot, vf_bacteria, vf_viruses):
    base = Path("./data")
    base.mkdir(parents=True, exist_ok=True)
    pickle.dump(uniprot, open(base / "victors_uniprot.pkl", "wb"))
    pickle.dump(vf_bacteria, open(base / "vf_bacteria.pkl", "wb"))
    pickle.dump(vf_viruses, open(base / "vf_viruses.pkl", "wb"))


def save_protegen(uniprot, vf_bacteria, vf_viruses):
    base = Path("./data")
    base.mkdir(parents=True, exist_ok=True)
    pickle.dump(uniprot, open(base / "protegen_uniprot.pkl", "wb"))
    pickle.dump(vf_bacteria, open(base / "protegen_bacteria.pkl", "wb"))
    pickle.dump(vf_viruses, open(base / "protegen_viruses.pkl", "wb"))


def train_test(input_file, sample_prob, num_samples, label):
    data = pickle.load(open(input_file, "rb"))
    x_train, x_test, y_train, y_test = [], [], [], []
    count = 0
    for _, features, _ in data:
        if random.random() < sample_prob:
            x_test += [features]
        else:
            x_train += [features]
        count += 1
        if count >= num_samples:
            break
    y_train = [label] * len(x_train)
    y_test = [label] * len(x_test)
    return {'x_train': x_train, 'x_test': x_test, 'y_train': y_train, 'y_test': y_test}


def train_test_bacteria(input_file, sample_prob, label):
    data = pickle.load(open(input_file, "rb"))
    x_train, x_test, y_train, y_test = [], [], [], []
    hits = get_bacterial_ids()
    for fastaID, features, _ in data:
        isMember = False
        for item in hits:
            if fastaID.find(item):
                isMember = True
                break
        if isMember:
            if random.random() < sample_prob:
                x_test += [features]
            else:
                x_train += [features]
    y_train = [label] * len(x_train)
    y_test = [label] * len(x_test)
    return {'x_train': x_train, 'x_test': x_test, 'y_train': y_train, 'y_test': y_test}


def get_bacterial_ids():
    """ Get a set of protein accession ids from input file. """
    col = "Protein Accession"
    path1 = "./data/victors_gneg.csv"
    path2 = "./data/victors_gpos.csv"
    s1 = [row[col] for row in csv.DictReader(open(path1)) if row[col]]
    s2 = [row[col] for row in csv.DictReader(open(path2)) if row[col]]
    return s1 + s2


def split_victors_train_test():
    b = Path("./feature_vectors")
    p1 = b / "victors_viruses.features.pkl"
    p2 = b / "victors_filtered_uniprot.features.pkl"
    p3 = b / "victors_all.features.pkl"
    vf_viruses = train_test(p1, 0.1, 5000, 1)
    vf_bacteria = train_test_bacteria(p3, 0.2, 1)
    length = len(vf_bacteria['x_train']) + len(vf_bacteria['x_test']) + \
        len(vf_viruses['x_train']) + len(vf_viruses['x_test'])
    uniprot = train_test(p2, 0.2, length, 0)
    save_victors(uniprot, vf_bacteria, vf_viruses)


def split_protegen_train_test():
    b = Path("./feature_vectors")
    p1 = b / "protegen_bacteria_proteins.features.pkl"
    p2 = b / "protegen_virus_proteins.features.pkl"
    p3 = b / "protegen_filtered_uniprot.features.pkl"
    vf_bacteria = train_test(p1, 0.1, 5000, 1)
    vf_viruses = train_test(p2, 0.1, 5000, 1)
    length = len(vf_bacteria['x_train']) + len(vf_bacteria['x_test']) + \
        len(vf_viruses['x_train']) + len(vf_viruses['x_test'])
    uniprot = train_test(p3, 0.2, length, 0)
    save_protegen(uniprot, vf_bacteria, vf_viruses)


if __name__ == "__main__":
    split_victors_train_test()
    split_protegen_train_test()
