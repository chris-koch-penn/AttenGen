# Chris Koch, 2020.
import os
from pathlib import Path


# Note that we use a really high E-value here. This is to catch any sequence that is even remotely similar - we will be filtering all sequences with greater than 30% percent identity and that will be our set of negative training examples. This is not supposed to be a typical blast search, it's more used as a way to calculate all of the possible percent identities in the dataset.
def run_blast(path, db, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / path.stem
    cmd = f"blastp -query {path.absolute()} -out {out} -db {db} -num_threads 16 -outfmt 7 -evalue 1000"
    os.system(cmd)


def make_blast_db(path):
    cmd = f"makeblastdb -in {path} -dbtype prot -out {path} -title {path.name}"
    os.system(cmd)


if __name__ == "__main__":
    # Make the database. Sometimes this will throw an error on Windows - Google it and follow steps to fix system settings.
    path_db = Path("./data/uniprot_sprot.fasta")
    make_blast_db(path_db)

    # Run Blast on Victors and Protegen proteins.
    db = "./data/uniprot_sprot.fasta"
    out = "./data/blast_output/"
    paths = ["./data/victors_all.faa", "./data/protegen_bacteria_proteins.faa",
             "./data/protegen_virus_proteins.faa"]
    for input_path in paths:
        run_blast(Path(input_path), db, Path(out))
