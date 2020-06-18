# Chris Koch, 2020.
import os
from uuid import uuid4
from pathlib import Path


def run_blast(path, db, output_dir):
    out = output_dir / str(uuid4())
    cmd = f"blastp -query {path.absolute()} -out {out} -db {db} -num_threads 16 -outfmt 7 -evalue 1000"
    os.system(cmd)


def make_blast_db(path):
    cmd = f"makeblastdb -in {path} -dbtype prot -out {path} -title {path.name}"
    os.system(cmd)


def blast_on_dir(input_dir, db, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    files = [input_dir / f for f in os.listdir(input_dir)
             if (input_dir / f).exists()]
    for query_file in files:
        run_blast(query_file, db, output_dir)


if __name__ == "__main__":
    # Make the database. 
    path_db = Path("./data/databases/uniprot_sprot.fasta")
    make_blast_db(path_db)

    # Run Blast on Victors and Protegen proteins.
    db = "./data/databases/uniprot_sprot.fasta"
    out = "./data/blast_output/"
    paths = [("./data/proteins/victors/", out + "victors"),
             ("./data/proteins/victors_viruses/", out + "victors_viruses"),
             ("./data/proteins/protegen/bacteria", out + "protegen_bacteria"),
             ("./data/proteins/protegen/viruses", out + "protegen_viruses")]
    for input_dir, output in paths:
        blast_on_dir(Path(input_dir), db, Path(output))
