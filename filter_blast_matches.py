# Chris Koch, 2020.
from pathlib import Path
from Bio import SearchIO, SeqIO
import os


def process_blast_output_dir(path):
    files = [path / f for f in os.listdir(path)
             if (path / f).exists()]
    hit_set = set()
    for query_file in files:
        try:
            qresult = SearchIO.read(query_file, 'blast-tab', comments=True)
            hit_set.update(
                [q.id for q in qresult.hits if q.hsps[0].ident_pct > 30])
        except Exception as e:
            pass
    return hit_set


def filter_fasta(input_file, output_file, exclusion_set):
    records = (r for r in SeqIO.parse(input_file, "fasta")
               if r.id not in exclusion_set)
    SeqIO.write(records, output_file, "fasta")


if __name__ == "__main__":
    path1 = Path("../data/proteins/blast_output/victors")
    path2 = Path("../data/proteins/blast_output/victors_viruses")
    path3 = Path("../data/proteins/blast_output/protegen_bacteria")
    path4 = Path("../data/proteins/blast_output/protegen_viruses")
    set1 = process_blast_output_dir(path1) | process_blast_output_dir(path2)
    set2 = process_blast_output_dir(path3) | process_blast_output_dir(path4)
    base = "../data/proteins/databases/"
    input1 = base + "uniprot_sprot.fasta"
    output1 = base + "victors_filtered_uniprot.faa"
    output2 = base + "protegen_filtered_uniprot.faa"
    filter_fasta(input1, output1, set1)
    filter_fasta(input1, output2, set2)
