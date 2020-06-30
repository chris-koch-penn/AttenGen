# Chris Koch, 2020.
from Bio import SearchIO, SeqIO


def process_blast_output_dir(path):
    hit_set = set()
    try:
        qresult = SearchIO.read(path, 'blast-tab', comments=True)
        print(qresult)
        print(qresult.hits)
        print(len(qresult.hits))
        print(len(qresult.hits.hsps))
        return
        for a in qresult:
            print(a)
            return
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
    path1 = "./data/blast_output/victors_all"
    path2 = "./data/blast_output/protegen_bacteria"
    path3 = "./data/blast_output/protegen_viruses"
    set1 = process_blast_output_dir(path1)
    set2 = process_blast_output_dir(path2) | process_blast_output_dir(path3)
    input1 = "./data/uniprot_sprot.fasta"
    output1 = "./data/filtered/victors_filtered_uniprot.faa"
    output2 = "./data/filtered/protegen_filtered_uniprot.faa"
    filter_fasta(input1, output1, set1)
    filter_fasta(input1, output2, set2)
