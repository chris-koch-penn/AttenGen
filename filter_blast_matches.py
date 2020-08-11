# Chris Koch, 2020.
from Bio import SearchIO, SeqIO
from io import StringIO
from pathlib import Path

PCT_IDENTITY_THRESHOLD = 30


def process_blast_output(query_results):
    hit_set = set()
    for query_res in query_results:
        try:
            f = StringIO(query_res)
            qresult = SearchIO.read(f, 'blast-tab', comments=True)
            hit_set.update(
                [hit.id for hit in qresult.hits if hit.hsps[0].ident_pct > PCT_IDENTITY_THRESHOLD])
        except Exception as e:
            print(e)
            pass
    return hit_set


def filter_fasta(input_file, output_file, exclusion_set):
    records = (r for r in SeqIO.parse(input_file, "fasta")
               if r.id not in exclusion_set)
    SeqIO.write(records, output_file, "fasta")


def extract_query_results(path):
    with open(path, "r") as f:
        text = f.read()
        sep = "# BLASTP"
        queries = [sep + query for query in text.split(sep)]
        return queries[1:]


if __name__ == "__main__":
    query_results1 = extract_query_results("./data/blast_output/victors_all")
    query_results2 = extract_query_results(
        "./data/blast_output/protegen_bacteria_proteins")
    query_results3 = extract_query_results(
        "./data/blast_output/protegen_virus_proteins")
    set1 = process_blast_output(query_results1)
    set2 = process_blast_output(
        query_results2) | process_blast_output(query_results3)
    input1 = "./data/uniprot_sprot.fasta"
    output1 = "./data/filtered/victors_filtered_uniprot.faa"
    output2 = "./data/filtered/protegen_filtered_uniprot.faa"
    Path("./data/filtered").mkdir(parents=True, exist_ok=True)
    filter_fasta(input1, output1, set1)
    filter_fasta(input1, output2, set2)
