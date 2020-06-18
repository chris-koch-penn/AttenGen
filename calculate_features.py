# Chris Koch, 2020.
# Feature vector calculation inspired by Vaxign-ML, by Edison Ong.
from propy.PyPro import GetProDes
import pickle
from multiprocessing import Pool, cpu_count
from pathlib import Path
from Bio import SeqIO


def get_descriptor(accession_no, seq):
    try:
        seq = seq.replace("X", "L")
        Des = GetProDes(seq)
        features = Des.GetAAComp()
        features.update(Des.GetMoreauBrotoAuto())
        features.update(Des.GetGearyAuto())
        features.update(Des.GetCTD())
        features.update(Des.GetQSO())
        return (accession_no, list(features.values()), seq)
    except Exception as e:
        print(e)
        return None


def run_descriptor(input_fasta, output_dir, max_num_samples=5000):
    output_dir.mkdir(parents=True, exist_ok=True)
    vals = SeqIO.parse(input_fasta, "fasta")
    print("done parsing fasta")
    vals = [(f.id, str(f.seq))
            for i, f in enumerate(vals) if i < max_num_samples]
    print(len(vals))
    output = []
    with Pool() as pool:
        output = pool.starmap(get_descriptor, vals)
    output_path = output_dir / (Path(input_fasta).stem + ".features.pkl")
    output = [item for item in output if item != None]
    pickle.dump(output, open(output_path, "wb"))


if __name__ == "__main__":
    output = Path("./feature_vectors")
    base = Path("./data")
    run_descriptor(base / "victors_viruses.faa", output)
    run_descriptor(base / "victors_filtered_uniprot.faa", output)
    run_descriptor(base / "victors_all.faa", output)
    run_descriptor(base / "protegen_filtered_uniprot.faa", output)
    run_descriptor(base / "protegen_bacteria_proteins.faa", output)
    run_descriptor(base / "protegen_virus_proteins.faa", output)
