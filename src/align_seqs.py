import csv
import logging
import os
import subprocess

from argparse import ArgumentParser, FileType

BUILD_DIR = "build"
OUTFMT = '"6 sseqid pident length mismatch gapopen qstart qend sstart send staxids sscinames"'

# Compound might be other way around for search (multiple genes in the given sequence)
# Flanked deletions could be combined into one sequence with Ns where the deletions happen?


def run_blast(database, output, short=False):
    if os.path.exists(output) and os.stat(output).st_size > 0:
        return True
    blast = f"./ncbi-blast/bin/blastn"
    if short:
        blast += " -task short"
    cmd = f"{blast} -db {database} -query query.fasta -outfmt {OUTFMT} -out {output}"
    os.chdir(BUILD_DIR)
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
        data, err = p.communicate()
        rc = p.returncode
    except subprocess.CalledProcessError as e:
        rc = 1
        err = e.output
    except ValueError as e:
        rc = 1
        err = str(e)
    finally:
        os.chdir("..")
    if rc > 0:
        logging.error("unable to run command: " + cmd)
        logging.error("cause: " + str(err))
        return False
    return True


def main():
    parser = ArgumentParser()
    parser.add_argument("sequences", type=FileType("r"))
    parser.add_argument("database")
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    reader = csv.DictReader(args.sequences, delimiter="\t")
    seqs = {}
    short_seqs = {}
    for row in reader:
        seq = row["sequence"]
        if len(seq) <= 50:
            short_seqs[row["seq_name"]] = seq
        else:
            seqs[row["seq_name"]] = seq

    missed = []

    logging.info(f"{len(seqs)} regular sequences found")
    for seq_name, seq in seqs.items():
        with open(f"{BUILD_DIR}/query.fasta", "w") as fq:
            fq.write(seq)
        logging.info("Running blast for " + seq_name)
        success = run_blast(args.database, f"results/{seq_name}.tsv")
        if not success:
            missed.append(seq_name)

    logging.info(f"{len(seqs)} short sequences found")
    for seq_name, seq in seqs.items():
        with open(f"{BUILD_DIR}/query.fasta", "w") as fq:
            fq.write(seq)
        logging.info("Running blast for " + seq_name)
        success = run_blast(args.database, f"results/{seq_name}.tsv", short=True)
        if not success:
            missed.append(seq_name)

    if missed:
        logging.error(f"Missed {len(missed)} sequences:")
        for m in missed:
            logging.error(m)


if __name__ == "__main__":
    main()
