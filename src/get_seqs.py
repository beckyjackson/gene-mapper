import csv
import logging
import psycopg2
import re

from argparse import ArgumentParser, FileType
from collections import defaultdict


def main():
    parser = ArgumentParser()
    parser.add_argument("user", help="database user name")
    parser.add_argument("seqs", type=FileType("w"))
    parser.add_argument("short_seqs", type=FileType("w"))
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    # This requires running a port forward to the remote postgresql database
    # ssh -N -L 55432:bicoid:5432 [USER]@merlot.lbl.gov &
    fasta = ""
    seq_count = 0
    short_fasta = ""
    short_count = 0
    with psycopg2.connect(host="localhost", port=55432, user=args.user, dbname="felix") as conn:
        cur = conn.cursor()
        # Get the sequences for each element
        cur.execute(
            """SELECT m.element_id, ps.seq_name, ps.sequence FROM parts_sequences ps
            JOIN modifications m ON m.parts_ptr_id = ps.part_id"""
        )
        for res in cur.fetchall():
            fasta_entry = f">{res[0]}|{res[1]}\n{res[2]}\n"
            if len(res[2]) < 50:
                short_fasta += fasta_entry
                short_count += 1
            else:
                fasta += fasta_entry
                seq_count += 1

    logging.info(f"Retrieved {seq_count} regular sequences and {short_count} short sequences")

    args.seqs.write(fasta)
    args.short_seqs.write(short_fasta)


if __name__ == '__main__':
    main()
