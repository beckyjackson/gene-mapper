import logging
import os
import pandas as pd
import subprocess

from argparse import ArgumentParser
from Bio import SeqIO
from collections import defaultdict

# Column names for BLAST output
ALIGN_HEADER = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "staxids",
]

# Build directory to download records to (BUILD/records)
BUILD_DIR = "build"

# NCBI API URL to fetch Genbank records
RECORD_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"


def add_proteins_to_output(output_df, output_row, sseqid, align_segs):
    """Add the protein identifiers and gene names to the output for this alignment."""
    proteins, comment = get_protein_ids(sseqid, align_segs)
    if not proteins:
        # logging.warning(f"For {element_id}: {comment}")
        return output_df, False
    genbank_ids = proteins["genbank"]
    uniprot_ids = proteins["uniprot"]
    if not genbank_ids and not uniprot_ids:
        # We found CDS in the record, but they did not have proteins
        return output_df, False
    output_row["genbank_ids"] = "|".join(genbank_ids)
    output_row["uniprot_ids"] = "|".join(uniprot_ids)
    output_row["genes_automatic"] = "|".join(proteins["genes"])
    output_row["aligned"] = sseqid
    return output_df.append(output_row, sort=False), True


def get_protein_ids(seqid, align_segs):
    """Get the Genbank and UniProt protein IDs for an aligned sequence from Genbank."""
    file = f"{BUILD_DIR}/records/{seqid}.gb"
    record = None
    if os.path.exists(file):
        try:
            record = SeqIO.read(open(file, "r"), "genbank")
        except ValueError:
            # File did not download correctly, we will try to get it again
            pass
    if not record:
        logging.debug("Downloading genbank record for " + seqid)
        query = f"{RECORD_URL}db=nuccore&id={seqid}&rettype=gb&retmode=text"
        cmd = f'curl -Lk "{query}" > {file}'
        success = run_cmd(cmd)
        if not success:
            if os.path.exists(file):
                os.remove(file)
            return None, "Could not download record " + seqid
        if os.stat(file).st_size == 0:
            logging.error("Could not download record for " + seqid)
            os.remove(file)
            return None, "Could not download record " + seqid
        try:
            record = SeqIO.read(open(file, "r"), "genbank")
        except ValueError:
            # File did not download correctly, we cannot parse it
            logging.error("Could not download record for " + seqid)
            os.remove(file)
            return None, "Could not download record " + seqid
    return parse_gb(record, seqid, align_segs)


def parse_gb(record, seqid, align_segs):
    """Parse the Genbank record to get the Genbank and UniProt protein IDs for the aligned sequence.
    """
    feature_pos = defaultdict(list)
    feature_idx = 0
    cds_count = 0

    # First iter through the features to get a map of position to feature(s)
    # Each pos is one nt in the sequence, and probably has one feature but we can handle overlaps
    for feature in record.features:
        if feature.type != "CDS":
            # Skip if not a coding seq, does not have protein ID
            feature_idx += 1
            continue
        cds_count += 1
        start = feature.location.nofuzzy_start
        end = feature.location.nofuzzy_end
        for i in range(start, end + 1):
            if i not in feature_pos:
                feature_pos[i] = list()
            feature_pos[i].append(feature_idx)
        feature_idx += 1

    # This record has no coding seq
    if not feature_pos:
        return None, "No coding sequences in genbank record " + seqid

    # This record has only one coding seq - this is our protein
    if cds_count == 1:
        genbank_ids = set()
        uniprot_ids = set()
        genes = set()
        # Get the index of the feature in the record
        feature_idx = feature_pos[list(feature_pos.keys())[0]][0]
        # Get the feature from the record
        feature = record.features[feature_idx]
        # Add the protein ID and UniProt db_xrefs, if they exist
        genbank_id = feature.qualifiers.get("protein_id")
        if not genbank_id:
            logging.warning(
                f"Genbank record {seqid} feature {feature_idx} does not have a protein_id"
            )
        else:
            genbank_ids.update(genbank_id)
        gene_name = feature.qualifiers.get("gene")
        if gene_name:
            genes.update(gene_name)
        for xref in feature.qualifiers.get("db_xref", []):
            if xref.startswith("UniProt"):
                uniprot_ids.add(xref.split(":")[1])
        return {"genbank": list(genbank_ids), "uniprot": list(uniprot_ids), "genes": genes}, None

    # For each position in the aligned sequence, get the feature at that position
    align_features = set()
    for align_start, align_end in align_segs:
        for i in range(align_start, align_end + 1):
            feature_idxs = feature_pos.get(i)
            if not feature_idxs:
                # No CDS for this part of the sequence
                continue
            align_features.update(feature_idxs)

    if not align_features:
        align_strs = []
        for align_start, align_end in align_segs:
            align_strs.append(f"{align_start}:{align_end}")
        align_str = ", ".join(align_strs)
        return None, f"No coding sequences in sequence alignment in {align_str} in " + seqid

    genbank_ids = set()
    uniprot_ids = set()
    genes = set()
    # For each feature in the aligned sequence, add the protein IDs
    for feature_idx in align_features:
        feature = record.features[feature_idx]
        genbank_id = feature.qualifiers.get("protein_id")
        if not genbank_id:
            logging.error(
                f"Genbank record {seqid} feature {feature_idx} does not have a protein_id"
            )
        else:
            genbank_ids.update(genbank_id)
        gene_name = feature.qualifiers.get("gene")
        if gene_name:
            genes.update(gene_name)
        for xref in feature.qualifiers.get("db_xref", []):
            if xref.startswith("UniProt"):
                uniprot_ids.add(xref.split(":")[1])
    return {"genbank": list(genbank_ids), "uniprot": list(uniprot_ids), "genes": genes}, None


def process_multiseq_alignments(element_id, output_df, output_row, element_align):
    # Filter where the sseqid is the same between all sequence regions that we aligned
    aligned_seqs = element_align[
        ["seqname", "sseqid", "pident", "length", "mismatch", "sstart", "send"]
    ].copy(deep=True)
    aligned_seqs["sstart"] = aligned_seqs["sstart"].apply(to_str)
    aligned_seqs["send"] = aligned_seqs["send"].apply(to_str)
    aligned_seqs = (
        aligned_seqs.groupby("sseqid")
        .agg(
            {
                "seqname": ", ".join,
                "pident": "mean",
                "length": "sum",
                "mismatch": "sum",
                "sstart": "|".join,
                "send": "|".join,
            }
        )
        .reset_index()
    )
    aligned_seqs = aligned_seqs.loc[aligned_seqs["seqname"].str.contains(", ")]
    if aligned_seqs.empty:
        # No overlap between regions
        output_row["comment"] = "No overlapping matches for multiple sequence regions"
        return output_df.append(output_row, sort=False)

    top_alignments = sort_and_filter(aligned_seqs)
    logging.debug(f"Found {len(top_alignments)} top matches")

    good_match = False
    for _, align_row in top_alignments.iterrows():
        # Get the genbank identifier
        sseqid = align_row["sseqid"].split("|")[1]
        align_start = align_row["sstart"].split("|")
        align_end = align_row["send"].split("|")
        # Create a list of tuples representing the aligned segments
        align_segs = []
        for x in range(0, len(align_start)):
            ast = int(align_start[x])
            ae = int(align_end[x])
            if ast > ae:
                # Reverse start & end so start is always less
                ast_copy = ast
                ast = ae
                ae = ast_copy
            align_segs.append((ast, ae))
        output_df, has_update = add_proteins_to_output(output_df, output_row, sseqid, align_segs)
        if has_update:
            good_match = True

    if not good_match:
        logging.warning(f"No protein IDs found in alignments for {element_id}")
        output_row["comment"] = "No protein IDs found in alignments"
        output_df = output_df.append(output_row, sort=False)

    return output_df


def process_singleseq_alignments(element_id, output_df, output_row, element_align):
    top_alignments = sort_and_filter(element_align)
    logging.debug(f"Found {len(top_alignments)} top matches")
    good_match = False
    for _, align_row in top_alignments.iterrows():
        # TODO: determine if this is a good match by finding protein acc in Genbank summary
        #       we may end up iterating through more than once (if there is more than one)
        sseqid = align_row["sseqid"].split("|")[1]
        align_start = align_row["sstart"]
        align_end = align_row["send"]
        if align_start > align_end:
            # Reverse start and end
            align_start_copy = align_start
            align_start = align_end
            align_end = align_start_copy
        output_df, has_update = add_proteins_to_output(
            output_df, output_row, sseqid, [(align_start, align_end)]
        )
        if has_update:
            good_match = True

    if not good_match:
        logging.warning(f"No protein IDs found in alignments for {element_id}")
        output_row["comment"] = "No protein IDs found in alignments"
        output_df = output_df.append(output_row, sort=False)
    return output_df


def to_int(x):
    """Convert a value to int."""
    if pd.isna(x):
        return None
    return int(x)


def to_str(x):
    """Convert a value to str."""
    if pd.isna(x):
        return None
    return str(x)


def run_cmd(cmd):
    """Run a command. Return True on success."""
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
    if rc > 0:
        logging.error("unable to run command: " + cmd)
        logging.error("cause: " + str(err))
        return False
    return True


def sort_and_filter(df):
    """Sort the alignments for the top % identity, longest length, and lowest mismatch.
    Then, filter for only those that have these same numbers."""
    df = (
        df.sort_values("mismatch")
        .sort_values("length", ascending=False)
        .sort_values("pident", ascending=False)
    )
    top_pident = df.iloc[0]["pident"]
    top_length = df.iloc[0]["length"]
    top_mismatch = df.iloc[0]["mismatch"]
    df = df.loc[df["pident"] >= top_pident]
    df = df.loc[df["length"] >= top_length]
    return df.loc[df["mismatch"] <= top_mismatch]


def main():
    parser = ArgumentParser()
    parser.add_argument("elements")
    parser.add_argument("seq_align")
    parser.add_argument("short_seq_align")
    parser.add_argument("output")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-vv", "--very-verbose", action="store_true")
    args = parser.parse_args()

    # Set up logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format="%(levelname)s: %(message)s")
    elif args.very_verbose:
        logging.basicConfig(level=logging.DEBUG, format="%(levelname)s: %(message)s")
    else:
        logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    elements_df = pd.read_csv(args.elements, sep="\t")
    elements_df["species id"] = elements_df["species id"].apply(to_int)

    # Combine regular & short alignments
    align_df = pd.read_csv(args.seq_align, sep="\t", header=None, names=ALIGN_HEADER)
    short_align_df = pd.read_csv(args.short_seq_align, sep="\t", header=None, names=ALIGN_HEADER)
    align_df = align_df.append(short_align_df, sort=False)

    # Sometimes entries have more than one tax ID
    # Split these into distinct rows, then apply int to all tax IDs
    align_df["staxids"] = align_df["staxids"].apply(to_str)
    align_df["staxids"] = align_df["staxids"].str.split(";")
    align_df = (
        align_df.set_index([x for x in ALIGN_HEADER if x != "staxids"])["staxids"]
        .apply(pd.Series)
        .stack()
        .reset_index()
        .drop("level_10", axis=1)
        .rename(columns={0: "staxids"})
    )
    align_df["staxids"] = align_df["staxids"].apply(to_int)
    align_df[["element", "seqname"]] = align_df["qseqid"].str.split("|", expand=True)

    output_df = pd.DataFrame()
    for _, row in elements_df.iterrows():
        element_id = row["id"]
        species_id = row["species id"]
        output_row = row.copy(deep=True)

        if not species_id or pd.isna(species_id):
            output_row["comment"] = "No species ID"
            output_df = output_df.append(output_row, sort=False)
            continue
        species_id = int(species_id)

        # Filter the alignments for this element & its species
        # Checking the sequence names that remain to make sure we don't lose any
        element_align = align_df.loc[align_df["element"] == element_id]
        seqnames = element_align["seqname"].unique().tolist()
        element_align = element_align.loc[align_df["staxids"].astype(int) == species_id]
        if element_align.empty:
            output_row["comment"] = "No hits for species"
            output_df = output_df.append(output_row, sort=False)
            continue

        seqnames_2 = element_align["seqname"].unique().tolist()
        if len(seqnames_2) < len(seqnames):
            output_row["comment"] = "Not all sequence parts have hits in species"
            output_df = output_df.append(output_row, sort=False)
            continue

        if len(seqnames_2) > 1:
            # Database provided multiple sequences for this one element
            logging.debug("Handling multiple sequence alignments for " + element_id)
            output_df = process_multiseq_alignments(
                element_id, output_df, output_row, element_align
            )
            continue

        logging.info("Handling one sequence alignment for " + element_id)
        output_df = process_singleseq_alignments(element_id, output_df, output_row, element_align)

    output_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
