import logging
import os
import pandas as pd
import urllib.error
import urllib.request
import xml.etree.ElementTree

from argparse import ArgumentParser
from Bio import SeqIO

BUILD_DIR = "build"
EVIDENCE = [
    "uncertain",
    "predicted",
    "inferred from homology",
    "evidence at transcript level",
    "evidence at protein level",
]
UNIPROT_URL = "https://www.uniprot.org/uniprot"


def get_evidence_score(evidence):
    """Convert the Uniprot evidence level to a numeric value."""
    if evidence in EVIDENCE:
        return EVIDENCE.index(evidence)
    return 0


def get_protein_details(uniprot_id, retry=True):
    """Get a protein record from UniProt and parse the details to get evidence level, gene name,
    gene synonyms, protein name, and protein synonyms."""
    uniprot_record = os.path.join(BUILD_DIR, "records", uniprot_id + ".xml")
    if not os.path.exists(uniprot_record):
        try:
            contents = urllib.request.urlopen(f"{UNIPROT_URL}/{uniprot_id}.xml").read()
        except (urllib.error.HTTPError, urllib.error.URLError):
            logging.error("Could not retrieve results for " + uniprot_id)
            return None
        with open(uniprot_record, "wb") as f:
            f.write(contents)
    try:
        record = SeqIO.read(open(uniprot_record, "r"), "uniprot-xml")
    except (ValueError, xml.etree.ElementTree.ParseError):
        if retry:
            # Remove record and retry downloading once
            os.remove(uniprot_record)
            return get_protein_details(uniprot_id, retry=False)
        logging.error("Could not read result for " + uniprot_id)
        return None

    return {
        "uniprot id": uniprot_id,
        "evidence": get_evidence_score(
            record.annotations.get("proteinExistence", ["uncertain"])[0]
        ),
        "gene": record.annotations.get("gene_name_primary"),
        "gene synonyms": "|".join(record.annotations.get("gene_name_synonym", [])),
        "protein": "|".join(record.annotations.get("recommendedName_fullName", [])),
        "protein synonyms": "|".join(record.annotations.get("alternativeName_fullName", [])),
    }


def process_multiple_proteins(element_id, output_df, output_row, uniprot_ids):
    """For alignments with more than one UniProt ID, select the 'best' protein based on the UniProt
    evidence level. Add this protein & its details to the output dataframe. If multiple entries
    share the same evidence level, include them all. If an element has aligned to more than one
    sequence (multiple features), take the average 'score' of the evidence levels."""
    ids_to_score = pd.DataFrame(columns=["id", "score"])
    id_to_details = {}
    # Iterate through uniprot IDs and get the protein details from UniProt
    for uniprot_id in uniprot_ids:
        evidence_score = 0
        ids_ok = True
        for uid in uniprot_id.split("|"):
            details = get_protein_details(uid)
            if not details:
                ids_ok = False
                break
            evidence_score += details["evidence"]
            id_to_details[uid] = details
        if not ids_ok:
            logging.error("Could not get all details for " + uniprot_id)
            continue
        evidence_score = evidence_score / len(uniprot_id.split("|"))
        ids_to_score = ids_to_score.append(
            {"id": uniprot_id, "score": evidence_score}, ignore_index=True
        )
    if ids_to_score.empty:
        # Uniprot records could not be returned - often due to an obsolete entry
        output_row["comment"] = "Could not retrieve UniProt protein details for " + ", ".join(
            uniprot_ids
        )
        output_df = output_df.append(output_row, sort=False, ignore_index=True)
        return output_df

    # Sort return proteins by their score and pick the top score
    ids_to_score = ids_to_score.sort_values("score", ascending=False)
    top_score = ids_to_score.iloc[0]["score"]
    top_uniprot_ids = ids_to_score.loc[ids_to_score["score"] == top_score]["id"].unique().tolist()
    if len(top_uniprot_ids) > 1:
        logging.warning("More than one set of top-scoring UniProt proteins for " + element_id)
    i = 0
    for uniprot_id in top_uniprot_ids:
        i += 1
        for uid in uniprot_id.split("|"):
            protein = id_to_details[uid].copy()
            del protein["evidence"]
            output_row.update(protein)
            if len(top_uniprot_ids) > 1:
                output_row["comment"] = f"Group {i}/{len(top_uniprot_ids)} of top-scoring proteins"
            output_df = output_df.append(output_row, sort=False, ignore_index=True)
    return output_df


def process_proteins(element_id, output_df, output_row, uniprot_ids):
    if len(uniprot_ids) == 1:
        logging.debug(element_id + " has one UniProt ID!")
        for uid in uniprot_ids[0].split("|"):
            protein = get_protein_details(uid)
            if not protein:
                output_row["comment"] = "Could not get protein details for " + uid
                output_df = output_df.append(output_row, sort=False, ignore_index=True)
                continue
            del protein["evidence"]
            output_row.update(protein)
            output_df = output_df.append(output_row, sort=False, ignore_index=True)
        return output_df
    # Multiple uniprot IDs - pick the one with the best evidence
    logging.debug(f"Selecting best from {len(uniprot_ids)} UniProt IDs for {element_id}")
    output_df = process_multiple_proteins(element_id, output_df, output_row, uniprot_ids)
    return output_df


def main():
    parser = ArgumentParser()
    parser.add_argument("alignments")
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

    sep = "\t"
    if args.alignments.endswith(".csv"):
        sep = ","
    align_df = pd.read_csv(args.alignments, sep=sep)

    output_df = pd.DataFrame(
        columns=[
            "id",
            "element",
            "species",
            "uniprot id",
            "gene",
            "gene synonyms",
            "protein name",
            "protein synonyms",
            "comment",
        ]
    )
    for element_id in align_df["id"].unique().tolist():
        element_align = align_df.loc[align_df["id"] == element_id]
        if element_align.empty:
            logging.warning("No alignments for " + element_id)
            continue
        output_row = {
            "id": element_id,
            "element": element_align.iloc[0]["element"],
            "species": element_align.iloc[0]["species"],
            "comment": element_align.iloc[0]["comment"],
        }
        if not element_align["aligned"].dropna().unique().tolist():
            # No alignments
            output_df = output_df.append(output_row, sort=False, ignore_index=True)
            continue

        # Check for UniProt IDs provided by Genbank record
        uniprot_ids = element_align["uniprot_ids"].dropna().unique().tolist()
        if uniprot_ids:
            output_df = process_proteins(element_id, output_df, output_row, uniprot_ids)
            continue

        # No Uniprot IDs provided by Genbank record
        # We will query UniProt for the Genbank protein ID
        genbank_ids = element_align["genbank_ids"].dropna().unique().tolist()
        uniprot_ids = []
        for genbank_id in genbank_ids:
            uniprot_id = []
            for gid in genbank_id.split("|"):
                query = f"{UNIPROT_URL}?query={gid}&format=tab&sort=score&columns=id"
                result = urllib.request.urlopen(query).read().decode("utf-8").split("\n")
                if len(result) > 1:
                    uniprot_id.append(result[1].strip())
            if uniprot_id:
                uniprot_ids.append("|".join(uniprot_id))

        if not uniprot_ids:
            # TODO: Maybe BLAST protein sequences from Genbank -> Uniprot and find best match?
            #       This won't work so well for synthetic, which most of the time this is
            logging.warning("Could not retrieve any UniProt IDs for " + element_id)
            output_row["comment"] = "Could not retrieve any UniProt IDs"
            output_df = output_df.append(output_row, sort=False, ignore_index=True)
            continue

        output_df = process_proteins(element_id, output_df, output_row, uniprot_ids)

    sep = "\t"
    if args.output.endswith(".csv"):
        sep = ","
    output_df.to_csv(args.output, sep=sep, index=False)


if __name__ == "__main__":
    main()
