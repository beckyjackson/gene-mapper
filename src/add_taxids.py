import csv
import yaml

from argparse import ArgumentParser, FileType

# Not mapped by LinkML
MANUAL = {
    "Klebsiella.pneumoniae": "573",
    "Lentivirus.human-immunodeficiency-virus1": "11676",
    "Nepovirus.Tobacco-ringspot-virus": "12282"
}


def main():
    parser = ArgumentParser()
    parser.add_argument("elements", type=FileType("r"))
    parser.add_argument("linkml", type=FileType("r"))
    parser.add_argument("output", type=FileType("w"))
    args = parser.parse_args()

    # Create a map of species value -> NCBITaxon ID
    data = yaml.safe_load(args.linkml)
    species_ids = {}
    for k, v in data["enums"]["species_enum"]["permissible_values"].items():
        if "meaning" in v:
            species_ids[k] = v["meaning"].split(":")[1]

    # Read in elements table and add species ID where we can
    reader = csv.DictReader(args.elements, delimiter="\t")
    fields = list(reader.fieldnames)
    fields.append("species id")
    out_rows = []
    for row in reader:
        row["species id"] = species_ids.get(row["species"])
        out_rows.append(row)

    writer = csv.DictWriter(args.output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(out_rows)


if __name__ == '__main__':
    main()
