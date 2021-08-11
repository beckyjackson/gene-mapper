## Synbio Element-Gene Mapper

This process attempts to map elements to their parent genes and provide UniProt identifiers for the proteins encoded by these genes. This is a work in progress.

This process requires a little bit of setup... follow the steps below to get started.

1. Install the Python requirements (a venv is not required):
```
$ python3 -m venv _venv
$ source _venv/bin/activate
$ python3 -m pip install -r requirements.txt
```
2. Download the elements table [here](https://docs.google.com/spreadsheets/d/1Ll3o7hbQywoHYy4OWW1qV_K10o9A3woUnLN-1AIuNiE/edit#gid=2022655928) and save this as `elements.tsv` in this directory.
3. Save the LinkML model that maps species names to NCBITaxonomy IDs as `linkml.yaml` in this directory.
4. Download the `nt` database from NCBI (warning - this takes a long time and may time out! You can just restart the process by running the same command again):
```
make download_nt
```
5. Begin a port forward to the remote PostgreSQL database (replace `[USER]` with your username)
```
ssh -N -L 55432:bicoid:5432 [USER]@merlot.lbl.gov &
```

Now you're ready to go! You can build the element-to-protein mapping file (`build/element_protein.tsv`) by just running the following, replacing `[USER]` with the *same* username you used in the last step:
```
make USER=[USER]
```

Note that the first time you run this will take a bit longer, as the task to create `build/elements_mapped.tsv` has to download sequence records from NCBI. Once the records are downloaded, they are stored for quicker re-runs in the future (`build/records`). The task for `build/element_protein.tsv` also downloads UniProt protein records and stores them in `build/records`.

The intermediate file, `build/elements_mapped.tsv`, contains all "best" nucleotide sequence alignments. These are the alignments with the highest percent identity, longest length, and lowest number of mismatches. Many elements will have more than one "best" alignment, each represented by a distinct record number. Most of these alignments will also have one or more Genbank IDs. These represent the proteins from the coding sequences of the aligned regions. Sometimes, especially when an element has more than one sequence, the aligned region will have more than one coding sequence and therefore more than one Genbank protein ID. In some cases, there are also UniProt proteins IDs provided by NCBI.

These "best" alignments are used to select the "best" UniProt protein in `build/element_protein.tsv`. If NCBI has provided UniProt proteins, we can just retrieve these records. If there is more than one option, we select the protein with the highest level of evidence. In cases where a single alignment has more than one proteins, an average score of the evidence level is used to select the "best" protein. Sometimes, UniProt will have more than one "best" protein - in which case, the output will include all the "best" matches with a comment about this case.
