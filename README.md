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

Now you're ready to go! You can build the mapped elements file by just running the following, replacing `[USER]` with the *same* username you used in the last step:
```
make USER=[USER]
```
