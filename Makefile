### Configuration
#
# These are standard options to make Make sane:
# <http://clarkgrubb.com/makefile-style-guide#toc2>

MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.SUFFIXES:
.SECONDARY:


UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	BLAST_URL := https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ncbi-blast-2.10.0+-x64-macosx.tar.gz
else
	BLAST_URL := https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ncbi-blast-2.10.0+-x64-linux.tar.gz
endif

all: build/element_protein.tsv


### SET UP

build build/records:
	mkdir -p $@

build/ncbi-blast.tar.gz: | build
	cd build && wget -O $(notdir $@) $(BLAST_URL)

build/ncbi-blast/bin: build/ncbi-blast.tar.gz
	mkdir -p build/ncbi-blast
	cd build && tar -xzf $(notdir $<) -C ncbi-blast --strip-components 1

# Download the RefSeq prokaryotic genomes
.PHONY: download_refseq
download_refseq: build/ncbi-blast/bin/update_blastdb.pl
	cd build && ../$< --decompress ref_prok_rep_genomes

# Download the NCBI nucleotide sequence database (excludes RefSeq)
.PHONY: download_nt
download_nt: build/ncbi-blast/bin/update_blastdb.pl
	cd build && ../$< --decompress nt

# Add species IDs from LinkML model to elements table
# TODO: get this from the database & automatically generate YAML
# Right now using the manual version from https://docs.google.com/spreadsheets/d/1Ll3o7hbQywoHYy4OWW1qV_K10o9A3woUnLN-1AIuNiE
# ... which includes some gene names for reference
build/elements.tsv: src/add_taxids.py elements.tsv linkml.yaml | build
	python3 $^ $@


### SEQUENCE ALIGNMENT

# This requires running a port forward to the remote postgresql database
# ssh -N -L 55432:bicoid:5432 bjackson@merlot.lbl.gov &
build/sequences.fasta build/short_sequences.fasta: src/get_seqs.py | build
	python3 $< $(USER) build/sequences.fasta build/short_sequences.fasta

# Filter the taxa to align against to just those in our elements file
# This should speed up the BLAST alignment
build/txids.txids: | build
	cut -f7 build/elements.tsv | sort | uniq | sed '1,1d' | sed '$d' > $@

# Align sequences > 50 nt
build/seq_align.tsv: build/sequences.fasta build/txids.txids | build/ncbi-blast/bin
	cd build && ./ncbi-blast/bin/blastn \
	-db nt \
	-query $(notdir $<) \
	-taxidlist $(notdir $(word 2,$^)) \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send staxids" \
	-num_threads 12 \
	-out $(notdir $@)

# Align sequences < 50 nt
build/short_seq_align.tsv: build/short_sequences.fasta build/txids.txids | build/ncbi-blast/bin
	cd build && ./ncbi-blast/bin/blastn \
	-db nt \
	-task blastn-short \
	-query $(notdir $<) \
	-taxidlist $(notdir $(word 2,$^)) \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send staxids" \
	-num_threads 12 \
	-out $(notdir $@)

.PHONY: seq_align
seq_align: build/short_seq_align.tsv build/seq_align.tsv


### RESULT PROCESSING

# All "top" mappings with their Genbank & Uniprot protein IDs
build/elements_mapped.tsv: src/process_alignments.py build/elements.tsv build/seq_align.tsv build/short_seq_align.tsv | build/records
	python3 $^ $@

# Only the "best" UniProt IDs
build/element_protein.tsv: src/select_protein.py build/elements_mapped.tsv
	python3 $^ $@
