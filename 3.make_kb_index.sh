#!/bin/bash

## create KB indexes suitable for expression quantification

mkdir /data/home/Stephen/oocyte_atlas/indexes/Homo_sapiens && cd $_

wget http://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
wget http://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.111.gtf

rm Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.111.gtf