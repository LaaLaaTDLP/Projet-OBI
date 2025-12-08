# Projet-OBI

Un fichier python nommé metagenomic avec toutes les fonctions. 
metagenomic.py est appeler après par import 

Rapport de 2 pages max
sur les différentes fonctions qu'on a fait 
préciser quelle version de python on a utilisée

### Some notions for Metagenomic Taxonomic Profiling
Metagenomic taxonomic profiling is the process of determining which microorganisms are present in an environmental sample (microbiome) and how abundant each is. In 16S rRNA studies, short DNA fragments that encode a conserved ribosomal gene are sequenced; because the gene contains both highly conserved and variable regions, the sequences are used as barcodes for bacterial taxonomy. By comparing these barcodes to a reference database such as Greengenes (exact match or partial alignment) we can assign each read to a taxonomic rank (e.g., genus). This project focuses on creating a reference-based classifier using a curated 16S database, specifically the Greengene database. The objective is to classify DNA reads based on the provided reference sequences and their corresponding taxonomic information.

### Input Data
gg_99.pds.ng.fasta Greengene database reference sequences
gg_99.pds.tax : Greengene database taxonomy. The sequence ID and taxonomy fields are separated by a
tab character. The taxonomy information field uses ‘;‘ as separator. ‘g__‘ corresponds to ‘genus’ information.
HQ_76bp_16SRNA.fa : synthetic high-quality 16S short reads (76bp)
SRR36139653.fasta : real dataset, gut microbiome 16S long reads (250pb)

