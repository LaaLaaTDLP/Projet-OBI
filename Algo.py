import metagenomic as meta

filename = "test.fasta"
taxofile = "gg_99.pds.tax"

meta.load_fasta(filename)
meta.load_tax_genus(taxofile)