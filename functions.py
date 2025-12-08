filename = "test.fasta"
taxofile = "gg_99.pds.tax"

def load_fasta(filename):
    """Lit le fichier fasta et retourne un dictionnaire de seqID et seqString"""
    seq_dict = {}
    with open(filename, "r") as f: #"r" for read #as f = la manière dont on va appeler le fichier 
        for line in f:
            line = line.strip() #les espaces au début et à la fin surtout le saut de ligne \n
            if line.startswith(">"): 
                cle=line[1:]
                seq_dict[cle]=""  #enlève la description de la séquence"
            else: 
                seq_dict[cle]+=line
    f.close()
    return seq_dict

def load_tax_genus(taxofile):
    """lis les taxo et retourne {seq_id: genus}."""
    dico_taxo = {}
    with open(taxofile, "r") as t: #"r" for read #as f = la manière dont on va appeler le fichier 
        for line in t:
            line = line.strip() #les espaces au début et à la fin surtout le saut de ligne \n

            seq_ID , reste = line.split()
            indicateur  = reste.split(";")

            dico_taxo[seq_ID] = None

            for i in indicateur :
                if i.startswith('g'):
                    dico_taxo[seq_ID] = i
                    break

    t.close()
    return dico_taxo

print(load_tax_genus(taxofile))

