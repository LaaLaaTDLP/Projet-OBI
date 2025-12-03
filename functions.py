def load_fasta(filename):
    """Reads a FASTA file and returns the sequence as a string."""
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

