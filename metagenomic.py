def load_fasta(filename):
    """la fonction filename prend en paramètre un fichier fasta.\n
    Elle lit le fichier fasta puis retourne un dictionnaire de seqID et seqString"""
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
    """La fonction load_tax_genus prend en paramètre un fichier tax.\n
     Elle lis les taxons et retourne son genre {seq_id: genus}."""
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


def dico_test(dico, n):
    """La fonction dico_test prend en paramètre un dictionnaire et un entier correspondant a la taille attendu\n
    la fonction permet de retourner un dictionnaire de moindre taille."""
    d = {}
    count = 0
    for ID, seq in dico.items():
        d[ID] = seq
        count += 1
        if count >= n:
            break
    return d


def exact_match_classifier(reads, ref_seq, ref_genre):
    """La fonction exact_match_classifier, prend en paramètre un dictionnaire """
    # résultat : {read_id: genre_ou_Unassigned}
    result = {}

    for reads_id, reads_seq in reads.items():
        genre_inconnu = "Unassigned"
        for ref_id, ref_sequence in ref_seq.items():
            if reads_seq in ref_sequence:
                genre = ref_genre.get(ref_id)
                if genre is None:
                    genre_inconnu = "Unassigned"
                else:
                    genre_inconnu = genre

        result[reads_id] = genre_inconnu
    return result