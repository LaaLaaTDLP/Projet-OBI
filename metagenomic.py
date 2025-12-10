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

def best_match(reads_seq, ref_seq, w): 
    """La fonction best_match prend en paramètre les séquences a comparer, 
    les séquences de références et la taille des morceaux de séquences a comparer."""
    # ref_seq : dict {ref_id: sequence}
    # ref_genre : dict {ref_id: genre}
    
    #reads_seq = reads.items()
    #reads_seq = reads.values()
    resultat=0 #longest  match
    L = len(reads_seq)
    l=len(ref_seq)
    i=0
    for i in range(L - w + 1,):  
        array = reads_seq[i:i + w]
        
        if array in ref_seq:
            for j in range(l-w+1):
                if ref_seq[j:j + w]==array:
                    res_temp = w
                    j_read = i+w
                    j_ref= j+w
                    while(j_read<L and j_ref<l and reads_seq[j_read] == ref_seq[j_ref] ):
                        res_temp += 1
                        j_read += 1
                        j_ref += 1


                        if res_temp> resultat:
                            resultat = res_temp 
    return resultat


def longest_substring_classifier(reads, ref_seq, ref_genre, w, result):
    """la fonction longest_substring_classifier prend comme paramètre :
    - reads : dict {read_id: séquence à comparer}
    - ref_seq : dict {ref_id: séquence de référence}
    - ref_genre : dict {ref_id: genre taxonomique}
    - w : longueur minimale du match
    - result : dict {read_id: genre ou 'Unassigned'} (résultat de l'exact match)
    """
    classification = {}

    for read_id, genre in result.items():
        # si déjà classé par exact_match_classifier, on recopie
        if genre != "Unassigned":
            classification[read_id] = genre
            continue

        read = reads[read_id]
        best_len = 0
        best_ref = None

        for ref_id, ref in ref_seq.items():
            L = best_match(read, ref, w)
            # print(f"Résultat de {read_id} et {ref_id}: {L}")
            if L > best_len:
                best_len = L
                best_ref = ref_id

        if best_len >= w and best_ref is not None:
            classification[read_id] = ref_genre.get(best_ref, "Unassigned")
        else:
            classification[read_id] = "Unassigned"

    return classification