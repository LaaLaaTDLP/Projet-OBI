import metagenomic as meta

reads = meta.load_fasta(filename='HQ_76bp_16SRNA.fa')
ref_seq = meta.load_fasta(filename='gg_99.pds.ng.fasta')
ref_genre = meta.load_tax_genus(taxofile='gg_99.pds.tax')

# print('test reads : ',test_reads)
# print('test ref seq ',test_ref_seq)

dico_test_reads = meta.dico_test(reads, 100)
dico_test_ref_seq = meta.dico_test(ref_seq, 500)


result = meta.exact_match_classifier(reads, ref_seq, ref_genre)
print(result)

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

# for reads_ID, reads_seq in reads.items() :
#     for ref_id, seq in ref_seq.items(): 
#         print(f"Résultat de {reads_ID} et {ref_id}: {best_match(reads_seq, seq, w=7)}") #affiche chaque lignes de comparaison.
   

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

        for ref_id, ref in dico_test_ref_seq.items():
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

print(longest_substring_classifier(dico_test_reads, dico_test_ref_seq, ref_genre, 7, result))