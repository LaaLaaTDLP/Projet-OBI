import metagenomic as meta

reads = meta.load_fasta(filename='HQ_76bp_16SRNA.fa')
ref_seq = meta.load_fasta(filename='gg_99.pds.ng.fasta')
ref_genre = meta.load_tax_genus(taxofile='gg_99.pds.tax')

def exact_match_classifier(reads, ref_seq, ref_genre):
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

result = exact_match_classifier(reads, ref_seq, ref_genre)


def best_match(reads_seq, ref_seq,w):
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



for reads_ID, reads_seq in reads.items() :
    for ref_id, seq in ref_seq.items(): 
        print(f"Résultat de {reads_ID} et {ref_id}: {best_match(reads_seq, seq, w=4)}")
   

