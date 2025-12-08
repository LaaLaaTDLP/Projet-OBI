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

def best_match(result):
    cut_reads={}
    for read_id, status in result.items():
        if status == "Unassigned":
            sequence = reads[read_id]   # On récupère la séquence associée à cet id
            cut_seq=sequence[0:len(sequence)//2]
            cut_reads[read_id] = cut_seq
    # return cut_reads
    result = exact_match_classifier(cut_reads, ref_seq, ref_genre)
    return result


cut_reads= best_match(result)
print(cut_reads)

# import metagenomic as meta


# reads = meta.load_fasta(filename='HQ_76bp_16SRNA.fa') # reads : dict {read_id: sequence}
# ref_seq = meta.load_fasta(filename='test.fasta') # ref_seq : dict {seq_ID : sequence}
# ref_genre = meta.load_tax_genus(taxofile='gg_99.pds.tax') # ref_genre : {seq_ID: genre_taxo}

# def classify_read(reads_seq, ref_seq, ref_genre):
#     # ref_seq : dict {ref_id: sequence}
#     # ref_genre : dict {ref_id: genre}
#     L = len(reads_seq)
#     for w in range(6, 0, -1):       # 6 -> 5 -> 4 -> 3

#         for it in range(0, L - w + 1):
#             motif = reads_seq[it:it + w]

#             matches = []
#             for ref_id, ref_s in ref_seq.items():
#                 if motif in ref_s:
#                     matches.append(ref_id)
#                     print(matches)
#             if len(matches) == 0:
#                 continue
#             elif len(matches) == 1:
#                 ref_id = matches[0]
#                 return ref_genre.get(ref_id, "Unassigned")
#             else:
#                 ref
#                 # plusieurs matches : tu peux décider de marquer "Ambiguous"
#                 # ou de continuer à glisser / réduire encore
#                 continue

#     return "Unassigned"

# def exact_match_classifier(reads, ref_seq, ref_genre):

#     result = {}
#     for reads_id, reads_seq in reads.items():
#         result[reads_id] = classify_read(reads_seq, ref_seq, ref_genre)
#     return result

# assignations = exact_match_classifier(reads, ref_seq, ref_genre)
# print(assignations)


# exact_match_classifier(reads, ref_seq, ref_genre)