import metagenomic as meta


reads = meta.load_fasta(filename='HQ_76bp_16SRNA.fa')
ref_seq = meta.load_fasta(filename='test.fasta')
ref_genre = meta.load_tax_genus(taxofile='gg_99.pds.tax')

def exact_match_classifier(reads, ref_seq, ref_genre):
    # résultat : {read_id: genre_ou_Unassigned}
    result = {}

    for read_id, seq in reads.items():
        motif = seq[:6]
        genre_inconnu = "Unassigned"

        for ref_id, ref_sequence in ref_seq.items():
            if motif in ref_sequence:
                genre = ref_genre.get(ref_id)
                if genre is None:
                    genre_inconnu = "Unassigned"
                else:
                    genre_inconnu = genre

        result[read_id] = genre_inconnu
    return result

print(exact_match_classifier(reads, ref_seq, ref_genre))
