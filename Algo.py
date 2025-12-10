import metagenomic as meta

def menu():
    '''Menu qui permet la configuration des fichiers et le lancement des différentes fonctions.'''
    # valeurs par défaut
    reads = meta.load_fasta('HQ_76bp_16SRNA.fa')
    ref_seq = meta.load_fasta('gg_99.pds.ng.fasta')
    ref_genre = meta.load_tax_genus('gg_99.pds.tax')
    #valeurs des dico de test | pour modifier la taille du test, changez les nombres en paramètre
    dico_test_reads = meta.dico_test(reads, 100)
    dico_test_ref_seq = meta.dico_test(ref_seq, 500)
    #Taille de la fenètre pour la fonction best_match
    w = 7
    #1er lancement de la fonction exact_match, on récup le dico de résultat pour les autres fonctions.
    result = meta.exact_match_classifier(reads, ref_seq, ref_genre)

    while True:
        print("\n================= Metagenomic =====================")
        print("1. Changer les fichiers de référence")
        print("2. Changer le fichier échantillon à comparer")
        print("3. Lancer les fonctions (sous-menu)")
        print("0. Quitter")
        rep = int(input("Votre choix : "))

        if rep == 1:
            # fichiers de référence
            filename = input("Chemin du fichier de séquences de référence (.fasta) : ")
            ref_seq = meta.load_fasta(filename)
            filename = input("Chemin du fichier de taxonomie (.tax) : ")
            ref_genre = meta.load_tax_genus(filename)
            # calculer exact match avec les nouvelles refs
            result = meta.exact_match_classifier(reads, ref_seq, ref_genre)

        elif rep == 2:
            # fichier échantillon (reads)
            filename = input("Chemin du fichier de l'échantillon (.fasta) : ")
            reads = meta.load_fasta(filename)
            # recalculer l'exact match avec les nouveaux reads
            result = meta.exact_match_classifier(reads, ref_seq, ref_genre)

        elif rep == 3:
            # sous-menu des fonctions
            print("\n--- Lancer un fonction ---")
            print("1. Exact match classifier")
            print("2. Longest substring classifier avec dico test (moins long)")
            print("3. Longest substring classifier (long)")
            print("0. Retour")
            choix = int(input("Votre choix : "))

            if choix == 1:
                result = meta.exact_match_classifier(reads, ref_seq, ref_genre)
                print("Résultat exact_match_classifier :")
                print(result)

            elif choix == 2:
                w = int(input("Choisissez une longueur minimale d'alignement w : "))
                result_longest = meta.longest_substring_classifier(dico_test_reads, dico_test_ref_seq, ref_genre, w, result)
                print("Résultat longest_substring_classifier :")
                print(result_longest)

            elif choix == 3:
                w = int(input("Choisissez une longueur minimale d'alignement w : "))
                result_longest = meta.longest_substring_classifier(reads, ref_seq, ref_genre, w, result)
                print("Résultat longest_substring_classifier :")
                print(result_longest)

            elif choix == 0:
                continue

        elif rep == 0:
            break

    return

menu()