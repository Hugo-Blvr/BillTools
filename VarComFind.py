#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Hugo
"""
import os
import sys
import venn 


def find_fichiers_vcf(dossier, listVcF):
    fichiers_vcf = []
    for racine, dossiers, fichiers in os.walk(dossier):
        for fichier in fichiers:
            if fichier.endswith('.vcf') and any(motif + '.' in fichier for motif in listVcF):
                chemin_complet = os.path.join(racine, fichier)
                fichiers_vcf.append(chemin_complet)
    fichiers_vcf = sorted(fichiers_vcf, key=lambda s: s.split('/')[-1])

    return fichiers_vcf


def generer_combinaisons(liste_elements, index=0, combinaison_actuelle=None, dictionnaire_combinaisons=None):
    if combinaison_actuelle is None:
        combinaison_actuelle = []
    if dictionnaire_combinaisons is None:
        dictionnaire_combinaisons = {}

    # Condition d'arrêt : si l'index est à la fin de la liste
    if index == len(liste_elements):
        if len(combinaison_actuelle) > 1:  # On exclut les combinaisons d'un seul élément
            cle = "_".join(combinaison_actuelle)
            dictionnaire_combinaisons[cle] = 0
        return dictionnaire_combinaisons

    # Ajouter l'élément courant à la combinaison actuelle
    generer_combinaisons(liste_elements, index + 1, combinaison_actuelle + [liste_elements[index]], dictionnaire_combinaisons)
    # Ne pas ajouter l'élément courant et passer au suivant
    generer_combinaisons(liste_elements, index + 1, combinaison_actuelle, dictionnaire_combinaisons)

    return dictionnaire_combinaisons


def recup_infos(input_files,legende=None):
    dicoPos_seq = {}
    dicoPos_noseq = {}
    nb_var = {}
    lrep = []
    for vcf_replicat in input_files:  # Boucle sur les fichiers vcf de l'échantillon
        id_rep = vcf_replicat
        rep = id_rep
        lrep.append(rep.split('/')[-1].split(".")[0])
        nbvar = 0
        with open(vcf_replicat, 'r') as vcf:
            for ligne in vcf:
                if not ligne.startswith('#'):
                    nbvar += 1
                    infos = ligne.strip().split("\t")
                    for i in infos[7].split(";"):
                        if i.startswith("SVLEN"):
                            _, long = i.split('=')
                        if i.startswith("SVTYPE"):
                            _, svtype = i.split('=')

                    pos = infos[1]
                    seqAlt = infos[4]

                    if ">" not in seqAlt:
                        if pos in dicoPos_seq:
                            dicoPos_seq[pos][f"{infos[2]}*_*{id_rep}"] = [id_rep, seqAlt, svtype, long]
                        else:
                            dicoPos_seq[pos] = {f"{infos[2]}*_*{id_rep}": [id_rep, seqAlt, svtype, long]}
                    else:
                        if pos in dicoPos_noseq:
                            dicoPos_noseq[pos][f"{infos[2]}*_*{id_rep}"] = [id_rep, svtype, long]
                        else:
                            dicoPos_noseq[pos] = {f"{infos[2]}*_*{id_rep}": [id_rep, svtype, long]}
            nb_var[rep] = nbvar


    if legende:
        dico_combin = generer_combinaisons(legende)
    else:
        dico_combin = generer_combinaisons(lrep)


    for dico in [dicoPos_seq, dicoPos_noseq]:
        soloPos = [pos for pos, var in dico.items() if len(var) == 1]
        for pos in soloPos:
            del dico[pos]

    return dicoPos_seq, dicoPos_noseq, dico_combin,nb_var


def calcul_id(seq1, seq2):
    # Identifie la séquence la plus longue et la plus courte
    longue, courte = seq1, seq2
    if len(seq1) < len(seq2):
        longue, courte = seq2, seq1

    #  GLISSEMENT DE LA SEQUENCE LA PLUS COURTE SUR LA PLUS LONGUE + CALCUL DU MAX DE SIMILARITE
    for i in range(len(longue)):
        similarite = sum(n1 == n2 for n1, n2 in zip(courte, longue[i: i + len(courte)]))
        #  ↑ calcule la similarité entre la courte le long de la longue avec un décalage de 1 par itération
        similarite /= len(longue)  # Normalisation par la longueur de la + longue séquence
        return similarite


def calcul_idNoSeq(nom, nom2):
    if nom[1] == nom2[1] and nom[2] == nom[2]:
        return 1
    return 0


def identifier_variants_communs(pos_variants, id_seq):
    # Étape 1: Organiser les séquences par réplicat
    dic_rep = {}
    for variant, details in pos_variants.items():
        dic_rep.setdefault(details[0], []).append((variant, details[1]))


    # Étape 2: Comparer les variants
    pairs = []
    traite = set()
    for replicat, variant_pairs in dic_rep.items():
        for other_replicat, other_variant_pairs in dic_rep.items():
            if replicat != other_replicat:
                for nom, seq in variant_pairs:
                    for nom2, seq2 in other_variant_pairs:
                        compname = tuple(sorted((nom, nom2)))
                        if compname not in traite:
                            if '<' not in seq:
                                idseq = calcul_id(seq, seq2)
                            else:
                                idseq = calcul_idNoSeq(pos_variants[nom], pos_variants[nom2])
                            if idseq >= id_seq:
                                pairs.append((compname, idseq))
                                traite.add(compname)
    # Trier les paires par score de similarité décroissant
    pairs.sort(key=lambda x: x[1], reverse=True)
    # Étape 3: Former des groupes de variants communs
    groupes = []
    variants_traites = set()
    for (variant1, variant2), score in pairs:
        replicatVar1 = pos_variants[variant1][0]
        replicatVar2 = pos_variants[variant2][0]

        if variant1 in variants_traites and variant2 in variants_traites:
            continue

        find_groupeVariant1 = next((g for g in groupes if variant1 in g), None)
        find_groupeVariant2 = next((g for g in groupes if variant2 in g), None)

        if find_groupeVariant1 is not None and find_groupeVariant2 is None:
            if replicatVar2 not in {pos_variants[variant][0] for variant in find_groupeVariant1}:
                find_groupeVariant1.add(variant2)
                variants_traites.add(variant2)
        elif find_groupeVariant2 is not None and find_groupeVariant1 is None:
            if replicatVar1 not in {pos_variants[variant][0] for variant in find_groupeVariant2}:
                find_groupeVariant2.add(variant1)
                variants_traites.add(variant1)
        elif find_groupeVariant1 is None and find_groupeVariant2 is None:
            nouveau_groupe = {variant1, variant2}
            groupes.append(nouveau_groupe)
            variants_traites.add(variant1)
            variants_traites.add(variant2)

    varcom=[]
    for groupe in groupes:
        varcom.append(list(groupe))

    return varcom


def diagVenn(dossier,inputFiles,outname,legende,id_seq):
    if outname == 'false':outname=None
    if legende == 'false':legende=None
    
    repfiles = find_fichiers_vcf(dossier, inputFiles)
    files = [f.split('/')[-1] for f in repfiles]
    print(f"{len(files)} fichiers trouver : {files}\n")
    if len(repfiles) < 2:
        print("Nombre de fichier trouver inférieure à 2")
        return 1
    elif len(repfiles) > 5:
        print("Le nombre de fichiers d'entrée maximum est de 5.")
        return 1
    dicoPos_seq, dicoPos_noseq, dico_combin, nb_var = recup_infos(repfiles,legende)
    listVar = []

    for pos in dicoPos_seq:
        vcs = identifier_variants_communs(dicoPos_seq[pos],id_seq)
        for vc in vcs:
            listVar.append(sorted(vc, key=lambda s: s.split('/')[-1]))
    for pos in dicoPos_noseq:
        vcs = identifier_variants_communs(dicoPos_noseq[pos], id_seq)
        for vc in vcs:
            listVar.append(sorted(vc, key=lambda s: s.split('/')[-1]))


    for i in range(len(listVar)):
        for j in range(len(listVar[i])):
            if legende:
                if len(set(legende)) == len(files):
                    listVar[i][j] = legende[files.index(listVar[i][j].split('/')[-1])]
                else:
                    print(f"Le nombre de légende doit être égale au fichier trouver ({len(files)}) : {files}")
                    return 1
            else:
                listVar[i][j] = listVar[i][j].split('/')[-1].split(".")[0]


    for combin in dico_combin.keys():
        vars = combin.split('_')
        compte = 0
        for varcom in listVar:
            if all(var in varcom for var in vars):
                compte += 1
        dico_combin[combin] = compte
    

    l = []
    for rep in nb_var:
        indice = files.index(rep.split('/')[-1])
        if legende:
            dico_combin[legende[indice]] = nb_var[rep]
            l.append(legende[indice])
        else:
            dico_combin[rep.split('/')[-1].split(".")[0]] = nb_var[rep]
            l.append(rep.split('/')[-1].split(".")[0])

    binary_map = {elem: format(1 << index, f'0{len(l)}b') for index, elem in enumerate(l)}
    # Fonction pour convertir une clé en sa représentation binaire
    def to_binary_key(key):
        binary_key = '0' * len(l)  # Commence avec une chaîne de zéros de la longueur de a
        for part in key.split('_'):
            binary_key = ''.join('1' if binary_key[i] == '1' or binary_map[part][i] == '1' else '0'
                                 for i in range(len(binary_key)))
        return binary_key

    # Créer le nouveau dictionnaire avec des clés binaires
    labels = {to_binary_key(key): value for key, value in dico_combin.items()}
    fig = None

    if len(files) == 2:
        fig, ax = venn.venn2(labels, names=l)
    elif len(files) == 3:
        fig, ax = venn.venn3(labels, names=l)
    elif len(files) == 4:
        fig, ax = venn.venn4(labels,names=l)
    elif len(files) == 5:
        fig, ax = venn.venn5(labels, names=l)

    if fig:
        if outname:
            fig.savefig(str(outname))
        else:
            fig.savefig("_".join(l))



if __name__ == "__main__":
    dossier = sys.argv[1] 
    inputFiles = sys.argv[2] 
    inputFiles = inputFiles.strip("[]")
    inputFiles = inputFiles.split(",")
    id_seq = float(sys.argv[3])
    outname = sys.argv[4]
    legende = sys.argv[5]
    if legende!='false':
        legende = legende.strip("[]")
        legende = legende.split(",")
 
    diagVenn(dossier,inputFiles,outname,legende,id_seq)
