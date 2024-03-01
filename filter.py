#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Hugo
"""
import os
import sys


def find_fichiers_vcf(dossier):
    fichiers_vcf = []
    for racine, dossiers, fichiers in os.walk(dossier):
        for fichier in fichiers:
            if fichier.endswith('.vcf'):
                chemin_complet = os.path.join(racine, fichier)
                fichiers_vcf.append(chemin_complet)
    return fichiers_vcf


def definir_echantillon(dossier):
    fichiers_vcf = find_fichiers_vcf(dossier)
    dico_fichier = {}
    for fichier_vcf in fichiers_vcf:
        id_echantillon = fichier_vcf.split('/')[-1].split('-')[0]
        if id_echantillon not in dico_fichier:
            dico_fichier[id_echantillon] = []
        dico_fichier[id_echantillon].append(fichier_vcf)

    for i in dico_fichier:
        dico_fichier[i] = sorted(dico_fichier[i])
    return dico_fichier


def recup_infos(input_files):
    dicoInfos_seq = {}
    dicoInfos_noseq = {}
    for id_echantillon in input_files:  # Boucle sur les tuples contenant les vcf d'un échantillon
        dicoPos_seq = {}
        dicoPos_noseq = {}
        for vcf_replicat in input_files[id_echantillon]:  # Boucle sur les fichiers vcf de l'échantillon
            id_rep = vcf_replicat
            with open(vcf_replicat, 'r') as vcf:
                for ligne in vcf:
                    if not ligne.startswith('#'):
                        infos = ligne.strip().split("\t")
                        inf = infos[9].split(":")

                        AF = 0.0
                        for i in infos[7].split(";"):
                            if i.startswith("SVLEN"):
                                _, long = i.split('=')
                            if i.startswith("SVTYPE"):
                                _, svtype = i.split('=')
                            if i.startswith("AF"):
                                _, AF = i.split('=')

                        if ">" not in infos[4]:
                            if infos[1] in dicoPos_seq:
                                dicoPos_seq[infos[1]][f"{infos[2]}*_*{id_rep}"] = [id_rep, svtype, long, int(infos[5]),
                                                                                   float(int(inf[-1]) + int(inf[-2])),
                                                                                   float(AF), infos[4]]
                                # 0:id_rep , 1:type, 2:len, 3:qual, 4:dp, 5:AF, 6:alt
                                # nom plusieurs fois dans diff replictas comme Sniffles2.INS.10S0 donc {infos[2]}_{id_rep}b1
                            else:
                                dicoPos_seq[infos[1]] = {}
                                dicoPos_seq[infos[1]][f"{infos[2]}*_*{id_rep}"] = {}
                                dicoPos_seq[infos[1]][f"{infos[2]}*_*{id_rep}"] = [id_rep, svtype, long, int(infos[5]),
                                                                                   float(int(inf[-1]) + int(inf[-2])),
                                                                                   float(AF), infos[4]]
                        else:
                            if infos[1] in dicoPos_noseq:
                                dicoPos_noseq[infos[1]][f"{infos[2]}*_*{id_rep}"] = [id_rep, svtype, long, int(infos[5]),
                                                                                     float(int(inf[-1]) + int(inf[-2])),
                                                                                     float(AF), infos[4]]
                                # 0:id_rep , 1:type, 2:len, 3:qual, 4:dp, 5:AF, 6:alt
                            else:
                                dicoPos_noseq[infos[1]] = {}
                                dicoPos_noseq[infos[1]][f"{infos[2]}*_*{id_rep}"] = {}
                                dicoPos_noseq[infos[1]][f"{infos[2]}*_*{id_rep}"] = [id_rep, svtype, long, int(infos[5]),
                                                                                     float(int(inf[-1]) + int(inf[-2])),
                                                                                     float(AF), infos[4]]
        dicoInfos_seq[id_echantillon] = dicoPos_seq
        dicoInfos_noseq[id_echantillon] = dicoPos_noseq

    return dicoInfos_seq, dicoInfos_noseq


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


def filtre_variant(DicoVar, Qual, DP, AF, FqRep=None):
    if FqRep is None:
        if DicoVar[0] >= Qual and DicoVar[1] >= DP and DicoVar[2] >= AF: return True
        else: return False
    else:
       #af = DicoVar[2]+(FqRep/2)
        af = DicoVar[2] + (FqRep*AF)
        dp = DicoVar[1] + (DicoVar[1]*FqRep)
        if DicoVar[0] >= Qual and dp >= DP and af >= AF: return True
        else: return False


def identifier_variants_communs(pos_variants):
    # Étape 1: Organiser les séquences par réplicat
    dic_rep = {}
    for variant, details in pos_variants.items():
        dic_rep.setdefault(details[0], []).append((variant, details[6]))

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
                            if idseq >= 0.75:
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

    vars = [variant for variant in set(pos_variants.keys())]
    variants_non_traites = set(vars) - variants_traites

    return groupes, variants_non_traites


def filteredVarSeq(dicoPosSeq, dicoVarSeq, inputs_files, Qual, DP, AF):
    variantTrie = {}
    suprime = {}
    for echantillon in dicoPosSeq:
        suprime[echantillon] = []
        variantTrie[echantillon] = {}
        variantTrie[echantillon]['VarCom'] = []
        variantTrie[echantillon]['VarSeul'] = []
        for pos in dicoPosSeq[echantillon]:
            groupe, nonTraite = identifier_variants_communs(dicoPosSeq[echantillon][pos])
            for g in list(groupe):
                if g: variantTrie[echantillon]['VarCom'].append(list(g))
            if nonTraite: variantTrie[echantillon]['VarSeul'].extend(list(nonTraite))

        for varCom in range(len(variantTrie[echantillon]['VarCom']) - 1, -1, -1):
            fqRep = len(variantTrie[echantillon]['VarCom'][varCom]) / len(inputs_files[echantillon])

            VarComTrie = []
            for var in variantTrie[echantillon]['VarCom'][varCom]:
                if filtre_variant(dicoVarSeq[echantillon][var], Qual=Qual, DP=DP, AF=AF,FqRep=fqRep):
                    VarComTrie.append(var)
                else:
                    suprime[echantillon].append(var)
            variantTrie[echantillon]['VarCom'][varCom] = VarComTrie
            if len(variantTrie[echantillon]['VarCom'][varCom]) == 0:
                # Supprimer la sous-liste si elle est vide
                del variantTrie[echantillon]['VarCom'][varCom]
            elif len(variantTrie[echantillon]['VarCom'][varCom]) == 1:
                # Déplacer l'élément dans liste_unique et supprimer la sous-liste
                variantTrie[echantillon]['VarSeul'].extend(variantTrie[echantillon]['VarCom'][varCom])
                del variantTrie[echantillon]['VarCom'][varCom]

        filtered = []
        for var in variantTrie[echantillon]['VarSeul']:
            if filtre_variant(dicoVarSeq[echantillon][var], Qual=Qual, DP=DP, AF=AF):
                filtered.append(var)
            else:
                suprime[echantillon].append(var)
        variantTrie[echantillon]['VarSeul'] = filtered

    return variantTrie, suprime


def filtered_vcf(intput_directory, Qual, DP, AF, output_directory):
    infos=False
    inputs_files = definir_echantillon(intput_directory)
    dicoPosSeq, dicoPosNoseq = recup_infos(inputs_files)

    dicoVarSeq = {
        key_p: {
            key_t: value_t[3:6]
            for key_s in value_p
            for key_t, value_t in value_p[key_s].items()
        }
        for key_p, value_p in dicoPosSeq.items()
    }
    VariantSeqTrie, VarSeqSup = filteredVarSeq(dicoPosSeq, dicoVarSeq, inputs_files,Qual, DP, AF)

    dicoVarNoSeq = {
        key_p: {
            key_t: value_t[3:6]
            for key_s in value_p
            for key_t, value_t in value_p[key_s].items()
        }
        for key_p, value_p in dicoPosNoseq.items()
    }
    VariantNoSeqTrie, VarNoSeqSup = filteredVarSeq(dicoPosNoseq, dicoVarNoSeq, inputs_files,Qual, DP, AF)

    VariantTrie = {}
    VariantsCommun = {} 
    for echantillon in VariantSeqTrie:
        VariantTrie[echantillon] = {}
        VariantsCommun[echantillon]  = {}
        VarComSeq = VariantSeqTrie[echantillon].get('VarCom', [])
        VarComNoSeq = VariantNoSeqTrie.get(echantillon, {}).get('VarCom', [])
        VariantTrie[echantillon]['VarCom'] = VarComSeq + VarComNoSeq


        b=[]
        for varcom in VariantTrie[echantillon]['VarCom']:
            a = [var.split('*')[0]+'_'+var.split('/')[-1].split('.')[0] for var in varcom]
            b.append(a)
        VariantsCommun[echantillon] = b
        
        
        VarSeulSeq = VariantSeqTrie[echantillon].get('VarSeul', [])
        VarSeulNoSeq = VariantNoSeqTrie.get(echantillon, {}).get('VarSeul', [])
        VariantTrie[echantillon]['VarSeul'] = VarSeulSeq + VarSeulNoSeq


    for echantillon in VariantNoSeqTrie:
        if echantillon not in VariantTrie:
            VariantTrie[echantillon] = {}
            VariantTrie[echantillon]['VarCom'] = VariantNoSeqTrie[echantillon].get('VarCom', [])
            VariantTrie[echantillon]['VarSeul'] = VariantNoSeqTrie[echantillon].get('VarSeul', [])


    VarSup = {}
    for dico in [VarSeqSup, VarNoSeqSup]:
        for echantillon, variants in dico.items():
            if echantillon not in VarSup:
                VarSup[echantillon] = {}
            for variant in variants:
                rep = variant.split('*_*')[-1]
                if rep not in VarSup[echantillon]:
                    VarSup[echantillon][rep] = []
                VarSup[echantillon][rep].append(variant.split('*_*')[0])


    os.makedirs(output_directory, exist_ok=True)
    print(f"dossier {output_directory} crée")
    nb_var = {}
    for echantillon in inputs_files:
        nb_var[echantillon] = {}
        chemin_sous_dossier = os.path.join(output_directory, echantillon)
        os.makedirs(chemin_sous_dossier, exist_ok=True)
        print(f"\tdossier {chemin_sous_dossier} crée")
        for vcf in inputs_files[echantillon]:
            nb_var[echantillon][vcf] = 0
            chemin_FilteredFile = os.path.join(chemin_sous_dossier, f"{vcf.split('/')[-1].split('.')[0]}.filtered.vcf")
           # print(f"\t\tfichier {chemin_FilteredFile} crée")
            with open(vcf, 'r') as intput_files, open(chemin_FilteredFile, 'w') as output_files:
                for ligne in intput_files:
                    if ligne.startswith('#'):
                        output_files.write(ligne)
                    else:
                        nb_var[echantillon][vcf] +=1
                        variant_id = ligne.split('\t')[2]
                        if vcf not in VarSup[echantillon]:
                            VarSup[echantillon][vcf] = []
                        if variant_id not in VarSup[echantillon][vcf]:
                            output_files.write(ligne)


    if True:
        print("\n")
        for echantillon, info in VariantTrie.items():
            print(f"\n--------------{echantillon}-------------")
            if 'VarCom' in info:
                nb_var_com = len(info['VarCom'])
                print(f"Il y a {nb_var_com} variant(s) commun(s)")
                compteur_replicats = {}
                for var_com in info['VarCom']:
                    nb_replicats = len(var_com)
                    if nb_replicats not in compteur_replicats:
                        compteur_replicats[nb_replicats] = 0
                    compteur_replicats[nb_replicats] += 1

                sorted_keys = sorted(compteur_replicats.keys())
                for sorted_key in sorted_keys:
                    print(f"\t{compteur_replicats[sorted_key]} sont retrouvés dans {sorted_key} réplicat(s)")

            if 'VarSeul' in info:
                nb_var_seul = len(info['VarSeul'])
                print(f"Il y a {nb_var_seul} variant(s) seul(s)\n")

            sorted_files = sorted(list(VarSup[echantillon].keys()),
                                  key=lambda x: int(x.split('/')[-1].split('-')[1].split('.')[0]))

            nb_echSup = 0
            for rep in VarSup[echantillon].values(): nb_echSup += len(rep)
            nb_echVar = 0
            for nb in nb_var[echantillon].values(): nb_echVar += nb

            x = 0
            if nb_echVar != 0: x = round((nb_echSup/nb_echVar)*100,2)
            print(f"\n{nb_echSup} sur {nb_echVar} variants suprimés dans l'échantillon {echantillon} ({x}%)")
            for sorted_file in sorted_files:
                x = 0
                if nb_var[echantillon][sorted_file] != 0: x = round((len(VarSup[echantillon][sorted_file]) / nb_var[echantillon][sorted_file]) * 100, 2)
                print(f"\t{len(VarSup[echantillon][sorted_file])} variants suprimés sur {nb_var[echantillon][sorted_file]}"
                      f" dans le replicat {sorted_file.split('/')[-1].split('.')[0]} "
                      f"({x}%)")

    return VariantsCommun

                
                
#a=filtered_vcf("VCF_sniffle",Qual=60, DP=300, AF=0.1, id_seq=0.2)
#print(a)
#./execute.sh -ipd Resultats -q 60 -dp 300 -id 0.8 -af 0.2
if __name__ == "__main__":
    intput_directory = sys.argv[1] 
    Qual = int(sys.argv[2])
    DP = int(sys.argv[3])
    AF = float(sys.argv[4])
    output_directory = sys.argv[5] 
    filtered_vcf(intput_directory, Qual, DP, AF, output_directory)
