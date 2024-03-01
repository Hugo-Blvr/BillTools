#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 19:35:44 2024

@author: vetea
"""

# résultats finaux dans le fichier "synthese_finale.csv"

# importation des packages de gestion de fichiers (système)
from os import listdir
from os.path import isfile
from os import remove
import os
import sys
import filter

in_dir = sys.argv[1]   
out_dir = sys.argv[2]  

# récupération des fichier vcf a partir du in_dir
dic_tri = filter.definir_echantillon(in_dir)
vcf_files = set()
for files_liste in dic_tri.values():
    for vcf_file in files_liste:
        vcf_files.add(vcf_file)


vcf_files = sorted(vcf_files)
print("Début du script d'extraction")

# ouverture en écriture d'un fichier de réception des informations extraites (synthèse.csv) et d'un fichier "rapport" pour les infos globales pour chaque réplicat biologique (nb variants)
os.makedirs(out_dir, exist_ok=True)
fws = open(out_dir + "/synthese_inter.csv", "w")
fwr = open(out_dir + "/rapport.txt", "w")
# écriture des entêtes
fws.write("Passage,Replicat_Biologique,Stress,Position,Reference,Alternative,Qualite")
fwr.write("Passage,Replicat_Biologique,Nb_Variants")
# parcours des fichiers du dossier courant pour trouver les vcf
for vcf_file in vcf_files :
        with open(vcf_file, "r") as fr:    
            nb_var = 0
            # l'information du passage, réplicat biologique et technique est récupérée depuis le nom du fichier
            # pour les passages précédents, pas d'informations au niveau du réplicat technique, donc besoi nde vérification de la présence d'un chiffre de réplicat, et sinon, utilisation de la valeur "NA" pour le traitement R
            pas = vcf_file.split("/")[-1].split("-")[0][1:3]
            rep_b = (".".join(vcf_file.split(".")[:2])).split("/")[-1].split("-")[1].split(".")[0]
            # pour chaque ligne du fichier, les informations pertinentes séparées par des tabulations sont récupérées et écrites dans le fichier de réception sous format csv (séparation par des virgules)
            for l in fr:
                # pour toute ligne ne faisant pas partie de l'entête du format vcf et dont l'identifiant est présent dans les fichiers filtrés du script d'Hugo...
                if (l[0] != "#"):
                    # récupération des informations tabulées
                    spl = l.split("\t")
                    if len(spl[4]) > 2:
                        spl[4] = spl[4].replace(",", "/") 
                    # écriture dans le fichier de réception, avec une valeur différente dans la colonne "Stress" selon le numéro du réplicat biologique
                    if int(rep_b) < 6 :
                        fws.write("\n{},{},{},{},{},{},{}".format(pas, rep_b, "cold",spl[1],spl[3],spl[4],spl[5]))
                    else :
                        fws.write("\n{},{},{},{},{},{},{}".format(pas, rep_b, "hot",spl[1],spl[3],spl[4],spl[5]))
                        
                    nb_var += 1
            # écriture du rapport global du réplicat biologique avec les informations de nombres de variants
            fwr.write(pas + "," + rep_b + "," + str(nb_var)+ "\n")   
                                                                  
fws.close()
fwr.close()

print("Début de la création du fichier final")

# Deuxième partie : créeation d'un fichier csv compilant tous les résultats

# ouverture en lecture des deux fichiers précédemment créés, et en écriture d'un fichier de réception final
fe = open(out_dir + "/synthese_inter.csv", "r")
ff = open(out_dir + "/synthese_finale.csv", "w")
# Création de l'entête du fichier de réception 
ff.write("Passage,Replicat_Biologique,Stress,Position,Reference,Alternative,Qualite\n")


for le in fe :
    with open(out_dir + "/rapport.txt", "r") as fg :
        for lg in fg : 
            sple = le.split(",")
            splg = lg.split(",")
            # si les numéros de réplicats correspondent entre les deux fichiers, écriture des lignes correspondantes concaténées (en enlevant un retour à la ligne et les numéros de réplicats de l'un) 
            if sple[0] == splg[0] and sple[1] == splg[1] and sple[1] != "Replicat_Biologique":    
               # print("Fusion des informations du réplicat P" + sple[0] + "-" + sple[1])
                ff.write(le[:-1] + "," + ",".join(lg.split(",")[2:]))



fe.close()
fg.close()
ff.close()

# nettoyage : suppression des fichiers intermédiaires
remove(out_dir + "/synthese_inter.csv")
remove(out_dir + "/rapport.txt")

print("Script terminé sans erreur\n")