#!/bin/bash

# Valeurs par défaut pour les paramètres de filtrage
id_default=0.75
name_fig_default=false
legende_default=false
quality_default=50
deep_default=300
allele_frequency_default=0.1
output_directory_default="vcf_filtre"


# Initialisation des variables avec les valeurs par défaut
id_sequence=$id_default
name_fig=$name_fig_default
legende=$legende_default
quality=$quality_default
deep=$deep_default
allele_frequency=$allele_frequency_default
output_directory=$output_directory_default
compare_mode=false
bcf_used=false
SNP_mode=false

# Fonction pour afficher l'aide
show_help() {
    echo "Usage: de $0:"
    echo ""
    echo "Options générales :"
    echo "  -h, --help                   Affiche ce message d'aide."
    echo "  -i, --input_directory        Chemin du dossier à parcourir pour les opérations générales. Obligatoire."
    echo "  -q, --quality                Qualité minimale requise. Par défaut $quality_default."
    echo "  -dp, --deep                  Profondeur minimale requise. Par défaut $deep_default."
    echo "  -af, --allele_frequency      Fréquence allélique minimale requise. Par défaut $allele_frequency_default."
    echo "  -o, --output_directory       Nom du répertoire de sortie pour les fichiers filtrés. Par défaut '$output_directory_default'."
    echo "  -bcf, --bcftools             Utilise bcftools pour le filtrage. Par défaut 'false'."
    echo "Exemple d'utilisation : $0 -i MesFichierVcf -dp 200 -o MesFichierVcfFiltre -bcf"
    echo ""
    echo "Options pour le mode SNP (--SNP) :"
    echo "  --SNP                    Active le mode SNP pour extraire les données de fichiers VCF SNP."
    echo "  -i, --input_directory        Chemin du dossier à parcourir pour l'extraction. Obligatoire."
    echo "Exemple d'utilisation : $0 --SNP -i SNPsFiles/"
    echo ""
    echo "Options pour le mode comparaison (--compare) :"
    echo "  --compare                    Active le mode comparaison pour analyser les différences entre les fichiers VCF spécifiés."
    echo "  -if, --inputFiles            Spécifie les fichiers d'entrée pour la comparaison. Requis en mode --compare."
    echo "                                        *A mettre entre crochet,[idvcf1,idvcf2..]; 5 maximum"
    echo "  -i, --input_directory        Chemin du dossier à parcourir pour la comparaison. Répertoire courant par default"
    echo "  -id, --id_sequence           Définit l'identité de séquences minimale requise. Par défaut $id_default."
    echo "  -nf, --name_fig              Nom de la figure générée pour les résultats de comparaison. Par default identifiant des fichier vcf joint par '_'."
    echo "  -l, --legende                Légende pour la figure générée dans le mode comparaison. Par défaut identifiant des fichier vcf."
    echo "                                        *A mettre entre crochet,[legende1,legende2..]"
    echo "Exemple d'utilisation : $0 --compare -i MesFichierVcfFiltre -if [idvcf1,idvcf2] -l[l1,l2] -nf diagVenn -"
    echo ""
    echo "Note: Le mode comparaison (--compare) crée un diagramme de Venn et requiert l'utilisation de l'option --inputFiles"
    echo "pour spécifier les fichiers à comparer. Les options comme --id_sequence, --name_fig et --legende permettent de personnaliser"
    echo "davantage l'analyse de comparaison. L'id d'un vcf contient la chaine de caractere avant le premier point dans le nom du fichier."
    echo ""
}


# Traitement des options de ligne de commande
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        --SNP)
            SNP_mode=true
            ;;
        --compare)
            compare_mode=true
            ;;
        -i|--input_directory)
            input_directory="$2"
            shift
            ;;
        -if|--inputFiles)
            inputFiles="$2"
            shift
            ;;
        -id|--id_sequence)
            id_sequence="$2"
            shift
            ;;
        -nf|--name_fig)
            name_fig="$2"
            shift
            ;;
        -l|--legende)
            legende="$2"
            shift
            ;;
        -q|--quality)
            quality="$2"
            shift
            ;;
        -dp|--deep)
            deep="$2"
            shift
            ;;
        -af|--allele_frequency)
            allele_frequency="$2"
            shift
            ;;
        -o|--output_directory)
            output_directory="$2"
            shift
            ;;
        -bcf|--bcftools)
            bcf_used=true
            ;;
        *)
            echo "Option invalide: $1"
            show_help
            exit 1
            ;;
    esac
    shift
done

# Fonction de comparaison
compare() {
    if [ -z "$input_directory" ]; then
        ipp='.'
    else
        ipp=$input_directory
    fi

    local inputFiles=$1  # Utilisez la valeur passée en argument
    python3 VarComFind.py "$ipp" "$inputFiles" "$id_sequence" "$name_fig" "$legende"
}



# Exécuter les opérations en fonction du mode
if [ "$compare_mode" = true ]; then
    if [ -z "$inputFiles"]; then
        echo "Les fichiers d'entrée sont requis pour l'option --compare."
        exit 1
    fi
    compare "$inputFiles"
    exit 0
elif [ "$SNP_mode" = true ]; then
    if [ -z "$input_directory" ]; then
        echo "Le répertoire d'entrée est requis."
        exit 1
    fi
    python3 extractSNP.py "$input_directory" "$input_directory"
    exit 0
else
    if [ -z "$input_directory" ]; then
        echo "Le répertoire d'entrée est requis."
        show_help
        exit 1
    fi
    python3 filter.py "$input_directory" "$quality" "$deep" "$allele_frequency" "$output_directory"
    python3 extract.py "$output_directory" "$output_directory"
fi


# Utilisation de bcftools si demandé
if [ "$bcf_used" = true ]; then
    echo "Utilisation de bcftools pour le filtrage."
    vcf_files=$(find "$input_directory" -type f -name "*.vcf")
    mkdir -p "vcfTools_filtres"
    for vcf_file in $vcf_files; do
        filename=$(basename "$vcf_file")
        name=$(echo "$filename" | cut -d'.' -f1)
        output_file="vcfTools_filtres/${name}.vcftools.filtered.vcf"
        bcftools filter -i "AF > $allele_frequency & (DR + DV > $deep)" "$vcf_file" -o "$output_file"
    done
    python3 extract.py "vcfTools_filtres" "vcfTools_filtres"
fi
