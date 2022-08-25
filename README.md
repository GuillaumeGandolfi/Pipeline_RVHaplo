# Workflow RVHaplo

[![forthebadge](http://forthebadge.com/images/badges/built-with-love.svg)](http://forthebadge.com)  [![forthebadge](https://forthebadge.com/images/badges/it-works-why.svg)](http://forthebadge.com)  [![forthebadge](https://forthebadge.com/images/badges/for-sharks.svg)](http://forthebadge.com)

Ce pipeline réalisé en Snakemake permet d'automatiser l'utilisation de l'outil bioinformatique RVHaplo (https://github.com/dhcai21/RVHaplo.git)

Le but de cet outil est de reconstruire les haplotypes viraux à partir de données de séquençage 'long reads'

## Pour commencer

Première étape : cloner le répertoire Git sur votre machine

```
git clone https://github.com/RainbowBishop/Pipeline_RVHaplo.git
```

Une fois le répertoire cloné, modifié les droits d'exécution du fichier rvhaplo.sh

```
chmod +x rvhaplo.sh
```


### Pré-requis

Pour faire fonctionner le pipeline, il est nécessaire de charger 2 modules sur le cluster
  - Miniconda3
  - python 3.8.12

```
module load system/Miniconda3/1.0

module load system/python/3.8.12
```

Enfin, il faut créer un dossier et y placer le ou les fichiers '.fastq' ou '.fasta' des échantillons séquencés et un dossier qui va contenir les résultats de RVHaplo et du pipeline en général. (Facultatif) : Créer un dossier contenant les fichiers 'sequencing_summary.txt' issus du basecalling

Exemple :

```
mkdir reads_files

mkdir results

(facultatif) mkdir sequencing_summary
```

### Configuration du pipeline

Avant de lancer le pipeline, il est nécessaire de le configurer. Direction le fichier config/config.yaml. 

Soit l'ouvrir directement sur le cluster avec un éditeur (vim, nano, ...) soit le transférer localement pour le modifier

Il faudra alors indiquer les différents chemins des dossiers / fichiers et les paramètres à appliquer


## Démarrage

La ligne de commande pour lancer le pipeline va nécessiter différents arguments

```
--cores [N] | nombre maximum de coeurs utilisables en parallèles pour les différents outils du pipeline

--snakefile [FILE] | fichier snakemake (si le fichier se nomme Snakemake ou snakemake, argument facultatif)

--configfile [FILE] | fichier de configuration

--use-envmodules | pour que le pipeline charge les modules nécessaires à l'utilisation des différents outils

--use-conda | pour que le pipeline installe les différents environnements conda

--conda-frontend conda | pour spécifier que l'on utilise conda et non mamba pour activer les environnements

(Facultatif) -np | lance en mode "dry-run", le pipeline va juste indiquer toutes les tâches qu'il va faire

```

Lors du premier lancement, je conseille le mode "dry-run" pour voir si le pipeline ne détecte pas d'erreur

Exemple : 

```
snakemake -np --configfile config/config.yaml --use-envmodules --use-conda --conda-frontend conda
```

Et pour un "vrai" lancement :

```
snakemake --cores 16 --configfile config/config.yaml --use-envmodules --use-conda --conda-frontend conda
```


## Résultats

Dans le dossier créer pour contenir les résultats, on trouvera un dossier par échantillon inséré dans le pipeline

A l'intérieur on retrouvera :
  - le fichier contenant les haplotypes reconstruits par RVHaplo avec leur abondance (rvhaplo_haplotypes.fasta)
  - les fichiers contenant ces haplotypes avec les séquences des souches / isolats (haplo_refs.fasta & haplo_refs_CP.fasta)
  - les fichiers d'alignement (alignment.fasta & alignment_CP.fasta) 
  - l'arbre phylogénétique construit (NJ_tree.nwk)

## Citation (RVHaplo)

Dehan Cai, Yanni Sun, Reconstructing viral haplotypes using long reads, Bioinformatics, Volume 38, Issue 8, 15 April 2022, Pages 2127–2134, https://doi.org/10.1093/bioinformatics/btac089

### E-mail: dhcai2-c@my.cityu.edu.hk
### Version: V3 (2022-08-02 updated)
