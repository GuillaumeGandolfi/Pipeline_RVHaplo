#!/bin/bash 

# Nom du job
#SBATCH --job-name=Workflow_RVHaplo

# Choix de la partition
#SBATCH -p highmemdell --nodelist=node30

# Déplacement vers répertoire de l'outil
cd $1

# Chargement des modules
module load system/python/3.8.12
module load system/Miniconda3/1.0

# lancement de snakemake

snakemake --cores $2 --configfile config/config.yaml --use-envmodules --use-conda --conda-frontend conda

#Infos sur le job
seff $SLURM_JOB_ID