from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


# Read the alignment file

alignment = AlignIO.read(snakemake.input[0], "fasta")
print(alignment)

# Calculare the distance matrix

calculator = DistanceCalculator('identity')
distance_Matrix = calculator.get_distance(alignment)
print(distance_Matrix)

# Create a DistanceTreeConstructor object

constructor = DistanceTreeConstructor()

# Construct the phlyogenetic tree using NJ algorithm

NJ_tree = constructor.nj(distance_Matrix)

# Draw the phlyogenetic tree using terminal

Phylo.draw_ascii(NJ_tree)

# Write tree in new file
Phylo.write(NJ_tree, snakemake.output[0], "newick")