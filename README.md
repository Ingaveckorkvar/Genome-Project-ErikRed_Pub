# Genome-Project-ErikRed_Pub

Undergrad Student First Computational Biology Project.

Genome project named after Erik Röde, the first explorer to have explored greenland. Similarly to how Erik explored greenland, this will be me exploring the world of bioinformatics, and considering this is my first time exploring computational biology, I believed the name to be fitting.

This repository will contain some solutions to rosalind problems, mostly done by me, though occasionally with the help of the rosalind forum (in cases where I may have been stumped or my code may be particularly chunky).

Repository will also contain code found with biopython documentation, mostly playing around. Datasets more likely than not either collected from rosalind or from NCBI library.

This repository is simply for me to both document and demonstrate in the future my skills within bioinformatics. Purely educational, no commercial use.

## What can ErikRed do?

ErikRed is a basic script that can read FASTA/CLUSTAL files and produce a "report" of sorts, which details multiple attributes about the protein it is given data about, including:
- GRAVY score
- Kyte & Doolittle hydrophobicity scores, position-by-position (helps visualise "inside" and "outside" parts of a protein)
- Shannon Entropy 
- Conservation of amino acids (at certain positions, helps determine which are essential for function)
- Phylogenetic trees (By comparing the same protein's AA sequence to other species provided, to produce a graph of both constant mutation rate & graph w/o mutation rate assumption)
- Pairwise Likeness % (To show, on a graph, which species are which most-like other species)
- AA Composition of the protein
- ErikRed can also save all terminal outputs to a .txt file for large datasets.
- All graphs produced by ErikRed are saved in the same folder.

### Credit where it's due
Big thank you to my friend, Matt Bergström, who helped me on this. I am nowhere near this level and could never accomplish this completely solo.
Likewise credit is due to biopython documentation, various forums and rosalind for helping me get somewhere.
