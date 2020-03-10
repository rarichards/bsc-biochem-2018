# bsc-biochem-2018

Here is the code used for my final year project for my BSc in Biochemistry. The project was set out to analyse the GC content at the third codon position (GC3) in 5’ and 3’ gene termini, and its effect on optimal growth temperature and gene expression levels.

These scripts were used to calculate the GC content at the third codon position (GC3) for the 5’ and 3’ termini of genes in 650 bacterial genomes and 79 archaeal genomes. Genomes were downloaded from the European Nucleotide Archive at https://www.ebi.ac.uk/genomes/. These were filtered down to one genome per genus and genomes of at least 500 kBP in length, resulting in 650 bacterial and 79 archaeal genomes for analysis. Optimal growth temperature was predicted from the GC content of structural RNAs, after accounting for underlying GC pressure. Gene expression levels were predicted using codon usage bias.

The results supported the hypothesis that 5’ GC content is under selection for mRNA structural stability, and identified an unexplored importance in 3’ nucleotide composition. The gene expression results were inconclusive.
