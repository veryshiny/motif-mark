# Python script using object-oriented code to visualize motifs on sequences

- Script is named motif-mark-oop.py and has the following file inputs: 

    -f: Input FASTA file (seqs ≤1000 bases)
    -m: Motifs file (≤10 bases each, one motif per line in a text file)

- Output file will have same prefix as input file (e.g. Figure_1.fa -> Figure_1.png) and have one single, well-labeled figure, per FASTA file
    - All features (motifs, introns, exons) are to scale
    - Shows different motifs in different colors identifiable with legend
    - Shows overlap of motifs 
    - Shows introns and exons
  

- Script is capable of handling
    - Up to 5 Motifs with ambiguous nucleotides (see https://en.wikipedia.org/wiki/Nucleic_acid_notation)
    - Multiple sequences 



