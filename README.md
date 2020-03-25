# motif-mark
Program capable of highlighting particular motifs in particular genes

The code in this script performs the following tasks: 

1. It iterates over the sequence lines and determines their lengths and makes the FASTA headers the labels for the genes.
2. The lines are also scaled by the sequence length. Draws capitalized sequences from the FASTA file as exons (boxes) based on the input FASTA file.
3. A motifmarks dictionary with key=motif and value=list of coordinates
4. Places lines on across exons and introns on the image with single lines scaled to the length of motif and per base distance apart from each other
5. Iterates over motifs in the motifmarks dictionary with a key=motif and value='list of coordinates designating where lines (designating a given motif) should be
6. It places a very thin rectangle (appears as line in picture) where that given motif is located.
7. A nested for loop ensures that the motifs that get marked on the output image are distinctly colored via the particular colormap (jet is the default internal to the script).
8. The counters ensure lines for each of the respective motifs and exons are properly spaced so there is no overlap across any lines drawn in the output image).


Likewise, it accounts for all overlapping sequences, repeat sequences, and accounts for base length of the motifs themselves.

