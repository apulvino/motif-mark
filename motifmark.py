#!/usr/bin/env python3

''' 
Anthony Pulvino
The code in this script performs the following tasks: 

1. It iterates over the sequence lines and determines their lengths and makes the FASTA headers the labels for the genes.
2. The lines are also scaled by the sequence length. Draws capitalized sequences from the FASTA file as exons (boxes) based on the input FASTA file.
3. A motifmarks dictionary with key=motif and value=list of coordinates
4. Places lines on across exons and introns on the image with single lines scaled to the length of motif and per base distance apart from each other
5. Iterates over motifs in the motifmarks dictionary with a key=motif and value='list of coordinates designating where lines (designating a given motif) should be
6. It places a very thin rectangle (appears as line in picture) where that given motif is located.
7. A nested for loop ensures that the motifs that get marked on the output image are distinctly colored via the particular colormap (jet is the default internal to the script).
8. The counters ensure lines for each of the respective motifs and exons are properly spaced so there is no overlap across any lines drawn in the output image).'''

## packages imported to run programs internal to this script, all of these will of course,
## need to be installed in the conda environment built to run cairo in
import sys
#print(version)
import argparse
import re
import math
import cairo
import numpy as np
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
## DO NOT FORGET TO CONDA ACTIVATE YOUR PYCAIRO ENV: "cairoenv"

## argparse arguments so the user can specify input motifs list and input FASTA file and opening of files
def get_args():
	parser = argparse.ArgumentParser(add_help=False, description="Parse file for particular motif sequence and output a visualization of where these motifs are along a particular gene")
	parser.add_argument("-f", "--file", help="FASTA input file", required=True)
	parser.add_argument("-m", "--motiflist", help="designates file containing list of motifs of interest inputted by user", default =True)
	parser.add_argument("-h", "--help", action = 'help', help="length of read; the index or biological read length", default=argparse.SUPPRESS)
	return parser.parse_args()

args = get_args()
f = open(args.file,'r')
m = open(args.motiflist, "r")

	 
## the base python code generating coordinate information to pass to cairo later 

## PLEASE NOTE:

##  THE FOLLOWING AWK COMMAND SHOULD BE USED FOR REMOVING NEWLINES FROM THE INPUT FASTA: 
#  awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' file

MOTIF_LIST = []
regex_list = []
regex_dict = {}
coords_dict = {}
exonpos_list = []
fasta_headers = []
fasta_seqs = []
motif_ct = 0

## appends motifs to the MOTIF_LIST
for motif in m:
	motif = motif.strip("\n")
	MOTIF_LIST.append(motif)
   		
    	

def convert_regex(MOTIF_LIST):
 	'''Takes a motif list and returns its complementary regex string to a list.'''
 	regex_str = "(?i)"
 	for motif in MOTIF_LIST:
 		regex_str = ''
 		#makes it so that previous regex strings for a motif aren't added to the next regex string
 		for char in motif:
 			#print(char)
 			if char == "T" or char == "t" or char == "u" or char == "U":
 				regex_str = regex_str + "[TtUu]"
 			elif char == "A" or char == "a":
 				regex_str = regex_str + "[Aa]"
 			elif char == "C" or char == "c":
 				regex_str = regex_str + "[Cc]"
 			elif char == "G" or char == "g":
 				regex_str = regex_str + "[Gg]"
 			elif char == "R" or char == "r":
 				regex_str = regex_str + "[AaGg]"
 			elif char == "Y" or char == "y":
 				regex_str = regex_str + "[TtCcUu]"
 			elif char == "K" or char == "k":
 				regex_str = regex_str + "[GgTtUu]"
 			elif char == "M" or char == "m":
 				regex_str = regex_str + "[AaCc]"
 			elif char == "S" or char == "s":
 				regex_str = regex_str + "[CcGg]"
 			elif char == "W" or char == "w":
 				regex_str = regex_str + "[AaTtUu]"
 			elif char == "B" or char == "b":
 				regex_str = regex_str + "[CcGgTtUu]"
 			elif char == "D" or char == "d":
 				regex_str = regex_str + "[AaGgTtUu]"
 			elif char == "H" or char == "h":
 				regex_str = regex_str + "[AaCcTtUu]"
 			elif char == "V" or char == "v":
 				regex_str = regex_str + "[AaCcGg]"
 			elif char == "N" or char == "n":
 				regex_str = regex_str + "[AaCcGgTtUu]"
 				my_str = ''.join(regex_str)
 		regex_dict[motif] = regex_str
 		#print(regex_dict)
 		## populates the fasta_dict with regex_str stored as a value and motif as the key
 		#print(regex_dict)
 	return regex_str
convert_regex(MOTIF_LIST)


def get_coords(seq):
	'''Takes an input FASTA file and returns, by referencing the motifs list, the coordinate positions per each motif in the FASTA file.'''
	for motif in MOTIF_LIST:
		# initiate a new key into your dict key being the current motif make the value a blank list
		coords_dict[motif] = []
		#print(seq)
		startpos = re.finditer(regex_dict[motif], seq)
		for match in startpos:
			coords_dict[motif].append(match.start())		
					#print(match.start(), end)
					#print(coords_dict)
	return(coords_dict)

f.close()
f = open(args.file,'r')

def get_exons(f):
	'''Takes an input FASTA file and returns, by referencing the input FASTA file, the coordinate positions 
	of all of the exons in the FASTA file.'''
	for line in f:
		if ">" not in line:
			#print(line)
			exon = re.search('[A-Z]+', line)
			# grabs number of exons, where exons are the capitalized strings of character
			# length 2 (min) to 2500(max)
			exonpos_list.append(exon.span())	
			#print(exon.span())	
			#print(exonpos_list)	
	return(exonpos_list)
exons = (exonpos_list)
get_exons(f)

def get_fasta(f):
	'''Takes in a FASTA file and returns a dictionary where the sequence is the value and the header is the key.'''
	f.close()
	f = open(args.file,'r')
	for line in f:
		if ">" in line:
			line = line.strip()
			fasta_headers.append(line)
		else:
			line = line.strip()
			fasta_seqs.append(line)
get_fasta(f)


## Drawing actual genes/marking motifs from the reference FASTA file info generated
## from the code above

#surface = cairo.ImageSurface(cairo.FORMAT_RGB24, 5000, 5000)
#ctx = cairo.Context(surface)

surface = cairo.SVGSurface("plot.svg", 1000, 1000)
ctx = cairo.Context(surface)


## reopening files
f.close()
f = open(args.file,'r')



## creates a list of tuples which are used to specify a colormap, here the "nipy_spectral" colormap
## from matplotlib was chosen
cmap = []
digit = (254/(len(MOTIF_LIST)+1))
for the_motif in MOTIF_LIST:
	motif_ct += 1
	nipy_spectral = cm.get_cmap('nipy_spectral')
	cmap.append(cm.jet(round(motif_ct * digit)-1))

## counters initialized for use below	
y = 50
i = 0
n = 0

## reopening files
f.close()
f = open(args.file,'r')



legend_y_offset = 30
for motif in MOTIF_LIST:
	ctx.set_font_size(10)
	ctx.set_source_rgba(0, 0, 0, 1)
	ctx.move_to(5, (y-15))
	#ctx.show_text(motif)
	ctx.set_line_width(5)
	ctx.move_to(0, 400)
	ctx.stroke()
for index,motif in enumerate(MOTIF_LIST):
	ctx.set_source_rgba(cmap[index][0],cmap[index][1],cmap[index][2],cmap[index][3])
	ctx.set_font_size(10)
	ctx.move_to(700, legend_y_offset)
	ctx.show_text(motif)
	ctx.set_line_width(5)
	ctx.move_to(0, 400)
	ctx.stroke()
	legend_y_offset += 10
for line in fasta_seqs:
	seqlength = len(line)    
	ctx.set_font_size(10)
	ctx.set_source_rgba(0, 0, 0, 1)
	ctx.move_to(5, (y-15))
	ctx.show_text(fasta_headers[i])
	ctx.set_line_width(5)
	ctx.move_to(0, y)
	ctx.line_to(seqlength, y)
	ctx.stroke()
	ctx.rectangle(exons[i][0], y-20, (exons[i][1]-exons[i][0]), 40)
	motifmarks = get_coords(fasta_seqs[i])
	motif_counter = 0
	for index,mark in enumerate(motifmarks):
		ctx.set_source_rgba(cmap[index][0], cmap[index][1], cmap[index][2], cmap[index][3])
		ctx.move_to(700, 90)
		ctx.set_line_width(2)
		ctx.stroke()
		ctx.fill
		for a_motif in motifmarks[mark]:
			ctx.rectangle(a_motif,y-10,3,20)
			ctx.fill()
			for color_val in cmap:
				motif_counter += 1
				ctx.set_source_rgba(cmap[0][0], cmap[0][1], cmap[1][2], cmap[2][3])
	y += 50
	i += 1
surface.finish()




