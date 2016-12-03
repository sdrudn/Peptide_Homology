#!/usr/bin/env python2.7
import sys
import re
import io
import argparse
import unicodedata

# Written by Shawn O'Neil (Center for Genome Research and Biocomputing)
# and 
# Soeren Nielsen
# Oregon State University
# Nov, 2016

#This software program and documentation are copyrighted by
#Oregon State University. The software program and
#documentation are supplied "as is", without any accompanying
#services from Oregon State University. OSU does not warrant
#that the operation of the program will be uninterrupted or
#error-free. The end-user understands that the program was
#developed for research purposes and is advised not to rely
#exclusively on the program for any reason.

#IN NO EVENT SHALL OREGON STATE UNIVERSITY BE LIABLE TO ANY
#PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
#CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
#THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
#OREGON STATE UNIVERSITYHAS BEEN ADVISED OF THE POSSIBILITY OF
#SUCH DAMAGE. OREGON STATE UNIVERSITY SPECIFICALLY DISCLAIMS
#ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#PURPOSE AND ANY STATUTORY WARRANTY OF NON-INFRINGEMENT. THE
#SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND
#OREGON STATE UNIVERSITY HAS NO OBLIGATIONS TO PROVIDE
#MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
#MODIFICATIONS.





############################
## Some handy functions
############################


### returns 1 if aa1 is equal to aa2, and 0 otherwise
def score_identical(aa1, aa2):
	if aa1 == aa2:
		return 1
	return 0

# use a global lookup so we don't recreate it all the time
# -> minor speed increase
map = dict()
# D=E, I=L=V, Q=N and S=T 
map["D"] = set(["E"])
map["E"] = set(["D"])
map["I"] = set(["L", "V"])
map["L"] = set(["I", "V"])
map["V"] = set(["L", "I"])
map["Q"] = set(["N"])
map["N"] = set(["Q"])
map["S"] = set(["T"])
map["T"] = set(["S"])

### returns 1 if aa1 is equal to aa2, or similar in chemistry, otherwise 0
def score_similar(aa1, aa2):
	global map
	if aa1 == aa2:
		return 1
	
	if map.has_key(aa1):
		if aa2 in map[aa1]:
			return 1
			
	return 0
	

### returns the percentage of an alignment
### e.g. ARNDC against ARQDV returns 0.8 (sim, sim, sim, sim, not sim) if use_identical is false 
### and 0.6 (sim, sim, not sim, sim, not sim) if true
### this assumes that the two sequences will be the same length! 
def score_align(seq1, seq2, use_identical):
	score = 0.0
	for i in range(0, len(seq1)):
		if use_identical:
			score = score + score_identical(seq1[i], seq2[i])
		else:
			score = score + score_similar(seq1[i], seq2[i])

	return score/len(seq2)



# +1hr to here
### returns a list of tuples (or lists)
### like [(start, stop, match, %homology), (start, stop, match, %homology)]
### (using base-1 indexing)
### for those that are more %homologous than the given threshold
### if use_identical = True, uses score_align_identical(), otherwise uses score_align_similar()
def find_matches(short, long, threshold, use_identical):
	matches = list()

	for i in range(0, len(long) - len(short) + 1):
		substr = long[i:(i + len(short))]
		align_score = score_align(short, substr, False)
		if use_identical:
			align_score = score_align(short, substr, True)
			
		if align_score >= threshold:
			match = (i + 1, i + len(short), substr, align_score)
			matches.append(match)

	return matches
	



############################
## Parse Arguments
############################

parser = argparse.ArgumentParser(description="Matches a set of peptides against a set of proteins using simple matching rules. Only reports matches with a percent identity above the pident threshold.")

parser.add_argument('--proteins', required = True, dest = "proteins_file", help = "The file to read protein sequences from. Assumes text and 4 tab-separated columns: protein name, species, sequence, uniprot id")
parser.add_argument('--peptides', required = True, dest = "peptides_file", help = "The file to read (shorter) peptide sequences from. Assumes text and 4 tab-separated columns: protein name, sequence, function, species")
parser.add_argument('--pident_threshold', required = True, dest = "pident_str", help = "Only report matches with this percent identity or higher, e.g. 80")
parser.add_argument('--use_identical', required = False, action = 'store_true', default = False, dest = "use_identical", help = "Use only exact matches for scoring. Otherwise, let D/E match, as well as I/L/V, Q/N, and S/T.")
args = parser.parse_args()
pident_thresh = float(args.pident_str)/100.0





############################
## Read protein file input (as list of list)
############################

proteins_list = list()
fhandle = io.open(args.proteins_file, "rU", errors = "ignore")
for line in fhandle:
	lineascii = unicodedata.normalize("NFKD", line).encode('ascii', 'ignore')
	line_list = lineascii.strip().split("\t")
	if len(line_list) != 4:
		print(line_list)
		raise ValueError("This line in the proteins file doesn't have 4 tab-separated columns (often caused by oddities in export from Excel): " + line)
	proteins_list.append(line_list)



############################
## Read peptide file input (as list of list)
############################

peptides_list = list()
fhandle = io.open(args.peptides_file, "rU", errors = "ignore")
for line in fhandle:
	lineascii = unicodedata.normalize("NFKD", line).encode('ascii', 'ignore')
	line_list = lineascii.strip().split("\t")
	if len(line_list) != 4:
		raise ValueError("This line in the proteins file doesn't have 4 tab-separated columns (often caused by oddities in export from Excel): " + line)
	peptides_list.append(line_list)
	





############################
## Doooo it
############################
print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("prot_name", "prot_species", "match_start", "match_end", "pep_sequence", "match_seq", "pep_name", "pep_species", "pep_function", "match_pident"))
for protein in proteins_list:
	sys.stderr.write("Processing protein: " + protein[0] + "\n")
	for peptide in peptides_list:
		prot_name, prot_species, prot_sequence, prot_uniprot = protein
		pep_name, pep_sequence, pep_function, pep_species = peptide

		matches = find_matches(pep_sequence, prot_sequence, pident_thresh, args.use_identical)
		for match in matches:
			match_start, match_end, match_seq, match_pident = match
			print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(prot_name, prot_species, match_start, match_end, pep_sequence, match_seq, pep_name, pep_species, pep_function, match_pident))

