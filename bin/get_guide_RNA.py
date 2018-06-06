#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  get_guide_RNA.py
#  
#  Copyright 2018 Junli Zhang <zhjl86@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
### Imported
#from subprocess import call
import getopt, sys, os, re
from functions import *



############
usage="""
get_guide_RNA.py
	-i <sequence.fa>
	-t <seqID1,seqID2>
	-o <output file name>
	-b <blast all produced primers against the genomes: 1 for YES and 0 for NO, default is NO>
	-l <sgRNA length, default is 20>
	-a <whether to target all homeologs: 1 for YES and 0 for NO. default is 1>
	-p <PAM sequence. default is NGG>
	-h <print the help>
"""

# parameters
pam = "NGG"
grna_length = 20
target_all = 1 # whether to design sg-RNA for all targets
seqfile = "sequence.fa" # input is a fasta file with homeologs
targets = [] # target sequence names
out = "" # output file
mainID = "" # the one used as template to search sgRNAs
max_price = 200 # maximum restriciton enzyme price
# steps

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:a:t:o:l:p:b:h", ["help"])
except getopt.GetoptError as err:
	# print help information and exit:
	print str(err)  # will print something like "option -a not recognized"
	print usage
	sys.exit(2)
for o, a in opts:
	if o == "-i":
		seqfile = a
	elif o in ("-h", "--help"):
		print usage
		sys.exit()
	elif o in ("-a"):
		target_all = int(a)
	elif o in ("-t"):
		targets = a.split(",")
		target_all = 0
		print targets
	elif o in ("-o"):
		out = a
	elif o in ("-l"):
		grna_length = int(a)
	elif o in ("-p"):
		pam = a
	elif o in ("-b"):
		blast = int(a)
	else:
		assert False, "unhandled option"
print "Options done"

## read fasta
# get target and ids and rename fasta seq names
fasta_raw = get_fasta(seqfile)
fasta_raw_RC = get_fasta(seqfile)
for k in fasta_raw_RC:
	fasta_raw_RC[k] = ReverseComplement(fasta_raw_RC[k])

seq_names = fasta_raw.keys()
if target_all == 1:
	targets = seq_names

mainID = targets[0]
non_targets = [i for i in seq_names if i not in targets] # other non-target sequence names
print "non_targets are ", non_targets

wild_seq = fasta_raw[mainID]
wild_seq_RC = ReverseComplement(wild_seq)
groupname = "-".join(targets)
if not out:
	out = 'selected_gRNAs_for_' + groupname + ".txt"

getcaps_path = os.path.dirname(os.path.realpath(__file__))


# software path
#primer3_path, muscle_path = get_software_path(getcaps_path)


## Get restriciton enzyme information
# step 1: read the enzyme file
RE_file = getcaps_path + "/NEB_parsed_REs.txt" # this one removed some duplicated cuttings
REs = parse_RE_file(RE_file) # get the list of restriction enzymes
# step 2: get the list of enzymes that can be used
caps_list_forward = []
caps_list_reverse = []
for k in REs: # k is enzyme name + price, such as BccI,66
	enzyme_seq = REs[k]
	price = int(k.split(',')[-1])
	if price > max_price:
		continue
	enzyme1 = test_enzyme(k, enzyme_seq, wild_seq)
	enzyme2 = test_enzyme(k, enzyme_seq, wild_seq_RC)
	if enzyme1.caps == "Yes":
		caps_list_forward.append(enzyme1)
	if enzyme2.caps == "Yes":
		caps_list_reverse.append(enzyme2)

#print "RE list is ", caps_list_forward + caps_list_reverse

## find all potential gRNAs in the sequences and RC sequence
forward_grnas = find_pam(wild_seq, pam, grna_length, "Forward")
reverse_grnas = find_pam(wild_seq_RC, pam, grna_length, "Reverse")

print "RAW forward gRNA number ", len(forward_grnas)
print "RAW reverse gRNA number ", len(reverse_grnas)

## filter gRNAs based on specificity
specific_forward_grnas = []
specific_reverse_grnas = []

for i in forward_grnas:
	seq = i.seq + i.pam
	if test_spec(seq, targets, non_targets, fasta_raw):
		specific_forward_grnas.append(i)

for i in reverse_grnas:
	seq = i.seq + i.pam
	if test_spec(seq, targets, non_targets, fasta_raw_RC):
		specific_reverse_grnas.append(i)

print "Specific forward gRNA number ", len(specific_forward_grnas)
print "Specific reverse gRNA number ", len(specific_reverse_grnas)

## find whether there are restrction enzymes (REs) on the selected sequences with PAM
## only REs with recognization site overlap with the 4th position from the end of sg-RNA
for i in specific_forward_grnas:
	cut_pos = i.end - 3
	#print "cut pos is ", cut_pos, i.seq
	for j in caps_list_forward:
		allpos = j.allpos
		for k in allpos:
			rr = range(k, k + j.length)
			if cut_pos in rr:
				#print j.name, j.seq, j.allpos
				i.REs.append(j.name + "," + j.seq)

for i in specific_reverse_grnas:
	cut_pos = i.end - 3
	for j in caps_list_reverse:
		allpos = j.allpos
		for k in allpos:
			rr = range(k, k + j.length)
			if cut_pos in rr:
				i.REs.append(j.name + "," + j.seq)


## Print output files

# write to file
outfile = open(out, 'w')
outfile.write("Start\tEnd\tStrand\tLength\tSequence\tGC_content\tReverseComplement\tPAM\tTemplate\tRestriction_Enzyme\n")
template_length = len(wild_seq)
for i in specific_forward_grnas:
	outfile.write("\t".join([str(i.start + 1), str(i.end + 1), i.direction, str(i.length), i.seq, str(i.gc), ReverseComplement(i.seq), i.pam, mainID, ";".join(i.REs)]) + "\n")

for i in specific_reverse_grnas:
	outfile.write("\t".join([str(template_length - i.start), str(template_length - i.end), i.direction, str(i.length), i.seq, str(i.gc), ReverseComplement(i.seq), i.pam, mainID, ";".join(i.REs)]) + "\n")


