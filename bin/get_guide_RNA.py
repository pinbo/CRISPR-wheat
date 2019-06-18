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
from subprocess import call
import getopt, sys, os, re
from functions import *



############
usage="""
get_guide_RNA.py
	-i <sequence.fa>
	-t <seqID1,seqID2>
	-o <output file name>
	-b <blast all produced gRNAs against the genomes: 1 for YES and 0 for NO, default is NO>
	-l <sgRNA length, default is 20>
	-a <whether to target all homeologs: 1 for YES and 0 for NO. default is 1>
	-p <PAM sequence. default is NGG>
	-q <PAM position: left or right>
	-h <print the help>
	-r <reference location>
	-c <Cas9 etc cut position, for testing restriction enzymes>
	-f <gff3 file>
"""

# parameters
pam = "NGG"
pam_pos = "right"
grna_length = 20
target_all = 1 # whether to design sg-RNA for all targets
seqfile = "sequence.fa" # input is a fasta file with homeologs
targets = [] # target sequence names
out = "" # output file
mainID = "" # the one used as template to search sgRNAs
max_price = 200 # maximum restriciton enzyme price
blast = 0 # whether to blast
chrom = "" # chromosome name if target only one specific chromosome; seprated with comma, such as 1A,1B
chr_group = 0 # chromsome group, such as 1, meaning 1A, 1B, 1D

#reference = "/Library/WebServer/Documents/blast/db/nucleotide/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
reference = "/Volumes/DATA3/users/junli/wheat_Refseqv1" # the batmis indexed chromosomes locations
gff = "/Volumes/DATA3/users/junli/wheat_Refseqv1/filtered_sorted_HC.gff3" # sorted

cut_pos = 17 # NGG is 17, TTN is about 18


# steps

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:a:t:o:l:p:q:b:r:c:f:h", ["help"])
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
		if a != "ALL":
			target_all = 0
		print targets
	elif o in ("-o"):
		out = a
	elif o in ("-l"):
		grna_length = int(a)
	elif o in ("-p"):
		pam = a
	elif o in ("-q"):
		pam_pos = a
	elif o in ("-b"):
		blast = int(a)
	elif o in ("-r"):
		blast = 1
		reference = a
	elif o in ("-c"):
		cut_pos = int(a)
	elif o in ("-f"):
		gff = a
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
print "target_all is", target_all
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
muscle_path = get_software_path(getcaps_path)
## STEP 0: create alignment file and primer3output file
RawAlignFile = "alignment_raw.fa"
alignmentcmd = muscle_path + " -in " + seqfile + " -out " + RawAlignFile + " -quiet"
print "Alignment command: ", alignmentcmd
call(alignmentcmd, shell=True)

## parse alignment file
#fasta = get_fasta(RawAlignFile)
#alignlen = len(fasta[mainID])
#print "Alignment length: ", alignlen

## get the target ID template base coordinate in the alignment
#t2a = {} # template to alignment
#a2t = {} # alignment to target
#ngap = 0 # gaps
#for i in range(alignlen):
	#if fasta[mainID][i] == "-":
		#ngap += 1
		#continue
	#t2a[i - ngap] = i
	#a2t[i] = i - ngap

#print "last key of t2a", i - ngap
#print "last key of a2t", i

#wild_seq = fasta[mainID].replace("-","") # remove all gaps
#wild_seq_RC = ReverseComplement(wild_seq)

## Get restriction enzyme information
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
forward_grnas = find_pam(wild_seq, pam, pam_pos, grna_length, "Forward")
reverse_grnas = find_pam(wild_seq_RC, pam, pam_pos, grna_length, "Reverse")

print "RAW forward gRNA number ", len(forward_grnas)
print "RAW reverse gRNA number ", len(reverse_grnas)

## filter gRNAs based on specificity
specific_forward_grnas = []
specific_reverse_grnas = []

for i in forward_grnas:
	#seq = i.seq + i.pam
	if test_spec(i.forblast, targets, non_targets, fasta_raw):
		specific_forward_grnas.append(i)

for i in reverse_grnas:
	#seq = i.seq + i.pam
	if test_spec(i.forblast, targets, non_targets, fasta_raw_RC):
		specific_reverse_grnas.append(i)

print "Specific forward gRNA number ", len(specific_forward_grnas)
print "Specific reverse gRNA number ", len(specific_reverse_grnas)

## find whether there are restrction enzymes (REs) on the selected sequences with PAM
## only REs with recognization site overlap with the 4th position from the end of sg-RNA
for i in specific_forward_grnas:
	cut_pos2 = i.start + cut_pos - 1
	if i.seq4score:
		i.on_target_score = get_on_target_score(i.seq4score) # just to add score here to save time
	#print "cut pos is ", cut_pos, i.seq
	for j in caps_list_forward:
		allpos = j.allpos
		for k in allpos:
			rr = range(k, k + j.length)
			if cut_pos2 in rr:
				#print j.name, j.seq, j.allpos
				i.REs.append(j.name + "," + j.seq)

for i in specific_reverse_grnas:
	cut_pos2 = i.start + cut_pos - 1
	#i.on_target_score = get_on_target_score(i.seq4score)
	for j in caps_list_reverse:
		allpos = j.allpos
		for k in allpos:
			rr = range(k, k + j.length)
			if cut_pos2 in rr:
				i.REs.append(j.name + "," + j.seq)

## blast against the genome
# create a blank blast output file for galaxy output
#call('echo -e "query id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tq. sequence\ts. sequence\tq. length\ts. length" > blast_out.txt', shell=True)
#if blast:
	#specific_forward_grnas = blast_check(specific_forward_grnas, reference, "blast_out_forward.txt")
	#specific_reverse_grnas = blast_check(specific_reverse_grnas, reference, "blast_out_reverse.txt")
	### merge two blast results
	#merge_blast = "cat blast_out_forward.txt blast_out_reverse.txt >> blast_out.txt"
	#call(merge_blast, shell=True)

call('echo "gRNA\tChromosome\tStrand\tPosition\tMismatches\tPotential_target" > blast_out.txt', shell=True)
if blast:
	grna_dict = prepare_blast_file(specific_forward_grnas + specific_reverse_grnas) # output is for_blast.fa
	dum1 = off_target_check("for_blast.fa", reference, getcaps_path + "/mybatmap") # output is out-test-whole.txt
	dum2 = parse_mismatches("out-test-whole.txt", pam, pam_pos, grna_dict, gff)
	call("sort -k1,1 -k5,5n out.temp.txt | awk 'a[$1]++ < 20' >> blast_out.txt", shell=True)

## Print output files

# write to file
outfile = open(out, 'w')
#outfile.write("ID\tStart\tEnd\tStrand\tLength\tSequence (5' -> 3')\tGC_content_All\tGC_content_first_10nt\tReverse Complement (5' -> 3')\tPAM\tTemplate used\tOn_target_score\tRestriction Enzyme\tOff target score of BLAST top 4 hits\n")
#template_length = len(wild_seq)
#for i in specific_forward_grnas:
	#outfile.write("\t".join([i.name, str(i.start + 1), str(i.end + 1), i.direction, str(i.length), i.seq, str(i.gc), str(i.gc10), ReverseComplement(i.seq), i.pam, mainID, "{0:.2f}".format(i.on_target_score), ";".join(i.REs), i.blast]) + "\n")
#for i in specific_reverse_grnas:
	#outfile.write("\t".join([i.name, str(template_length - i.start), str(template_length - i.end), i.direction, str(i.length), i.seq, str(i.gc), str(i.gc10), ReverseComplement(i.seq), i.pam, mainID, "{0:.2f}".format(i.on_target_score), ";".join(i.REs), i.blast]) + "\n")

outfile.write("ID\tStart\tEnd\tStrand\tLength\tOn_target_score\tSequence (5' -> 3')\tGC_content_All\tGC_content_first_10nt\tReverse Complement (5' -> 3')\tPAM\tTemplate used\tRestriction Enzyme\t" + "\t".join(["gRNA", "Chromosome", "Strand", "Position", "Mismatches", "Potential_target"]) + "\n")
template_length = len(wild_seq)
for i in specific_forward_grnas + specific_reverse_grnas:
	outfile.write("\t".join([i.name, str(i.start + 1), str(i.end + 1), i.direction, str(i.length), str(i.on_target_score), i.seq, str(i.gc), str(i.gc10), ReverseComplement(i.seq), i.pam, mainID, ";".join(i.REs), i.blast.strip()]) + "\n")



#print "example score for on target ", get_on_target_score("ATGGGGAACAGAATAGGAGGGAGGAGGAAG")
#print "example score for off target ", get_off_target_score("CCAGGATGGGGCATTTCGAGAGG", "CCAGGATGGGGCATTTCTAAAGG")
