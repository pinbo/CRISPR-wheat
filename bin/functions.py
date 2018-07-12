#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  functions.py
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
from subprocess import call, check_output
import getopt, sys, os, re, copy, pickle
from Rule_Set_2_scoring_v1 import model_comparison
from CFD_Score import cfd_score_calculator as cfd


###############
def main(args):
	return 0

if __name__ == '__main__':
	import sys
	sys.exit(main(sys.argv))

##################### CLASSES   ######################
# guide RNA object
class gRNA(object):
	"""Information of guide RNA"""
	def __init__(self):
		self.name = ""
		self.start = 0 # 0 based
		self.end = 0 # 0 based
		self.length = 0
		self.gc = 0.0
		self.gc10 = 0.0
		self.seq = "" # sequence without PAM
		self.pam = "" # PAM
		self.pam_pos = "" # left or right
		self.forblast = "" # seq + pam or pam + seq for blast
		self.direction = ""
		self.template = ""
		self.REs = [] # list of restriction enzymes
		self.seq4score = "" # 30 nt for calculating on target score
		self.on_target_score = 0.0
		self.off_target_score = ""
		self.blast = "" # blast hits
		#self.on_target_score = get_on_target_score(self.seq4score)
	#@property
	#def on_target_score(self):
	#	return get_on_target_score(self.seq4score)


class Restriction_Enzyme(object):
	def __init__(self, name, seq):
		self.name = name
		self.seq = seq.lower()
		self.length = len(seq)
		self.template_seq = ""
		self.caps = "No"
		self.allpos = [] # all the match positions in the template
		self.price = int(name.split(',')[-1])
		self.primer_direction = "" # to use the end positions as left or right primer


##################### FUNCTIONS ######################
# function to extract sequences from a fasta file 
def get_fasta(infile):
	fasta = {} # dictionary for alignment
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line.startswith(">"):
				sequence_name = line.split()[0].lstrip(">")
			else:
				fasta.setdefault(sequence_name, "")
				fasta[sequence_name] += line.strip().replace(" ", "") # remove spaces in case
	return fasta

# get reverse complement sequence
def ReverseComplement(seq):
	s1 = "BDHKMNRSVWYATGCbdhkmnrsvwyatgc"
	s2 = "VHDMKNYSBWRTACGvhdmknysbwrtacg"
	seq_dict = {s1[i]:s2[i] for i in range(len(s1))}
	return "".join([seq_dict[base] for base in reversed(seq)])

#from sys import platform
def get_software_path(base_path):
	if sys.platform.startswith('linux'): # linux
		primer3_path = base_path + "/primer3_core"
		muscle_path = base_path + "/muscle"
	elif sys.platform == "win32" or sys.platform == "cygwin": # Windows...
		primer3_path = base_path + "/primer3_core.exe"
		muscle_path = base_path + "/muscle.exe"
	elif sys.platform == "darwin": # MacOSX
		primer3_path = base_path + "/primer3_core_darwin64"
		muscle_path = base_path + "/muscle3.8.31_i86darwin64"
	return primer3_path, muscle_path

# parse the restriciton enyzme file
def parse_RE_file(RE_file):
	REs = {}
	with open(RE_file) as file_one:
		for line in file_one:
			enzyme, seq = line.rstrip().split("\t")
			#REs[enzyme] = Restriction_Enzyme(enzyme, seq.strip("N")) # remove leading and ending Ns
			REs[enzyme] = seq.strip("N")
	return REs

def find_substring(substring, string): # find all the starting index of a substring
	substring = substring.lower()
	string = string.lower()
	return [m.start() for m in re.finditer(substring, string)]

def seq2pattern(seq):
	iupac = {
		"B": "[CGT]",
		"D": "[AGT]",
		"H": "[ACT]",
		"K": "[GT]",
		"M": "[AC]",
		"N": "[ACGT]",
		"R": "[AG]",
		"S": "[CG]",
		"V": "[ACG]",
		"W": "[AT]",
		"Y": "[CT]"
	}
	seq = seq.upper()
	seq2 = ""
	for i in seq:
		if i in "ATGC":
			seq2 += i
		else:
			seq2 += iupac[i]
	# r'(?=([atgc]gg))' for searching overlappig patterns
	return '(?=(' +  seq2 + '))'

# test whether an enzyme can be modifed to fit a dCAPS
def test_enzyme(enzyme_name, enzyme_seq, seq):
	enzyme = Restriction_Enzyme(enzyme_name, enzyme_seq)
	enzyme_seq_RC = ReverseComplement(enzyme_seq) # 1
	wild_seq = seq.lower()
	wild_allpos = find_substring(seq2pattern(enzyme_seq), wild_seq)
	wild_allpos += find_substring(seq2pattern(enzyme_seq_RC), wild_seq) # also check reverse complement sequences of enzyme
	enzyme.allpos = list(set(wild_allpos))
	if enzyme.allpos:
		enzyme.caps = "Yes"
	return enzyme

# calculate GC content of a sequence
def Calc_GC(seq):
	t = 0.0 # float
	for a in seq:
		if a in 'GCgc':
			t += 1
	return t / len(seq) * 100

# find PAM positions
def find_pam(seq, pam, pam_pos, grna_length, direction):
	# seq is the template, direction is the seq orientation
	# pam is the pam seqeuence
	# pam_pos is the postion of the pam, left or right
	allpos = find_substring(seq2pattern(pam), seq)
	grna_list = []
	for i in allpos:
		grna = gRNA()
		grna.name = "gRNA_" + str(i) + "_" + direction
		grna.direction = direction
		grna.length = grna_length
		if pam_pos == "right":
			grna.start = i - grna_length
			grna.end = i - 1
			#seq4score_start = i - grna_length - 4
			#seq4score_end = i + 5
		else:
			grna.start = i + len(pam)
			grna.end = i + len(pam) + grna_length - 1
		if grna.start < 0 or grna.end > len(seq):
		#if seq4score_start < 0 or seq4score_end >= len(seq):
			continue
		grna.seq = seq[grna.start:(grna.end + 1)].upper()
		grna.gc = Calc_GC(grna.seq)
		grna.gc10 = Calc_GC(grna.seq[-10:])
		grna.pam = seq[i:(i+len(pam))].upper()
		grna.pam_pos = pam_pos
		if pam_pos == "left":
			grna.forblast = grna.pam + grna.seq
		else:
			grna.forblast = grna.seq + grna.pam
		#grna.seq4score = seq[seq4score_start:(seq4score_end + 1)].upper()
		grna_list.append(grna)
	return grna_list

# test specificity of a gRNA
def test_spec(seq, targets, non_targets, fasta_raw):
	seq = seq.lower()
	spec = 1 # whether specific
	# seq should be in all targets
	for j in targets:
		tt = fasta_raw[j].lower()
		#if seq not in tt:
		if not find_substring(seq2pattern(seq), tt):
			spec = 0
			break
	if not spec:
		return spec
	# seq should NOT be in any targets
	# for non-target, I think I only need to test the last 12 bp
	# with confirm with Yanpeng
	for k in non_targets:
		tt2 = fasta_raw[k].lower()
		#if seq in tt2:
		if find_substring(seq2pattern(seq), tt2):
			spec = 0
			break
	return spec


# test specificity of a gRNA
# based on alignment
def test_spec2(seq, targets, non_targets, fasta, template, a2t, t2a): # template is the one that was used to search potential gRNAs
	seq = seq.lower()
	template_start_pos = find_substring(seq2pattern(seq), template)[0]
	template_end_pos = template_start_pos + len(seq) - 1
	alignment_end_pos = t2a[template_end_pos]
	spec = 1 # whether specific
	# seq should be in all targets
	for j in targets:
		tt = fasta[j].lower().replace("-","")
		#if seq not in tt:
		if not find_substring(seq2pattern(seq), tt):
			spec = 0
			break
	if not spec:
		return spec
	# seq should NOT be in any non-targets
	# for non-target, I think I only need to test the last 12 bp
	# with confirm with Yanpeng
	for k in non_targets:
		tt2 = fasta[k].lower().replace("-","")
		#if seq in tt2:
		if find_substring(seq2pattern(seq), tt2):
			spec = 0
			break
	return spec

# on target score calculation based on machine learning trained model
# Doench et al. 2006. Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9
def get_on_target_score(seq): #30mer nucleotide sequence, which should be of the form NNNN20merNGGNNN
	current_path = os.path.dirname(os.path.realpath(__file__))
	if not seq:
		return 0.0
	model_file = current_path + '/saved_models/V3_model_nopos.pickle'
	seq = seq.upper()
	if len(seq)!=30: 
		print "Please enter a 30mer sequence."
		sys.exit(1)
	aa_cut = -1
	per_peptide = -1
	with open(model_file, 'rb') as f:
		model = pickle.load(f)
	if seq[25:27] == 'GG':
		score = model_comparison.predict(seq, aa_cut, per_peptide, model=model)
	return score

# off target score
def get_off_target_score(wt, off):
	# 20mer sgRNA + PAM sequence for both wt and off sequences
	wt = wt.upper()
	off = off.upper()
	m_wt = re.search('[^ATCG]',wt)
	m_off = re.search('[^ATCG]',off)
	if (m_wt is None) and (m_off is None):
		pam = off[-2:]
		sg = off[:-3]
		cfd_score = cfd.calc_cfd(wt,sg,pam)
		return cfd_score

#print get_off_target_score("CCAGGATGGGGCATTTCTAGAGG", "CCAGGATGGGGCATTTCTAAAGG")
#print get_off_target_score("CCAGGATGGGGCATTTCTAGAGG", "CCAGGATGGGGCATTTCTAAAGT")
#print get_off_target_score("CCAGGATGGCGCATTTCTAGAGG", "CCAGGATGGGGCATTTCTAAAGG")

# Get software path
#from sys import platform
def get_software_path(base_path):
	if sys.platform.startswith('linux'): # linux
		#primer3_path = base_path + "/primer3_core"
		muscle_path = base_path + "/muscle"
	elif sys.platform == "win32" or sys.platform == "cygwin": # Windows...
		#primer3_path = base_path + "/primer3_core.exe"
		muscle_path = base_path + "/muscle.exe"
	elif sys.platform == "darwin": # MacOSX
		#primer3_path = base_path + "/primer3_core_darwin64"
		muscle_path = base_path + "/muscle3.8.31_i86darwin64"
	#return primer3_path, muscle_path
	return muscle_path

# function to blast and parse the output of each primer in the wheat genome
def blast_check(primer_list, reference, blastoutfile, pam):
	# primer_list is a list of primer objects
	# reference is the blastdb
	forblast = open("for_blast.fa", 'w') # for blast against the gnome
	primer_for_blast = {} # dict of all primers with name and seq
	for pp in primer_list:
		#if pp.REs: # only blast those with retriction enzymes
		if pp.pam_pos == "left":
			forblast.write(">" + pp.name + "\n" + pp.seq + pp.pam + "\n")
		else:
			forblast.write(">" + pp.name + "\n"  + pp.pam + pp.seq + "\n")
		primer_for_blast[pp.name] = pp
	forblast.close()
	### for blast
	num_threads = 4
	#cmd2 = 'blastn -task blastn -db ' + reference + ' -query for_blast.fa -outfmt "6 std qseq sseq qlen slen" -num_threads 3 -word_size 7 -out blast_out.txt'
	#cmd2 = "blastn -task blastn-short -db " + reference + " -query for_blast.fa -evalue 30000 -word_size 7 -gapopen 2 -gapextend 1 -reward 1 -penalty -1 -perc_identity 70 -max_target_seqs 13 -max_hsps 20 -num_threads " + str(num_threads) + " -dust no  -outfmt '6 std qseq sseq qlen slen' -out blast_out.txt"
	cmd2 = "blastn -task blastn-short -db " + reference + " -query for_blast.fa -ungapped -perc_identity 70 -word_size 8 -max_target_seqs 10 -max_hsps 2 -num_threads " + str(num_threads) + " -outfmt '6 std qseq sseq qlen slen' -out " + blastoutfile
	print "Step 2: Blast command:\n", cmd2
	call(cmd2, shell=True)
	# blast fields
	# 1: query id, subject id, % identity, alignment length, mismatches, gap opens, 
	# 7: q. start, q. end, s. start, s. end, evalue, bit score
	# 13: q. sequence, s. sequence, q. length s. length
	sortcmd = "sort -k1,1 -k12,12nr " + blastoutfile + " | awk 'a[$1]++ < 4' > sorted_blast_out_top4.txt"
	call(sortcmd, shell=True)
	for line in open("sorted_blast_out_top4.txt"):
		if line.startswith("#"):
			continue
		fields = line.split("\t")
		#query, subject, pct_identity, align_length, mismatches = fields[:5]
		#qstart, qstop, sstart, sstop = fields[6:10]
		query, subject = fields[:2]
		if query in primer_for_blast:
			primer = primer_for_blast[query]
		else:
			continue
		align_length, mismatches, ngap, qstart, qstop, sstart, sstop = [int(x) for x in fields[3:10]]
		qlen = int(fields[14]) # query length
		if sstart > sstop: # plus:minus
			sleft = sstart + (qstart - 1)
			sright = sstop - (qlen - qstop)
		else: # plus:plus
			sleft = sstart - (qstart - 1)
			sright = sstop + (qlen - qstop)
		if align_length == qlen and mismatches == 0:
			off_target_score = 1
		else:
			extended_subject_algn = extractseq(reference, subject, sleft, sright)
			#print query
			#print primer.seq + primer.pam
			#print extended_subject_algn
			off_target_score = get_off_target_score(primer.seq + primer.pam, extended_subject_algn)
		#primer.blast += subject + ": " + qstart + "-" + qstop + ", " + mismatches + "; "
		primer.blast += subject + ": " + "{0:.1f}".format(off_target_score) + "; "

	return(primer_list) # I think the primers should be changed

## function to extract sequences from the refernce
def extractseq(reference, chromosome, start, end):
	strand = "plus"
	if start > end:
		strand = "plus"
		start, end = end, start
	cmd = "blastdbcmd -db " + reference + " -entry " + chromosome + " -range " + str(start) + "-" + str(end) + " -strand " + strand
	output = check_output(cmd, shell=True)
	#>chr7D 
	#AGGGTTTAGGG
	seq = output.splitlines()[1] # only the 2nd line
	return seq

# prepare blast file
# function to blast and parse the output of each primer in the wheat genome
def prepare_blast_file(primer_list):
	# primer_list is a list of primer objects
	# reference is the blastdb
	forblast = open("for_blast.fa", 'w') # for blast against the gnome
	primer_for_blast = {} 
	for pp in primer_list:
		forblast.write(">" + pp.name + "\n" + pp.forblast + "\n")
		primer_for_blast[pp.name] = pp
	forblast.close()
	return primer_for_blast

## find mismatches
def off_target_check(query, reference, mybatmap):
	# default is 4 mismatches
	cmd = mybatmap + " -q " + query + " -d " + reference + " -o out-test-whole.txt -n4 -mall"
	call(cmd, shell=True)
	return 0

## check gene names
def check_gene(infile, gff):
	outfile="for.gene.check.bed"
	out = open(outfile, 'w')
	chrloc = {} # dictionary of chromosome locations
	with open(infile) as file_one:
		for line in file_one:
			li=line.strip()
			ll = li.split("\t")
			#chrloc[":".join([ll[1], ll[3]])] = "" # chr1A:12345 as key
			end = int(ll[3]) + len(ll[5])
			out.write('\t'.join([ll[1], ll[3], str(end)]) + "\n")
	out.close()
	# bedtools interact
	bedfile = "sorted." + outfile
	cmd1 = "bedtools sort -i " + outfile + " > " + bedfile
	print cmd1
	call(cmd1, shell=True)
	checkfile = "check.result.txt"
	cmd2 = "bedtools intersect  -sorted -wao -a " + bedfile + " -b " + gff + " > " + checkfile
	print cmd2
	call(cmd2, shell=True)
	# process the check results
	with open(checkfile) as file_two:
		for line in file_two:
			li=line.strip()
			ll = li.split("\t")
			kk = ":".join(ll[0:2]) # key
			if kk not in chrloc or ll[5] != "gene":
				chrloc[kk] = "; ".join([ll[5], ll[11], ll[12]])
	# return check information
	return chrloc

## parse batmis output for mismatch search
def parse_mismatches(infile, pam, pam_pos, grna_dict, gff):
	# 1. reformat the output
	# 2. filter out these without intact pam
	grna_name = ""
	outfile="out.temp.txt"
	out = open(outfile, 'w')
	with open(infile) as file_one:
		for line in file_one:
			li=line.strip()
			if li.startswith("@"):
				continue
			if li.startswith(">"):
				li = li.lstrip(">")
				grna_name, grna_seq = li.split("\t")[0:2]
				grna_seq = grna_seq.lower()
				continue
			cols = li.split("\t")
			chrom, strand, pos, nmis = cols[1:5]
			#out.write("\t".join([grna_name] + cols[1:5]))
			if int(nmis) == 0:
				out.write("\t".join([grna_name] + cols[1:5]))
				out.write("\t" + grna_seq + "\n")
			else:
				ll = list(grna_seq) # conver to list for replace nth nt
				isindel = 0
				for i in cols[7:]:
					pos, nt = i.split(">")
					if nt == nt.upper() or nt == "D": # insertion or deletion
						break
						isindel = 1
					ll[int(pos)] = nt.upper()
				if isindel:
					continue
				mism = "".join(ll) # convert to string
				# check whether intact pam
				if pam_pos == "left":
					current_pam = mism[0:len(pam)]
				else:
					current_pam = mism[-len(pam):]
				current_pam
				allpos = find_substring(seq2pattern(pam), current_pam)
				if allpos: # if intact pam
					out.write("\t".join([grna_name] + cols[1:5]))
					out.write("\t" + mism + "\n")
	out.close()
	# extract top 10 mismatches of each gRNAs
	outfile2 = "sorted." + outfile
	#out2 = open(outfile2, 'w')
	#out2.write("\t".join(["gRNA", "Chromosome", "Strand", "Position", "Mismatches", "Potential_target"]) + "\n")
	#out2.close()
	cmd1 = "sort -k1,1 -k5,5n " + outfile + " | awk 'a[$1]++ < 10' > " + outfile2
	call(cmd1, shell=True)
	
	# get locations on genes
	genecheck = check_gene(outfile2, gff)
	
	# get the list of off targets for gRNA objects
	with open(outfile2) as file_two:
		for line in file_two:
			line = line.strip()
			ll = line.split("\t")
			grna_name = ll[0]
			genename = genecheck[ll[1] + ":" + ll[3]]
			if genename.startswith("gene"):
				genename = genename.replace("gene", "intron")
			elif genename.startswith("."):
				genename = "intergenic"
			pp = grna_dict[grna_name]
			pp.blast +=  line + "\t" + genename + "\n" + "\t" * 12
	return 0

## test parse_mismatches
#temp = parse_mismatches("../mismatches/out-test-whole.txt", "NGG", "right")
