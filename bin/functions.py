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

import getopt, sys, os, re, copy

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
		self.seq = "" # sequence without PAM
		self.pam = "" # PAM
		self.direction = ""
		self.template = ""
		self.REs = [] # list of restriction enzymes

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
	return seq2.lower()

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
def find_pam(seq, pam, grna_length, direction):
	allpos = find_substring(seq2pattern(pam), seq)
	grna_list = []
	for i in allpos:
		grna = gRNA()
		grna.direction = direction
		grna.length = grna_length
		grna.start = i - grna_length
		grna.end = i - 1
		if grna.start < 0:
			continue
		grna.seq = seq[grna.start:i].upper()
		grna.gc = Calc_GC(grna.seq)
		grna.pam = seq[i:(i+len(pam))].upper()
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
