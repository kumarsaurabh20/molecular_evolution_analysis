#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO

##grep -f Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > single_copy_orthologues.list

def is_tool(name):
	"""Check whether `name` is on PATH and marked as executable."""
	# from whichcraft import which
	from shutil import which
	return which(name) is not None



def prep_files(orthologs_list, combined_fastas, columns=[]):

	sco = pd.read_csv(orthologs_list, sep="\t", header=None)
	sco.columns = columns
	sco.set_index("OGs", inplace=True)
	#
	for index, row in sco.iterrows():
		print("Writing file for the OG {}".format(index))
		outfile = index + ".fasta"
		records = (r for r in SeqIO.parse(combined_fastas, "fasta") if r.id in row.tolist())
		count = SeqIO.write(records, outfile, "fasta")
		#temp = subprocess.Popen([cmd, 'combined.p.fasta', '' server], stdout = subprocess.PIPE)


def perform_alignments():

	cwd = os.getcwd()
	#script = os.path.realpath(__file__)
	path = cwd + "/alignments"
	os.mkdir(path)
	for filename in glob.glob('*.fasta'):
		print("Aligning {}".format(filename))
		temp=os.path.basename(filename)
		base = os.path.splitext(temp)[0]
		outfilename = path + "/" + base + ".aln" + ".fasta"
		outfile = open(outfilename, 'w')
		subprocess.call(["mafft", "--auto", filename], stdout = outfile)
#temp = subprocess.Popen([cmd, '-c 1', server], stdout = subprocess.PIPE) 
## get the output as a string
#output = str(temp.communicate()) 
## store the output in the list
#outputlist.append(output)

def trim_alignments():
	script = os.path.realpath(__file__)
	base_dir_path = os.path.dirname(script)
	print(base_dir_path)
	work_dir = base_dir_path + "/alignments"
	print(work_dir)
	if os.path.isdir(work_dir):
		pass
	else:
		print("Alignments folder does not exist!")
		sys.exit()
	#
	if is_tool('bmge'):
		pass
	else:
		print("BMGE could not be found!")
		sys.exit()
	
	for filename in glob.glob(work_dir + "/*.aln.fasta"):
		index = os.path.basename(filename).split('.', 1)[0]
		outfile = work_dir + "/" + index + ".bmge.fasta"
		subprocess.call(["bmge", "-i", filename, "-t", "AA" , "-of", outfile])
	

def codon_alignments():
	script = os.path.realpath(__file__)
	base_dir_path = os.path.dirname(script)
	pal2nal_path = base_dir_path + "/pal2nal.v14"
	pal2nal_script = pal2nal_path + "/pal2nal.pl"
	if os.path.isdir(pal2nal_path):
		if os.path.isfile(pal2nal_script):
			pass
		else:
			print("PAL2NAL perl script does not exist!")
			sys.exit()
	else:
		print("PAL2NAL folder does not exist!")
		sys.exit()
	fasta_file_path = base_dir_path + "/alignments"
	if os.path.isdir(fasta_file_path):
		pass
	else:
		print("Alignments folder does not exist!")
		sys.exit()
	

	for filename in glob.glob(fasta_file_path + "/*.fasta"):
		index = os.path.basename(filename).split('.', 1)[0]
		outfile = fasta_file_path + "/" + index + ".paml.fasta"
		outfile2 = fasta_file_path + "/" + index + ".codon.fasta"
		infile1 = fasta_file_path + "/" + index + ".aln.fasta"
		infile2 = fasta_file_path + "/" + index + ".fasta"
		with open(outfile, "w+") as o:
			subprocess.call(["pal2nal.pl", infile1, infile2, "-output", "paml", "-nogap", "-nomismatch"], stdout = o)
		with open(outfile2, "w+") as o:
			subprocess.call(["pal2nal.pl", infile1, infile2, "-output", "fasta", "-nogap", "-nomismatch"], stdout = o)
		#
		finaloutfile = fasta_file_path + "/combined.paml.fasta"
		with open(finaloutfile, 'w+') as finalout:
			for filename in sorted(glob.glob(fasta_file_path + "/*.paml.fasta")):
				with open(filename) as infile:
					for line in infile:
						finalout.write(line)

def concatenate_file(filetype):
	script = os.path.realpath(__file__)
	base_dir_path = os.path.dirname(script)
	file_path = ""
	filename = ""
	pattern = ""
	if filetype == "paml":
		file_path = base_dir_path + "/alignments"
		filename = file_path + "/combined.paml.fasta"
		pattern = file_path + "/*.paml.fasta"
	elif filetype == "tree":
		file_path = base_dir_path + "/tree"
		filename = file_path + "/combined.tree.phy"
		pattern = file_path + "/*.treefile"
	else:
		print("Filetype is not recognized")
	#
	files = sorted(glob.glob(pattern))
	with open(filename, 'w') as finalout:
		for f in files:
			with open(f) as infile:
				for line in infile:
					finalout.write(line)

def generate_trees():
	script = os.path.realpath(__file__)
	base_dir_path = os.path.dirname(script)
	work_dir = base_dir_path + "/alignments"
	#
	if os.path.isdir(work_dir):
		pass
	else:
		print("Alignments folder does not exist!")
		sys.exit()
        #
	tree_dir = base_dir_path + "/trees"
	if os.path.isdir(tree_dir):
		pass
	else:
		os.makedirs(tree_dir)

	for filename in sorted(glob.glob(work_dir + "/*.codon.fasta")):
		index = os.path.basename(filename).split('.', 1)[0]
		outfile = tree_dir + "/" + index + ".phy"
		alignment = AlignIO.read(open(filename), "fasta")
		print("Alignment length %i" % alignment.get_alignment_length())
		count = 1
		for record in alignment:
			record.id = str(count)
			count += 1
		print(alignment)	
		AlignIO.write(alignment, outfile, "phylip")
		#
		if is_tool('iqtree'):
			pass
		else:
			print("IQtree executable was not found!")
			sys.exit()
				
		subprocess.call(["iqtree", "-s", outfile])
	#
	finaloutfile = tree_dir + "/combined.tree.phy"
	with open(finaloutfile, 'w+') as finalout:
		for filename in sorted(glob.glob(tree_dir + "/*.treefile")):
			with open(filename) as infile:
				for line in infile:
					finalout.write(line)

def main():
	#prep_files("single_copy_orthologues.list", "combined.m.fasta", ["OGs", "Ameliferra", "Lmalachurum", "Mgenalis", "Nvitripennis", "Obicornis"])
	
	concatenate_file("paml")
	return 

if __name__ == "__main__":
	main()
