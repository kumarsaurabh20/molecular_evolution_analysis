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

def check_folder(folder, remove=False):
	
	print("Checking the folder {}...".format(folder))
	if os.path.isdir(folder):
		pass
	elif os.path.isdir(folder) == False and remove == False:
		print("{} folder does not exist!".format(folder))
		os.mkdir(folder)
		print("{} is created!".format(folder))
	else:
		print("{} folder does not exist!".format(folder))
		sys.exit()


def prep_scp_file(full_ortho_file, scp_file):
	full_OG_df = pd.read_csv(full_ortho_file, sep="\t")
	columns = full_OG_df.columns
	scp_list_df = pd.read_csv(scp_file, sep="\t", names=["Orthogroup"], header=None)
	merged_df = pd.merge(scp_list_df, full_OG_df, on="Orthogroup", how="inner")
	return merged_df, columns


def prep_files(sco, combined_fastas, dna=False):
	cwd = os.getcwd()
	path = cwd + "/alignments"
	check_folder(path)

	#sco.columns = columns
	sco.set_index("Orthogroup", inplace=True)
	#
	if dna == True:
		for index, row in sco.iterrows():
			print("Writing file for the OG {}".format(index))
			outfile = path + "/" + index + ".m"+ ".fasta"
			records = (r for r in SeqIO.parse(combined_fastas, "fasta") if r.id in row.tolist())
			SeqIO.write(records, outfile, "fasta")	
	elif dna == False:
		for index, row in sco.iterrows():
			print("Writing file for the OG {}".format(index))
			outfile = path + "/"  + index + ".fasta"
			records = (r for r in SeqIO.parse(combined_fastas, "fasta") if r.id in row.tolist())
			SeqIO.write(records, outfile, "fasta")
			#temp = subprocess.Popen([cmd, 'combined.p.fasta', '' server], stdout = subprocess.PIPE)


def perform_alignments():

	cwd = os.getcwd()
	#script = os.path.realpath(__file__)
	path = cwd + "/alignments"
	check_folder(path)
	for filename in glob.glob(path + '/*.fasta'):
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
	work_dir = base_dir_path + "/alignments"
	check_folder(work_dir, remove=True)

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

	print("Performing codon alignments ...")
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
	

	for filename in glob.glob(fasta_file_path + "/*.m.fasta"):	
		index = os.path.basename(filename).split('.', 1)[0]
		outfile = fasta_file_path + "/" + index + ".paml.fasta"
		outfile2 = fasta_file_path + "/" + index + ".codon.fasta"
		infile1 = fasta_file_path + "/" + index + ".aln.fasta"
		infile2 = filename
		with open(outfile, "w+") as o:
			subprocess.call(["pal2nal.pl", infile1, infile2, "-output", "paml", "-nogap", "-nomismatch"], stdout = o)
		with open(outfile2, "w+") as o:
			subprocess.call(["pal2nal.pl", infile1, infile2, "-output", "fasta", "-nogap", "-nomismatch"], stdout = o)



def concatenate_file(destination_folder="", outfilename="", pattern="", path=None):
	
	print("Concatenating {}s ...".format(filetype))
	script = os.path.realpath(__file__)
	base_dir_path = os.path.dirname(script)
	#
	file_path = ""
	if path != None:
		file_path = path
	else:
		file_path = base_dir_path + "/" + destination_folder
	filename = file_path + "/" + outfilename
	pattern = file_path + "/" + pattern
	#
	files = sorted(glob.glob(pattern))
	with open(filename, 'w') as finalout:
		for f in files:
			with open(f) as infile:
				for line in infile:
					finalout.write(line)

def generate_trees():

	print("Generating trees from the alignments ...")
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
		if os.stat(filename).st_size == 0:
			pass
		else:
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

def main():
	scp_df = prep_scp_file("Orthogroups.tsv", "Orthogroups_SingleCopyOrthologues.txt")	
	prep_files(scp_df, "combined.p.fasta", dna=False)
	perform_alignments()
	trim_alignments()
	prep_files(scp_df, "combined.m.fasta", dna=True)
	codon_alignments()	
	concatenate_file(destination_folder="alignments", outfilename="combined.paml.fasta", pattern="*.paml.fasta", path=None)
	generate_trees()
	concatenate_file(destination_folder="trees", outfilename="combined.tree.nwk", pattern="*.treefile", path=None)
	return 

if __name__ == "__main__":

	#cwd = os.getcwd()
        #script = os.path.realpath(__file__)
        #path = cwd + "/" + "Data"
	#OF_dir = "OrthoFinder/Results_Jan13/Orthogroups"	
	#file_path = cwd + path + "/" OF_dir

	main()
