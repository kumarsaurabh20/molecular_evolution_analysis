#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import pandas as pd
import numpy as np
import os
import glob
import urllib
import gzip
import subprocess

'''
This script uses a genbank file in a gbff format and generates an output file with matching protein_id and mRNA_id in two columns.
Using this script it is also possible to create two fasta files with peptides and mRNAs separately with a common header.
This format is useful for PAML analysis.    
'''

__author__="Kumar (ks575@exeter.ac.uk)"


def get_assembly_summary(id):
	"""Get esummary for an entrez id"""
	esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
	esummary_record = Entrez.read(esummary_handle)
	return esummary_record

def get_assemblies(term="", download=True, path='assemblies', mymail=""):
	"""Download genbank assemblies for a given search term.
	Args:
	term: search term, usually organism name
	download: whether to download the results
	path: folder to save to
	"""
	#provide your own mail here
	Entrez.email = mymail
	handle = Entrez.esearch(db="assembly", term=term, retmax='200')
	
	outfile_terms =  "_".join(term.split( )) 
	record = Entrez.read(handle)
	ids = record['IdList']
	print (f'found {len(ids)} ids')
	links = []
	label = ""
	for id in ids:
		#get summary
		summary = get_assembly_summary(id)
		#get ftp link
		url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
		if url == '':
			continue
		label = os.path.basename(url)
		#get the fasta link - change this to get other formats
		link = os.path.join(url,label+'_genomic.gbff.gz')
		links.append(link)
	for each in links:
		count = links.index(each) + 1
		print("{}: {}".format(count, each))
	if download == True:
		#download link
		print("Please select a file number from the list above that you would like to download: ")
		filenumber = int(input())
		filename = links[filenumber - 1]
		label = os.path.basename(filename)
		urllib.request.urlretrieve(filename, f'{label}')
	
	return outfile_terms


def retrieve_paired_cds_and_translation(term):
	#
	combined_precords = []
	combined_mrecords = []
	script = os.path.realpath(__file__)
	base_dir_path = os.path.dirname(script)
	filename = ""
	files = sorted(glob.glob(base_dir_path + "/*.gbff.gz"))
	if len(files) > 1:
		print("There are multiple zipped file available in the folder!")
		sys.exit()
	elif len(files) == 1:
		tmpname = files[0]
		subprocess.call(["gunzip", tmpname])
		filename = os.path.splitext(tmpname)[0]
	else:
		print("No zipped files found in the folder!")
		sys.exit()
	#
	for record in SeqIO.parse(filename, "genbank"):
		for f in record.features:
			dna=""
			dnaseq = ""
			protein=""
			annotation=""
			product=""
			translation=""
			if f.type == "CDS":
				dnaseq = str(f.extract(record.seq))
				protein = f.qualifiers["protein_id"][0]
				annotation = f.qualifiers["product"][0]
				translation = f.qualifiers["translation"][0]
		#
				combined_precords.append(SeqRecord(Seq(translation),id=protein,description=annotation,))
				combined_mrecords.append(SeqRecord(Seq(dnaseq),id=protein,description=annotation,))
	#
	poutfilename = term + ".p.fasta"
	moutfilename = term + ".m.fasta"
	SeqIO.write(combined_precords, poutfilename, "fasta")
	SeqIO.write(combined_mrecords, moutfilename, "fasta")

def main():
	outfilename = get_assemblies(term="osmia bicornis", download=True, mymail="ks575@exeter.ac.uk")
	retrieve_paired_cds_and_translation(outfilename)

	return

if __name__ == "__main__":
	main()
