#!/usr/bin/env python

import re, os, sys, logging, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from mimetypes import guess_type
from functools import partial
import gzip
from mimetypes import guess_type

#checking if the transcript pattern is supplied by the user and making
#regular expression objects of the pattern
def get_longest_transcript(input_file,output,gene_start,search_pattern="\.[0-9]+",replace_pattern="",compressed=False):
    search_pattern = re.compile(search_pattern)
    
    if replace_pattern == None:
      replace_pattern=""

    '''
    Parsing the gene start pattern so that we can filter out unwanted genes
    '''
    gene_pattern = re.compile(gene_start)

    #object to store the gene sequences
    seqs = {}

    '''
    Looping through to read and store the longest transscript sequence for each gene
    for each iteration regex replace the trans id to get gene id
    if length is longer than existing sequence replace or add the gene id sequence to dict
    '''
    encoding = guess_type(input_file)[1]  # uses file extension
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    
    with _open(input_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            gene_id = re.sub(search_pattern,replace_pattern,record.id)
            if gene_pattern.match(gene_id) is None:
                continue
            if gene_id in seqs:
                if len(seqs[gene_id].seq) < len(record.seq):
                    seqs[gene_id] = record
            else:
                seqs[gene_id] = record

    '''
    This creates a list of sequences which can be saved into a file
    '''
    out_seqs = []
    for key in seqs.keys():
        curr_seq = seqs[key]
        curr_seq.id = key +" "+curr_seq.id
        curr_seq.seq = Seq(re.sub(r"[^A-Za-z]","",str(curr_seq.seq)))
        out_seqs.append(curr_seq)

    #Write the output file as fasta
    if (os.path.splitext(output)[1] != ".gz") & compressed:
        output = output + ".gz"
    
    SeqIO.write(out_seqs,output,"fasta")


'''
Parse arguments using argumentparser
'''

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",help="Input sequence fasta file. Can be DNA or AA sequences",required=True)
    parser.add_argument("-o","--output",help="Output file. the file will be saved as fasta",required=True)
    parser.add_argument("-z","--compressed",help="Compress the output file",action='store_true')
    parser.add_argument("-s","--gene_start",help="Pattern seen for gene IDs from the start",required=True)
    parser.add_argument("-p","--search-pattern",help="String/RegEx Pattern that would recognize individual transcripts of a gene")
    parser.add_argument("-r","--replace-pattern",help="String/RegEx Pattern that would replace the individual transcripts of a gene")
    trans_args = parser.parse_args()
    print(trans_args)
    get_longest_transcript(trans_args.input,trans_args.output,trans_args.gene_start,trans_args.search_pattern,trans_args.replace_pattern,  compressed=trans_args.compressed)
