#!/usr/bin/python3
"""calculates the cytosine and quanine content of a sequence"""


import sys
from Bio import SeqRecord
from Bio import Seq


def GC_content(input_file): 
    """calculates the cytosine and quanine content of a sequence"""
    try:
        with open(input_file, "r") as in_handle:
            # Parse input FASTA file
            sequences = SeqIO.parse(in_handle, "fasta")
            GCcount = 0
            total_bases = 0


        for seq_record in sequences:
            sequence = seq_record.seq
            total_bases += len(sequence)
                # Count the number of 'G' and 'C' nucleotides
            GCcount += sequence.count('G') + sequence.count('C')

        return GCcount, total_bases
    except Exception as e:
        print("Error:", str(e))
        return None

if __name__ == "__main":
    if len(sys.argv) != 2:
        print("Usage: python GC_content.py input.fasta")
    else:
        input_file = sys.argv[1]
        GCcount, total_bases = GC_content(input_file)
        if GCcount is not None:
            GC_percent = (GCcount / total_bases) * 100
            print(f"GC content: {GC_percent:.2f}%")
            