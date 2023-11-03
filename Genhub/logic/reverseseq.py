#!/usr/bin/python3
"""reverse complement of a DNA sequence"""

from Bio.Seq import Seq
import sys
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def reverse_complement_sequence(input_file, output_file):
    try:
        with open(input_file, "r") as in_handle:
            # Parse input FASTA file
            sequences = SeqIO.parse(in_handle, "fasta")
            records = []

            for seq_record in sequences:
                seq = seq_record.seq
                # Calculate the reverse complement of the sequence
                reversed_seq = seq.reverse_complement()
                reversed_record = SeqRecord(reversed_seq, id=seq_record.id, description="")
                records.append(reversed_record)

            with open(output_file, "w") as out_handle:
                # Write the reversed and complemented sequences to the output FASTA file
                SeqIO.write(records, out_handle, "fasta")

        print("Reverse complement completed and saved to", output_file)
    except Exception as e:
        print("Error:", str(e))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        # Check if the correct number of command-line arguments is provided
        print("Usage: python reverse_complement.py input.fasta output.fasta")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        reverse_complement_sequence(input_file, output_file)
