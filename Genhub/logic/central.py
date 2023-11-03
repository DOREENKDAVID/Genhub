#!/usr/bin/python3
""" translates a DNA strand to mRNA and mRNA to protein"""
"""
1.users can upload  fasta file or genbank
    if file has one seq use  SeqIO.read_sequence_from_file
    if file has many sequnces use SeqIO.parse
2.second method 
users can use databases to get sequence
"""


from Bio.Seq import Seq
import sys
from Bio import SeqIO
from Bio import ExPASy
from Bio import Entrez
import os


def read_sequence_from_file(file_path):
    if os.path.isfile(file_path):
            try:
                with open(file_path, "r") as file:
                    # Check if the file is in FASTA format
                    if file_path.endswith(".fasta") or file_path.endswith(".fa"):
                    # Parse FASTA file
                        sequences = list(SeqIO.parse(file, "fasta"))
                        if len(sequences) == 1:
                            return sequences[0].seq
                    # Assuming there's many sequences in the file
                        else:
                            for seq_record in sequences:
                                #(seq_record.id,seq_record.name, seq_record.description, seq_record.seq)
                                 return seq_record
                                 
                    else:
                        # Read the file as plain text of one sequence
                        return Seq(file.read())
            except Exception as e:
                print("Error:", str(e))
                #return None
def read_sequence_from_database(db,  accession_code):
            """read sequences from a database
                args:
                        db: database name
                        acession_code: id of the sequence"""
            #default email and search
            Entrez.email = "A.N.Other@example.com"
            #if user requests data from swisprot db
            if db.lower() == 'swissprot':
                try:
                    with ExPASy.get_sprot_raw(accession_code) as handle:
                        seq_record = SeqIO.read(handle, "swiss")
                    or_sequence = print(seq_record.id,seq_record.name, seq_record.description, seq_record.seq)
                    return str(seq_record.seq)
                except Exception as e:
                    print("Error during online query:", str(e))
            else:
                #if user request data from other databasesa
                try:
       
                    with Entrez.efetch(db=db, rettype="fasta", retmode="text", id=accession_code) as handle:
                        seq_record = SeqIO.read(handle, "fasta")
                    or_Sequence = print(seq_record.id, seq_record.name, seq_record.description, seq_record.seq)
                    return str(seq_record.seq)
                except Exception as e:
                    print("Error during online query:", str(e))
                    return None                     
def complement(sequence):
    seq_complement = sequence.complement()
    #print("sequence complement:", seq_complement)
    return seq_complement

def reverse_complement(sequence):
    # Reverse complement the DNA sequence
    template_dna = sequence.reverse_complement()
    #print("Reverse Complement:", template_dna)
    return template_dna

def transcribe(sequence):
    # Transcribe the DNA sequence to mRNA
    mRNA = sequence.transcribe()
    #print("mRNA:", mRNA)
    return mRNA

def reverse_transcribe(mRNA):
    # Back-transcribe the mRNA to DNA
    coding_dna = mRNA.back_transcribe()
    #print("Back-transcribed DNA:", coding_dna)
    return coding_dna

def translate_mRNA(mRNA):
    # Translate the mRNA to a protein sequence
    #standard table
    table = 1
    # Extract the sequence data from mRNA (remove any metadata)
    sequence_data = mRNA[mRNA.find("\n") + 1:].replace("\n", "")
    codon_remainder = len(mRNA) % 3
    if codon_remainder > 0:
        mRNA += "N" * (3 - codon_remainder)
    protein_seq_mRNA = mRNA.translate(table=table)
    #print("Protein Sequence (from mRNA):", protein_seq_mRNA)
    return protein_seq_mRNA

def translate_dna(sequence):
    # Translate the DNA sequence to a protein sequence
    table = 1
    codon_remainder = len(sequence) % 3
    if codon_remainder > 0:
        # Add "N" nucleotides to the end to make it a multiple of three
        sequence += "N" * (3 - codon_remainder)

    protein_seq_dna = sequence.translate(table=table)
    #print("Protein Sequence (from DNA):", protein_seq_dna)
    return protein_seq_dna

def save_sequence_to_file(sequence, accession_code, filename):
    """save sequence ouput to files with
    the same accession codes as the original file"""
    header = f">{accession_code}\n"
    sequence = str(sequence)
    with open(filename, "w") as file:
        file.write(header)
        file.write(sequence)




if __name__ == "__main__":

    # Check if a file is provided as a command-line argument
    file_path = input("Enter the path to a local file (or press Enter to skip this option): ")
    db = None

    if file_path and (file_path.lower().endswith(('fasta', 'fa', 'genbank', 'gbk')) or os.path.isfile(file_path)):
        # Read sequence from the provided file
        sequence = read_sequence_from_file(file_path)
    else:
        sequence = None

    if sequence is not None:
        # Perform central dogma operations on the provided sequence
        #print("Original Sequence:", sequence)
        # Perform the rest of the operations
        seq_complement = complement(sequence)
        template_dna = reverse_complement(sequence)
        mRNA = transcribe(sequence)
        coding_dna = reverse_transcribe(mRNA)
        protein_seq_mRNA = translate_mRNA(mRNA)
        protein_seq_dna = translate_dna(sequence)

        # Save the sequences to files
        save_sequence_to_file(seq_record, accession_code, "original.fasta")
        save_sequence_to_file(seq_complement, accession_code, "complement.fasta")
        save_sequence_to_file(template_dna, accession_code, "template_dna.fasta")
        save_sequence_to_file(mRNA, accession_code, "mRNA.fasta")
        save_sequence_to_file(coding_dna, accession_code, "coding_dna.fasta")
        save_sequence_to_file(protein_seq_mRNA, accession_code, "protein_seq_mRNA.fasta")
        save_sequence_to_file(protein_seq_dna, accession_code, "protein_seq_dna.fasta")

    else:
        # No valid file provided, ask the user for database and ID
        db = input("Enter the database (e.g., 'nucleotide'): ")
        accession_code = input("Enter the ID of the record: ")

        if db and accession_code:
            sequence = read_sequence_from_database(db, accession_code)
            if sequence:
                # Convert the retrieved sequence to a BioPython Seq object
                sequence = Seq(sequence)
                # Perform central dogma operations on the retrieved sequence
                print("Original Sequence:", sequence)
                # Perform the rest of the operations
                seq_complement = complement(sequence)
                template_dna = reverse_complement(sequence)
                mRNA = transcribe(sequence)
                coding_dna = reverse_transcribe(mRNA)
                protein_seq_mRNA = translate_mRNA(mRNA)
                protein_seq_dna = translate_dna(sequence)

                # Save the sequences to files
                save_sequence_to_file(seq_complement, accession_code, "complement.fasta")
                save_sequence_to_file(template_dna, accession_code, "template_dna.fasta")
                save_sequence_to_file(mRNA, accession_code, "mRNA.fasta")
                save_sequence_to_file(coding_dna, accession_code, "coding_dna.fasta")
                save_sequence_to_file(protein_seq_mRNA, accession_code, "protein_seq_mRNA.fasta")
                save_sequence_to_file(protein_seq_dna, accession_code, "protein_seq_dna.fasta")

            else:
                print("Error retrieving sequence from the database.")
        else:
            print("Invalid input. Database query requires both 'db' and 'id'.")

    