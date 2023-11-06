#!/usr/bin/python3
""" translates a DNA strand to mRNA and mRNA to protein"""


from Bio.Seq import Seq
import sys
from Bio import SeqIO
from Bio import ExPASy
from Bio import Entrez
import os
import sqlite3


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
                    return str(or_sequence)
                except Exception as e:
                    print("Error during online query:", str(e))
            else:
                #if user request data from other databasesa
                try:
       
                    with Entrez.efetch(db=db, rettype="fasta", retmode="text", id=accession_code) as handle:
                        seq_record = SeqIO.read(handle, "fasta")
                        or_sequence = print(seq_record.id, seq_record.name, seq_record.description, seq_record.seq)
                    return str(or_sequence)
                except Exception as e:
                    print("Error during online query:", str(e))
                    return None                     
def complement(sequence):
    seq_complement = sequence.complement()
    #print("sequence complement:", seq_complement)
    return seq_complement

def reverse_complement(sequence):
    # Reverse complement the DNA sequence
    reverse_complement_dna = sequence.reverse_complement()
    #print("Reverse Complement:", template_dna)
    return reverse_complement_dna

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
    # Ensure the mRNA sequence is a multiple of three
    codon_remainder = len(mRNA) % 3
    if codon_remainder > 0:
        # Pad the sequence with 'N' nucleotides
        mRNA += "N" * (3 - codon_remainder)

    try:
        # Translate the mRNA to a protein sequence
        protein_seq_mRNA = mRNA.translate(table=1, to_stop=True)
        return str(protein_seq_mRNA)
    except Exception as e:
        print(f"Translation Error: {str(e)}")
        return "Translation Error"

def translate_dna(sequence):
    # Ensure the mRNA sequence is a multiple of three
    codon_remainder = len(sequence) % 3
    if codon_remainder > 0:
        # Pad the sequence with 'N' nucleotides
        sequence += "N" * (3 - codon_remainder)

    try:
        # Translate the mRNA to a protein sequence
        protein_seq_dna = sequence.translate(table=1, to_stop=True)
        return str(protein_seq_dna)
    except Exception as e:
        print(f"Translation Error: {str(e)}")
        return "Translation Error"



def save_sequence_to_file(sequence, accession_code, filename):
    """save sequence ouput to files with
    the same accession codes as the original file"""
    #>accession_code
    header = f">{accession_code}\n"
    #sequence
    sequence = str(sequence)
    #write them in file
    with open(filename, "w") as file:
        file.write(header)
        file.write(sequence)
# database logic
def create_database(database_name):
    """function that creates a database in sqlite"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()

    conn.commit()
    conn.close()



def create_tables(database_name):
    """function that creates tables in database"""

    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()

    #table that stores original sequences
    cursor.execute('''
        
        CREATE TABLE IF NOT EXISTS sequences (
            accession_code TEXT PRIMARY KEY NOT NULL,
            seq_record TEXT NOT NULL
            )
    ''')


     #create table to store the modified sequences
    cursor.execute ('''
    CREATE TABLE IF NOT EXISTS modifed_sequences (
            accession_code TEXT PRIMARY KEY NOT NULL,
            seq_record TEXT NOT NULL,
            seq_compliment   TEXT NOT NULL,
            seq_reverse_compliment TEXT NOT NULL,
            mRNA_seq TEXT NOT NULL,
            Protein_seq NOT NULL
    )
    ''')


    #table that stores user logins and will be usrd to reset passwords
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS users (
            user_email TEXT PRIMARY KEY,
            user_password TEXT NOT NULL,
            user_name TEXT
    )
    ''')
    conn.commit()
    conn.close()


def insert_or_update_sequences(database_name, accession_code, sequence, table_name):
    """function to query databases"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()

    cursor.execute( f'''
    INSERT OR REPLACE INTO sequences (accession_code, seq_record)
    VALUES (?, ?)
    ''',(accession_code, sequence)
    )

def insert_or_update_modified_sequences(database_name, accession_code, seq_complement, reverse_complement_dna, mRNA, Protein_seq_mRNA ):
    """function to query databases"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()

    cursor.execute(f'''
    INSERT OR REPLACE INTO modifed_sequences (database_name, accession_code, seq_record, seq_compliment, seq_reverse_compliment, mRNA_seq, protein_seq)
    VALUES (?, ?, ?, ?, ?, ?)
    ''', ('sequence_data.db',accession_code, str(seq_complement), sequence, str(reverse_complement_dna), str(mRNA), str(Protein_seq_mRNA)
    )
    )
    conn.commit()
    conn.close()

def query_database(database_name, accession_code,table_name):
    """querys the database to get sequences bassed on the accession_code"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()
    query = f'SELECT seq_record FROM {table_name} WHERE accession_code = ?'
    cursor.execute (query, (accession_code,))
    results = cursor.fetchone()
    conn.close()
    return results[0] if results else None

if __name__ == "__main__":
    #create the database
    create_database('sequences_data.db')
    create_tables('sequences_data.db')

    # Check if a file is provided as a command-line argument
    file_path = input("Enter the path to a local file (or press Enter to skip this option): ")
    db = None

    if file_path and (file_path.lower().endswith(('fasta', 'fa', 'genbank', 'gbk')) or os.path.isfile(file_path)):
        # Read sequence from the provided file
        sequence = read_sequence_from_file(file_path)
        accession_code = sequence.id

        #insert the sequence into database
        insert_or_update_sequences('sequences_data.db', accession_code, str(sequence), 'sequences')
    else:
        sequence = None

    if sequence is not None:
        # Perform central dogma operations on the provided sequence
        #print("Original Sequence:", sequence)
        # Perform the rest of the operations
        seq_complement = complement(sequence)
        reverse_complement_dna = reverse_complement(sequence)
        mRNA = transcribe(sequence)
        coding_dna = reverse_transcribe(mRNA)
        protein_seq_mRNA = translate_mRNA(mRNA)
        protein_seq_dna = translate_dna(sequence)
        

        # Save the sequences to files
        save_sequence_to_file(sequence, accession_code, "original.fasta")
        save_sequence_to_file(seq_complement, accession_code, "complement.fasta")
        save_sequence_to_file(reverse_complement_dna, accession_code, "reverse_complement_dna.fasta")
        save_sequence_to_file(mRNA, accession_code, "mRNA.fasta")
        save_sequence_to_file(coding_dna, accession_code, "coding_dna.fasta")
        save_sequence_to_file(protein_seq_mRNA, accession_code, "protein_seq_mRNA.fasta")
        save_sequence_to_file(protein_seq_dna, accession_code, "protein_seq_dna.fasta")
        insert_or_update_modified_sequences('sequences_data.db', accession_code, sequence, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)

    else:
        # No valid file provided, ask the user for database and ID
        db = input("Enter the database (e.g., 'nucleotide'): ")
        accession_code = input("Enter the ID of the record: ")

        if db and accession_code:
            sequence = read_sequence_from_database(db, accession_code)
            if sequence:
                # Convert the retrieved sequence to a BioPython Seq object
                sequence = Seq(sequence)
                insert_or_update_sequences('sequences_data.db', accession_code, str(sequence), 'sequences')
                # Perform central dogma operations on the retrieved sequence
                print("Original Sequence:", sequence)
                # Perform the rest of the operations
                seq_complement = complement(sequence)
                reverse_complement_dna = reverse_complement(sequence)
                mRNA = transcribe(sequence)
                coding_dna = reverse_transcribe(mRNA)
                protein_seq_mRNA = translate_mRNA(mRNA)
                #protein_seq_dna = translate_dna(DNA)

                # Save the sequences to files for user download
                save_sequence_to_file(sequence,accession_code, "sequence")
                save_sequence_to_file(seq_complement, accession_code, "complement.fasta")
                save_sequence_to_file(reverse_complement_dna, accession_code, "reverse_complement_dna.fasta")
                save_sequence_to_file(mRNA, accession_code, "mRNA.fasta")
                save_sequence_to_file(coding_dna, accession_code, "coding_dna.fasta")
                save_sequence_to_file(protein_seq_mRNA, accession_code, "protein_seq_mRNA.fasta")
                #save_sequence_to_file(protein_seq_dna, accession_code, "protein_seq_dna.fasta")
                insert_or_update_modified_sequences('sequences_data.db', accession_code, sequence, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)

            else:
                print("Error retrieving sequence from the database.")
        else:
            print("Invalid input. Database query requires both 'db' and 'id'.")

    

