#!/usr/bin/python3
"""Translates a DNA strand to mRNA and mRNA to protein"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import ExPASy
from Bio import Entrez
import os
import sqlite3
from tabulate import tabulate
# Define functions for reading sequences from files and databases
def read_sequence_from_file(file_path):
    if os.path.isfile(file_path):
        try:
            with open(file_path, "r") as file:
                if file_path.endswith((".fasta", ".fa", ".genbank", ".gbk")):
                    sequences = list(SeqIO.parse(file, "fasta"))
                    if len(sequences) == 1:
                        return sequences[0].seq
                else:
                    return Seq(file.read())
        except Exception as e:
            print("Error:", str(e))
            return None

def read_sequence_from_database(db, accession_code):
    Entrez.email = "A.N.Other@example.com"
    if db.lower() == 'swissprot':
        try:
            with ExPASy.get_sprot_raw(accession_code) as handle:
                seq_record = SeqIO.read(handle, "swiss")
                return seq_record.seq
        except Exception as e:
            print("Error during online query:", str(e))
    else:
        try:
            with Entrez.efetch(db=db, rettype="fasta", retmode="text", id=accession_code) as handle:
                seq_record = SeqIO.read(handle, "fasta")
                return seq_record.seq
        except Exception as e:
            print("Error during online query:", str(e))
            return None

# Define functions for DNA/RNA/protein operations
def complement(sequence):
    seq_complement = sequence.complement()
    return seq_complement

def reverse_complement(sequence):
    reverse_complement_dna = sequence.reverse_complement()
    return reverse_complement_dna

def transcribe(sequence):
    mRNA = sequence.transcribe()
    return mRNA

def reverse_transcribe(mRNA):
    coding_dna = mRNA.back_transcribe()
    return coding_dna

def translate_mRNA(mRNA):
    codon_remainder = len(mRNA) % 3
    if codon_remainder > 0:
        mRNA += "N" * (3 - codon_remainder)

    try:
        protein_seq_mRNA = mRNA.translate(table=1, to_stop=True)
        return str(protein_seq_mRNA)
    except Exception as e:
        print(f"Translation Error: {str(e)}")
        return "Translation Error"

def translate_dna(sequence):
    codon_remainder = len(sequence) % 3
    if codon_remainder > 0:
        sequence += "N" * (3 - codon_remainder)

    try:
        protein_seq_dna = sequence.translate(table=1, to_stop=True)
        return str(protein_seq_dna)
    except Exception as e:
        print(f"Translation Error: {str(e)}")
        return "Translation Error"

# Define a function to save a sequence to a file
def save_sequence_to_file(sequence, accession_code, filename):
    header = f">{accession_code}\n"

    sequence = str(sequence)
    with open(filename, "w") as file:
        file.write(header)
        file.write(sequence)

# Database logic
def create_database(database_name):
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()
    # Commit and close the connection to create the database file
    conn.commit()
    conn.close()

def create_tables(database_name):
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS sequences (
            accession_code TEXT PRIMARY KEY NOT NULL,
            seq_record TEXT NOT NULL
        )
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS modifed_sequences (
            accession_code TEXT PRIMARY KEY NOT NULL,
            seq_complement TEXT NOT NULL,
            seq_reverse_complement TEXT NOT NULL,
            mRNA_seq TEXT NOT NULL,
            Protein_seq TEXT NOT NULL
        )
    ''')

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
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO sequences (accession_code, seq_record)
    VALUES (?, ?)
    ''', (accession_code, sequence))
    conn.commit()
    conn.close()

def insert_or_update_modified_sequences(database_name, accession_code, seq_complement, reverse_complement_dna, mRNA, Protein_seq_mRNA):
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO modifed_sequences (accession_code, seq_complement, seq_reverse_complement, mRNA_seq, protein_seq)
    VALUES (?, ?, ?, ?, ?)
    ''', (accession_code, str(seq_complement), str(reverse_complement_dna), str(mRNA), str(Protein_seq_mRNA))
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

def query_modified_seq_table(database_name, accession_code):
    """querys the database to get sequences bassed on the accession_code"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()

    # Query the modified_sequences table
    query = f'SELECT seq_complement, seq_reverse_complement, mRNA_seq, protein_seq FROM modifed_sequences WHERE accession_code = ?'
    cursor.execute (query, (accession_code,))
    modified_results = cursor.fetchone()

    # Query the sequences table to get the original sequence
    query = f'SELECT seq_record FROM sequences WHERE accession_code = ?'
    cursor.execute (query, (accession_code,))
    original_results = cursor.fetchone()
    conn.close()
    # Combine the results, the original  and modified sequence
    if modified_results and original_results:
        seq_record = original_results[0]
        seq_complement, seq_reverse_complement, mRNA_seq, protein_seq = modified_results

        return ( seq_record, seq_complement, seq_reverse_complement, mRNA_seq, protein_seq)
    else:
        None

if __name__ == "__main__":
    create_database('sequences_data.db')
    create_tables('sequences_data.db')

    file_path = input("Enter the path to a local file (or press Enter to skip this option): ")
    db = None

    if file_path and (file_path.lower().endswith(('.fasta', '.fa', '.genbank', '.gbk')) or os.path.isfile(file_path)):
        sequence = read_sequence_from_file(file_path)
        accession_code = "local_file"
        insert_or_update_sequences('sequences_data.db', accession_code, str(sequence), 'sequences')
    else:
        sequence = None

    if sequence is not None:
        seq_complement = complement(sequence)
        reverse_complement_dna = reverse_complement(sequence)
        mRNA = transcribe(sequence)
        coding_dna = reverse_transcribe(mRNA)
        protein_seq_mRNA = translate_mRNA(mRNA)
        protein_seq_dna = translate_dna(sequence)

        save_sequence_to_file(sequence, accession_code, "original.fasta")
        save_sequence_to_file(seq_complement, accession_code, "complement.fasta")
        save_sequence_to_file(reverse_complement_dna, accession_code, "reverse_complement_dna.fasta")
        save_sequence_to_file(mRNA, accession_code, "mRNA.fasta")
        save_sequence_to_file(coding_dna, accession_code, "coding_dna.fasta")
        save_sequence_to_file(protein_seq_mRNA, accession_code, "protein_seq_mRNA.fasta")
        save_sequence_to_file(protein_seq_dna, accession_code, "protein_seq_dna.fasta")
        insert_or_update_modified_sequences('sequences_data.db', accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)
    else:
        db = input("Enter the database (e.g., 'nucleotide'): ")
        accession_code = input("Enter the ID of the record: ")

        if db and accession_code:
            sequence = read_sequence_from_database(db, accession_code)
            if sequence:
                sequence = Seq(sequence)
                insert_or_update_sequences('sequences_data.db', accession_code, str(sequence), 'sequences')
                accession_code = input("Enter the ID of the query record: ")
                result = query_database('sequences_data.db', accession_code, 'sequences')
                if result:
                    table = [["Record_ID", "sequence"], [accession_code, result]]
                    print("Result:")
                    print(tabulate(table, headers="firstrow", tablefmt="grid"))
                else:
                    print("No result found.")
                seq_complement = complement(sequence)
                reverse_complement_dna = reverse_complement(sequence)
                mRNA = transcribe(sequence)
                coding_dna = reverse_transcribe(mRNA)
                protein_seq_mRNA = translate_mRNA(mRNA)
                save_sequence_to_file(sequence, accession_code, "sequence.fasta")
                save_sequence_to_file(seq_complement, accession_code, "complement.fasta")
                save_sequence_to_file(reverse_complement_dna, accession_code, "reverse_complement_dna.fasta")
                save_sequence_to_file(mRNA, accession_code, "mRNA.fasta")
                save_sequence_to_file(coding_dna, accession_code, "coding_dna.fasta")
                save_sequence_to_file(protein_seq_mRNA, accession_code, "protein_seq_mRNA.fasta")
                # Corrected function call
                insert_or_update_modified_sequences('sequences_data.db', accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)
                accession_code = input("Enter the ID of the query record: ")
                result = query_modified_seq_table('sequences_data.db', accession_code)

                if result:
                    seq_record, seq_complement, seq_reverse_complement, mRNA_seq, protein_seq = result
                    table = [["Record_ID", "sequence"], [accession_code, seq_record]]
                    
                    print("Original Sequence:")
                    print(tabulate(table, headers="firstrow", tablefmt="grid"))
                    
                    print("Modified Sequences:")
                    table = [
                        ["Sequence Type", "Sequence"],
                        ["Complement", seq_complement],
                        ["Reverse Complement", seq_reverse_complement],
                        ["mRNA", mRNA_seq],
                        ["Protein", protein_seq],
                    ]
                    
                    print(tabulate(table, headers="firstrow", tablefmt="grid"))
                else:
                    print("No result found.")


            else:
                print("Error retrieving sequence from the database.")
        else:
            print("Invalid input. Database query requires both 'db' and 'id'.")
