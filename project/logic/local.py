from Bio import SeqIO
import os

def read_sequences_from_file(file_path):
    """function that reads a file from user input gets the sequence ans stores in db"""
    if os.path.isfile(file_path):
        #try extracting from file and getting the accession code and sequence
        try:
            sequences = list(SeqIO.parse(file_path, "fasta"))
            if sequences:
                return sequences
            else:
                print("No sequences found in the file.")
        except Exception as e:
            print("Error:", str(e))
            return None
    else:
        print(f"File not found: {file_path}")
        return None

if __name__ == "__main__":
    file_path = input("Enter the path to a local file (or press Enter to skip this option): ")

    if file_path:
        sequences = read_sequences_from_file(file_path)
        if sequences:
            for sequence in sequences:
                accession_code = sequence.id
                seq = sequence.seq
                print(f"Accession Code: {accession_code}")
                print(seq)
