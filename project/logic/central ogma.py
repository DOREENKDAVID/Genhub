from Bio.Seq import Seq

sequence = Seq("AGTACACTGGT")
print(sequence)

coding_dna = Seq(sequence)
template_dna = coding_dna.reverse_seq()
#seq i read from 5' to 3' so reverse strand will be 3'to 5'
return template_dna


#transcribe
mRNA = coding_dna.transcribe()
return mRNA


#back transcription
coding_dna = mRNA.back_transcribe()
return mRNA


#translate
protein_seq = mRNA.translate()
return protein_seq

#if no file is uploaded search on daabase
        if db is None:
            print("No local file found. Please specify an online database query:")
            #ask for db and accession code
            db = input("Enter the database (e.g., 'nucleotide'): ")
            accession_code = input("Enter the ID of the record: ")
            if not db or not accession_code:
                print("Invalid input. Database query requires both 'db' and 'id'.")
                #return None
#direct translate from dna
protein_seq = coding_dna.translate()
return protein_seq
