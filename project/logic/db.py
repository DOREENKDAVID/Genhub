import sqlite3


def create_database(database_name):
    """function that creates a database in sqlite"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()


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


def insert_or_update_database(database_name, accession_code, sequence, table_name):
    """function to query databases"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()

    cursor.execute(f'''
    INSERT OR REPLACE INTO {table_name} (accession_code, seq_record)
    VALUES (?, ?)
    ''',(accession_code, sequence)
            
    )
    conn.commit()
    conn.close()

def query_database(database_name, accession_code,table_name):

    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()
    query = f'SELECT seq_record FROM {table_name} WHERE accession_code = ?'
    cursor.execute (query, (accession_code,))
    results = cursor.fetchone()
    conn.close()
    return results[0] if results else None

# Example usage:
create_database('mydatabase.db')
create_tables('mydatabase.db')

# Insert a sequence into the 'sequences' table
insert_or_update_database('mydatabase.db', 'ACCESSION123', 'ATCGATCGATCG', 'sequences')

# Query the sequence from the 'sequences' table
sequence = query_database('mydatabase.db', 'ACCESSION123', 'sequences')
print(sequence)

    