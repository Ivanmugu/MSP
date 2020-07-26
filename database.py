# Scripts to be used in fetch_features.py and biosample.py
# Created by: Ivan Munoz-Gutierrez
# Date: 2020-07-25

# Importing relevant modules
from Bio import Entrez
from Bio import SeqIO
from datetime import datetime
import csv
import sqlite3
import sys


#############################################################################
#         Function to make a list of the BioSample accession numbers        #
#############################################################################
# It takes as argument a list of accession numbers and email address
def BioSample_list(accession_number_list, email_address):
    # Provide email address to GeneBank
    Entrez.email = email_address

    # Using '.efetch' to retrieve the information (BioSample number) of the
    # requested accesion number.
    # db -> database, nuccore -> nucleotide, id -> id number of the requested
    # information, in this case the accesion number provided in argv[1],
    # rettype -> retrieval type, retmode -> determines the format of the return
    # output
    handle = Entrez.efetch(db='nuccore', id=accession_number_list,
                           rettype='gb', retmode='text')

    # Copying the information in computer memory
    record_acc = SeqIO.read(handle, 'gb')

    handle.close()

    # Convert the list dbxrefs into dictionary to get the BioSample number
    dictionary_dbxrefs = {}
    for index in range(len(record_acc.dbxrefs)):
        list_dbxrefs = record_acc.dbxrefs[index].split(':')
        dictionary_dbxrefs[list_dbxrefs[0]] = list_dbxrefs[1]

    # Getting the BioSample number
    biosample_number = dictionary_dbxrefs['BioSample']

    print(f'\nRequested accesion number: {record_acc.id}')
    print(f'Corresponding BioSample number: {biosample_number}')

    return biosample_number


#############################################################################
#           Function to parse information fetched from nuccore NCBI         #
#############################################################################
# It takes as arguments the handle to fetch information, the set or batch
# batch number, and a conunter to keep track of the sequences retrieved.
def parser(fetch_handle, set, seq_counter):
    # Creating a list to save records in dictionaries
    records = []
    info = {}

    counter = seq_counter

    # Parsing throw the fetched information
    for seq_record in SeqIO.parse(fetch_handle, "gb"):
        # keeping track of the set (or batch) analyzed
        info['set_batch'] = set

        # Keeping track of the number of sequence saved
        info['counter'] = counter

        # Extracting description
        info['description'] = seq_record.description

        # Extracting sequence id, i.e. accession number
        info['accession'] = seq_record.id

        # '.seq' is an object with the sequence itself
        info['size'] = len(seq_record.seq)

        # '.annotations' is a dictionary of aditional information about the
        # sequence as last modification date, topology, sequence_version,
        # organims, references, etc.
        if 'date' in seq_record.annotations:
            mod_date = seq_record.annotations['date']
            mod_date = datetime.strptime(mod_date, '%d-%b-%Y')
            mod_date = datetime.strftime(mod_date, "%Y-%m-%d")
            info['mod_date'] = mod_date
        else:
            info['mod_date'] = 'missing'

        if "topology" in seq_record.annotations:
            info['topology'] = seq_record.annotations["topology"]
        else:
            info['topology'] = "missing"

        # Checking whether it is chromosome or plasmid
        if 'chromosome' in info['description']:
            info['molecule'] = 'chromosome'
        elif info['size'] >= 4000000 and info['topology'] == 'circular':
            info['molecule'] = 'chromosome'
        elif 'plasmid' in info['description']:
            info['molecule'] = 'plasmid'
        else:
            info['molecule'] = 'missing'

        # '.features' is a list of SeqFeatures objects with more structured
        # information about the features on a sequence
        feature = seq_record.features

        # Looping throw list feature
        for index in feature:

            # '.type' is only a description of the type of feature
            # that could be source, CDS, gene, etc.
            # In source we can find organism, strain, host, country, etc.
            if index.type == "source":

                # Creating a dictionary of the qualifiers from source
                dictionary = dict(index.qualifiers)

                # '.get' gives a list
                if "mol_type" in dictionary:
                    info['mol_type'] = dictionary.get("mol_type")[0]
                else:
                    info['mol_type'] = 'missing'

                if 'organism' in dictionary:
                    info['organism'] = dictionary.get('organism')[0]
                else:
                    info['organism'] = 'missing'

                if 'strain' in dictionary:
                    info['strain'] = dictionary.get('strain')[0]
                else:
                    info['strain'] = 'missing'

                if 'isolation_source' in dictionary:
                    info['isolation_source'] =\
                        dictionary.get('isolation_source')[0]
                else:
                    info['isolation_source'] = 'missing'

                if 'host' in dictionary:
                    info['host'] = dictionary.get('host')[0]
                else:
                    info['host'] = 'missing'

                if 'plasmid' in dictionary:
                    info['plasmid'] = dictionary.get('plasmid')[0]
                else:
                    info['plasmid'] = 'missing'

                if 'country' in dictionary:
                    info['country'] = dictionary.get('country')[0]
                else:
                    info['country'] = 'missing'

                if "lat_lon" in dictionary:
                    info['lat_lon'] = dictionary.get("lat_lon")[0]
                else:
                    info['lat_lon'] = 'missing'

                if "collection_date" in dictionary:
                    info['collection_date'] =\
                        dictionary.get("collection_date")[0]
                else:
                    info['collection_date'] = 'missing'

                if "note" in dictionary:
                    info['note'] = dictionary.get("note")[0]
                else:
                    info['note'] = 'missing'

                if "serovar" in dictionary:
                    info['serovar'] = dictionary.get("serovar")[0]
                else:
                    info['serovar'] = 'missing'

                if "collected_by" in dictionary:
                    info['collected_by'] = dictionary.get("collected_by")[0]
                else:
                    info['collected_by'] = 'missing'

                if "genotype" in dictionary:
                    info['genotype'] = dictionary.get("genotype")[0]
                else:
                    info['genotype'] = 'missing'

                break

        # '.dbxrefs' is a list populated from any PROJECT or DBLINK
        # Checking if .dbxrefs has any information
        if len(seq_record.dbxrefs) == 0:
            info['BioProject'] = 'missing'
            info['BioSample'] = 'missing'

        # Converting the list dbxrefs into dictionary
        dictionary_dbxrefs = {}
        for i in range(len(seq_record.dbxrefs)):
            s = seq_record.dbxrefs[i].split(":")
            dictionary_dbxrefs[s[0]] = s[1]

        # Saving BioProject and BioSample
        if "BioProject" in dictionary_dbxrefs:
            info['BioProject'] = dictionary_dbxrefs.get("BioProject")
        else:
            info['BioProject'] = 'missing'

        if "BioSample" in dictionary_dbxrefs:
            info['BioSample'] = dictionary_dbxrefs.get("BioSample")
        else:
            info['BioSample'] = 'missing'

        # Getting the Genome-Assembly-Data
        # Checking if the sequence has structured_comment
        if "structured_comment" in seq_record.annotations and\
                "Genome-Assembly-Data" in\
                seq_record.annotations["structured_comment"]:
            if "Assembly Method" in\
                    seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]:
                info['Assem_Method'] =\
                    seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]["Assembly Method"]
            else:
                info['Assem_Method'] = "missing"

            if "Genome Coverage" in\
                    seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]:
                info['Gen_Coverage'] =\
                    seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]["Genome Coverage"]
            else:
                info['Gen_Coverage'] = "missing"

            if "Sequencing Technology" in\
                    seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]:
                info['Seq_Technol'] =\
                    seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]["Sequencing Technology"]
            else:
                info['Seq_Technol'] = "missing"

        elif "structured_comment" in seq_record.annotations and\
                "Assembly-Data" in\
                seq_record.annotations["structured_comment"]:
            if "Assembly Method" in\
                    seq_record.annotations["structured_comment"]["Assembly-Data"]:
                info['Assem_Method'] =\
                    seq_record.annotations["structured_comment"]["Assembly-Data"]["Assembly Method"]
            else:
                info['Assem_Method'] = "missing"

            if "Genome Coverage" in\
                    seq_record.annotations["structured_comment"]["Assembly-Data"]:
                info['Gen_Coverage'] =\
                    seq_record.annotations["structured_comment"]["Assembly-Data"]["Genome Coverage"]
            else:
                info['Gen_Coverage'] = "missing"

            if "Sequencing Technology" in\
                    seq_record.annotations["structured_comment"]["Assembly-Data"]:
                info['Seq_Technol'] =\
                    seq_record.annotations["structured_comment"]["Assembly-Data"]["Sequencing Technology"]
            else:
                info['Seq_Technol'] = "missing"

        else:
            info['Assem_Method'] = "missing"
            info['Gen_Coverage'] = "missing"
            info['Seq_Technol'] = "missing"

        # Copying info (dictionary) into records (list)
        records.append(info.copy())

        counter += 1

    return (records, counter)


#############################################################################
#         Function to clean up features obtained from a BioSample           #
#         and get the most uptaded information                              #
#############################################################################
# It takes as argument raw_data that is a list of dictionaries holding all the
# data downloaded from an specific BioSample associated accession numbers
def clean_features(raw_data):

    # Creating empty features.db
    open('features.db', 'w').close()

    # Opening features.db for sqlite3
    conn = sqlite3.connect('features.db')

    c = conn.cursor()

    # Creating table for features raw data
    c.execute("""CREATE TABLE features_raw (
                set_batch int,
                counter integer,
                description text,
                accession text,
                size integer,
                molecule text,
                mod_date text,
                topology text,
                mol_type text,
                organism text,
                strain text,
                isolation_source text,
                host text,
                plasmid text,
                country text,
                lat_lon text,
                collection_date text,
                note text,
                serovar text,
                collected_by text,
                genotype text,
                BioProject text,
                BioSample text,
                Assem_Method text,
                Gen_Coverage text,
                Seq_Technol text
                )""")
    conn.commit()

    # Transfering results obtained from NCBI into table features_raw
    for i in range(len(raw_data)):
        c.execute("""INSERT INTO features_raw
                     VALUES(:set_batch, :counter, :description, :accession,
                            :size, :molecule, :mod_date, :topology, :mol_type,
                            :organism, :strain, :isolation_source, :host,
                            :plasmid, :country, :lat_lon, :collection_date,
                            :note, :serovar, :collected_by, :genotype,
                            :BioProject, :BioSample, :Assem_Method,
                            :Gen_Coverage, :Seq_Technol)""",
                  {'set_batch': int(raw_data[i]['set_batch']),
                   'counter': int(raw_data[i]['counter']),
                   'description': raw_data[i]['description'],
                   'accession': raw_data[i]['accession'],
                   'size': int(raw_data[i]['size']),
                   'molecule': raw_data[i]['molecule'],
                   'mod_date': raw_data[i]['mod_date'],
                   'topology': raw_data[i]['topology'],
                   'mol_type': raw_data[i]['mol_type'],
                   'organism': raw_data[i]['organism'],
                   'strain': raw_data[i]['strain'],
                   'isolation_source': raw_data[i]['isolation_source'],
                   'host': raw_data[i]['host'],
                   'plasmid': raw_data[i]['plasmid'],
                   'country': raw_data[i]['country'],
                   'lat_lon': raw_data[i]['lat_lon'],
                   'collection_date': raw_data[i]['collection_date'],
                   'note': raw_data[i]['note'],
                   'serovar': raw_data[i]['serovar'],
                   'collected_by': raw_data[i]['collected_by'],
                   'genotype': raw_data[i]['genotype'],
                   'BioProject': raw_data[i]['BioProject'],
                   'BioSample': raw_data[i]['BioSample'],
                   'Assem_Method': raw_data[i]['Assem_Method'],
                   'Gen_Coverage': raw_data[i]['Gen_Coverage'],
                   'Seq_Technol': raw_data[i]['Seq_Technol']})
    conn.commit()

    # Getting the most updated files and ordering by size
    c.execute("""SELECT * FROM (
                 SELECT table_one.*
                 FROM features_raw table_one
                 WHERE table_one.mod_date = (SELECT MAX(table_two.mod_date)
                                             FROM features_raw table_two
                                             WHERE table_two.size = table_one.size)
                                             )
                 WHERE molecule != "missing"
                 ORDER BY size DESC
                 """)

    results = c.fetchall()

    # Converting results (list of tuples) into a list of dictionaries
    records = []
    info = {}
    for i in range(len(results)):
        info['set_batch'] = int(results[i][0])
        info['counter'] = int(results[i][1])
        info['description'] = results[i][2]
        info['accession'] = results[i][3]
        info['size'] = int(results[i][4])
        info['molecule'] = results[i][5]
        info['mod_date'] = results[i][6]
        info['topology'] = results[i][7]
        info['mol_type'] = results[i][8]
        info['organism'] = results[i][9]
        info['strain'] = results[i][10]
        info['isolation_source'] = results[i][11]
        info['host'] = results[i][12]
        info['plasmid'] = results[i][13]
        info['country'] = results[i][14]
        info['lat_lon'] = results[i][15]
        info['collection_date'] = results[i][16]
        info['note'] = results[i][17]
        info['serovar'] = results[i][18]
        info['collected_by'] = results[i][19]
        info['genotype'] = results[i][20]
        info['BioProject'] = results[i][21]
        info['BioSample'] = results[i][22]
        info['Assem_Method'] = results[i][23]
        info['Gen_Coverage'] = results[i][24]
        info['Seq_Technol'] = results[i][25]
        records.append(info.copy())

    return records

    conn.close()
