# biosamble.py version 1.0
# Created by: Ivan Munoz-Gutierrez
# Date: July 05, 2020
# Function: the program fetches information from a collection of sequences from
# GeneBank. By providing an accession number, the program looks for the
# BioSample number and fetches features from all the sequences associated to
# that BioSample number.
# Notes: this program uses the history feature and the WebEnv session cookie to
# download large data in batches.
# Usage: python biosamble.py accession_number

# Importing Biopython modules
from Bio import Entrez
from Bio import SeqIO
# Importing csv and sys modules
import csv
import sys

# Checking the correct usage of the program
if len(sys.argv) != 2:
    sys.exit('usage: python biosamble.py accession_number')

# Provide email address to GeneBank
Entrez.email = "ivan.munoz.guterrez@gmail.com"

###########################################################################
#      Getting the BioSample number of the requested accession number     #
###########################################################################

# Using '.efetch' to retrieve the information (BioSample number) of the
# requested accesion number.
# db -> database, nuccore -> nucleotide, id -> id number of the requested
# information, in this case the accesion number provided in argv[1],
# rettype -> retrieval type, retmode -> determines the format of the return output
handle = Entrez.efetch(db='nuccore', id=sys.argv[1], rettype='gb',
                       retmode='text')

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

print(f'Requested accesion number: {record_acc.id}')
print(f'Corresponding BioSample number: {biosample_number}')

##########################################################################
#        Obtaining the features of all the BioSample sequences           #
##########################################################################

# Using ".esearch" to find the information.
# Also we have to implement "usehistory" to get the cookie and query key.
# db -> database to search, term -> Entrez text query
search_handle = Entrez.esearch(db="nucleotide", term=biosample_number,
                               usehistory="y")

# Copying the information in computer memory
search_results = Entrez.read(search_handle)
search_handle.close()

# Counting the number of results (number of sequences)
count = int(search_results["Count"])
print(f"Number of requested sequences from BioSample: {count}")

# Copying cookie and query from history to keep track of our batch fetching
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

# Number of sequences to be requested by batch.
# A batch of 500 is the max that we can request.
batch_size = 500

# Number to keep track of sequences, it is important in case the connection
# to NCBI is interrupted so we can know where to continue downloading.
seq_counter = 1

# Opening our results file to write the fetched data in csv format
with open("results.csv", "w") as results:
    writer = csv.writer(results)

    # Field names or headers in the csv table
    fields = ["counter", "description", "accession", "size", "molecule",
              "mod_date", "topology", "mol_type", "organism", "strain",
              "isolation_source", "host", "plasmid", "country", "lat_lon",
              "collection_date", "BioProject", "BioSample", "Assem_Method",
              "Gen_Coverage", "Seq_Technol"]

    # Writing headers
    writer.writerow(fields)

    # Fetching the information from GenBank by batches
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)

        # Printing download batch record
        print("Going to download record %i to %i" % (start + 1, end))

        # Getting the batch information
        # db -> database, nuccore -> nuleotide, rettype -> retrieval type
        # retmode -> determines the format of the return output
        # retstart -> sequential index of the first UID in the retrieved set to be shown in the XML output
        # retmax -> total number of UIDs from the retrieved set to be shown in the XML output
        # idtype-> specifies the type of identifier to return for sequence databases, acc -> accesion number
        fetch_handle = Entrez.efetch(
            db="nuccore",
            rettype="gb",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
            idtype="acc"
        )

        # Parsing throw the fetched information
        for seq_record in SeqIO.parse(fetch_handle, "gb"):
            # Keeping track of the number of sequences saved
            field0 = seq_counter
            seq_counter += 1

            # Extracting description
            field1 = seq_record.description

            # Extracting sequence id, i.e. accession number
            field2 = seq_record.id

            # '.seq' is an object with the sequence itself
            field3 = len(seq_record.seq)

            # Checking whether it is chromosome or plasmid
            if 'chromosome' in field1:
                field4 = 'chromosome'
            elif 'plasmid' in field1:
                field4 = 'plasmid'
            else:
                field4 = 'missing'

            # '.annotations' is a dictionary of aditional information about the sequence as
            # last modification date, topology, sequence_version, organims, references, etc.
            if 'date' in seq_record.annotations:
                field5 = seq_record.annotations['date']
            else:
                field5 = 'missing'

            if "topology" in seq_record.annotations:
                field6 = seq_record.annotations["topology"]
            else:
                field6 = "missing"

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
                        mol_type = dictionary.get("mol_type")
                        field7 = mol_type[0]
                    else:
                        field7 = "missing"

                    if 'organism' in dictionary:
                        bacterium = dictionary.get('organism')
                        field8 = bacterium[0]
                    else:
                        field8 = "missing"

                    if 'strain' in dictionary:
                        strain = dictionary.get('strain')
                        field9 = strain[0]
                    else:
                        field9 = "missing"

                    if 'isolation_source' in dictionary:
                        isolation_source = dictionary.get('isolation_source')
                        field10 = isolation_source[0]
                    else:
                        field10 = "missing"

                    if 'host' in dictionary:
                        host = dictionary.get('host')
                        field11 = host[0]
                    else:
                        field11 = "missing"

                    if 'plasmid' in dictionary:
                        plasmid = dictionary.get('plasmid')
                        field12 = plasmid[0]
                    else:
                        field12 = "missing"

                    if 'country' in dictionary:
                        country = dictionary.get('country')
                        field13 = country[0]
                    else:
                        field13 = "missing"

                    if "lat_lon" in dictionary:
                        lat_lon = dictionary.get("lat_lon")
                        field14 = lat_lon[0]
                    else:
                        field14 = "missing"

                    if "collection_date" in dictionary:
                        collection_date = dictionary.get("collection_date")
                        field15 = collection_date[0]
                    else:
                        field15 = "missing"

            # '.dbxrefs' is a list populated from any PROJECT or DBLINK
            # Checking if .dbxrefs has any information
            if len(seq_record.dbxrefs) == 0:
                field16 = "missing"
                field17 = "missing"

            # Converting the list dbxrefs into dictionary
            dictionary_dbxrefs = {}
            for i in range(len(seq_record.dbxrefs)):
                s = seq_record.dbxrefs[i].split(":")
                dictionary_dbxrefs[s[0]] = s[1]

            # Saving BioProject and BioSample
            if "BioProject" in dictionary_dbxrefs:
                field16 = dictionary_dbxrefs.get("BioProject")
            else:
                field16 = "missing"

            if "BioSample" in dictionary_dbxrefs:
                field17 = dictionary_dbxrefs.get("BioSample")
            else:
                field17 = "missing"

            # Getting the Genome-Assembly-Data
            # Checking if the sequence has structured_comment
            if "structured_comment" in seq_record.annotations and "Genome-Assembly-Data" in seq_record.annotations["structured_comment"]:
                if "Assembly Method" in seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]:
                    field18 = seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]["Assembly Method"]
                else:
                    field18 = "missing"

                if "Genome Coverage" in seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]:
                    field19 = seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]["Genome Coverage"]
                else:
                    field19 = "missing"

                if "Sequencing Technology" in seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]:
                    field20 = seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]["Sequencing Technology"]
                else:
                    field20 = "missing"
            else:
                field18 = "missing"
                field19 = "missing"
                field20 = "missing"

            # Copying all the obtained results in the list fields
            fields = [field0, field1, field2, field3, field4,
                      field5, field6, field7, field8, field9,
                      field10, field11, field12, field13, field14,
                      field15, field16, field17, field18, field19,
                      field20]

            # Saving the retrived data in the csv file
            writer.writerow(fields)
        fetch_handle.close()

# If everything was done OK print Done and exit the program
print("Done")
sys.exit(0)
