# fetch_features.py version 1.0
# Created by Ivan Munoz-Gutierrez
# Date July 09, 2020
# Function: the program fetches information from a list of Genebank accession
# numbers. The program enters the nulecotide database and collect all the
# features of the list.
# Notes: this program uses the history feature and the WebEnv session cookie to
# download large data in batches. Don't forget to type your email address in
# line 25.
# Usage: python fetch_features.py list_of_accession_numbers.txt


# Importing Biopython modules
from Bio import Entrez
from Bio import SeqIO
# Importing csv module and sys
import csv
import sys

# Checking the correct useage of the program
if len(sys.argv) != 2:
    sys.exit("usage: python fetch_features.py accession_list.txt")

#############################################################################
#        Making a list of the accession numbers to be analyzed              #
#############################################################################
with open(sys.argv[1], 'r') as reader:

    # Skip the header
    next(reader)

    list_accessions = []

    # Creating a list of accession numbers
    for row in reader:
        list_accessions.append(row.replace('\n', ''))

    # Counting the number of results (number of sequences)
    count = len(list_accessions)
    print(f"Number of requested sequences: {count}")

# Converting the list into string
list_accessions = ','.join(list_accessions)

#############################################################################
#                      Working with GenBank                                 #
#############################################################################

# IMPORTANT: always provide your email address to GenBank
Entrez.email = "ivan.munoz.gutierrez@gmail.com"

# Because we are requesting information from a huge list of accession numbers
# we have to use the ".epost" function which uploads a list of UIs (accession
# numbers) for use in subsequent searches.
# From .epost we can get the QueryKey and the WebEnv which define our history
# session and can be used to performe searches of data.
posting = Entrez.epost('nuccore', id=list_accessions)
search_results = Entrez.read(posting)

# Copying cookie "WebEnv" and query "QueryKey" from our history session to keep
# track of our batch fetching. WevEnv -> Web environment string returned from a
# previous ESearch, EPost or ELink call; QueryKey -> Integer query key returned
# by a previous ESearch, EPost or ELink call
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

# Number of sequences to be requested by batch.
# IMPORTANT: a batch of 500 is the max that we can request from NCBI
batch_size = 500

# Number to keep track of sequences, it is important in case the connection
# to NCBI is interrupted so we can know where to continue downloading
seq_counter = 1

# Opening our results file to write the fetched data in csv format
with open("results.csv", "w") as results:
    writer = csv.writer(results)

    # Field names or headers in the csv table
    fields = ["counter", "description", "accession", "size", "molecule",
              "mod_date", "topology", "mol_type", "organism", "strain",
              "isolation_source", "host", "plasmid", "country", "lat_lon",
              "collection_date", "note", "serovar", "collected_by", "genotype",
              "BioProject", "BioSample", "Assem_Method", "Gen_Coverage",
              "Seq_Technol"]

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

                    if "note" in dictionary:
                        note = dictionary.get("note")
                        field16 = note[0]
                    else:
                        field16 = "missing"

                    if "serovar" in dictionary:
                        serovar = dictionary.get("serovar")
                        field17 = serovar[0]
                    else:
                        field17 = "missing"

                    if "collected_by" in dictionary:
                        collected_by = dictionary.get("collected_by")
                        field18 = collected_by[0]
                    else:
                        field18 = "missing"

                    if "genotype" in dictionary:
                        genotype = dictionary.get("genotype")
                        field19 = genotype[0]
                    else:
                        field19 = "missing"

                    break

            # '.dbxrefs' is a list populated from any PROJECT or DBLINK
            # Checking if .dbxrefs has any information
            if len(seq_record.dbxrefs) == 0:
                field20 = "missing"
                field21 = "missing"

            # Converting the list dbxrefs into dictionary
            dictionary_dbxrefs = {}
            for i in range(len(seq_record.dbxrefs)):
                s = seq_record.dbxrefs[i].split(":")
                dictionary_dbxrefs[s[0]] = s[1]

            # Saving BioProject and BioSample
            if "BioProject" in dictionary_dbxrefs:
                field20 = dictionary_dbxrefs.get("BioProject")
            else:
                field20 = "missing"

            if "BioSample" in dictionary_dbxrefs:
                field21 = dictionary_dbxrefs.get("BioSample")
            else:
                field21 = "missing"

            # Getting the Genome-Assembly-Data
            # Checking if the sequence has structured_comment
            if "structured_comment" in seq_record.annotations and "Genome-Assembly-Data" in seq_record.annotations["structured_comment"]:
                if "Assembly Method" in seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]:
                    field22 = seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]["Assembly Method"]
                else:
                    field22 = "missing"

                if "Genome Coverage" in seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]:
                    field23 = seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]["Genome Coverage"]
                else:
                    field23 = "missing"

                if "Sequencing Technology" in seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]:
                    field24 = seq_record.annotations["structured_comment"]["Genome-Assembly-Data"]["Sequencing Technology"]
                else:
                    field24 = "missing"
            else:
                field22 = "missing"
                field23 = "missing"
                field24 = "missing"

            # Copying all the obtained results in the list fields
            fields = [field0, field1, field2, field3, field4,
                      field5, field6, field7, field8, field9,
                      field10, field11, field12, field13, field14,
                      field15, field16, field17, field18, field19,
                      field20, field21, field22, field23, field24]

            # Saving the retrived data in the csv file
            writer.writerow(fields)
        fetch_handle.close()

# If everything was done OK print Done and exit the program
print("Done")
sys.exit(0)
