# fetch_features.py version 1.0
# Created by Ivan Munoz-Gutierrez
# Date 2020-07-26
#
# Function:
# This program needs the Biophyton and cs50 modules. If you don't have those
# modules, please visit https://biopython.org/ and
# https://github.com/cs50/python-cs50 for more information.
#
# The program fetches information from a list of Genebank accession numbers or
# a list of BioSample numbers. The program enters the nuccore database and
# collects all the features of the corresponding list of accession numbers.
# When you enter a list of accession numbers you have two options. The first
# option is to get the features of the provided accession list. The second one
# is to get the features of all the accession numbers associated with an
# specific BioSample number. In this second option, the program gets the
# BioSample number of every accession number in the list, accesses the nuccore
# database with the BioSample number and selects the most updated information
# of every molecule (chromosome and/or plasmid(s)).
#
# You can create one list of accession numbers in Excel by saving the file as
# txt. The list needs a header, if it doesn't have one the first accession
# numbers is not going to be included.
#
# Usage: python fetch_features.py

#############################################################################
#                       Importing relevant modules                          #
#############################################################################
from Bio import Entrez
from Bio import SeqIO
from database import *
import csv
import sys
import cs50

#############################################################################
#                       Getting input from the user                         #
#############################################################################
# Checking the correct useage of the program
if len(sys.argv) != 1:
    sys.exit("usage: python fetch_features.py")

# Getting type of input data.
while True:
    type_list = cs50.get_string(
        "Does your list have accession numbers or "
        "biosample numbers (accession or biosample)? ")
    type_list = type_list.lower()
    if type_list == 'accession' or type_list == 'biosample':
        break

# Asking about getting data from all BioSample related acc numbers
while True:
    get_biosample = cs50.get_string(
        "If you have a list of accession numbers, do you want to get the "
        "most updated features of \nall the related accession numbers that "
        "belong to the same BioSample (yes or no)? ")
    get_biosample = get_biosample.lower()
    if get_biosample == 'yes' or get_biosample == 'no':
        break

# Gettin name of infile.txt
infile = cs50.get_string('Provide the name of your infile.txt: ')

# IMPORTANT: always provide your email address to GenBank
email_address = cs50.get_string("Provide your email address to the NCBI: ")
Entrez.email = email_address

# Opening infile.txt
with open(infile, 'r') as reader:

    # Skip the header
    next(reader)

    # Creating a list of accession numbers
    list_accessions = reader.readlines()

# Counting the number of results (number of sequences)
count = len(list_accessions)
print(f"Number of requested sequences: {count}")

# Creating batches of acc numbers for the specified case in the if statement
if type_list == 'accession' and get_biosample == 'no':
    # Number of sequences to be requested by batch.
    # A batch of 500 is the max that we can request.
    batch_size = 100

    # This is going to be a list of strings containg the batches of requested
    # accession numbers, i.e. every string in the list is going to have in this
    # case 500 accesion numbers separaded by comas.
    submission_list = []

    # Counter to access the list_accessions
    counter_accessions = 0

    # Loop to create the list of accession numbers by batches of 500
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        # This list is going to save temporarily the batch of accession numbers
        # that are goingo to be converted into a string separed by commas
        set_list = []
        for set in range(start, end):
            set_list.append(list_accessions[counter_accessions].
                            replace('\n', ''))
            counter_accessions += 1
        # Converting the list into string
        set_list = ','.join(set_list)
        submission_list.append(set_list)

#############################################################################
#                      Working with GenBank                                 #
#############################################################################
# Number to keep track set of sequences (or batches) and the sequences,
# it is important in case the connection to NCBI is interrupted so we can
# know where to continue downloading
set, seq_counter = 1, 1

# Opening our results file to write the fetched data in csv format
with open("results.csv", "w") as results:
    # Field names or headers in the csv table
    fields = ["set_batch", "counter", "description", "accession", "size",
              "molecule", "mod_date", "topology", "mol_type", "organism",
              "strain", "isolation_source", "host", "plasmid", "country",
              "lat_lon", "collection_date", "note", "serovar", "collected_by",
              "genotype", "BioProject", "BioSample", "Assem_Method",
              "Gen_Coverage", "Seq_Technol"]

    # Create DictWriter
    writer = csv.DictWriter(results, fields)

    # Writing headers
    writer.writeheader()

    ###########################################
    # Workig with a list of accession numbers #
    ###########################################
    if type_list == 'accession' and get_biosample == 'no':
        # Declaring end
        end = 0

        # Fetching the information from GenBank by batches
        for submission in range(len(submission_list)):
            start = end
            # submission_list is a list of accession numbers separated by
            # commas. Therefore, the number of commas indicate the number of
            # accession numbers.
            end = end + submission_list[submission].count(',') + 1

            # Printing download batch record
            print("Going to download record %i to %i" % (start + 1, end))

            # Posting the submission_list.
            # Because we are requesting information from a huge list of acc
            # numbers, we have to use the ".epost" function which uploads a
            # list of UIs (acc numbers) for use in subsequent searches.
            # From .epost we can get the QueryKey and the WebEnv which define
            # our history session and can be used to performe searches of data.
            posting = Entrez.epost('nuccore', id=submission_list[submission])
            search_results = Entrez.read(posting)

            # Copying cookie "WebEnv" and query "QueryKey" from our history
            # session to keep track of our batch fetching.
            # WevEnv -> Web environment string returned from a previous
            # ESearch, EPost or ELink call; QueryKey -> Integer query key
            # returned by a previous ESearch, EPost or ELink call
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            # Getting the batch information
            # db -> database, nuccore -> nuleotide, rettype -> retrieval type
            # retmode -> determines the format of the return output
            # retstart -> sequential index of the first UID in the retrieved
            # set to be shown in the XML output
            # retmax -> total number of UIDs from the retrieved set to be shown
            # in the XML output
            # idtype-> specifies the type of identifier to return for sequence
            # databases, acc -> accesion number
            fetch_handle = Entrez.efetch(
                db="nuccore",
                rettype="gb",
                retmode="text",
                retstart=0,
                retmax=batch_size,
                webenv=webenv,
                query_key=query_key,
                idtype="acc"
            )

            # Parsing the data fetched from NCBI
            records = parser(fetch_handle, set, seq_counter)

            # Recording the number set and sequences downloded
            set += 1
            seq_counter = records[1]

            # Saving the retrived data in the csv file
            for i in range(len(records[0])):
                writer.writerow(records[0][i])

            # Closing fetch_handle
            fetch_handle.close()

    ###################################################################
    # Workig with a list of accession numbers and getting all related #
    # BioSample associated accession numbers and their features       #
    ###################################################################
    elif type_list == 'accession' and get_biosample == 'yes':
        seq_counter = 1

        # Iterating over the list of accession numbers
        for query in range(len(list_accessions)):
            # Number to keep track set of sequences (or batches) and the
            # sequences, it is important in case the connection to NCBI is
            # interrupted so we can know where to continue downloading
            set = query + 1

            # Getting the BioSample number of the requested accession number.
            # BioSample() gets two arguments, a list of accession numbers and
            # the email address of the user
            biosample_number = BioSample_list(list_accessions[query],
                                              email_address)

            # Using ".esearch" to find the information.
            # Also we have to implement "usehistory" to get the cookie and
            # query key.
            # db -> database to search, term -> Entrez text query
            search_handle = Entrez.esearch(db="nucleotide",
                                           term=biosample_number,
                                           usehistory="y")

            # Copying the information in computer memory
            search_results = Entrez.read(search_handle)

            # Closing the handle
            search_handle.close()

            # Counting the number of results (number of sequences)
            count = int(search_results["Count"])
            print(f"Number of requested sequences from BioSample: {count}")

            # Copying cookie "WebEnv" and query "QueryKey" from history to keep
            # track of our batch fetching.
            # WevEnv -> Web environment string returned from a previous
            # ESearch, EPost or ELink call.
            # QueryKey -> Integer query key returned by a previous ESearch,
            # EPost or ELink call
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            # Number of sequences to be requested by batch.
            # A batch of 500 is the max that we can request.
            batch_size = 500

            # I need to think about how to clean features of a BioSample number
            # that has more than 500 records
            # TODO

            # Fetching the information from GenBank by batches
            for start in range(0, count, batch_size):
                end = min(count, start + batch_size)

                # Printing download batch record
                print(f"Going to download record {start + 1} to {end} "
                      f"from set {query + 1}")

                # Getting the batch information
                # db -> database, nuccore -> nuleotide, rettype -> retrieval
                # type, retmode -> determines the format of the return output
                # retstart -> sequential index of the first UID in the
                # retrieved set to be shown in the XML output, retmax -> total
                # number of UIDs from the retrieved set to be shown in the
                # XML output, idtype-> specifies the type of identifier to
                # return for sequence databases, acc -> accesion number
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

                # Parsing the data fetched from NCBI
                records = parser(fetch_handle, set, seq_counter)

                # Recording the number of sequences downloded
                seq_counter = records[1]

                # Using dabase.py to clean data and obtain the most updated
                # information.
                # clean_features() returns a list of dictionaries
                uptaded_features = clean_features(records[0])

                # Saving the updated retrived data in the csv file
                for i in range(len(uptaded_features)):
                    writer.writerow(uptaded_features[i])

                print(f"Number of sequences saved after processing: "
                      f"{len(uptaded_features)}")

                # Closing handle
                fetch_handle.close()

    ############################################
    # Working with a list of BioSample numbers #
    ############################################
    else:
        lenght_acc_list = len(list_accessions)

        # Fetching the information from GenBank
        for submission in range(lenght_acc_list):

            # Printing download record
            print(f"Going to download record {submission} "
                  f"of {lenght_acc_list}")

            # Searching for the BioSample accession number. We need usehistory
            # to get the QueryKey and the WebEnv which define our history
            # session and can be used to performe searches of data.
            search_handle = Entrez.esearch(db="nuccore",
                                           term=list_accessions[submission],
                                           usehistory="y")
            search_results = Entrez.read(search_handle)

            # Copying cookie "WebEnv" and query "QueryKey" from our history
            # session.
            # WevEnv -> Web environment string returned from a previous
            # ESearch, EPost or ELink call; QueryKey -> Integer query key
            # returned by a previous ESearch, EPost or ELink call
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            # Getting information
            # db -> database, nuccore -> nuleotide, rettype -> retrieval type,
            # retmode -> determines the format of the return output,
            # retstart -> sequential index of the first UID in the retrieved
            # set to be shown in the XML output, retmax -> total number of
            # UIDs from the retrieved set to be shown in the XML output,
            # idtype-> specifies the type of identifier to return for sequence
            # databases, acc -> accesion number
            fetch_handle = Entrez.efetch(
                db="nuccore",
                rettype="gb",
                retmode="text",
                retstart=0,
                retmax=1,
                webenv=webenv,
                query_key=query_key,
                idtype="acc"
            )

            # Parsing the data fetched from NCBI
            records = parser(fetch_handle, submission + 1, submission + 1)

            # Saving the retrived data in the csv file
            for i in range(len(records[0])):
                writer.writerow(records[0][i])

            # Closing fetch_handle
            fetch_handle.close()

# If everything was done OK print Done and exit the program
print("""Done!
You should have a results.csv file in your folder""")

sys.exit(0)
