# program to convert a fasta file into a csv file
# created by Ivan Munoz-Gutierrez
# June 24, 2020
# usage: python fasta_to_csv.py infile.fasta

import sys

# checking correct usage of the program
if len(sys.argv) != 2:
    sys.exit("usage: python fasta_to_csv.py infile.fasta")

# declaring a dictionary to save the infile information in memory
dictionary = {}

# opening the infile
with open(sys.argv[1], "r") as reader:
    # counter to know the number sequences
    counter = 0

    # variable to save the header of the fasta sequence
    head = ""

    for row in reader:
        # finding the header of the fasta sequence
        if row[0] == '>':
            # creating a nested dictionary to save in numerical order all the fasta sequences
            dictionary[counter] = {}
            # removing the ">" and "\n" characters
            head = row.replace('>', '', 1)
            head = head.replace('\n', '')
            # saving the header in the nested dictionary
            dictionary[counter]["header"] = head
            # initializing the seq vatiable to store the sequence
            seq = ""
            # the counter increments when detecting a '>' character in the header
            counter += 1
        # saving the sequence in the nested array
        else:
            # replacing the '\n' character
            row = row.replace('\n', '')
            # concatenated the sequence
            seq = seq + row
            # saving the sequence in dictionary
            # counter is -1 to save the sequence with the same key number as the header
            dictionary[counter - 1]["sequence"] = seq

# printing the number of sequences to be converted
print(f"converting {counter} sequences")

# opening the results file
with open("results.csv", 'w') as writer:
    # writing the headers
    writer.write("Hits in genome seqs,Allele\n")

    # looping into the dictionary to extract and write the information
    for i in range(counter):
        head = dictionary[i]['header']
        seq = dictionary[i]["sequence"]
        # concatenating the header and the sequence to save them as csv
        concatenated = head + ',' + seq + '\n'
        writer.writelines(concatenated)

print("done")

sys.exit(0)
