# MSP
Helpful scripts to process data in our lab

In order to run these scripts you need python3.
Additionaly, some script use the Biopython and cs50 modules.
Therefore, you have to install the Biopython and cs50 modules in your computer.
You can get more information about those modules by visiting:
https://biopython.org/
https://github.com/cs50/python-cs50

# About fetch_features.py
Function:
This program needs Biophyton and cs50 modules. If you don't have those
modules, please visit https://biopython.org/ and
https://github.com/cs50/python-cs50 for more information.

The program fetches information from a list of Genebank accession numbers or
a list of BioSample numbers. The program enters the nuccore database and
collects all the features of the corresponding list of accession numbers.
When you enter a list of accession numbers you have two options. The first
option is to get the features of the provided accession list. The second one
is to get the features of all the accession numbers associated with a
specific BioSample number. In this second option, the program gets the
BioSample number of every accession number in the list, accesses the nuccore
database with the BioSample number and selects the most updated information
of every molecule (chromosome and/or plasmid(s)).

You can create one list of accession numbers in Excel by saving the file as
txt. The list needs a header if it doesn't have one the first accession
numbers is not going to be included.

Usage: python fetch_features.py
