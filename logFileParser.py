#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 08:42:24 2020

@author: ivanmugu
"""

import os
import re
import csv
from itertools import zip_longest

def dirNamesAndAddresses(DIRECTORY_NAME, FILE_NAME):
    """
    Gives the path to all the files that have the same FILE_NAME that are
    contained in the all the first subfolders (subdirectories) of a given
    folder (DIRECTORY_NAME). Additionaly, provides the names af all the
    first subfolders cantaining FILE_NAME.
    
    Parameters
    ----------
    DIRECTORY_NAME : str
        Name of the directory to analyze.
    FILE_NAME : str
        Name of the file which path is needed
    
    Returns
    -------
    tuple
        Index 0 has a list of the addresses.
        Index 1 has a list of the directories' names.
    """
    # list of files' addresses
    file_addresses = []
    # list of directory's names
    dir_names = []
    # get all files' and folders' names in the indicated directory
    filesAndDirectoryNames = os.listdir(DIRECTORY_NAME)
    # iterate over all the files and folders contained in DIRECTORY_NAME
    for filename in filesAndDirectoryNames:
        # check if the currect object is a folder or not
        if os.path.isdir(os.path.join(DIRECTORY_NAME, filename)):
            # getting folder's name
            dir_names.append(filename)
            # getting path of 'unicycler.log'
            file_addresses.append(os.path.join(DIRECTORY_NAME, filename, FILE_NAME))
    return (file_addresses, dir_names)

def extractor(row, table, headers):
    """
    Reads a row (line) of a table (infile) and convert it into a dictionary.
    The dictionary's keys are the headers. Then, this dictionary is appended
    into a list.
    
    Parameters
    ----------
    row : str
        First row of the table to be process.
    table : iterable object
        File being process.
    headers : list object
        Contains a list of the table's headers
    
    Returns
    -------
    dictionary
        A dictionary of dictionaries. The key of the main dictionary is the
        Lenght of the molecule that correspond to every row the analyzed table
    """
    extracted_table = []
    # iterate over file (table) to extract data
    for row in table:
        # break when reach the end of the table
        if row == '\n':
            break
        # if 'none found' in row replace with 'none_found'
        if 'none found' in row:
            row = re.sub('none found', 'none_found', row)
        # replace line's spaces with tabs and convert line into a list
        line_list = re.sub('\s+', '\t', row.strip()).split('\t')
        # if data in first column is not numberic don't get info
        if not(line_list[0].isnumeric()):
            continue
        # convert list into dictionary using 'headers' as key and append the
        # dictionary to the extracted_table list
        extracted_table.append(dict(zip_longest(headers, line_list)))
    # convert list of dictionaries into dictionary of dictionaries
    final = {}
    for index in extracted_table:
        final[index.get('Length')] = index
    return final

def molecules(file_addresses, dir_names):
    """
    Creates a csv file with relevant information retrieved from all the unicycler.log
    files contined in all primary subfolders of a given directory
    """
    # opening the outfile to save the summary
    outfile = open('molecules_summary.csv', 'a')
    # headers of table molecules summary 
    fieldnames = ['Folder_name','Component', 'Segments', 'Links', 'Length', 'N50',
                  'Longest_segment', 'Status', 'Depth', 'Starting_gene', 'Position',
                  'Strand', 'Identity', 'Coverage']
    # creting an object to write the csv file
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    # write header of table in outfile
    writer.writeheader()
    # iterate over each directory
    for i in range(len(file_addresses)):
        # opening log file
        with open(file_addresses[i], 'r') as log_file:
            # iterate over log file
            for line in log_file:
                # if 'Component' and 'Status' are found in line extract table status
                if re.search('^Component.*Status', line):
                    # convert header 'Longest segment' into 'Longest_segment'
                    headers = re.sub('Longest segment', 'Longest_segment', line)
                    # replace line's spaces with tabs and convert headers into a list
                    headers = re.sub('\s+', '\t', headers.strip()).split('\t')
                    status = extractor(line, log_file, headers)
                # if 'Segment' and 'Depth' are found in line extract table depth
                if re.search('^Segment.*Depth', line):
                    # convert header 'Starting gene' into 'Starting_gene'
                    headers = re.sub('Starting gene', 'Starting_gene', line)
                    # replace line's spaces with tabs and convert headers into a list
                    headers = re.sub('\s+', '\t', headers.strip()).split('\t')
                    depth = extractor(line, log_file, headers)
            # saving relevant info from status and depth lists into a dictionary
            for key in status:
                # if Length from status table is in depth table get Depth
                if key in depth:
                    # saving relevant information in dictionary
                    relevant = {'Folder_name': dir_names[i],
                                'Component': status.get(key).get('Component'),
                                'Segments': status.get(key).get('Segments'),
                                'Links': status.get(key).get('Links'),
                                'Length': status.get(key).get('Length'),
                                'N50': status.get(key).get('N50'),
                                'Longest_segment': status.get(key).get('Longest_segment'),
                                'Status': status.get(key).get('Status'),
                                'Depth': depth.get(key).get('Depth'),
                                'Starting_gene': depth.get(key).get('Starting_gene'),
                                'Position': depth.get(key).get('Position'),
                                'Strand': depth.get(key).get('Strand'),
                                'Identity': depth.get(key).get('Identity'),
                                'Coverage': depth.get(key).get('Coverage')}
                # otherwhise put 'None' in Depth
                else:
                    # saving relevant information in dictionary
                    relevant = {'Folder_name': dir_names[i],
                                'Component': status.get(key).get('Component'),
                                'Segments': status.get(key).get('Segments'),
                                'Links': status.get(key).get('Links'),
                                'Length': status.get(key).get('Length'),
                                'N50': status.get(key).get('N50'),
                                'Longest_segment': status.get(key).get('Longest_segment'),
                                'Status': status.get(key).get('Status'),
                                'Depth': None,
                                'Starting_gene': None,
                                'Position': None,
                                'Strand': None,
                                'Identity': None,
                                'Coverage': None}
                # saving relevant information in outfile
                writer.writerow(relevant)
    # closing outfile file
    outfile.close()

def assemblies(file_addresses, dir_names):
    """
    Explain fuction
    """
    # opening the outfile to save the summary
    outfile = open('assemblies_summary.csv', 'a')
    # headers of table assemblies summary
    fieldnames = ['Folder_name', 'K-mer_best','Contigs_best', 'Dead_ends_best',
                  'Score_best', 'Total_read_count', 'Fully_aligned_reads',
                  'Partially_aligned_reads', 'Unaligned_reads',
                  'Total_bases_aligned', 'Mean_alignment_identity']
    # creting an object to write the csv file
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    # write header of table in outfile
    writer.writeheader()
    # iterate over each directory
    for i in range(len(file_addresses)):
        # opening log file
        with open(file_addresses[i], 'r') as log_file:
            # iterate over log file
            for line in log_file:
                # if 'K-mer', 'Contigs', 'Dead ends' and 'Score' are found in line extract table
                if re.search('^K-mer.*Contigs.*Dead ends.*Score', line):
                    # convert header 'Dead ends' into 'Dead_ends'
                    headers = re.sub('Dead ends', 'Dead_ends', line)
                    # replace line's spaces with tabs and convert headers into a list
                    headers = re.sub('\s+', '\t', headers.strip()).split('\t')
                    # Looking for the best in table
                    for best in log_file:
                        if 'best' in best:
                            # replace line's spaces with tabs and convert line into a list
                            best = re.sub('\s+', '\t', best.strip()).split('\t')
                            break
                # if 'Read alignment summary' in line extract table
                if re.search('Read alignment summary', line):
                    # list to save info
                    alignment_summary_list = []
                    for alignment_summary in log_file:
                        if alignment_summary == '\n':
                            break
                        if '--' in alignment_summary:
                            continue
                        # replace single line's spaces with '_'
                        alignment_summary = re.sub(r'([^\s])(\s)([^\s])', r'\1_\3', alignment_summary)
                        # replace multiple line's spaces with '\t' and conver line in list
                        alignment_summary = (re.sub(r'\s+', r'\t', alignment_summary)).split('\t')
                        # extract relevant data
                        alignment_summary_list.append(alignment_summary[1])
        # write relevant info in outfile
        writer.writerow({'Folder_name': dir_names[i],
                         'K-mer_best': best[0],
                         'Contigs_best': best[1],
                         'Dead_ends_best': best[2],
                         'Score_best': best[3],
                         'Total_read_count': alignment_summary_list[0],
                         'Fully_aligned_reads': alignment_summary_list[1],
                         'Partially_aligned_reads': alignment_summary_list[2],
                         'Unaligned_reads': alignment_summary_list[3],
                         'Total_bases_aligned': alignment_summary_list[4],
                         'Mean_alignment_identity': alignment_summary_list[5]})
    # closing outfile file
    outfile.close()

# the '.' means current directory
DIRECTORY_NAME = '.'
# getting a list of the file addresses and directory names
file_addresses = dirNamesAndAddresses(DIRECTORY_NAME, 'unicycler.log')[0]
dir_names = dirNamesAndAddresses(DIRECTORY_NAME, 'unicycler.log')[1]
# making new files with extracted information from unicycler.log
molecules(file_addresses, dir_names)
assemblies(file_addresses, dir_names)
