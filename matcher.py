import pandas as pd
import os
#This short script is reading a file in fasta format which contains multiple sequences
#the longest common substring is searched for within all of the sequences
#a table is exported in the end which contains the longest common substrings
#multiple of the longest common substrings are allowed

def readFile(file_path):
    file_content = []
    gene_names = []
    sequences = []
    current_gene_name = None
    current_sequence = []
    with open(file_path,"r") as f:
        content = f.read().splitlines()
        file_content = content
    # Iterate through the lines
    for line in file_content:
        if line.startswith(">"):
            # If the line starts with '>', it's a gene name
            # Save the previous gene name and sequence (if any)
            if current_gene_name is not None:
                gene_names.append(current_gene_name)
                sequences.append("".join(current_sequence))
            # Extract the gene name (remove the '>' character)
            current_gene_name = line[1:]
            # Initialize an empty list for the new sequence
            current_sequence = []
        else:
            # If the line doesn't start with '>', it's part of the sequence
            # Append it to the current sequence
            current_sequence.append(line)

    # Add the last gene name and sequence to the lists
    if current_gene_name is not None:
        gene_names.append(current_gene_name)
        sequences.append("".join(current_sequence))

    # Create a Pandas DataFrame
    df = pd.DataFrame({"gene": gene_names, "sequence": sequences})
    return df


def makeSubstrings(seq_table):
    #just take the first sequence in the table and make substrings of it
    #as the substring needs to appear in every sequence anyway, it doesnt matter from sequence i derive the substring
    template_string = seq_table.loc[0]["sequence"]
    #buffer size means how long the substring is
    buffer_size = 2
    #get the length of the template string
    template_string_len = len(template_string)
    #the maximal length of a substring can be the length of this template
    max_buffer_size = int(template_string_len)
    #variable to collect all the buffer strings
    buffers = []

    #loop until the end of the max buffer size and increment buffer size every loop
    for k in range(0, max_buffer_size-buffer_size):
        buffer_size += 1
        #loop until the end of the template string
        #every position in the template string + buffer size is selected
        for i in range(0, template_string_len-buffer_size+1):
            buffer = template_string[i:i+buffer_size]
            buffers.append(buffer)

    #this table contains all possible substrings
    substring_table = pd.DataFrame({"buffers": buffers})
    substring_table = substring_table.drop_duplicates(subset='buffers', keep='last')
    return substring_table

def searchSubstrings(seq_table, sub_table):
    #remove duplicates from substring table first!!


    matches = []
    seq_table = seq_table.drop(0).reset_index(drop=True)
    #loop through all substring sequences
    for sub in sub_table['buffers']:

        sub_counter = 0
        #loop through all sequences in the sequence table and look up if the current substring is in one of the sequences
        #increment the sub_counter if the substring was found in any of the sequences
        for seq in seq_table['sequence']:
            is_match = seq.find(sub)
            if is_match != -1:
                sub_counter += 1
            else:
                break
        matches.append({"match": sub, "length": len(sub), "counter": sub_counter})

    #the match table contains which substring was matched how often as well as the substrings length
    match_table = pd.DataFrame(matches)
    #count how many rows there are in the sequence table
    max_rows = seq_table.shape[0]
    #keep only the rows where the number of substrings matched equals the number of given sequences
    #this means that the substring is found in every given sequence, others are "dropped"
    filtered_matches_by_max_counter = match_table[match_table['counter']==max_rows]
    #filter by length to get the longest substring which was found
    filtered_matches_by_max_counter = filtered_matches_by_max_counter.sort_values(by=["length"], ascending=False)
    return filtered_matches_by_max_counter

def exportFile(file_name, table):
    dir = "results/"
    if os.path.exists(dir):
        table.to_csv("{0}{1}.csv".format(dir, file_name), sep=";")
    else:
        os.makedirs(dir)
    print("Results exported to {0}{1}.csv".format(dir, file_name))
seq_table = readFile("sequences.txt")
print("Read file")
sub_table = makeSubstrings(seq_table)
print("Substrings created")
print("Searching for substrings...")
matches = searchSubstrings(seq_table, sub_table)
file_name = input("Please enter file name for export: ")
exportFile(file_name, matches)
