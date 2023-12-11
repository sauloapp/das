import glob
import csv
import json
import os
from precomputed_tables import PrecomputedTables
import re

class TablesHeaderCollector:

    def clear_table_name(self, file_name):
        if ".tsv" in file_name:
            return re.sub(r"_fb_\d{4}_\d\d\.tsv", "", file_name)
        else:
            return re.sub(r".fb", "", file_name)  # added by saulo to handle gene_association.fb table file

    def generate_pairs_input(self, output_file):
        PRECOMPUTED_DIR = "/home/saulo/snet/hyperon/das/das/flybase2metta/fb_data/2023_04/precomputed"
        precomputed = PrecomputedTables(PRECOMPUTED_DIR, "2023_04") if PRECOMPUTED_DIR else None
        input_file = open(output_file, "w")
        for table in precomputed.all_tables:
            table_data = self.clear_table_name(table.name) + '\t'
            for col_num in range(len(table.header)):
                if col_num == len(table.header) - 1:
                    table_data += table.header[col_num] + '\n'
                else:
                    table_data += table.header[col_num] + '\t'
            input_file.write(table_data)


def process_file(file_name):
    relevant_pairs_file = open('relevant_pairs.tsv', "w")
    non_relevant_pairs_file = open('relevant_pairs_not.tsv', "w")
    with open(file_name, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        #next(reader)  # Skip header row

        for row in reader:
            relevant_pairs = []
            non_relevant_pairs = []
            table_name = row[0]
            column_names = row[1:]
            pairs = [(column1, column2) for i, column1 in enumerate(column_names) for column2 in column_names[i + 1:]]
            for pair in pairs:
                print(f"Table: {table_name}, Columns: {pair[0]}, {pair[1]}")
                user_input = input("Is this pair relevant? (y/n): ")

                if user_input.lower() == 'y':
                    relevant_pairs.append(pair)
                else:
                    non_relevant_pairs.append(pair)
            save_pairs_to_file(relevant_pairs, relevant_pairs_file , table_name)
            save_pairs_to_file(non_relevant_pairs, non_relevant_pairs_file, table_name)


def save_pairs_to_file(pairs, pairs_file, table_name):
    pairs_list = table_name + "\t"
    for pair in pairs:
        if pair == pairs[-1]:
            pairs_list += str(pair[0]) + "\t" + str(pair[1]) + "\n"
        else:
            pairs_list += str(pair[0]) + "\t" + str(pair[1]) + "\t"

    print(pairs_list)
    pairs_file.write(pairs_list)
    pairs_file.flush()


pairs_input = "/home/saulo/snet/hyperon/github/das/flybase2metta/input_to_pairing.tsv"
pairs_input = "/home/saulo/snet/hyperon/github/das/flybase2metta/input_to_pairing-lasts.tsv"
# these two lines are to get table names
#collector = TablesHeaderCollector()
#collector.generate_pairs_input(pairs_input)

#file_name = '/home/saulo/snet/hyperon/github/das/flybase2metta/pairs_input.tsv'
process_file(pairs_input)
