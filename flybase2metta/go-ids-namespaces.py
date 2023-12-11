#!/usr/bin/env python3
"""
 2023-07-24
 Script to generate a dictionary of (GO_id, GO_namespace) from the go.obo file. 
 * The dictionary is checked in "ids perspective": the ids list from go.obo should be 
 the same as the ids list from go-plus.csv. If they are not, the list from
 the later is kept and the dictionary is built from it.
 The dictionary is stored in the 'go-namespace.json' file and is used by the 
 'go-plus.py' script.

 Requires: file go.obo from https://bioportal.bioontology.org/ontologies/GO
"""


import wget
import os
import json
import pandas as pd
import sys


def equal_dicts(dict1, dict2):
	#if len(dict1) != len(dict2):
	print(f"len obo: {len(dict1)}: len plus: {len(dict2)}", file=sys.stderr)
#    	return False

	diff_pairs = []
	dif = 0
	for key, value in dict1.items():
		if key not in dict2 or dict2[key] != value:
			dif += 1
			diff_pairs.append(("dict1_" + str(dif), key, value))

	for key, value in dict2.items():
		if key not in dict1 or dict1[key] != value:
			dif += 1
			diff_pairs.append(("dict2_" + str(dif), key, value))

	if diff_pairs:
		print("Differing key-value pairs:")
		for dict_num, key, value in diff_pairs:
			print(f"{dict_num} --> {key}: {value}", file=sys.stderr)
		return False

	return True



obo_source = "https://data.bioontology.org/ontologies/GO/submissions/1815/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb"
output = 'dict_data/go-namespace.json'


if not os.path.exists('dict_data/go.obo'):
	wget.download(obo_source,"dict_data/")
go_file = open('dict_data/go.obo')
go_lines = go_file.readlines()


# Get id / namespace from go.obo file.
id_prefix = "id: "
ns_prefix = "namespace: "
ids = []
namespaces = []
i = 0
while i < len(go_lines):
	if go_lines[i].startswith('[Term]'):
		term_info = {}
		i += 1
		while i < len(go_lines) and go_lines[i].strip():
			key, value = go_lines[i].strip().split(': ', 1)
			term_info[key] = value
			i += 1
		#if 'is_obsolete' in term_info:
#	        print("is obsolete")
			#continue
		if 'id' in term_info and 'namespace' in term_info:
			if term_info['namespace'] != "":
				ids.append(term_info['id'])
				if term_info['namespace'] == "biological_process":
					term_info['namespace'] = "BiologicalProcess"
				elif term_info['namespace'] == "cellular_component":
					term_info['namespace'] = "CellularComponent"
				else:
					term_info['namespace'] = MOLECULAR_FUNCTION = "MolecularFunction"
				namespaces.append(term_info['namespace'])
	else:
		i += 1

obo_dict = {}
for i in range( len(ids) ):
	obo_dict[ids[i]] = namespaces[i]

##############################################################################
#source_csv = "https://gitlab.com/opencog-bio/pln_mozi/blob/master/raw_data/GO-PLUS.csv.gz"
source_csv_latest = "http://data.bioontology.org/ontologies/GO-PLUS/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv"

if not os.path.exists("dict_data/GO-PLUS.csv.gz"):
	dataset = wget.download(source_csv_latest, "dict_data")
df = pd.read_csv("dict_data/GO-PLUS.csv.gz", dtype=str)

#df = remove_empty_columns(df)
# Count elements in the 'Class ID' column that contain the substring "PR_"
#count_with_substring = df["Class ID"].str.contains("PR_").sum()
#print(f"Number of elements with 'PR_': {count_with_substring}")  # 482 protein ids to date (2023-06-25)
#exit()

def get_go_id(class_id):
	if str(class_id) == "nan":
		term = False
	else:
		term = class_id.split("/")[-1]
		if term.startswith("GO_"):
			term = term.replace("_",":")
		else:
			term = False
	return term

##############################################################################

df = df[["Class ID", "has_obo_namespace", "Obsolete"]]        
plus_dict = {}

for i in range(len(df)):
	term = get_go_id(df.iloc[i]["Class ID"])
	obsolete = df.iloc[i]["Obsolete"]
	if term != False:
		nspace = df.iloc[i]["has_obo_namespace"]
		if str(nspace) != "nan":
			plus_dict[term] = df.iloc[i]["has_obo_namespace"]

plus_dict.update(obo_dict)

if equal_dicts(obo_dict, plus_dict):
	print("OBO dict == PLUS dict")
	result_dict = obo_dict
else:
	if len(obo_dict) > len(plus_dict):	
		print("OBO dict > PLUS dict")
		result_dict = obo_dict
	else:
		if len(obo_dict) < len(plus_dict):	
			print("OBO dict < PLUS dict")
			result_dict = plus_dict
		else:
			print("OBO dict == PLUS dict")
			result_dict = plus_dict

		
with open(output, "w") as json_file:
	json.dump(plus_dict, json_file, indent=2)


print(f"\nData saved to {output}")

