from numpy import full
from simple_ddl_parser import parse_from_file, DDLParser
from pathlib import Path
import os, shutil
import subprocess
from enum import Enum, auto
from netAct_table_types import *
from netAct_precomputed_tables import PrecomputedTables
import sqlparse
import re
from datetime import datetime

EXPRESSIONS_PER_CHUNK = 100000000

# MeTTa structures are written to disk. Look at emit_precomputed_tables
# and checkpoint methods

SQL_FILE = "/home/saulo/snet/hyperon/das/das/flybase2metta/fb_data/2023_04/FB2023_04.sql"
SQL_FILE = "/home/saulo/snet/hyperon/das/das/flybase2metta/fb_data/2023_05/FB2023_05.sql"
# SQL_FILE = "/mnt/hdd_2/saulo/snet/hyperon/das/data/flybase/input/2023_04/FB2023_04.sql"

PRECOMPUTED_DIR = "/home/saulo/snet/hyperon/das/das/flybase2metta/fb_data/2023_04/precomputed"
PRECOMPUTED_DIR = "/home/saulo/snet/hyperon/das/das/flybase2metta/fb_data/2023_05/precomputed"
#PRECOMPUTED_DIR = "/mnt/hdd_2/saulo/snet/hyperon/das/data/flybase/input/2023_05/precomputed"

OUTPUT_DIR = "/home/saulo/snet/hyperon/das/das/flybase2metta/fb_data/2023_04/flybase_metta"
OUTPUT_DIR = "/home/saulo/snet/hyperon/das/das/flybase2metta/fb_data/2023_05/flybase_metta"
# OUTPUT_DIR = "/mnt/hdd_2/saulo/snet/hyperon/das/data/flybase/output/2023_04/metta"

SHOW_PROGRESS = True

FILE_SIZE = 657572301  # _file_line_count(SQL_FILE)  # for 2023_04 release sql has  lines
# FILE_SIZE = _file_line_count(SQL_FILE)  # for 2023_04 release sql has  lines
FILE_SIZE = 659729808  # _file_line_count(SQL_FILE)  # for 2023_05 release sql has those number of lines


SCHEMA_ONLY = False
SKIP_PRECOMPUTED_MATCH_BUILD = True
USE_PRECOMPUTED_NEAR_MATCHES = False
PRINT_PRECOMPUTED_NEAR_MATCHES = False
SKIP_SQL_JOIN = False
SKIP_FKEY_FOLLOWING = SKIP_SQL_JOIN or False


def _file_line_count(file_name):
    output = subprocess.run(["wc", "-l", file_name], stdout=subprocess.PIPE)
    return int(output.stdout.split()[0])


if SHOW_PROGRESS:
    print("Checking SQL file size...")


class AtomTypes(str, Enum):
    BIOLOGICAL_PROCESS = "BiologicalProcess"
    CELLTYPE = "Cell"
    CELLULAR_COMPONENT = "CellularComponent"
    CHEBI = "Chebi"
    CHEBI_ONTOLOGY = "ChebiOntology"
    CONCEPT = "Concept"
    DISEASE_ONTOLOGY = "DiseaseOntology"
    ECO_ONTOLOGY = "EcoOntology"
    ENZYME = "Enzyme"
    ENZYME_ONTOLOGY = "EnzymeOntology"
    EVALUATION = "Evaluation"
    EXECUTION = "Execution"
    FB_ANATOMY_ONTOLOGY = "FbAnatomyOntology"
    FB_CONTROLLED_VOCABULARY_ONTOLOGY = "FbControlledVocabularyOntology"
    FB_DEVELOPMENT_ONTOLOGY = "FbDevelopmentOntology"
    INHERITANCE = "Inheritance"
    LIST = "List"
    MOLECULAR_FUNCTION = "MolecularFunction"
    MOLECULAR_INTERACTION_ONTOLOGY = "MolecularInteractionOntology"
    NUMBER = "Number"
    PREDICATE = "Predicate"
    SCHEMA = "Schema"  # p execution link
    SEQUENCE_ONTOLOGY = "SequenceOntology"
    UBERON = "Uberon"
    VERBATIM = "Verbatim"
    CVTERM = "cvterm"
    DATABASE = "db"
    DBXREF = "feature"
    FEATURELOC = "featureloc"
    FEATUREPROP = "featureprop"
    FEATURE_SYNONYM = "feature_synonym"
    GROUP = "grp"
    GROUP_SYNONYM = "grp_synonym"
    LIBRARY = "library"
    ORGANISM = "organism"
    PUB = "pub"
    PUBPROP = "pubprop"
    SYNONYM = "synonym"


# Types that need different formatting: see _build_node()
TYPED_NAME = [AtomTypes.CONCEPT, AtomTypes.SCHEMA]

CREATE_TABLE_PREFIX = "CREATE TABLE "
CREATE_TABLE_SUFFIX = ");"
ADD_CONSTRAINT_PREFIX = "ADD CONSTRAINT "
PRIMARY_KEY = " PRIMARY KEY "
FOREIGN_KEY = " FOREIGN KEY "
COPY_PREFIX = "COPY "
COPY_SUFFIX = "\."


# remove  "fb_2023_04.tsv", for instance.

def _clear_table_name(file_name):
    if ".tsv" in file_name:
        return re.sub(r"_fb_\d{4}_\d\d\.tsv", "", file_name)
    else:
        return re.sub(r".fb", "", file_name)  # added by saulo to handle gene_association.fb table file


class State(int, Enum):
    WAIT_KNOWN_COMMAND = auto()
    READING_CREATE_TABLE = auto()
    READING_COPY = auto()


def non_mapped_column(column):
    # saulo
    # flybase db columns to not retrieve data from:
    # "residues" (public.residues), (public.featureloc) holds biological sequences
    non_mapped_columns = ["md5checksum", "Confidence Value(s)", "Expansion Method(s)",
                          "residues", "residue_info"]
    return column.startswith("time") or "timestamp" in column or (column in non_mapped_columns)


def filter_field(line):
    return \
            "timestamp" in line or \
            "CONSTRAINT" in line


def _compose_name(str_list):
    return "_".join(str_list).replace(" ", "_")

# removes "public"
def short_name(long_table_name):
    return long_table_name.split(".")[1] if long_table_name is not None else None


class LazyParser():

    def __init__(self, sql_file_name, precomputed=None):
        self.sql_file_name = sql_file_name
        self.table_schema = {}
        self.current_table = None
        self.current_table_header = None
        self.current_output_file_number = 1
        base_name = sql_file_name.split("/")[-1].split(".")[0]
        self.target_dir = f"{OUTPUT_DIR}/{base_name}"
        self.current_output_file = None
        self.error_file_name = f"{OUTPUT_DIR}/{base_name}_errors.txt"
        self.error_file = None
        self.schema_file_name = f"{OUTPUT_DIR}/{base_name}_schema.txt"
        self.precomputed_mapping_file_name = f"{OUTPUT_DIR}/{base_name}_precomputed_tables_mapping.txt"
        self.errors = False
        self.current_node_set = set()
        self.current_typedef_set = set()
        self.current_link_list = []
        self.all_types = set()
        self.current_field_types = {}
        self.discarded_tables = []
        self.line_count = None
        self.precomputed = precomputed
        self.relevant_tables = None
        self.expression_chunk_count = 0
        self.relevant_fkeys = {}
        # clear target dir
        Path(self.target_dir).mkdir(parents=True, exist_ok=True)
        for filename in os.listdir(self.target_dir):
            file_path = os.path.join(self.target_dir, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(e)

    def _print_progress_bar(self, iteration, total, length, step, max_step):
        filled_length = int(length * iteration // total)
        previous = int(length * (iteration - 1) // total)
        if iteration == 1 or filled_length > previous or iteration >= total:
            percent = ("{0:.0f}").format(100 * (iteration / float(total)))
            fill = '█'
            bar = fill * filled_length + '-' * (length - filled_length)
            print(f'\r STEP {step}/{max_step} Progress: |{bar}| {percent}% complete ({iteration}/{total})', end='\r')
            if iteration >= total:
                print()

    def _table_info(self, table_name):
        answer = [table_name]
        table = self.table_schema[table_name]
        for column in table['columns']:
            prefix = "  "
            suffix = ""
            if column['name'] == table['primary_key']:
                prefix = "PK"
            elif column['name'] in table['foreign_keys']:
                prefix = "FK"
                referenced_table, referenced_field = table['foreign_key'][column['name']]
                suffix = f"-> {referenced_table} {referenced_field}"
            answer.append(f"    {prefix} {column['type']} {column['name']} {suffix}")
        return "\n".join(answer)

    def _error(self, message):
        self.error_file.write(message)
        self.error_file.write("\n")
        self.errors = True

    def _emit_file_header(self):
        # metta
        for t in AtomTypes:
            self.current_output_file.write(f"(: {t.value} Type)\n")
        self.current_output_file.flush()

    def _open_new_output_file(self):
        if self.current_output_file_number > 1:
            self.current_output_file.close()
        fname = f"{self.target_dir}/file_{str(self.current_output_file_number).zfill(3)}.metta"
        self.current_output_file_number += 1
        self.current_output_file = open(fname, "w")
        self._emit_file_header()

    def _list_values(self, value, separator='|'):
        return value.split(separator)

    def _emit_precomputed_tables(self):
        # def _emit_precomputed_tables(self, output_file):
        # aqui...
        table_count = 0
        processed_rows = 0
        for table in self.precomputed.all_tables:
            if SHOW_PROGRESS:
                self._print_progress_bar(table_count, len(self.precomputed.all_tables), 50, 3, 5)
            for row in table.rows:
                # print(row)
                for key1, value1 in zip(table.header, row):  # handle mapped columns of precomputed tables
                    if key1 not in table.mapped_fields:
                        # print(f"::::::::::::::::key1: {key1} not in table.mapped_fields. Value of " + str(value1))
                        continue
                    if value1 == "":
                        # print("--------------------------> Lambda found in column " + key1 + ".  Nothing to do!")
                        continue
                    # print(f"key1: {key1} in table.mapped_fields. Value of " + str(value1))
                    # access mapping from precomputed table column to sql
                    # print("1111 -- " + sql_table1 + " 1 --  " + sql_field1)
                    if table._is_list_column(key1):  # no necessary if using lists but..
                        if isinstance(table, Dmel_unique_protein_isoforms_table):
                            value1_list = self._list_values(value1, separator=',')
                        # this is to create a list with phenotype symbols plus each allele/
                        # transposable element insertion (transposon) that compose the phenotype
                        elif isinstance(table, Genotype_phenotype_data_table):
                            #print("Geno: ------------------------------------key1----------------> " + key1 + "(" + value1)
                            # list of ontologies' ids
                            if key1 == "genotype_FBids" or key1 == "genotype_symbols":
                                value1_list = table._columns_list(value1)
                                #print("GENO_IDs/symbs:values1  "+str(value1_list))
                            else:
                                value1_list = self._list_values(value1)
                        else:
                            value1_list = self._list_values(value1)
                    else:
                        value1_list = self._list_values(value1)
                    for key2, value2 in zip(table.header, row):
                        if value2 == "":
                            # print("--------------------------> Lambda found in column " + key2 + ".  Nothing to do!")
                            continue
                        if key2 != key1:
                            sql_table1, sql_field1 = table.mapping[key1]
                            sql_table2, sql_field2 = table.mapping[key2] if key2 in table.mapping else (None, None)
                            # constructs a group of genes hierarchy: only group tables have these column names
                            if key1 == "FB_group_id" and key2 == "Parent_FB_group_id":
                                node1 = self._add_value_node(short_name(sql_table1),
                                                             self._get_type(sql_table1, sql_field1), value1)
                                node2 = self._add_value_node(short_name(sql_table2),
                                                             self._get_type(sql_table2, sql_field2), value2)
                                self._add_inheritance(node1, node2)
                                continue
                            # List of values: _values_list returns a list of values (that could contain only one value) if the value
                            # parameter is a concatenation of symbols/names/etc which are synonyms (in general). FB uses
                            # pipes ('|') to separate symbols. One exception is the precomputed table
                            # [Dmel_unique_protein_isoforms_fb_table (only releases 2023_0(2|3|4) were verified]
                            # If value1 or value2 represent a list of values then they are strings like "val1|val2|val3|..."
                            # sql_table2, sql_field2 = table.mapping[key2] if key2 in table.mapping else (None, None)
                            if table._is_list_column(key2):  # not needed if using lists but...
                                if isinstance(table, Dmel_unique_protein_isoforms_table):
                                    value2_list = self._list_values(value2, separator=',')
                                # this is to create a list with phenotype symbols plus each allele/
                                # transposable element insertion (transposon)
                                elif isinstance(table, Genotype_phenotype_data_table):
                                    #print("Geno: --------------------------key2--------------------------> " + key2 + "(" + value2)
                                    # THIS WOULD BE SOLVED AFTER...
                                    # this list can be built from 4 ontologies. It is handled as linked to non-mapped columns only which includes "qualifier_ids" nodes
                                    #if key2 == "qualifier_names":
                                     #   print("Geno qualifier: ----------------------------------------------------> " + key2 + "(" + value2 + " --key1:  "+ key1)
                                      #  continue
                                    if key2 == "genotype_FBids" or key2 == "genotype_symbols":
                                        value2_list = table._columns_list(value2)  # IN FACT, every table should have this method to simplify this code...
                                        #print("GENO_IDs:values2  " + str(value2_list))
                                    else:
                                        value2_list = self._list_values(value2)
                                else:
                                    value2_list = self._list_values(value2)
                            # print("2:", node2)
                            else:
                                value2_list = self._list_values(value2)
                            for l_value1 in value1_list:  # most commonly these two lists will hold only one element
                                # to link ontology nodes
                                if table._is_ontology_column(key1):
                                    node_type = table._ontology_node_type(l_value1)
                                    if node_type != None:
                                        node1 = self._add_node(node_type, l_value1)
                                    else:
                                        node1 = self._add_value_node(short_name(sql_table1),
                                                                     self._get_type(sql_table1, sql_field1), l_value1)
                                else:
                                    node1 = self._add_value_node(short_name(sql_table1),
                                                                 self._get_type(sql_table1, sql_field1), l_value1)
                                for l_value2 in value2_list:
                                    # I need to check this:
                                    # this test is not totally correct: consider the case (I DON'T KNOW IF THIS HAPPENS
                                    # IN FLYBASE DATA...) in which value1 and value2 repreent genes and key1 and/or
                                    # key2 a type of "regulatory" relation. So it's possible that a gene regulates itself...
                                    # But this is not the general case --> post-processing (?)
                                    # if value1 != l_value2 and value2 != l_value1:
                                    # if node1 == None:
                                    #    node1 = self._add_value_node(short_name(sql_table1), self._get_type(sql_table1, sql_field1), l_value1)
                                    #print(f"sql table1: {sql_table1}, sql field1: {sql_field1}. Value:::::::::{l_value1}")
                                    #print(f"sql table222222: {sql_table2}, sql field22222: {sql_field2}. Value:::::::::{l_value2}")
                                    # to link ontology nodes
                                    if table._is_ontology_column(key2):
                                        node_type = table._ontology_node_type(l_value2)
                                        if node_type != None:
                                            node2 = self._add_node(node_type, l_value2)
                                        else:
                                            node2 = self._add_value_node(short_name(sql_table2),
                                                                         self._get_type(sql_table2, sql_field2),
                                                                         l_value2)
                                    else:
                                        node2 = self._add_value_node(short_name(sql_table2), self._get_type(sql_table2, sql_field2), l_value2)
                                    schema = self._add_node(AtomTypes.SCHEMA, _compose_name([_clear_table_name(table.name), key2]))
                                    self._add_execution(schema, node1, node2)

                # non-mapped columns linked here
                for key1, value1 in zip(table.header, row):
                    if value1 == "":
                        # print("--------------------------> Lambda found in column " + key1 + ".  Nothing to do!")
                        continue
                    # List of values: _values_list returns a list of values (that could contain only one value) if the value
                    # parameter is a concatenation of symbols/names/etc which are synonyms (in general). FB uses
                    # pipes ('|') to separate symbols. One exception is the precomputed table
                    # Dmel_unique_protein_isoforms_fb_table (only releases 2023_0(2|3|4) were verified) that uses comma
                    # (',') to separate symbols
                    # If value1 or value2, if they represent a list, are strings like "val1|val2|val3|..."
                    if table._is_list_column(key1):  # no necessary if using lists but..
                        if isinstance(table, Dmel_unique_protein_isoforms_table):
                            value1_list = self._list_values(value1, separator=',')
                        # this is to create a list with phenotype symbols plus each allele/
                        # transposable element insertion (transposon) that compose the phenotype
                        elif isinstance(table, Genotype_phenotype_data_table):
                            #print("Geno: ----------------------------------------------------> " + key1 + "(" + value1)
                            if key1 == "genotype_FBids" or key1 == "genotype_symbols":
                                value1_list = table._columns_list(value1)
                            else:
                                value1_list = self._list_values(value1)
                        else:
                            value1_list = self._list_values(value1)
                    else:
                        value1_list = self._list_values(value1)
                    for key2, value2 in zip(table.header, row):
                        if value2 == "":
                            # print("--------------------------> Lambda found in column " + key2 + ".  Nothing to do!")
                            continue
                        if key2 != key1:
                            if table._is_list_column(key2):
                                if isinstance(table, Dmel_unique_protein_isoforms_table):
                                    value2_list = self._list_values(value2, separator=',')
                                # this is to create a list with phenotype symbols plus each allele/
                                # transposable element insertion (transposon) that compose the phenotype
                                elif isinstance(table, Genotype_phenotype_data_table):
                                    #print("Geno: ------------------------------------> " + key2 + "(" + value2)
                                    if key2 == "genotype_FBids" or key2 == "genotype_symbols":
                                        value2_list = table._columns_list(value2)
                                    else:
                                        value2_list = self._list_values(value2)
                                else:
                                    value2_list = self._list_values(value2)
                            else:
                                value2_list = self._list_values(value2)
                            for l_value1 in value1_list:  # most commonly these two lists will hold only one element
                                for l_value2 in value2_list:
                                    # I need to check this:
                                    # this test is not totally correct: consider the case (I DON'T KNOW IF THIS HAPPENS
                                    # IN FLYBASE DATA...) in which value1 and value2 represent genes and key1 and/or
                                    # key2 a type of "regulatory" relation. So it's possible that a gene regulates itself...
                                    # But this is not the general case --> post-processing (?)
                                    # if value1 != l_value2 and value2 != l_value1:
                                    schema = self._add_node(AtomTypes.SCHEMA,
                                                            _compose_name([_clear_table_name(table.name), key2]))
                                    # se a coluna key1 ou key2 for de hierarquia então chama _add_node com o tipo apropriado...
                                    if table._is_ontology_column(key1):
                                        node_type = table._ontology_node_type(l_value1)
                                        if node_type != None:
                                            node1 = self._add_node(node_type, l_value1)
                                        else:
                                            node1 = self._add_node(AtomTypes.VERBATIM, l_value1)
                                    else:
                                        node1 = self._add_node(AtomTypes.VERBATIM, l_value1)
                                    if table._is_ontology_column(key2):
                                        node_type = table._ontology_node_type(l_value2)
                                        if node_type != None:
                                            node2 = self._add_node(node_type, l_value2)
                                        else:
                                            node2 = self._add_node(AtomTypes.VERBATIM, l_value2)
                                    else:
                                        node2 = self._add_node(AtomTypes.VERBATIM, l_value2)
                                    self._add_execution(schema, node1, node2)


            table_count += 1
            if SHOW_PROGRESS:
                self._print_progress_bar(table_count, len(self.precomputed.all_tables), 50, 3, 5)
        # exit(9)
        #self._checkpoint(True)



    def _checkpoint(self, create_new):
        if SCHEMA_ONLY:
            return
        for metta_string in self.current_typedef_set:
            self.current_output_file.write(metta_string)
            self.current_output_file.write("\n")
        self.current_output_file.flush()
        for metta_string in self.current_node_set:
            self.current_output_file.write(metta_string)
            self.current_output_file.write("\n")
        self.current_output_file.flush()
        for metta_string in self.current_link_list:
            self.current_output_file.write(metta_string)
            self.current_output_file.write("\n")
        self.current_output_file.flush()
        self.current_node_set = set()
        self.current_typedef_set = set()
        self.current_link_list = []
        self.expression_chunk_count = 0
        if create_new:
            self._open_new_output_file()

    def _setup(self):
        self._open_new_output_file()
        self.error_file = open(self.error_file_name, "w")

    def _tear_down(self):
        self.current_output_file.close()
        self.error_file.close()

    def _create_table(self, text):
        parsed = DDLParser(text).run()
        full_name = f"{parsed[0]['schema']}.{parsed[0]['table_name']}"

        # added by saulo
        # tables like gene.gene, gene.allele, gene_group.pathway exist in the sql data file but are not populated and
        # they don't exist for querying!
        if full_name.startswith("public."):
            self.table_schema[full_name] = parsed[0]
            # saulo
            # print("table inserted in schema : " + full_name)# +"\n"+str(parsed[0]))
            assert len(parsed[0]['primary_key']) <= 1
            '''
            parsed[0]['primary_key'] = None
            parsed[0]['foreign_key'] = {}
            parsed[0]['foreign_keys'] = []
            parsed[0]['fields'] = [column['name'] for column in parsed[0]['columns']]
            parsed[0]['types'] = [column['type'] for column in parsed[0]['columns']]
            for column in parsed[0]['columns']:
                self.all_types.add(f"{column['type']} {column['size']}")
            '''
            self.table_schema[full_name]['primary_key'] = None
            self.table_schema[full_name]['foreign_key'] = {}
            self.table_schema[full_name]['foreign_keys'] = []
            self.table_schema[full_name]['fields'] = [column['name'] for column in parsed[0]['columns']]
            self.table_schema[full_name]['types'] = [column['type'] for column in parsed[0]['columns']]
            for column in parsed[0]['columns']:
                self.all_types.add(f"{column['type']} {column['size']}")

    def _start_copy(self, line):
        self.current_table = line.split(" ")[1]
        if SCHEMA_ONLY or self.current_table in self.discarded_tables or \
                (self.relevant_tables is not None and self.current_table not in self.relevant_tables):
            return False
        columns = line.split("(")[1].split(")")[0].split(",")
        columns = [s.strip() for s in columns]
        schema_columns = [column['name'] for column in self.table_schema[self.current_table]['columns']]
        assert all(column in schema_columns or non_mapped_column(column) for column in columns)
        self.current_table_header = columns
        self.current_field_types = {}
        table = self.table_schema[self.current_table]
        for name, ctype in zip(table['fields'], table['types']):
            self.current_field_types[name] = ctype
        return True

    def _get_type(self, table_name, field):
        if table_name is not None:
            table = self.table_schema[table_name]
            for name, ctype in zip(table['fields'], table['types']):
                if name == field:
                    if name == table['primary_key']:
                        return "pk"
                    else:
                        return ctype
        return "text"

    def _build_node(self, node_type, node_name):
        node_name = node_name.replace("(", "[")
        node_name = node_name.replace(")", "]")
        node_name = node_name.replace('"', "")
        if node_type in TYPED_NAME:
            quoted_node_name = f'"{node_type}:{node_name}"'
            quoted_canonical_node_name = f'"{node_type} {node_type}:{node_name}"'
        else:
            quoted_node_name = f'"{node_name}"'
            quoted_canonical_node_name = f'"{node_type} {node_name}"'
        return tuple([quoted_canonical_node_name, f"(: {quoted_node_name} {node_type})"])

    def _add_node_to_internal_sets(self, quoted_canonical_node_name, node):
        if node not in self.current_node_set:
            self.current_node_set.add(node)
        node_type = quoted_canonical_node_name.strip('"').split()[0]
        self.current_typedef_set.add(f"(: {node_type} Type)")
        self.expression_chunk_count += 1

    def _add_node(self, node_type, node_name):
        # metta
        # print(f"add_node {node_type} {node_name}")
        quoted_canonical_node_name, node = self._build_node(node_type, node_name)
        self._add_node_to_internal_sets(quoted_canonical_node_name, node)
        return quoted_canonical_node_name

    def _add_inheritance(self, node1, node2):
        # metta
        # print(f"add_inheritance {node1} {node2}")
        if node1 and node2:
            # saulo
            link = f"({AtomTypes.INHERITANCE} {node1} {node2})"
            if link not in self.current_link_list:
                self.current_link_list.append(link)
            # self.current_link_list.append(f"({AtomTypes.INHERITANCE} {node1} {node2})")
        self.expression_chunk_count += 1

    # def _add_evaluation(self, predicate, node1, node2):
    #    # metta
    #    #print(f"add_evaluation {predicate} {node1} {node2}")
    #    if predicate and node1 and node2:
    #        self.current_link_list.append(f"({AtomTypes.EVALUATION} {predicate} ({AtomTypes.LIST} {node1} {node2}))")
    #    self.expression_chunk_count += 1

    def _add_execution(self, schema, node1, node2):
        # metta
        # print(f"add_execution {schema} {node1} {node2}")
        if schema and node1 and node2:
            self.current_link_list.append(f"({AtomTypes.EXECUTION} {schema} {node1} {node2})")
        self.expression_chunk_count += 1

    def _add_value_node(self, table_short_name, field_type, value, build_only=False):
        if value == "\\N":  # null for PostgreSQL
            return None
        if field_type == "pk":
            assert table_short_name is not None
            if build_only:
                return self._add_node(table_short_name, value)
            else:
                return self._add_node(table_short_name, value)
        elif field_type == "boolean":
            if build_only:
                return self._add_node(AtomTypes.CONCEPT, "True" if value.lower() == "t" else "False")
            else:
                return self._add_node(AtomTypes.CONCEPT, "True" if value.lower() == "t" else "False")
        elif field_type in ["bigint", "integer", "smallint", "double precision"]:
            if build_only:
                return self._add_node(AtomTypes.NUMBER, value)
            else:
                return self._add_node(AtomTypes.NUMBER, value)
        elif "character" in field_type or field_type in ["date", "text"]:
            if build_only:
                return self._add_node(AtomTypes.VERBATIM, value)
            else:
                return self._add_node(AtomTypes.VERBATIM, value)
        elif field_type in ["jsonb"]:
            return None
        else:
            assert False

    def _new_row_precomputed(self, line):
        if SCHEMA_ONLY:
            return
        table = self.table_schema[self.current_table]
        fkeys = table['foreign_keys']
        data = line.split("\t")
        if len(self.current_table_header) != len(data):
            self._error(
                f"Invalid row at line {self.line_count} Table: {self.current_table} Header: {self.current_table_header} Raw line: <{line}>")
            return
        for name, value in zip(self.current_table_header, data):
            if (not non_mapped_column(name)) and (name not in fkeys):
                self.precomputed.check_field_value(self.current_table, name, value)

    def _new_row_relevant_fkeys(self, line):
        if SCHEMA_ONLY or (self.current_table not in self.relevant_fkeys):
            return
        table = self.table_schema[self.current_table]
        table_short_name = short_name(self.current_table)
        pkey = table['primary_key']
        fkeys = table['foreign_keys']
        assert pkey, f"self.current_table = {self.current_table} pkey = {pkey} \n{table}"
        data = line.split("\t")
        if len(self.current_table_header) != len(data):
            self._error(
                f"Invalid row at line {self.line_count} Table: {self.current_table} Header: {self.current_table_header} Raw line: <{line}>")
            return
        pkey_node = None
        for name, value in zip(self.current_table_header, data):
            if name == pkey:
                if value not in self.relevant_fkeys[self.current_table]:
                    return
                pkey_node = self._add_node(table_short_name, value)
                break
        assert pkey_node is not None
        for name, value in zip(self.current_table_header, data):
            if non_mapped_column(name):
                continue
            if name in fkeys:
                referenced_table, referenced_field = table['foreign_key'][name]
                schema_node = self._add_node(AtomTypes.SCHEMA, referenced_table)
                fkey_node = self._add_node(AtomTypes.CONCEPT, _compose_name([referenced_table, value]))
                self._add_execution(schema_node, pkey_node, fkey_node)
            elif name != pkey:
                ftype = self.current_field_types.get(name, None)
                if not ftype:
                    continue
                value_node = self._add_value_node(table_short_name, ftype, value)
                if not value_node:
                    continue
                schema_node = self._add_node(AtomTypes.SCHEMA, _compose_name([table_short_name, name]))
                self._add_execution(schema_node, pkey_node, value_node)

    def _new_row(self, line):
        if SCHEMA_ONLY or (self.relevant_tables is not None and self.current_table not in self.relevant_tables):
            return
        table = self.table_schema[self.current_table]  # sql table
        table_short_name = short_name(self.current_table)
        pkey = table['primary_key']
        fkeys = table['foreign_keys']
        assert pkey, f"self.current_table = {self.current_table} pkey = {pkey} \n{table}"
        data = line.split("\t")
        if len(self.current_table_header) != len(data):
            self._error(
                f"Invalid row at line {self.line_count} Table: {self.current_table} Header: {self.current_table_header} Raw line: <{line}>")
            return
        relevant_row = False
        for name, value in zip(self.current_table_header, data):
            for precomputed_table in self.precomputed.all_tables:
                for key in precomputed_table.mapping:
                    sql_table, sql_field = precomputed_table.mapping[key]
                    if sql_table == self.current_table and sql_field == name and value in precomputed_table.values[key]:
                        relevant_row = True
                        break
                if relevant_row:
                    break
            if relevant_row:
                break
        if not relevant_row:
            return
        pkey_node = None
        for name, value in zip(self.current_table_header, data):
            if name == pkey:
                pkey_node = self._add_node(table_short_name, value)
                break
        assert pkey_node is not None
        for name, value in zip(self.current_table_header, data):
            if non_mapped_column(name):
                continue
            if name in fkeys:
                referenced_table, referenced_field = table['foreign_key'][name]
                if not SKIP_FKEY_FOLLOWING:
                    if referenced_table not in self.relevant_fkeys:  # stores relevant foreign keys
                        self.relevant_fkeys[referenced_table] = set()
                    self.relevant_fkeys[referenced_table].add(value)
                schema_node = self._add_node(AtomTypes.SCHEMA, _compose_name([referenced_table]))
                fkey_node = self._add_node(AtomTypes.CONCEPT, _compose_name([referenced_table, value]))
                self._add_execution(schema_node, pkey_node, fkey_node)
            elif name != pkey:
                ftype = self.current_field_types.get(name, None)
                if not ftype:
                    continue
                value_node = self._add_value_node(table_short_name, ftype, value)
                if not value_node:
                    continue
                schema_node = self._add_node(AtomTypes.SCHEMA, _compose_name([table_short_name, name]))
                self._add_execution(schema_node, pkey_node, value_node)

    def _primary_key(self, first_line, second_line):
        try:
            line = first_line.split()
            table = line[2] if line[2] != "ONLY" else line[3]
            line = second_line.split()
            field = line[-1][1:-2]
            # saulo: only "public" tables are available for querying (2023_04)
            if table.startswith("public."):
                assert not self.table_schema[table]['primary_key']
                assert field in self.table_schema[table]['fields']
                self.table_schema[table]['primary_key'] = field
                if self.precomputed:
                    self.precomputed.set_sql_primary_key(table, field)
        except Exception as e:
            print("\nPK-->" + str(e) + "PPK -- table: " + table + " PK: " + str(self.table_schema[table]))

    def _foreign_key(self, first_line, second_line):
        try:
            line = first_line.split()
            table = line[2] if line[2] != "ONLY" else line[3]
            line = second_line.split()
            field = line[5][1:-1]
            # saulo: only "public" tables are available for querying (2023_04)
            if table.startswith("public."):
                reference = line[7].split("(")
                referenced_table = reference[0]
                referenced_field = reference[1].split(")")[0]
                assert field in self.table_schema[table]['fields']
                assert referenced_field in self.table_schema[referenced_table]['fields']
                self.table_schema[table]['foreign_key'][field] = tuple([referenced_table, referenced_field])
                self.table_schema[table]['foreign_keys'].append(field)
        except Exception as e:
            print("\nFK-->" + str(e) + "     FK -- table: " + table + " FK: " + str(self.table_schema[table]))

    def _parse_step_1(self):
        print("\nStep 1: " + str(datetime.now()))
        text = ""
        self.line_count = 0
        file_size = FILE_SIZE

        # Searches for SQL commands and takes the adequate action
        state = State.WAIT_KNOWN_COMMAND
        with open(self.sql_file_name, 'r') as file:
            line = file.readline()
            previous_line = None
            while line:
                self.line_count += 1
                if SHOW_PROGRESS:
                    self._print_progress_bar(self.line_count, file_size, 50, 1, 5)
                line = line.replace('\n', '').strip()
                if state == State.WAIT_KNOWN_COMMAND:
                    if line.startswith(CREATE_TABLE_PREFIX):
                        text = line
                        state = State.READING_CREATE_TABLE
                    elif line.startswith(ADD_CONSTRAINT_PREFIX) and PRIMARY_KEY in line:
                        self._primary_key(previous_line, line)
                    elif line.startswith(ADD_CONSTRAINT_PREFIX) and FOREIGN_KEY in line:
                        self._foreign_key(previous_line, line)
                elif state == State.READING_CREATE_TABLE:
                    if not filter_field(line):
                        text = f"{text}\n{line}"
                    if line.startswith(CREATE_TABLE_SUFFIX):
                        self._create_table(text)

                        state = State.WAIT_KNOWN_COMMAND
                        text = ""
                else:
                    print(f"Invalid state {state}")
                    assert False
                previous_line = line
                line = file.readline()

    def _parse_step_2(self):
        print("\nStep 2: " + str(datetime.now()))
        text = ""
        self.line_count = 0
        file_size = FILE_SIZE

        # tables without primary key are discarded.
        for key, table in self.table_schema.items():
            if not table['primary_key']:
                self.discarded_tables.append(key)
                self._error(f"Discarded table {key}. No PRIMARY KEY defined.")

        state = State.WAIT_KNOWN_COMMAND
        with open(self.sql_file_name, 'r') as file:
            line = file.readline()
            if SKIP_PRECOMPUTED_MATCH_BUILD:
                if SHOW_PROGRESS:
                    self._print_progress_bar(file_size, file_size, 50, 2, 5)
            else:
                while line:
                    self.line_count += 1
                    if SHOW_PROGRESS:
                        self._print_progress_bar(self.line_count, file_size, 50, 2, 5)
                    if not self.precomputed.all_tables_mapped():
                        line = line.replace('\n', '').strip()
                        if state == State.WAIT_KNOWN_COMMAND:
                            if line.startswith(COPY_PREFIX):
                                if self._start_copy(line):
                                    state = State.READING_COPY
                        elif state == State.READING_COPY:
                            if line.startswith(COPY_SUFFIX):
                                state = State.WAIT_KNOWN_COMMAND
                            else:
                                self._new_row_precomputed(line)
                        else:
                            print(f"Invalid state {state}")
                            assert False
                    line = file.readline()
            if USE_PRECOMPUTED_NEAR_MATCHES:
                self.precomputed.check_nearly_matched_tables()
            if PRINT_PRECOMPUTED_NEAR_MATCHES:
                self.precomputed.print_matched_tables()
            #  saulo
            # self._emit_precomputed_tables(self.current_output_file)
            self._emit_precomputed_tables()
            # self._checkpoint(True)
        self.relevant_tables = self.precomputed.get_relevant_sql_tables()

    def _parse_step_3(self):
        print("\nStep 3: " + str(datetime.now()))
        text = ""
        self.line_count = 0
        file_size = FILE_SIZE

        if not self.precomputed:
            for key, table in self.table_schema.items():
                if not table['primary_key']:
                    self.discarded_tables.append(key)
                    self._error(f"Discarded table {key}. No PRIMARY KEY defined.")

        state = State.WAIT_KNOWN_COMMAND
        with open(self.sql_file_name, 'r') as file:
            line = file.readline()
            while line:
                self.line_count += 1
                # if self.expression_chunk_count >= EXPRESSIONS_PER_CHUNK:
                # self._checkpoint(True)
                if SHOW_PROGRESS:
                    self._print_progress_bar(self.line_count, file_size, 50, 4, 5)
                line = line.replace('\n', '').strip()
                if state == State.WAIT_KNOWN_COMMAND:
                    if line.startswith(COPY_PREFIX):
                        if self._start_copy(line):
                            state = State.READING_COPY
                elif state == State.READING_COPY:
                    if line.startswith(COPY_SUFFIX):
                        state = State.WAIT_KNOWN_COMMAND
                    else:
                        if not SKIP_SQL_JOIN:
                            self._new_row(line)
                else:
                    print(f"Invalid state {state}")
                    assert False
                line = file.readline()
            # self._checkpoint(False)

    def _parse_step_4(self):
        print("\nStep 4: " + str(datetime.now()))
        text = ""
        self.line_count = 0
        file_size = FILE_SIZE

        if not self.precomputed:
            for key, table in self.table_schema.items():
                if not table['primary_key']:
                    self.discarded_tables.append(key)
                    self._error(f"Discarded table {key}. No PRIMARY KEY defined.")

        state = State.WAIT_KNOWN_COMMAND
        with open(self.sql_file_name, 'r') as file:
            line = file.readline()
            while line:
                self.line_count += 1
                # if self.expression_chunk_count >= EXPRESSIONS_PER_CHUNK:
                # self._checkpoint(True)
                if SHOW_PROGRESS:
                    self._print_progress_bar(self.line_count, file_size, 50, 5, 5)
                line = line.replace('\n', '').strip()
                if state == State.WAIT_KNOWN_COMMAND:
                    if line.startswith(COPY_PREFIX):
                        if self._start_copy(line):
                            state = State.READING_COPY
                elif state == State.READING_COPY:
                    if line.startswith(COPY_SUFFIX):
                        state = State.WAIT_KNOWN_COMMAND
                    else:
                        if not SKIP_SQL_JOIN:
                            self._new_row_relevant_fkeys(line)
                else:
                    print(f"Invalid state {state}")
                    assert False
                line = file.readline()
            # self._checkpoint(False)
            print("\nStep 4 finished: " + str(datetime.now()))

    def parse(self):
        self._setup()
        self._parse_step_1()
        if self.precomputed:
            self._parse_step_2()
            f = open(self.precomputed_mapping_file_name, "w")
            f.write(self.precomputed.mappings_str())
            f.close()
        f = open(self.schema_file_name, "w")
        for table in self.table_schema:
            if self.relevant_tables is None or table in self.relevant_tables:
                print("Schema f:  " + str(table))
                f.write(self._table_info(table))
                f.write("\n\n")
        f.close()
        self._parse_step_3()
        if not SKIP_FKEY_FOLLOWING:
            self._parse_step_4()
        if self.errors:
            print(f"Errors occurred while processing this SQL file. See them in {self.error_file_name}")
        self._checkpoint(False)
        self._tear_down()


def main():
    precomputed = PrecomputedTables(PRECOMPUTED_DIR, "2023_05") if PRECOMPUTED_DIR else None
    parser = LazyParser(SQL_FILE, precomputed)
    parser.parse()


if __name__ == "__main__":
    main()