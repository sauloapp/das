from datetime import date
from enum import Enum, auto
import pandas
import json
import re
from datetime import datetime
#from sql_reader import AtomTypes


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

SKIP_FULL_TABLE_COVERAGE_CHECK = True

# saulo: GO dictionary for creating correct GO nodes type
with open("dict_data/go-namespace.json", "r") as ns:
    go_plus_dict = json.load(ns)

class Table:

    def __init__(self, name):
        self.header = None
        self.rows = []
        #self.dataframe = None
        self.values = {}
        self.name = name
        print("Table {} created at {}".format(name, str(datetime.now())))
        self.covered_by = {}
        self.mapped_fields = set()
        self.unmapped_fields = set()
        self.mapping = {}
        #print("\nTable::::CONSTRUCTOR:::::::::INSTANTIATIING "+name)
        # Flybase ids ("uniquenames"): FBxxDDDDDDD  x E [a-z], D E [0-9]
        #self.flybase_id_re = re.compile("^(\S+:)?(FB[a-zA-Z]{2}[0-9]{5,10})$")
        self.flybase_id_re = re.compile("^(\S+:)?(FB[a-z]{2}[0-9]{5,10})$")


    # saulo
    def _dataframe_from_header_rows(self):
        dataframe = pandas.DataFrame(self.rows, columns=self.header)
        return dataframe

    # saulo
    def _preprocess(self):
        return

        # saulo
    def _is_list_column(self, column):
        return False

    def _is_ontology_column(self, column):
        return False

    def _ontology_node_type(self, term_id):
        return None #

    def set_header(self, header):
        self.header = [h.strip() for h in header]
        print(self.header)
        for key in self.header:
            assert key
            self.values[key] = set()
            self.covered_by[key] = {}
            self.unmapped_fields.add(key)
            #print("table:set_header_key: "+key)
        assert len(self.unmapped_fields) == len(self.header)


    def process_row_value(self, v):
        v = v.strip()
        m = self.flybase_id_re.search(v)
        if m is not None:
            v = m.group(2)
        return v


    def add_row(self, pre_row):
        row = [self.process_row_value(value) for value in pre_row]
        assert len(self.header) == len(row), f"header = {self.header} row = {row}"
        self.rows.append(row)
        for key, value in zip(self.header, row):
            if value:
                self.values[key].add(value)
                if value not in self.covered_by[key]:
                    self.covered_by[key][value] = set()


    def print_values(self):
        for key in self.values.keys():
            print(f"{key}: {self.values[key]}")


    def get_relevant_sql_tables(self):
        return set([sql_table for sql_table, _ in self.mapping.values()]
                      + ["public.grp_synonym", "public.feature_synonym"])


    def check_field_value(self, sql_table, sql_field, value):
        # saulo 2023/09/03
        if sql_table == None or sql_field == None or value == None:
            print(f"table {str(sql_table)}, column {str(sql_field)}, value {value}")
        for key, values in self.values.items():
            #saulo 2023/09/03
            if key == None or  values == None or value == None:
                print(f"column {str(key)}, values {str(values)}, value {value}")
                continue
            if key in self.unmapped_fields and value in values:
                tag = tuple([key, value])
                sql_tag = tuple([sql_table, sql_field])
                self.covered_by[key][value].add(sql_tag)
                if not SKIP_FULL_TABLE_COVERAGE_CHECK:
                    if all(sql_tag in s for s in self.covered_by[key].values()):
                        self.unmapped_fields.remove(key)
                        self.mapped_fields.add(key)
                        self.mapping[key] = sql_tag




    def print_near_match(self):
        for key in self.unmapped_fields:
            tag_count = {}
            for value in self.covered_by[key]:
                for sql_tag in self.covered_by[key][value]:
                    if sql_tag not in tag_count:
                        tag_count[sql_tag] = 0
                    tag_count[sql_tag] += 1
            for tag in tag_count:
                if (tag_count[tag] / len(self.values[key])) >= 0.8:     #?????????????????????????
                    table, field = tag
                    print(f"{(tag_count[tag] / len(self.values[key]))};{self.name};{key};{table};{field}")


    def check_near_match(self):
        finished = []
        for key in self.unmapped_fields:
            tag_count = {}
            max_count = 0
            max_tag = None
            for value in self.covered_by[key]:
                for sql_tag in self.covered_by[key][value]:
                    if sql_tag not in tag_count:
                        tag_count[sql_tag] = 0
                    tag_count[sql_tag] += 1
                    if tag_count[sql_tag] > max_count:
                        max_count = tag_count[sql_tag]
                        max_tag = sql_tag
            if max_count > 0 and max_count >= (0.9 * len(self.values[key])):        #?????????????????????????
                finished.append(tuple([key, max_tag]))
        for key, tag in finished:
            self.unmapped_fields.remove(key)
            self.mapped_fields.add(key)
            self.mapping[key] = tag

    def all_fields_mapped(self):
        return len(self.unmapped_fields) == 0


########################################################################################################################
class Dmel_enzyme_data_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _convert_ec_format(self, ec_number):
        '''
        From "2.2.-.-" to "2.2", "3.-.-.-" to "3." "5.1.2.4"  to "EC 5.1.2.4", for example
        :param ec_number:
        :return: string in correct ExplorEnz format.
        '''
        parts = ec_number.split('.')
        non_dash_parts = [part for part in parts if part != '-']
        if len(non_dash_parts) == 1:
            non_dash_parts.append('')
        new_ec_number = '.'.join(non_dash_parts)
        if new_ec_number.count('.') == 3 and len(new_ec_number) >= 7:
            return "EC " + new_ec_number
        return new_ec_number

    def _process_piped_string(self, input_string):
        parts = input_string.split("|")
        processed_parts = [self._convert_ec_format(part) for part in parts]
        result = '|'.join(processed_parts)
        return result

    def _preprocess(self):
        df = self._dataframe_from_header_rows()
        df["gene_group_EC_number(s)"] = df["gene_group_EC_number(s)"].apply(lambda x: self._process_piped_string(x) if x else x)
        df["gene_EC_number(s)"] = df["gene_EC_number(s)"].apply(lambda x: self._process_piped_string(x) if x else x)

        df = remove_empty_columns(df)
        # these columns are already represented in the ontologies
        df = df.drop(["gene_group_GO_name(s)", "gene_group_EC_name(s)", "gene_EC_name(s)"], axis=1)
        self.header = list(df.columns)
        self.rows = df.values.tolist()
        #print(df)

    def _is_list_column(self, column):
        if column == "gene_group_EC_number(s)" or\
            column == "gene_group_EC_name(s)" or \
            column == "gene_EC_number(s)" or\
            column == "gene_EC_name(s)":
            return True
        return False

    def _is_ontology_column(self, column):
        return column in ['gene_group_GO_id(s)', 'gene_group_EC_number(s)', 'gene_EC_number(s)']

    def _ontology_node_type(self, term_id):
        if term_id.startswith('GO'):        # Gene Ontology ID
            go_namespace = go_plus_dict[term_id]
            #print(go_namespace)
            return go_namespace
        elif len(term_id) >= 7:
            return AtomTypes.ENZYME
        else:
            return AtomTypes.ENZYME_ONTOLOGY

########################################################################################################################
class Fbrf_pmid_pmcid_doi_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _add_prefix(self, pmid_number):
        return 'PMID:' + pmid_number

    def _insert_colon(self, pmid):
        pmid = pmid.replace("PMC", "PMC:")
        return pmid

    def _preprocess(self):
        df = self._dataframe_from_header_rows()
        df["PMID"] = df["PMID"].apply(lambda x: self._add_prefix(x) if x else x)
        df["PMCID"] = df["PMCID"].apply(lambda x: self._insert_colon(x) if x else x)

        df = remove_empty_columns(df)
        self.header = list(df.columns)
        self.rows = df.values.tolist()

########################################################################################################################
class Gene_association_table(Table):

#remove Aspect column after linking GO
    def __init__(self, name):
        print("CONSTRUCTOR:::::::::INSTANTIATIING " + name)
        super().__init__(name)

    def _get_FB_reference(selfself, db_reference_string):
        return db_reference_string.split('|')[0].split(':')[-1]

    def _preprocess(self):
        df = self._dataframe_from_header_rows()
        df["DB:Reference"] = df["DB:Reference"].apply(lambda x: self._get_FB_reference(x) if x else x)

        df = df.drop(["Aspect"], axis=1)
        df = remove_empty_columns(df)
        self.header = list(df.columns)
        self.rows = df.values.tolist()
        #print(df)

    def _is_list_column(self, column):
        if column == "DB_Object_Synonym":
            return True
        return False

    def _is_ontology_column(self, column):
        return column in ['GO ID']

    def _ontology_node_type(self, term_id):
        go_namespace = go_plus_dict[term_id]
        #print(go_namespace)
        return go_namespace

########################################################################################################################

class Fb_synonym_table(Table):

    def __init__(self, name):
        print("CONSTRUCTOR:::::::::INSTANTIATIING " + name)
        super().__init__(name)

    def _get_FB_reference(selfself, db_reference_string):
        return db_reference_string.split('|')[0].split(':')[-1]

    def _preprocess(self):
        return

    def _is_list_column(self, column):
        if column == "symbol_synonym(s)" or column == "fullname_synonym(s)":
            return True
        return False

########################################################################################################################

class Genotype_phenotype_data_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        return

    def _is_list_column(self, column):
        if column in ["qualifier_ids", "qualifier_names", "genotype_FBids", "genotype_symbols"]:
            return True
        return False

    '''
    Split column value in a list of values [allele (FBal#) or "transposable element insertion" (FBti#)]
    in which the first element is the very column value (the entire string given by genotypes parameter).
    '''

    def _columns_list(self, genotypes):
        split_parts = re.split(r'[ /+]', genotypes)
        col_list = [genotypes] + [part for part in split_parts if part]
        #print("List of alleles::::::::======> " + str(col_list))
        return col_list

    def _is_ontology_column(self, column):
        return column in ['phenotype_id', 'genotype_FBids', "qualifier_ids"]

    def _ontology_node_type(self, term_id):
        if 'FBbt:' in term_id:
            return AtomTypes.FB_ANATOMY_ONTOLOGY
        if 'FBdv:' in term_id:
            return AtomTypes.FB_DEVELOPMENT_ONTOLOGY
        if 'FBcv:' in term_id:
            return AtomTypes.FB_CONTROLLED_VOCABULARY_ONTOLOGY
        if 'GO:' in term_id:
            go_namespace = go_plus_dict[term_id]
            #print(go_namespace)
            return go_namespace
        return None

########################################################################################################################

class Physical_interactions_mitab_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _remove_flybase_str(self, fbgn_string):
        return fbgn_string.split('|')[0].split(':')[-1]

    # Examples of alias:    flybase:ATPsynO(gene name)
    #                       flybase:"l(2)gl"(gene name)
    def _get_symbol_str(self, alias_string):
        #flybase:ATPsynO(gene name)
        return alias_string.split(':')[-1].split('(')[0].strip('"')
        #match = re.search(r'(?<=:)([^"]*|\([^)]*\))(?=\()', alias_string)
        #return match.group(1)

    def _get_FB_reference(self, db_reference_string):
        return db_reference_string.split('|')[0].split(':')[-1]

    def _get_MI_code(self, mi_string):
        try:
            parts = mi_string.split('"')
            return parts[1]
        except Exception as e:
            print(f"mi_string:: {mi_string}.  {e}")
            return ''

    def _extract_ids(self, input_string):
        ids = []
        parts = input_string.split("|")
        for part in parts:
            id_part = part.split(":")[-1]
            ids.append(id_part)
        return "|".join(ids)

    def is_invalid_string(self, string):
        return not ("FBig" in string or "FBlc" in string)

    def _strip_comments(self, string):
        output_string = re.sub(r'comment:', '', string)
        output_string = re.sub(r'\|', '"|"', output_string)
        return output_string.strip('"')

    def _preprocess(self):
        df = self._dataframe_from_header_rows()
        # replace cells equal to "-" by ""
        df = df.applymap(lambda x: '' if str(x) == '-' else x)
        df["ID(s) Interactor A"] = df["ID(s) Interactor A"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        df["ID(s) Interactor B"] = df["ID(s) Interactor B"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        df["Alt ID(s) Interactor A"] = df["Alt ID(s) Interactor A"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        df["Alt ID(s) Interactor B"] = df["Alt ID(s) Interactor B"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        df["Alias(es) Interactor A"] = df["Alias(es) Interactor A"].apply(lambda x: self._get_symbol_str(x) if x else x)
        df["Alias(es) Interactor B"] = df["Alias(es) Interactor B"].apply(lambda x: self._get_symbol_str(x) if x else x)
        df["Publication ID(s)"] = df["Publication ID(s)"].apply(lambda x: self._get_FB_reference(x) if x else x)
        df["Interaction Xref(s)"] = df["Interaction Xref(s)"].apply(lambda x: self._extract_ids(x) if x else x)

        df["Annotation(s) Interactor A"] = df["Annotation(s) Interactor A"].apply(lambda x: self._strip_comments(x) if x else x)
        df["Annotation(s) Interactor B"] = df["Annotation(s) Interactor B"].apply(lambda x: self._strip_comments(x) if x else x)
        df["Interaction Annotation(s)"] = df["Interaction Annotation(s)"].apply(lambda x: self._strip_comments(x) if x else x)

        # These columns will be linked to PSI-MI Ontology
        mi_columns = [
            "Interaction Detection Method(s)",
            "Interaction Type(s)",
            "Source Database(s)",
            "Experimental Role(s) Interactor A",
            "Experimental Role(s) Interactor B",
            "Type(s) Interactor A",
            "Type(s) Interactor B"
        ]
        for mi_column in mi_columns:
            df[mi_column] = df[mi_column].apply(lambda x: self._get_MI_code(x) if x else x)

        df = remove_empty_columns(df)
        self.header = list(df.columns)
        self.rows = df.values.tolist()


    def _is_ontology_column(self, column):
        mi_columns = [
            "Interaction Detection Method(s)",
            "Interaction Type(s)",
            "Source Database(s)",
            "Experimental Role(s) Interactor A",
            "Experimental Role(s) Interactor B",
            "Type(s) Interactor A",
            "Type(s) Interactor B"
        ]
        return column in mi_columns


    def _ontology_node_type(self, term_id):
        return AtomTypes.MOLECULAR_INTERACTION_ONTOLOGY

########################################################################################################################

class Pathway_group_data_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        return

########################################################################################################################

class Gene_group_data_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        return


########################################################################################################################

class Gene_groups_HGNC_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        return


########################################################################################################################

class Gene_map_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        return


########################################################################################################################

class Disease_model_annotations_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        df = self._dataframe_from_header_rows()

        df = df.drop(["DO term"], axis=1) # this column is in the ontologies data as "DiseaseOntology"
        df = remove_empty_columns(df)
        self.header = list(df.columns)
        self.rows = df.values.tolist()
        # print(df)

    def _is_ontology_column(self, column):
        return column in ['DO ID']

    def _ontology_node_type(self, term_id):
        return AtomTypes.DISEASE_ONTOLOGY



########################################################################################################################

class Gene_genetic_interactions_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        return


########################################################################################################################

class Dmel_unique_protein_isoforms_table(Table):
    def __init__(self, name):
        super().__init__(name)

    def _is_list_column(self, column):
        if column == "identical_protein(s)":
            return True
        return False

########################################################################################################################

class Dmel_gene_sequence_ontology_annotations_table(Table):

    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        df = self._dataframe_from_header_rows()
        df = df.drop(["so_term_name"], axis=1)
        df = remove_empty_columns(df)
        self.header = list(df.columns)
        self.rows = df.values.tolist()
        # print(df)

    def _is_ontology_column(self, column):
        return column in ['so_term_id']

    def _ontology_node_type(self, term_id):
        return AtomTypes.SEQUENCE_ONTOLOGY

########################################################################################################################

class Dmel_human_orthologs_disease_table(Table):
    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        df = self._dataframe_from_header_rows()
        df = df.drop(["OMIM_Phenotype_IDs"], axis=1)
        df = remove_empty_columns(df)
        self.header = list(df.columns)
        self.rows = df.values.tolist()

########################################################################################################################

class ScRNA_Seq_gene_expression_table(Table):
    def __init__(self, name):
        super().__init__(name)

    def _preprocess(self):
        df = self._dataframe_from_header_rows()
        #num_unique_values = df['uniquename'].nunique()
        #print("Número de valores únicos na coluna 'uniquename':", num_unique_values)

        #df = df.drop(["Cluster_Cell_Type_Name"], axis=1)
        df = remove_empty_columns(df)
        self.header = list(df.columns)
        self.rows = df.values.tolist()

    def _is_ontology_column(self, column):
        return column in ['Cluster_Cell_Type_ID', 'Cluster_Cell_Type_Name']

    def _ontology_node_type(self, term):
        return AtomTypes.FB_ANATOMY_ONTOLOGY


########################################################################################################################

class NcRNA_genes_main_table(Table):

    def __init__(self, name):
        super().__init__(name)

    '''
    def _preprocess(self):
        df = self._dataframe_from_header_rows()
        df = remove_empty_columns(df)
        self.header = list(df.columns)
        self.rows = df.values.tolist()
        print(df)
    '''

    def _is_ontology_column(self, column):
        return column in ["soTermId"]

    def _ontology_node_type(self, term_id):
        return AtomTypes.SEQUENCE_ONTOLOGY



########################################################################################################################
########################################################################################################################



def instantiate_table(data_file_name):
    if "Dmel_enzyme_data_fb" in data_file_name:
        return Dmel_enzyme_data_table(data_file_name)
    elif "dmel_gene_sequence_ontology_annotations_fb" in data_file_name:
        return Dmel_gene_sequence_ontology_annotations_table(data_file_name)
    elif "dmel_human_orthologs_disease_fb" in data_file_name:
        return Dmel_human_orthologs_disease_table(data_file_name)
    elif "dmel_unique_protein_isoforms_fb" in data_file_name:
        return Dmel_unique_protein_isoforms_table(data_file_name)
    elif "disease_model_annotations_fb" in data_file_name:
        return Disease_model_annotations_table(data_file_name)
    elif "fbrf_pmid_pmcid_doi_fb" in data_file_name:
        return Fbrf_pmid_pmcid_doi_table(data_file_name)
    elif "fb_synonym_fb" in data_file_name:
        return Fb_synonym_table(data_file_name)
    elif "gene_association.fb" in data_file_name:
        return Gene_association_table(data_file_name)
    elif "gene_genetic_interactions_fb" in data_file_name:
        return Gene_genetic_interactions_table(data_file_name)
    elif "gene_group_data_fb" in data_file_name:
        return Gene_group_data_table(data_file_name)
    elif "gene_map_table_fb" in data_file_name:
        return Gene_map_table(data_file_name)
    elif "genotype_phenotype_data_fb" in data_file_name:
        return Genotype_phenotype_data_table(data_file_name)
    elif "gene_groups_HGNC" in data_file_name:
        return Gene_groups_HGNC_table(data_file_name)
    elif "pathway_group_data_fb" in data_file_name:
        return Pathway_group_data_table(data_file_name)
    elif "physical_interactions_mitab_fb" in data_file_name:
        return Physical_interactions_mitab_table(data_file_name)
    elif "scRNA-Seq_gene_expression_fb" in data_file_name:
        return ScRNA_Seq_gene_expression_table(data_file_name)
    elif "ncRNA_genes" == data_file_name:
        return NcRNA_genes_main_table(data_file_name)
    elif data_file_name in ["ncRNA_genes_synonyms", "ncRNA_genes_cross_references",
                            "ncRNA_genes_related_sequences", "ncRNA_genes_gene_synonyms",
                            "ncRNA_genes_publications", "ncRNA_genes_genome_locations"]:
        #print(f"J------------------------SON table: No preprocessing for  table {data_file_name}...")
        return Table(data_file_name)
    else:
        print(f"No preprocessing for  table {data_file_name}...")
        return Table(data_file_name)

def remove_empty_columns(df):
    # Print the names of empty columns
    empty_columns = [column for column in df.columns if df[column].dropna().empty]
    #data_columns = [column for column in df.columns if not df[column].dropna().empty]
    df = df.drop(empty_columns, axis=1)

    return df