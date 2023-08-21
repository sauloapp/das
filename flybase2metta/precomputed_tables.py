from datetime import date
import os
import glob
import csv
import json
import re
import pandas
import table_types





class PrecomputedTables:
    '''
    2023 August: changed to handle more information from "mapping(_2023_04).txt" file:
    1. Data type column
    2. PK (Primary Key) or FK (Foreign Key) columns

    '''
    def __init__(self, dir_name, release_id):
        '''
        release_id is in format yyyy_nn  where:
        yyyy is the year
        nn is the release number in the year
        For example: 2022_06, 2023_2, 2023_04

        :param dir_name:
        :param release:
        '''
        print(dir_name)
        self.all_tables = []
        self.unmapped_tables = {}
        self.mapped_tables = {}
        self.sql_primary_key = {}
        self.sql_tables = None
        self.preloaded_mapping = False
        os.chdir(dir_name)
        # This is to output a tsv given an original mappings file (generated by sql_reader)
        '''
        if os.path.exists(f"{dir_name}/mapping.txt"):
            with open(f"{dir_name}/mapping.txt", "r") as f:
                for line in f:
                    line = line.strip("\n")
                    if not line.startswith("\t"):
                       fname = line
                    else:
                        line = line.strip("\t")
                        n = line.find(" ->")
                        pre = line[:n]
                        pos = line[n+4:]
                        if pos == "???":
                            continue
                        table, field = tuple(pos.split())
                        print("\t".join([fname, pre, table, field]))
        '''
        for file_name in glob.glob("ncRNA_genes_*.json"):
            with open(file_name) as f:
                json_dict = json.load(f)
            self._process_ncrna(json_dict)
        #  + glob.glob("*.fb") added by saulo to handle the "gene_association.fb" table
        for file_name in glob.glob("*.tsv") + glob.glob("*.fb"):
            table = table_types.instantiate_table(file_name)
            self.unmapped_tables[file_name] = table
            self.all_tables.append(table)
            #self._process_tsv(file_name)
            # saulo  --> It should be necessary change constructor parameters to get a dict of (fb_file, column names)
            # gene_association.fb is precomputed since 2006_01 Flybase release
            if file_name.endswith(".fb"):
                fb_column_names = "DB	DB_Object_ID	DB_Object_Symbol	Qualifier	GO ID	DB:Reference	Evidence	With (or) From	Aspect	DB_Object_Name	DB_Object_Synonym	DB_Object_Type	Taxon	Date	Assigned_by"
                self._process_tsv_fb(file_name, start_symbol='!', fb_column_names=fb_column_names)
            else:
                self._process_tsv_fb(file_name, start_symbol='#')
            print(table.name)
            print(str(table.header))
            print(str(table.rows))
            table._preprocess()
            print("\n\nDF:\n"+str(table.header))
            print(str(table.rows))
            #print("\n\nDF:\n"+str(table.dataframe))
#        exit(9)

        if os.path.exists(f"{dir_name}/mapping_{release_id}.txt"):
            print(f"{dir_name}/mapping_{release_id}.txt")
            self.preloaded_mapping = True
            mappings = {}
            '''
            adicionar colunas extras aqui...
            '''
            with open(f"{dir_name}/mapping_{release_id}.txt", "r") as f:
                for line in f:
                    line = line.strip("\n").split("\t")
                    fname, column, referenced_table, referenced_column = tuple(line)
                    # fname, column, referenced_table, referenced_column, data_type, key_data = tuple(line)
                    if not fname.startswith("ncRNA_genes"):
                        if not fname.endswith(".fb"):
                            fname = fname + '_' + release_id + ".tsv"
                    if fname not in mappings:
                        mappings[fname] = []
                    mappings[fname].append(tuple([column, referenced_table, referenced_column]))
                    #mappings[fname].append(tuple([column, referenced_table, referenced_column, data_type, key_data]))

            finished = []
            for key, table in self.unmapped_tables.items():
                if key not in mappings:     # to not process unmapped files right now
                    print(f"Table {key} not in mappings...")
                    continue
                print(f"Table --> {key} <-- in mappings...")
                print(str(mappings[key]))
                for column, referenced_table, referenced_colum in mappings[key]:    # separates mapped from non-mapped columns (fields)
                #for column, referenced_table, referenced_colum, data_type, key_data in mappings[key]:
                    table.unmapped_fields.remove(column)
                    table.mapped_fields.add(column)
                    table.mapping[column] = tuple([referenced_table, referenced_column])    # inserts columns data into the current table's map
                    #table.mapping[column] = tuple([referenced_table, referenced_column, data_type, key_data])
                if table.all_fields_mapped():
                    finished.append(key)
                #print(str(table.mapping))
            for key in finished:    # move tables with all columns mapped to the mapped_tables structure
                print(f"Table {key} totally in mappings...")
                self.mapped_tables[key] = self.unmapped_tables.pop(key)


    # saulo

    # added by saulo 2023/08/10 to process "fb" file that has no column names
    def _process_tsv_fb(self, file_name, start_symbol='#', fb_column_names=None):
        header = None
        with open(file_name) as f:
            rows = csv.reader(f, delimiter="\t", quotechar='"')
            for row in rows:
                # strip() added by saulo (2023/08/09) because of some "blank" lines in tsvs
                if not row or not row[0].strip():
                    continue
                # print(row[0])
                if not row[0].startswith(start_symbol):
                    if header is None:
                        if start_symbol != None and start_symbol == "!":
                            header = fb_column_names.split('\t')
                            #print(header)
                        else:
                            header = [previous[0].lstrip(start_symbol), *previous[1:]]
                        print("File " + file_name + " header: " + str(header))
                        self._set_header(file_name, header)
                    if start_symbol != None and start_symbol == "!":    # file gene_association.fb (2023_03) has two blank columns more than header :-(
                        row = row[:len(header)]
                        #self._add_row(file_name, row)
                    self._add_row(file_name, row)
                    #print(row)
                if not row[0].startswith(start_symbol + "-----"):
                    previous = row

    def _process_tsv(self, file_name):
        header = None
        with open(file_name) as f:
            rows = csv.reader(f, delimiter="\t", quotechar='"')
            for row in rows:
                # strip() added by saulo (2023/08/09) because of some "blank" lines in tsvs
                if not row or not row[0].strip():
                    continue
                #print(row[0])
                if not row[0].startswith("#"):
                    if header is None:
                        header = [previous[0].lstrip("#"), *previous[1:]]
                        self._set_header(file_name, header)
                        #print(header)
                    self._add_row(file_name, row)
                    #print(row)
                if not row[0].startswith("#-----"):
                    previous = row
        #self.unmapped_tables[file_name].print_values()


    def mappings_str(self):
        output = []
        output.append(f"Fully mapped tables: {len(self.mapped_tables)}\n")
        for key, table in self.mapped_tables.items():
            output.append(f"{key}")
            for key in table.mapped_fields:
                sql_table, sql_field = table.mapping[key]
                output.append(f"\t{key} -> {sql_table} {sql_field}")
        if len(self.unmapped_tables) == 0:
            return "\n".join(output)
        output.append(f"Non (or partially) mapped tables: {len(self.unmapped_tables)}\n")
        for key, table in self.unmapped_tables.items():
            output.append(f"{key}")
            for key in table.mapped_fields:
                sql_table, sql_field = table.mapping[key]
                output.append(f"\t{key} -> {sql_table} {sql_field}")
            for key in table.unmapped_fields:
                output.append(f"\t{key} -> ???")
        return "\n".join(output) + "\n"



    def _add_row(self, file_name, row):
        self.unmapped_tables[file_name].add_row(row)



    def _set_header(self, file_name, header):
        self.unmapped_tables[file_name].set_header(header)


    def _process_ncrna(self, json_dict):
        known_keys = [
            "primaryId",
            "symbol",
            "sequence",
            "taxonId",
            "soTermId",
            "gene",
            "symbolSynonyms",       # symbol synonyms are a list: 1, 2, 3,...
            "publications",
            "genomeLocations",
            "url",
            "crossReferenceIds",
            "relatedSequences",
        ]
        main_table_header = [
            "primaryId",
            "symbol",
            "sequence",
            "taxonId",
            "soTermId",
            "gene_geneId",
            "gene_symbol",
            "gene_locusTag"
        ]
        main_table_rows = []
        synonyms_table_header = ["symbol1", "symbol2"]
        synomyms_table_rows = []
        cross_reference_table_header = ["symbol1", "symbol2"]
        cross_reference_table_rows = []
        related_sequences_table_header = ["primaryId", "sequenceId", "relationship"]
        related_sequences_table_rows = []
        gene_synonyms_table_header = ["symbol1", "symbol2"]
        gene_synomyms_table_rows = []
        publications_table_header = ["primaryId", "publication"]
        publications_table_rows = []
        genome_locations_table_header = [
            "primaryId", 
            "assembly", 
            "gca_accession", 
            "INSDC_accession", 
            "chromosome", 
            "strand", 
            "startPosition", 
            "endPosition"
        ]
        genome_locations_table_rows = []
        for row in json_dict["data"]:
            for key in row:
                assert key in known_keys, f"Invalid key: {key}"
                    
            #fbid = row["primaryId"].split(":")[1]
            fbid = row["primaryId"]
            symbol = row["symbol"]
            sequence = row["sequence"]
            taxonid = row["taxonId"]
            sotermid = row["soTermId"]
            gene_geneid = row["gene"]["geneId"]
            gene_symbol = row["gene"]["symbol"]
            gene_name = row["gene"]["name"]
            gene_locustag = row["gene"]["locusTag"]
            main_table_rows.append([
                fbid, symbol, sequence, taxonid, sotermid, 
                gene_geneid, gene_symbol, gene_locustag])
            if "symbolSynonyms" in row:
                for synonym in row["symbolSynonyms"]:
                    synomyms_table_rows.append([symbol, synonym])
                    synomyms_table_rows.append([synonym, symbol])
            if "crossReferenceIds" in row:
                for cross_reference in row["crossReferenceIds"]:
                    cross_reference_table_rows.append([symbol, cross_reference])
                    cross_reference_table_rows.append([cross_reference, symbol])
            if "relatedSequences" in row:
                for related_sequence in row["relatedSequences"]:
                    related_sequences_table_rows.append([
                        fbid, 
                        related_sequence["sequenceId"], 
                        related_sequence["relationship"]])
            if "synonyms" in row["gene"]:
                for synonym in row["gene"]["synonyms"]:
                    gene_synomyms_table_rows.append([gene_symbol, synonym])
                    gene_synomyms_table_rows.append([synonym, gene_symbol])
            if "publications" in row:
                for publication in row["publications"]:
                    publications_table_rows.append([fbid, publication])
            for genome_location in row["genomeLocations"]:
                for exon in genome_location["exons"]:
                    genome_locations_table_rows.append([
                        fbid,
                        genome_location["assembly"],
                        genome_location["gca_accession"],
                        exon["INSDC_accession"],
                        exon["chromosome"],
                        exon["strand"],
                        str(exon["startPosition"]),
                        str(exon["endPosition"])])
        table_list = [
            ("ncRNA_genes", main_table_header, main_table_rows),
            ("ncRNA_genes_synonyms", synonyms_table_header, synomyms_table_rows),
            ("ncRNA_genes_cross_references", cross_reference_table_header, cross_reference_table_rows),
            ("ncRNA_genes_related_sequences", related_sequences_table_header, related_sequences_table_rows),
            ("ncRNA_genes_gene_synonyms", gene_synonyms_table_header, gene_synomyms_table_rows),
            ("ncRNA_genes_publications", publications_table_header, publications_table_rows),
            ("ncRNA_genes_genome_locations", genome_locations_table_header, genome_locations_table_rows)
        ]
        for table_name, header, rows in table_list:
            table = table_types.Table(table_name)
            table.set_header(header)
            for row in rows:
                table.add_row(row)
            self.unmapped_tables[table_name] = table
            self.all_tables.append(table)

    def set_sql_primary_key(self, sql_table, field):
        self.sql_primary_key[sql_table] = field

    def all_tables_mapped(self):
        return self.preloaded_mapping or len(self.unmapped_tables) == 0

    def check_field_value(self, sql_table, sql_field, value):
        finished = []
        for key, table in self.unmapped_tables.items():
            table.check_field_value(sql_table, sql_field, value)
            if table.all_fields_mapped():
                finished.append(key)
        for key in finished:
            self.mapped_tables[key] = self.unmapped_tables.pop(key)

    def get_relevant_sql_tables(self):
        answer = set()
        for table in self.all_tables:
            answer = answer.union(table.get_relevant_sql_tables())
        return answer

    def print_matched_tables(self):
        for table in self.mapped_tables.values():
            table.print_near_match()
        for table in self.unmapped_tables.values():
            table.print_near_match()

    def check_nearly_matched_tables(self):
        finished = []
        for key, table in self.unmapped_tables.items():
            table.check_near_match()
            if table.all_fields_mapped():
                finished.append(key)
        for key in finished:
            self.mapped_tables[key] = self.unmapped_tables.pop(key)

    def get_table(self, table_name):
        if table_name in self.mapped_tables:
            return self.mapped_tables[table_name]
        elif table_name in self.unmapped_tables:
            return self.unmapped_tables[table_name]
        else:
            return None

