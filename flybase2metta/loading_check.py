from hyperon_das_atomdb import UNORDERED_LINK_TYPES, WILDCARD

from hyperon_das.api import DistributedAtomSpace
from hyperon_das.exceptions import QueryParametersException
from hyperon_das.pattern_matcher.pattern_matcher import And, Link, Variable, Node
#import hyperon_das
import csv

class DASLoadingVerifier:

    '''
        Generates a list like this:
        table_strings = [
            "organism_abbreviation",
            "synonym_name",
            "gene_group_data_FB_group_id",
            "gene_group_data_FB_group_symbol" ,
            "gene_group_data_FB_group_name"
        ]

        This is the name of "Schema nodes" in DAS lingo. A "Schema node" represents a column in
        a Flybase precomputed table and it links all data in that column. The name of a Schema
        node is composed of  table name underscore column name.

        Parameter:  file_absolute_name
                    Look the "precomputed" directory for the "essential_pairs.txt" file. It should
                    be passed to this function.
    '''
    def read_precomputed_table_columns(self, file_absolute_name):
        #with open(PRECOMPUTED_DIR + "/essential_pairs.txt", 'r') as file:
        with open(file_absolute_name, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            schema_nodes_list = []

            for row in reader:
                # print(row)
                table_name = row[0].strip()
                table_name = table_name.replace(' ', '_')
                column_names = row[1:]      #  repetitions can be stored in this list, so
                column_names = list( set( column_names ) )  # remove duplicates

                for i in range(0, len(column_names)):
                    if column_names[i] != '':
                        schema_nodes_list.append( table_name + "_" + column_names[i].strip().replace(' ', '_') )

        return schema_nodes_list


    def print_schema_nodes_data(self, file_absolute_name, print_full_data=False):
        schema_names = self.read_precomputed_table_columns(file_absolute_name)

        for column in schema_names:
            v1 = Variable("v1")
            v2 = Variable("v2")
            s =  Node("Schema", "Schema:" + column)
            # s1 = Node("Schema", "Schema:gene_residues")
            q1 = Link("Execution", ordered=True, targets=[s, v1, v2])  # linka schema=s à v1 e v2. v1 = pk e v2 = column
            # q2 = Link("Execution", ordered=True, targets=[s1, v1, v2])
            assignments = query(q1, True)
            # ass2 = query(q2, True)

            cont = 0
            print(f"Number of relationships in column {column}: {len(assignments)}")
            if print_full_data:
                for assignment in assignments:
                    # pkey_handle = assignment.mapping["v1"]
                    # pkey = db.get_node_name(pkey_handle)
                    un_handle = assignment.mapping["v2"]
                    unique_name = db.get_node_name(un_handle)
                    print(str(cont) + ": " + unique_name)
                    cont += 1
                    if cont > 100:
                        break
                    # if unique_name == verificar_prefixo(unique_name):
                    #    print("v2: " + unique_name)
            print("FINISHED for column: " + column + "\n")


    '''
        Returns a list of pairs 
        
        [(FBgg, [(FBgn id, gene_symbol), (FBgn id, gene_symbol),...]),
         (FBgg, [(FBgn id, gene_symbol), (FBgn id, gene_symbol)])     ]
          
          of all genes belonging to the same groups that the gene passed as prmeter belongs to.
    '''
    def get_gene_groups(self, gene_symbol):
        pass

    '''
        Returns the hierarchy of the groups that the gene belongs to
        
        @parameter:  the gene_symbol
        @return:    a list: if the gene belongs only to one group the list has one element
                    that is itself a list of group FBgg ids: the first element is the group that the gene 
                    belongs to; the second is its parent group; the third is its grandparent;...
                    
                    if the gene belongs to more than one group the list contains as many elements as the number of groups. 
    '''
    def get_group_hierarchy(self, gene_symbol):
        pass

    """
    This is an expanded version of the above tailored to performing queries using only the  "essential pairing" of columns.

    In fact, the three functions use the same code structure. So, the three could be merged into only one adding a "schema parameter"
    (and changing identifiers, of course
    """

    # get FB id (FBgn#) and group id(s) (FBgg#) of gene designated by gene_symbol
    def get_gene_group_ids(gene_symbol):
        fb_id = get_feature_fb_id(gene_symbol)
        n1 = Node("Verbatim", fb_id)
        v1 = Variable("v1")
        # groups
        s = Node("Schema", "Schema:gene_group_data_FB_group_id")

        q1 = Link("Execution", ordered=True, targets=[s, n1, v1])
        return fb_id, get_mappings(q1, "v1")

    # get FB id (FBgn#) and PATHWAY group id(s) (FBgg#) of gene designated by gene_symbol
    def get_pathway_gene_group_ids(gene_symbol):
        fb_id = get_feature_fb_id(gene_symbol)
        n1 = Node("Verbatim", fb_id)
        v1 = Variable("v1")
        # groups
        s = Node("Schema", "Schema:pathway_group_data_FB_group_id")

        q1 = Link("Execution", ordered=True, targets=[s, n1, v1])
        return fb_id, get_mappings(q1, "v1")

    def get_gene_group_symbols(gene_symbol):
        gene_fb_id, gene_groups_ids = get_gene_group_ids(gene_symbol)

        gg_symbols = []
        for gg_id in gene_groups_ids:
            n1 = Node("Verbatim", gg_id)
            v1 = Variable("v1")
            # groups
            s = Node("Schema", "Schema:gene_group_data_FB_group_symbol")

            q1 = Link("Execution", ordered=True, targets=[s, n1, v1])
            for symb in get_mappings(q1, "v1"):
                gg_symbols.append(symb)
        return gg_symbols

    # get group name(s) of gene designated by gene_symbol
    def get_gene_group_names(gene_symbol):
        gene_fb_id, gene_groups_ids = get_gene_group_ids(gene_symbol)

        gg_names = []
        for gg_id in gene_groups_ids:
            n1 = Node("Verbatim", gg_id)
            v1 = Variable("v1")
            # groups
            s = Node("Schema", "Schema:gene_group_data_FB_group_name")

            q1 = Link("Execution", ordered=True, targets=[s, n1, v1])
            for gg_name in get_mappings(q1, "v1"):
                gg_names.append(gg_name)
        return gg_names

    # Returns list(s) of gene symbols that are members of the same group(s) of
    # the gene desinated by gene_symbol
    def get_gene_group_members(gene_symbol):
        gene_fb_id, gene_groups_ids = get_gene_group_ids(gene_symbol)

        gene_group_symbols = []
        for gg_id in gene_groups_ids:
            a_gg_symbols = []
            n1 = Node("Verbatim", gg_id)
            v1 = Variable("v1")
            # gets all gene ids (FBgn#) of group given by gg_id
            s = Node("Schema", "Schema:gene_group_data_Group_member_FB_gene_id")
            q1 = Link("Execution", ordered=True, targets=[s, n1, v1])
            gene_ids = get_mappings(q1, "v1")
            for gene_id in gene_ids:  # gets the symbol for each gene
                n2 = Node("Verbatim", gene_id)
                v2 = Variable("v2")
                sf_name = Node("Schema", "Schema:gene_group_data_Group_member_FB_gene_symbol")
                q1 = Link("Execution", ordered=True, targets=[sf_name, n2, v2])
                a_gg_symbols.append(get_mappings(q1, "v2")[0].replace("[", "(").replace("]", ")"))
            gene_group_symbols.append((a_gg_symbols, len(a_gg_symbols)))

        return gene_fb_id, gene_group_symbols

    # same as get_gene_group_members() for pathway groups
    def get_pathway_gene_group_members(gene_symbol):
        gene_fb_id, gene_groups_ids = get_pathway_gene_group_ids(gene_symbol)

        gene_group_symbols = []
        for gg_id in gene_groups_ids:
            a_gg_symbols = []
            n1 = Node("Verbatim", gg_id)
            v1 = Variable("v1")
            # gets all gene ids (FBgn#) of group given by gg_id
            sp = Node("Schema", "Schema:pathway_group_data_Group_member_FB_gene_id")
            q1 = Link("Execution", ordered=True, targets=[sp, n1, v1])
            gene_ids = get_mappings(q1, "v1")
            for gene_id in gene_ids:
                n2 = Node("Verbatim", gene_id)
                v2 = Variable("v2")
                ss = Node("Schema", "Schema:pathway_group_data_Group_member_FB_gene_symbol")
                q1 = Link("Execution", ordered=True, targets=[ss, n2, v2])
                a_gg_symbols.append(get_mappings(q1, "v2")[0].replace("[", "(").replace("]", ")"))
            gene_group_symbols.append((a_gg_symbols, len(a_gg_symbols)))

        return gene_fb_id, gene_group_symbols

    # In fact, a group could have more than one parent! E.g: FBgg0000275
    def get_parent_group_ids(group_id):
        group_node = Node("Verbatim", group_id)
        parent_var = Variable("v1")
        parent_query = Link("Inheritance", ordered=True, targets=[group_node, parent_var])

        return get_mappings(parent_query, "v1")

    # In fact, a group could have more than one parent! E.g: FBgg0000275
    def get_parent_group_symbols(group_id):
        parent_ids = get_parent_group_ids(group_id)
        parents_symbols = []
        for parent_id in parent_ids:
            n1 = Node("Verbatim", parent_id)
            v1 = Variable("v1")
            s = Node("Schema", "Schema:gene_group_data_Parent_FB_group_symbol")
            q1 = Link("Execution", ordered=True, targets=[s, n1, v1])
            parents_symbols.extend(get_mappings(q1, "v1"))
        return parents_symbols
        '''
        n1 = Node("Verbatim", group_id)
        v1 = Variable("v1")
        s = Node("Schema", "Schema:gene_group_data_Parent_FB_group_symbol")        
        q1 = Link("Execution", ordered=True, targets=[s, n1, v1])

        return get_mappings(q1, "v1")
        '''

    # In fact, a group could have more than one parent! E.g: FBgg0000275
    # In the precomputed tables there are only group PARENT id and symbol...
    def get_parent_group_names(group_id):
        parent_ids = get_parent_group_ids(group_id)
        parents_names = []
        for parent_id in parent_ids:
            n1 = Node("Verbatim", parent_id)
            v1 = Variable("v1")
            s = Node("Schema", "Schema:gene_group_data_FB_group_name")
            q1 = Link("Execution", ordered=True, targets=[s, n1, v1])
            parents_names.extend(get_mappings(q1, "v1"))
        return parents_names


