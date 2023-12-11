from hyperon_das import DistributedAtomSpace
from hyperon_das.utils import QueryOutputFormat
from hyperon_das.factory import DatabaseFactory, database_factory
from hyperon_das.logger import logger
from hyperon_das.pattern_matcher import (
    LogicalExpression,
    PatternMatchingAnswer,
    Variable, Node, Link, Variable, Not, And, Or
)



def _name(link, index, typed = False):
    named_type = f"{link['targets'][index]['type']}:" if typed else ''
    return f"{named_type}{link['targets'][index]['name']}"

def _print_query_answer(query_answer, typed = False):
    if query_answer:
        for link in query_answer:
            if len(link['targets']) == 2:
                print(f"{link['type']}: {_name(link, 0)} -> {_name(link, 1)}")
            elif len(link['targets']) == 3:
                print(f"{link['type']}: {_name(link, 0)}({_name(link, 1, typed)}) -> {_name(link, 2, typed)}")
            else:
                assert False
        print(f'Number of nodes linked: {len(query_answer)}')

host = '104.238.183.115'
port = '8081'
query_params = {
    "toplevel_only": False,
    "return_type": QueryOutputFormat.ATOM_INFO,
}

das = DistributedAtomSpace()
das.attach_remote(host=host, port=port)
print(f"Connected to DAS at {host}:{port}")
#print("(nodes, links) =", server.count_atoms())

def _filter(query_answer, index, value):
    filtered = []
    for link in query_answer:
        if link['targets'][index]['type'] == value:
            filtered.append(link)
    return filtered

def _fbgns(das, symbol, handles=False):
    server = das.remote_das[0]
    answer = server.query({
        "atom_type": "link",
        "type": "Execution",
        "targets": [
            #{"atom_type": "node", "type": "Schema", "name": "Schema:fb_synonym_primary_FBid"},
            {"atom_type": "node", "type": "Schema", "name": "Schema:gene_association_DB_Object_ID"},
            #{"atom_type": "node", "type": "Schema", "name": "Schema:feature_feature_id"},
            {"atom_type": "node", "type": "Verbatim", "name": symbol},
            {"atom_type": "variable", "name": "v1"},
        ]
    })
    print(answer)
    if handles:
        return [link['targets'][2]['handle'] for link in answer]
    else:
        return [link['targets'][2]['name'] for link in answer]

def query1(das, symbol, node_type = None):
    server = das.remote_das[0]
    fbgns = _fbgns(das, symbol)
    print(f"FBgn: {fbgns}")
    answer = []
    for fbgn in fbgns:
        query_answer = server.query({
            "atom_type": "link",
            "type": "Execution",
            "targets": [
                {"atom_type": "variable", "name": "v0"},
                {"atom_type": "variable", "name": "v1"},
                {"atom_type": "node", "type": "Verbatim", "name": fbgn},
            ]
        })
        if node_type:
            query_answer = _filter(query_answer, 1, node_type)
        answer.extend(query_answer)
        query_answer = server.query({
            "atom_type": "link",
            "type": "Execution",
            "targets": [
                {"atom_type": "variable", "name": "v0"},
                {"atom_type": "node", "type": "Verbatim", "name": fbgn},
                {"atom_type": "variable", "name": "v1"},
            ]
        })
        if node_type:
            query_answer = _filter(query_answer, 0, node_type)
        answer.extend(query_answer)
    return answer

_print_query_answer(query1(das, "Myc", "BiologicalProcess"))
_print_query_answer(query1(das, "Myc") )
#exit(9)
#_print_query_answer(query1(das, "Top3beta", "BiologicalProcess"))
    #_print_query_answer(query1(das, "Top3beta", "CellularComponent"))
#_print_query_answer(query1(das, "Top3beta", "MolecularFunction"))

#_print_query_answer(_get_FBids(das, "Top3beta", "MolecularFunction"))
#print(das.remote_das[0].get_nodes("Schema", output_format=QueryOutputFormat.ATOM_INFO))

result = das.remote_das[0].get_nodes(
    node_type='Schema',
    output_format=QueryOutputFormat.ATOM_INFO
)

print(result)

result = das.remote_das[0].get_nodes(
    node_type='Schema',
    node_name='Schema:gene_association_DB_Object_ID',
    output_format=QueryOutputFormat.ATOM_INFO
)

print(result)

result = das.remote_das[0].get_link(
    link_type='Inheritance',
    targets=['FBgg0001581', 'FBgg0001782'],
    output_format=QueryOutputFormat.HANDLE
)
print(f'Link:\n{result}')
exit(9)
result = das.remote_das[0].get_links(
    link_type='Inheritance',

    target_types=['Verbatim', 'Verbatim'],
    #output_format=QueryOutputFormat.HANDLE
    output_format=QueryOutputFormat.ATOM_INFO
)
#print(f'Linksssssss:\n{result}')
print(f"len get_links Inheritance type(result) = {len(result)}")

V1 = Variable("V1")
V2 = Variable("V2")
V3 = Variable("V3")

logical_expression = And([
    Link("Inheritance", ordered=True, targets=[V1, V2]),
    Link("Inheritance", ordered=True, targets=[V2, V3])
])

query_params = {
                "toplevel_only": False,
                "return_type": QueryOutputFormat.ATOM_INFO,
            }
#result = das.remote_das[0].pattern_matcher_query(query=logical_expression, {'return_type': QueryOutputFormat.HANDLE})
# deu pau
#result = das.remote_das[0].query(query=logical_expression, extra_parameters=query_params)

print("And Link")
print(result)
#print( das.remote_das[0].get_node("Verbatim", "Myc", output_format=QueryOutputFormat.ATOM_INFO) )

print( das.remote_das[0].get_node("Verbatim", "Myc", output_format=QueryOutputFormat.ATOM_INFO) )
print( das.remote_das[0].get_node("Schema", "Schema:gene_association_DB_Object_ID", output_format=QueryOutputFormat.ATOM_INFO) )
#{"atom_type": "node", "type": "", "name": "Schema:gene_association_DB_Object_ID"},

'''
symbol = "Top3beta"

server = das.remote_das[0]
#print("(nodes, links) =", server.count_atoms())
answer = server.query({
    "atom_type": "link",
    "type": "Execution",
    "targets": [
        # {"atom_type": "node", "type": "Schema", "name": "Schema:fb_synonym_primary_FBid"},
        #{"atom_type": "node", "type": "Schema", "name": "Schema:gene_association_DB_Object_ID"},
        {"atom_type": "node", "type": "Schema", "name": "Schema:feature_name"},
        {"atom_type": "variable", "name": "v1"},
        {"atom_type": "node", "type": "Verbatim", "name": symbol}
    ]
})
print(answer)
print([link['targets'][1]['type'] for link in answer])
print([link['targets'][1]['name'] for link in answer])
#print(das.remote_das[0].get_node("Verbatim", symbol, output_format=QueryOutputFormat.ATOM_INFO))
#das.get_node("Verbatim", symbol)
'''