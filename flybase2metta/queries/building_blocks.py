from typing import Any, Dict, List, Optional, Set, Tuple, Union
from hyperon_das import DistributedAtomSpace
from hyperon_das.utils import QueryOutputFormat
from hyperon_das.factory import DatabaseFactory, database_factory
from hyperon_das.logger import logger
from hyperon_das.pattern_matcher import (
    LogicalExpression,
    PatternMatchingAnswer,
    Variable, Node, Link, Variable, Not, And, Or
)


host = '104.238.183.115'
port = '8081'
query_params = {
    "toplevel_only": False,
    "return_type": QueryOutputFormat.ATOM_INFO,
}


class FBQueryBuildingBlocks():

    def __init__(self, das):
        self.server = das.remote_das[0]

    def getFBid(self, entity_symbol: str, entity_type = "FBgn") -> str:
        """
        Returns the current non-obsolete FBid for the entity referred to by its symbol (name) (Myc, Clk, GAPDH,...)
        if there are many "uniquenames" (or obsoletes or for internal FB stuff) search for the one that is both:
            1) not obsolete
            2) an FB id
        """
        schema = "feature_name"
        if entity_type == "FBgg":
            schema = "grp_name"
        answer = self._primary_key(entity_symbol.replace('(', '[').replace(')', ']'), schema)

        # if there are many "uniquenames" (or obsoletes or for internal FB stuff) search for the one that is both:
        # 1) not obsolete
        # 2) an FB id
        for link in answer:
            pk = link['targets'][1]['name']
            type_pk = link['targets'][1]['type']
            fb_id = self._getFBid(pk, type_pk)
            if not fb_id.startswith("FB"):
                #print(fb_id)
                continue
            if not self._is_obsolete(pk, type_pk):
                return fb_id

    # In fact, a group could have more than one parent! E.g: FBgg0000275
    def get_parent_groups_ids(self, group_id):
        answer = self.server.query({
            "atom_type": "link",
            "type": "Inheritance",
            "targets": [
                {"atom_type": "node", "type": "Verbatim", "name": group_id},
                {"atom_type": "variable", "name": "v1"}
            ]
        })

        print(answer)#        return get_mappings(parent_query, "v1")
        return [link['targets'][1]['name'] for link in answer]

#######################################################################################################################
    def _primary_key(self, entity_symbol: str, schema: str) -> List[Dict[str, Any]]:
        answer = self.server.query({
            "atom_type": "link",
            "type": "Execution",
            "targets": [
                {"atom_type": "node", "type": "Schema", "name": f"Schema:{schema}"},
                {"atom_type": "variable", "name": "v1"},
                {"atom_type": "node", "type": "Verbatim", "name": entity_symbol}
            ]
        })
        # print("Primary_key answer:")
        # print(answer)
        # print("\nPrimary_keys:")
        # print([link['targets'][1]['name'] for link in answer])
        return answer


    def _getFBid(self, pk: str, type_pk: str) -> str:
        answer = self.server.query({
            "atom_type": "link",
            "type": "Execution",
            "targets": [
                {"atom_type": "node", "type": "Schema", "name": f"Schema:{type_pk}_uniquename"},
                {"atom_type": "node", "type": type_pk, "name": pk},
                {"atom_type": "variable", "name": "v1"},
            ]
        })
        # print(f'___getFBid::\n{answer}\n{answer[0]["targets"][2]["name"]}')
        return answer[0]['targets'][2]['name']


    def _is_obsolete(self, pk, type_pk):
        answer = self.server.query({
            "atom_type": "link",
            "type": "Execution",
            "targets": [
                {"atom_type": "node", "type": "Schema", "name": f"Schema:{type_pk}_is_obsolete"},
                {"atom_type": "node", "type": type_pk, "name": pk},
                {"atom_type": "variable", "name": "v1"},
            ]
        })
        # print(f'is_obsolete::\n{answer}\n{answer[0]["targets"][2]["name"]}')
        for link in answer:
            if link['targets'][2]['name'] == "Concept:False":
                return False
        return True


    def _name(self, link, index, typed = False):
        named_type = f"{link['targets'][index]['type']}:" if typed else ''
        return f"{named_type}{link['targets'][index]['name']}"

    def _print_query_answer(self, query_answer, typed = False):
        if query_answer:
            for link in query_answer:
                if len(link['targets']) == 2:
                    print(f"{link['type']}: {self._name(link, 0)} -> {self._name(link, 1)}")
                elif len(link['targets']) == 3:
                    print(f"{link['type']}: {self._name(link, 0)}({self._name(link, 1, typed)}) -> {self._name(link, 2, typed)}")
                else:
                    assert False
            print(f'Number of nodes linked: {len(query_answer)}')



das = DistributedAtomSpace()
das.attach_remote(host=host, port=port)
print(f"Connected to DAS at {host}:{port}")
#print("(nodes, links) =", server.count_atoms())

primers = FBQueryBuildingBlocks(das)
##print(primers.getFBid("AGO2"))
#print(primers.getFBid("Diedel"))
'''
# the FBgn# for the 20 NetAct TFs.
# In this small "GeneGroup DAS" "Nipped-B" and "Rbf" are not in the "feature_uniquename" schema
for name in ["Su(var)205", "Top3beta", "Mef2", "Clk", "Dref", "TfIIB", "Myc", "AGO2", "Nipped-B",
             "Cp190", "TfIIA-L", "Trl", "ash1", "Raf", "Abd-B", "Orc2", "Rbf", "mof", "msl-1", "Hmr"]:
    fb_id = primers.getFBid(name)
    print(str(fb_id) + "  ::---->    " + name)

print("\n\nFBgn# for gene Argonaut:")
fb_ids = primers.getFBid("AGO2")
print (fb_ids)
'''
print("\n\nFBgn# for gene cap-n-collar-RG:")
fb_ids = primers.getFBid("cnc-RG")
print (fb_ids)
print("\nFBgn# for gene Diedel:")
fb_ids = primers.getFBid("Diedel")
print (fb_ids)

fb_ids = primers.get_parent_groups_ids("FBgg0000275")  #FBgg0001581
print (fb_ids)
fb_ids = primers.getFBid("HP1", "FBgg")
print (fb_ids)