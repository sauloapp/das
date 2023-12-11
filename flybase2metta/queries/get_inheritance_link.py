from hyperon_das import DistributedAtomSpace
from hyperon_das.utils import QueryOutputFormat

host = '104.238.183.115'
port = '8081'
das = DistributedAtomSpace()
das.attach_remote(host=host, port=port)
print(f"Connected to DAS at {host}:{port}")
#print("(nodes, links) =", server.count_atoms())

result = das.remote_das[0].get_link(
    link_type='Inheritance',
    targets=['FBgg0001581', 'FBgg0001782'],
    output_format=QueryOutputFormat.HANDLE
)
print(f'Link:\n{result}')