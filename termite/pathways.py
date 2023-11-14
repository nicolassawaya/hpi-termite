# Parsing and analyzing pathways



from Bio.KEGG import REST
# human_pathways = REST.kegg_list("pathway", "hsa").read()
from xml.etree import ElementTree


import networkx as nx
import itertools

# import matplotlib.pyplot as plt



def get_pathway_graph(kegg_pathway_id):
    '''Get networkx graph of KEGG pathway'''

    kgml_pathway = REST.kegg_get("hsa05168", "kgml").read()

    pathroot = ElementTree.fromstring(kgml_pathway)

    G = nx.MultiDiGraph()

    for child in pathroot:
        # print()
    #     print("------")
    #     print(dir(child))
    #     print(child.tag) # says if entry, relation, etc. [if relation or reaction, add edge]
    #     print(child.attrib)
        print(child.tag,child.attrib['type'])
        
        if child.tag=='entry':
            # 'ortholog'
            # 'enzyme'
            # 'reaction'
            # 'gene'
            # 'group'
            # 'compound'
            # 'map'
            # 'brite'
            # 'other'
            G.add_node(child.attrib['id'],**child.attrib)

            if child.attrib['type']=='group':
                components = []
                for subgroup in child:
                    if subgroup.tag=="graphics":
                        continue
                    # print('**',subgroup.tag,subgroup.attrib['id'])
                    if subgroup.tag=='component':
                        components.append(subgroup.attrib['id'])
                G.nodes[child.attrib['id']]['components'] = tuple(components)
                


        elif child.tag=='relation':
            rel_attr = {}
            for subtype in child:
                rel_attr[subtype.attrib['name']] = subtype.attrib['value']
            
            G.add_edge(child.attrib['entry1'],child.attrib['entry2'],**rel_attr)
            
        for subchild in child:
            if subchild.tag=='graphics':
                continue
            if child.tag=='relation':  # This is if child=='relation'
                print(subchild.tag,subchild.attrib['name'])
            else:
                print(subchild.tag, subchild.tag)


    # # Additional processing
    # for node in G.nodes():
    #     # Split up the "name" node, since sometimes it's i.e. "hsa:5818" "hsa:5819"
    #     names_split = node['name'].split(' ')
    #     node['name'] = tuple(names_split)

    return G



def save_pathway_nxgraph():
    '''Save pathway nxgraph to file'''
    pass


def load_pathway_nxgraph():
    '''Load pathway nxgraph from file'''
    pass

"""
    # https://www.kegg.jp/kegg/xml/docs/

    # 'entry'-->'type' has to be in 
    'ortholog'
    'enzyme'
    'reaction'
    'gene'
    'group'
    'compound'
    'map'
    'brite'
    'other'

    # 'relation'-->'type' has to be in
    'ECrel'
    'PPrel'
    'GErel'
    'PCrel'
    'maplink'

    # 'reaction'-->'type' has to be in
    'reversible'
    'irreversible'


    # 'relation' entries have 'subtype'-->'name' equal one of many 
    # things, incl e.g. 'activation,' 'phosphorylation,' 'glycosylation,' etc.
"""







