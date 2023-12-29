# Parsing and analyzing pathways



from Bio.KEGG import REST
# human_pathways = REST.kegg_list("pathway", "hsa").read()
from xml.etree import ElementTree

from copy import deepcopy

import networkx as nx
import itertools
import json

import re

# import matplotlib.pyplot as plt



def get_pathway_graph(kegg_pathway_id, verbose=5):
    '''Get networkx graph of KEGG pathway'''
    
    assert verbose in list(range(11)), "Verbose must be in [0,10]"

    if verbose>=3:
        print("** TODO: Modify to have reactions be *nodes* **")
        print("** TODO: Grab all genes for orthologs, not just 'SEQUENCE' entry **")


    kgml_pathway = REST.kegg_get(kegg_pathway_id, "kgml").read()

    pathroot = ElementTree.fromstring(kgml_pathway)

    G = nx.MultiDiGraph()

    for child in pathroot:
        # print()
    #     print("------")
    #     print(dir(child))
    #     print(child.tag) # says if entry, relation, etc. [if relation or reaction, add edge]
    #     print(child.attrib)
        # print(child.tag,child.attrib['type'])
        
        if child.tag=='entry':
            # 'ortholog'
            # 'enzyme'
            # 'reaction'
            # 'gene'
            # 'group' **
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
            attr = deepcopy(child.attrib)
            attr['edgetype'] = 'relation'
            attr['subtypes'] = []
            for subtype in child:
                attr['subtypes'].append( (subtype.attrib['name'],subtype.attrib['value']) )
            
            G.add_edge(child.attrib['entry1'],child.attrib['entry2'],**attr)
        
        elif child.tag=='reaction':
            attr = deepcopy(child.attrib)
            attr['edgetype'] = 'reaction'
            # substrates as list of pairs. products as list of pairs. 
            substrates = []
            products = []
            for subtype in child:
                if subtype.tag=='substrate':
                    substrates.append( (subtype.attrib['id'],subtype.attrib['name']) )
                elif subtype.tag=='product':
                    products.append( (subtype.attrib['id'],subtype.attrib['name']) )
                else:
                    raise NotImplementedError("Reaction encountered, not yet implemented")
            attr['substrates'] = substrates
            attr['products'] = products
            
            for substrate,product in itertools.product(substrates,products):
                G.add_edge(substrate[0],product[0],**attr)



            # G.add_edge(child.attrib['entry1'],child.attrib['entry2'],**attr)


    # print(G)

    # Additional processing
    for node,attr in G.nodes(data=True):
        # Split up the "name" node, since sometimes it's i.e. "hsa:5818" "hsa:5819"
        if 'name' in attr:
            names_split = attr['name'].split(' ')
            attr['names'] = tuple(names_split)

            del attr['name']

    return G


def get_ortholog_genes_from_keggtext(inpstr):
    '''Get ortholog genes from KEGG text
    
    Note: Current version of code just grabs from "SEQUENCE" line,
    and for now takes only first regex match.
    '''

    match = re.search(r"SEQUENCE\s+\[(.*?)\]", inpstr)
    strgene = ""
    if match:
        strgene = match.group(1)
        if len(strgene)>50:
            strgene = "error_ortholog_texttoolong"
    else:
        strgene = "no_match_for_ortholog"

    return strgene




def convert_orthologs_to_genes(G,verbose=5):
    '''Convert orthologs to gene list in a pathway graph
    
    Passing by reference (i.e. modifies in-place).
    '''

    for node,attr in G.nodes(data=True):

        # print('###',node,attr)
        if 'type' not in attr:
            continue

        if attr['type']=='ortholog':
            
            new_names = []

            for name in attr['names']:
                if verbose>=6:
                    print("Getting ortholog genes for",name)
                ko = REST.kegg_get( name ).read()
                strgene = get_ortholog_genes_from_keggtext(ko)
                new_names.append(strgene)
            
            # new_names[node] = tuple(new_names)
            attr['names'] = tuple(new_names)

    return True


def get_pathway_genes(pathways_graph):
    '''Get all genes in a given pathway's node's "names" field
    
    Returns a list of kegg id numbers for genes.
    '''

    genes = set()

    for node,attr in pathways_graph.nodes(data=True):

        if 'names' not in attr:
            continue

        for name in attr['names']:
            genes.add(name)

    return genes








def save_pathway_graph(G,fname):
    '''Save pathway nxgraph to file'''

    # Convert the graph to JSON-compatible dictionary
    data = nx.node_link_data(G)

    # Dump the dictionary to a JSON file
    with open(fname, 'w') as f:
        json.dump(data, f)


def load_pathway_graph(fname):
    '''Load pathway nxgraph from file'''
    
    # Load the JSON data from the file
    with open(fname, 'r') as f:
        data = json.load(f)

    # Create a NetworkX graph from the JSON data
    G = nx.node_link_graph(data)

    return G
    


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







