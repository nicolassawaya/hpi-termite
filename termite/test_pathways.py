# Tests for pathways.py

import networkx as nx
from Bio.KEGG import REST

from . import pathways

import os


def test_hsv1():
    '''Herpes Simplex Virus 1 (HSV-1) pathway test'''

    # hsa05168: Herpes simplex virus 1 infection
    str_hsv1_kegg = 'hsa05168'

    G = pathways.get_pathway_graph(str_hsv1_kegg)

    # print(G)
    # print("*** nodes ***")
    # for node,attr in G.nodes(data=True):
    #     print(node)
    #     print(attr)
    # print("*** edges ***")
    # for u,v,attr in G.edges(data=True):
    #     print(u,v)
    #     print(attr)


def test_fatty_acid_metabolism():
    '''Fatty acid metabolism - Homo sapiens (human)'''

    # hsa:01212 Fatty acid metabolism - Homo sapiens (human)
    str_fatty = 'hsa01212'

    G = pathways.get_pathway_graph(str_fatty)

    # print(G)
    # print("*** nodes ***")
    # for node,attr in G.nodes(data=True):
    #     print(node)
    #     print(attr)
    # print("*** edges ***")
    # for u,v,attr in G.edges(data=True):
    #     print(u,v)
    #     print(attr)


def test_ortholog_regex():
    '''Test regular expression for grepping ortholog's gene(s)'''

    # Quick test of the REST (result not used)
    ko = REST.kegg_get("ko:K19264").read()

    # Snippet from "ko:K19264" https://www.kegg.jp/entry/K19264
    raw_ortholog_text = '''
             Cytokine-cytokine receptor interaction
              K19264
GENES       VG: 1487358(US6) 1487459(US6) 18533958(US6) 19621712(US6) 2703444(US6) 3190329(US6) 3850196(US6)
REFERENCE   PMID:27723487
  AUTHORS   Bhargava AK, Rothlauf PW, Krummenacher C
  TITLE     Herpes simplex virus glycoprotein D relocates nectin-1 from intercellular contacts.
  JOURNAL   Virology 499:267-277 (2016)
            DOI:10.1016/j.virol.2016.09.019
  SEQUENCE  [vg:2703444]
REFERENCE   PMID:17295428
  AUTHORS   Reske A, Pollara G, Krummenacher C, Chain BM, Katz DR
  '''

    gene = pathways.get_ortholog_genes_from_keggtext(raw_ortholog_text)
    assert gene=="vg:2703444"



def test_ortholog_to_genes():
    '''Test turning e.g. ko:***** into hsa:*****'''


    ko = REST.kegg_get("ko:K19264").read()
    # print(ko)


    '''
    EXAMPLE FROM path:hsa05168 ("Herpes simplex virus 1 infection")

    <entry id="20" name="ko:K19264" type="ortholog"
        link="https://www.kegg.jp/dbget-bin/www_bget?K19264">
        <graphics name="K19264" fgcolor="#000000" bgcolor="#FFFFFF"
             type="rectangle" x="183" y="162" width="46" height="17"/>
    </entry>
    '''

    G = nx.Graph()
    G.add_node(1)
    G.add_node(2)
    node_dict = {'names':("ko:K19264",), 'type':"ortholog"}
    G.add_node(20, **node_dict)

    pathways.convert_orthologs_to_genes(G)
    assert G.nodes[20]['names']==('vg:2703444',)



    '''
    EXAMPLE FROM path:hsa01212 ("Fatty acid metabolism")

    <entry id="169" name="ko:K10256 ko:K22993" type="ortholog" reaction="rn:R11043"
        link="https://www.kegg.jp/dbget-bin/www_bget?K10256+K22993">
        <graphics name="K10256..." fgcolor="#000000" bgcolor="#FFFFFF"
             type="line" coords="571,1580,880,1580"/>
    </entry>
    '''


    G = nx.Graph()
    G.add_node(1)
    G.add_node(2)
    node_dict = { 'names':("ko:K10256","ko:K22993",), 'type':"ortholog"}
    G.add_node(169, **node_dict)

    pathways.convert_orthologs_to_genes(G)
    assert G.nodes[169]['names']==('ath:AT3G12120', 'fvr:FVEG_12065')


    
def test_get_pathway_genes():
    '''Test getting genes from pathway graph'''

    G = nx.Graph()
    G.add_node(1)
    G.add_node(2)
    node_dict = {'names':("ko:K19264",), 'type':"ortholog"}
    G.add_node(20, **node_dict)
    node_dict = { 'names':('ath:AT3G12120', 'fvr:FVEG_12065'), 'type':"ortholog"}
    G.add_node(169, **node_dict)

    genes = pathways.get_pathway_genes(G)
    assert genes=={'ath:AT3G12120', 'fvr:FVEG_12065', 'ko:K19264'}
    



def test_saving_loading():
    '''Test saving pathway graph'''

    G = nx.Graph()
    G.add_node(1)
    G.add_node(2)
    node_dict = {'names':("ko:K19264",), 'type':"ortholog"}
    G.add_node(20, **node_dict)
    node_dict = { 'names':('ath:AT3G12120', 'fvr:FVEG_12065'), 'type':"ortholog"}
    G.add_node(169, **node_dict)

    fname = "test.json"
    
    pathways.save_pathway_graph(G, fname)


    Gloaded = pathways.load_pathway_graph(fname)
    assert Gloaded.nodes[169]['names']==['ath:AT3G12120', 'fvr:FVEG_12065']
    # print("Gloaded: ",Gloaded)
    # for node in Gloaded.nodes():
    #     print(node)
    #     print(Gloaded.nodes[node])


    # # Remove test file
    # os.remove(fname)


    





