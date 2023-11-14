# Tests for pathways.py



from . import pathways




def test_hsv1():
    '''Herpes Simplex Virus 1 (HSV-1) pathway test'''

    # hsa05168: Herpes simplex virus 1 infection
    str_hsv1_kegg = 'hsa05168'

    G = pathways.get_pathway_graph(str_hsv1_kegg)

    print(G)

    for node,attr in G.nodes(data=True):
        print(node)
        print(attr)















