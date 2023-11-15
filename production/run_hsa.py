# Run termite for all kegg's hsa

from termite import pathways
from Bio.KEGG import REST

# Read in HSA pathway list from KEGG (is missing some)
human_pathways = REST.kegg_list("pathway", "hsa").read()


hsa_list = []
for line in human_pathways.splitlines():
    # print(line)
    if line[:3]=='hsa':
        hsa_list.append(line[:8])


print(hsa_list)
print("Number of pathways:",len(hsa_list))

for idhsa in hsa_list:
    print(idhsa)
    G = pathways.get_pathway_graph(idhsa,verbose=0)
    print(G)










