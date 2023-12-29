# Run termite for all kegg's hsa

import time
import math
from termite import pathways
from Bio.KEGG import REST

start_time = time.time()

# Read in HSA pathway list from KEGG (is missing some)
human_pathways = REST.kegg_list("pathway", "hsa").read()


hsa_list = []
for line in human_pathways.splitlines():
    # print(line)
    if line[:3]=='hsa':
        hsa_list.append(line[:8])
'''
- ko:** xx:** dict
- enumerate pathways
- timing
- *save* the pathways
'''

# dict for converting ko:*** to xx:***
# ko_to_xx = {}

print(hsa_list)
print("Number of pathways:",len(hsa_list))

print("Removing hsa01100 (Which is a 'super pathway')")
hsa_list.remove('hsa01100')
print("Removing hsa01200")
hsa_list.remove('hsa01200')


# Advice from Aldo:
aldo_koen_set = {'hsa00982','hsa03022','hsa03008','hsa03060',
'hsa04141','hsa04120','hsa03050','hsa03018',
'hsa03030','hsa03410','hsa03420','hsa03430',
'hsa03440','hsa03450','hsa03082','hsa04010',
'hsa04014','hsa04310','hsa04330','hsa04340',
'hsa04350','hsa04630','hsa04064','hsa04668',
'hsa04020','hsa04151','hsa04150','hsa04144',
'hsa04145','hsa04142','hsa04110','hsa04210',
'hsa04215','hsa04217','hsa04115','hsa04610',
'hsa04611','hsa04620','hsa04622','hsa04623',
'hsa04612','hsa04614','hsa04721','hsa05200',
'hsa05202','hsa05206','hsa05203','hsa05235',
'hsa05224','hsa05332','hsa05010','hsa05012',
'hsa05020','hsa05022','hsa01521'}

print("Taking union with Aldo/Koen set")
hsa_list = set(hsa_list).intersection( aldo_koen_set )


all_kegg_genes = set()

ctr = 0
max_ctr = math.inf

print("len(hsa_list):",len(hsa_list))

for i,idhsa in enumerate(hsa_list):
    print(f"=========== {i} ===========")
    print(idhsa)

    # Create pathway graph
    G = pathways.get_pathway_graph(idhsa,verbose=0)
    print(G)

    # Convert orthologs to genes
    pathways.convert_orthologs_to_genes(G,verbose=9)

    # Output genes
    genes = pathways.get_pathway_genes(G)
    all_kegg_genes = all_kegg_genes.union( genes )

    # Print genes and ko-to-xx up to this point
    print("** All genes so far **")
    print(all_kegg_genes)
    # print("** ko-to-xx **")
    print()

    # Timing
    elapsed_time = time.time() - start_time
    print("Elapsed time:", elapsed_time, "seconds")

    # Save nx graph
    fname = f"{idhsa}.json"
    pathways.save_pathway_graph(G,fname)

    ctr += 1
    if ctr>=max_ctr:
        print(f"NOTE: Breaking early after max_ctr={max_ctr} pathways.")
        break


print(all_kegg_genes)

# Save all kegg genes to file
with open('all_kegg_genes.txt','w') as f:
    for gene in all_kegg_genes:
        f.write(gene+'\n')





