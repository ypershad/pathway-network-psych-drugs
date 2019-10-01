import gensim
import networkx as nx
from node2vec import Node2Vec
import pandas as pd
import numpy as np
import pickle
import collections
import statistics
from sklearn.metrics.pairwise import cosine_similarity

embed = 'string_embeddings_64'

# Load pre-computed gene sets for drugs
with open('/scratch/PI/rbaltman/yash_margaret/randomDrugDictionary_string.pickle', 'rb') as handle:
    randomDrugDictionary = pickle.load(handle)

# Load pre-computed disease gene sets
with open('diseaseSig_psygenet_GWAS.pickle', 'rb') as handle:
    diseaseSig_psygenet_GWAS = pickle.load(handle)

with open('/scratch/PI/rbaltman/yash_margaret/cosine_distance/' + embed + '.pickle', 'rb') as handle:
    cos_sim_matrix_genes = pickle.load(handle)
    
print("Loaded")

for disease, diseaseSig in diseaseSig_psygenet_GWAS.items():
    cos = collections.defaultdict(list)
    for randindex, genes in randomDrugDictionary.items():
        drug = randindex
        drugSet = genes
        for x in diseaseSig:
            for y in drugSet:
                if (str(y) in cos_sim_matrix_genes.columns) and (str(x) in cos_sim_matrix_genes.columns):
                    cos[drug].append(cos_sim_matrix_genes.get_value(str(x), str(y)))
    mean_cos_distance = {}
    for key, value in cos.items():
        new = list(map(int, value))  
        mean_cos_distance[key] = statistics.mean(new)
    mean_cos = pd.DataFrame(list(mean_cos_distance.items()))
    mean_cos.to_csv('/scratch/PI/rbaltman/yash_margaret/cosine_distance/random/' + disease + '_' + embed + '.csv')
    print("Finished " + disease)