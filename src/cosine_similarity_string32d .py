import gensim
import networkx as nx
from node2vec import Node2Vec
import pandas as pd
import numpy as np
import pickle
import collections
import statistics
from sklearn.metrics.pairwise import cosine_similarity


embed = 'string_embeddings_32'

# Load pre-trained Word2Vec model.
model = gensim.models.Word2Vec.load("/scratch/PI/rbaltman/yash_margaret/" + embed + ".model")

# Load pre-computed gene sets for drugs
with open('/scratch/PI/rbaltman/yash_margaret/cmap_data/expert_drug_gene_list.pickle', 'rb') as handle:
    expert_drug_gene_list = pickle.load(handle)

# Load pre-computed disease gene sets
with open('/scratch/PI/rbaltman/yash_margaret/new_classifier_genes/phenotype_refinement_dict.pickle', 'rb') as handle:
    diseaseSig_psygenet_GWAS = pickle.load(handle)

nodes = list(model.wv.vocab)
embeddings = np.array([model.wv[x] for x in nodes])
cos_sim = cosine_similarity(embeddings)

cos_sim_matrix_genes = pd.DataFrame(data=cos_sim[0:,0:], index=nodes, columns=nodes)
cos_sim_matrix_genes.to_pickle('/scratch/PI/rbaltman/yash_margaret/cosine_distance/' + embed + '.pickle')

for disease, diseaseSig in diseaseSig_psygenet_GWAS.items():
    cos = collections.defaultdict(list)
    for index, row in expert_drug_gene_list.iterrows():
        drug = row['drug_name']
        drugSet = row['union']
        for x in diseaseSig:
            for y in drugSet:
                if (str(y) in cos_sim_matrix_genes.columns) and (str(x) in cos_sim_matrix_genes.columns):
                    cos[drug].append(cos_sim_matrix_genes.get_value(str(x), str(y)))
    mean_cos_distance = {}
    for key, value in cos.items():
        new = list(map(int, value))  
        mean_cos_distance[key] = statistics.mean(new)
    mean_cos = pd.DataFrame(list(mean_cos_distance.items()))
    mean_cos.to_csv('/scratch/PI/rbaltman/yash_margaret/cosine_distance/pheno_enriched/' + disease + '_' + embed + '.csv')
    print("Finished " + disease)