import pandas as pd
import pickle
import networkx as nx
import itertools
import collections

expert_drug_gene_list = pd.read_pickle('cmap_data/expert_drug_gene_list.pickle')

edgelist_file = '/scratch/groups/rbaltman/yash_margaret/db_string/string.txt'
string_edgelist = pd.read_csv(edgelist_file, sep = '\t', header = None)

mesh_atc_dictionary = pickle.load( open('cmap_data/mesh_atc_dictionary.pickle', "rb" ) )
mesh_drug = pd.read_csv('/scratch/PI/rbaltman/yash_margaret/gnbr_chemical_nodes_vfinal.csv')
drug_mesh_dictionary = pd.Series(mesh_drug.name.values,index=mesh_drug.uri).to_dict()

with open('/scratch/groups/rbaltman/yash_margaret/diseaseSig_psygenet_GWAS.pickle', 'rb') as handle:
    psygenet_dictionary = pickle.load(handle)

def makeDisNetwork(disease_signature):
    x = set()
    for gene in disease_signature:
        x.add(gene)

    edgelist = string_edgelist[string_edgelist[0].isin(x) | string_edgelist[1].isin(x)]
    network = nx.from_pandas_edgelist(edgelist, source = 0, target = 1, edge_attr=2)
    return network, x

def returnShortestPaths(disSet, disNetwork):
    subgraph = nx.Graph()
    #missing connections
    for x in disSet:
        for y in disSet:
            if x in disNetwork.nodes and y in disNetwork.nodes and nx.has_path(disNetwork, source=x, target=y):
                path = nx.shortest_path(disNetwork, source=x, target=y)
                subgraph.add_path(path)
    return subgraph

def rankDrugsName(disease, disSet):
    disease_drug = {}
    total = len(disSet)
    for index, row in expert_drug_gene_list.iterrows():
        intersection = disSet.intersection(row['union'])
        intersecting = len(intersection)
        drugSigLength = len(row['union'])
        frac_intersecting = intersecting/(total*drugSigLength)
        disease_drug[(disease, row['drug_name'])] = frac_intersecting
    df = pd.DataFrame(list(disease_drug.items()))
    disease_drug = pd.DataFrame(df[0].tolist(), index=df.index)
    disease_drug[2] = df[1]
    disease_drug[3] = disease_drug.groupby(0)[2].rank(ascending=False)
    disease_drug.columns = ['disease', 'drug', 'fraction_intersecting', 'rank']
    return disease_drug

def rankDrugsMESH(disease, disSet):
    disease_drug = {}
    total = len(disSet)
    for index, row in expert_drug_gene_list.iterrows():
        intersection = disSet.intersection(row['union'])
        intersecting = len(intersection)
        drugSigLength = len(row['union'])
        frac_intersecting = intersecting/(total*drugSigLength)
        disease_drug[(disease, row['drug'])] = frac_intersecting
    df = pd.DataFrame(list(disease_drug.items()))
    disease_drug = pd.DataFrame(df[0].tolist(), index=df.index)
    disease_drug[2] = df[1]
    disease_drug[3] = disease_drug.groupby(0)[2].rank(ascending=False)
    disease_drug.columns = ['disease', 'drug', 'fraction_intersecting', 'rank']
    disease_drug['drug_name'] = disease_drug['drug'].map(drug_mesh_dictionary)
    disease_drug['atc'] = disease_drug['drug'].map(mesh_atc_dictionary)
    return disease_drug

#for all diseases
alcohol_network, alcohol_set = makeDisNetwork(psygenet_dictionary['Alcohol_use_disorders'])
bpd_network, bpd_set = makeDisNetwork(psygenet_dictionary['Bipolar_disorders_and_related_disorders'])
cannabis_network, cannabis_set = makeDisNetwork(psygenet_dictionary['Cannabis_use_disorders'])
cocaine_network, cocaine_set = makeDisNetwork(psygenet_dictionary['Cocaine_use_disorders'])
depress_network, depress_set = makeDisNetwork(psygenet_dictionary['Depressive_disorders'])
drug_ind_psych_network, drug_ind_psych_set = makeDisNetwork(psygenet_dictionary['Drug-induced_psychosis'])
scz_network, scz_set = makeDisNetwork(psygenet_dictionary['Schizophrenia_spectrum_and_other_psychotic_disorders'])
sub_ind_depress_network, sub_ind_depress_set = makeDisNetwork(psygenet_dictionary['Substance_drug_induced_depressive_disorder'])

#create submodules for all diseases
alcohol_subgraph = returnShortestPaths(alcohol_set, alcohol_network)
print("alcohol")
bpd_subgraph = returnShortestPaths(bpd_set, bpd_network)
print("bpd")
cannabis_subgraph = returnShortestPaths(cannabis_set, cannabis_network)
print("cannabis")
cocaine_subgraph = returnShortestPaths(cocaine_set, cocaine_network)
print("cocaine")
depress_subgraph = returnShortestPaths(depress_set, depress_network)
print("depression")
drug_ind_psych_subgraph = returnShortestPaths(drug_ind_psych_set, drug_ind_psych_network)
print("dip")
scz_subgraph = returnShortestPaths(scz_set, scz_network)
print("scz")
sub_ind_depress_subgraph = returnShortestPaths(sub_ind_depress_set, sub_ind_depress_network)
print("sid")


alcohol_df = rankDrugsMESH('Alcohol_use_disorders', set(alcohol_subgraph.nodes))
bpd_df = rankDrugsMESH('Bipolar_disorders_and_related_disorders', set(bpd_subgraph.nodes))
cannabis_df = rankDrugsMESH('Cannabis_use_disorders', set(cannabis_subgraph.nodes))
cocaine_df = rankDrugsMESH('Cocaine_use_disorders', set(cocaine_subgraph.nodes))
depress_df = rankDrugsMESH('Depressive_disorders', set(depress_subgraph.nodes))
drug_ind_psych_df = rankDrugsMESH('Drug-induced_psychosis', set(drug_ind_psych_subgraph.nodes))
scz_df = rankDrugsMESH('Schizophrenia_spectrum_and_other_psychotic_disorders', set(scz_subgraph.nodes))
sub_ind_depress_df = rankDrugsMESH('Substance_drug_induced_depressive_disorder', set(sub_ind_depress_subgraph.nodes))


alcohol_df.to_csv('/scratch/groups/rbaltman/yash_margaret/updatedSignatureResults/alcohol_drugrank.csv')
bpd_df.to_csv('/scratch/groups/rbaltman/yash_margaret/updatedSignatureResults/bpd_drugrank.csv')
cannabis_df.to_csv('/scratch/groups/rbaltman/yash_margaret/updatedSignatureResults/cannabis_drugrank.csv')
cocaine_df.to_csv('/scratch/groups/rbaltman/yash_margaret/updatedSignatureResults/cocaine_drugrank.csv')
depress_df.to_csv('/scratch/groups/rbaltman/yash_margaret/updatedSignatureResults/depress_drugrank.csv')
drug_ind_psych_df.to_csv('/scratch/groups/rbaltman/yash_margaret/updatedSignatureResults/drug_ind_psych_drugrank.csv')
scz_df.to_csv('/scratch/groups/rbaltman/yash_margaret/updatedSignatureResults/scz_drugrank.csv')
sub_ind_depress_df.to_csv('/scratch/groups/rbaltman/yash_margaret/updatedSignatureResults/sub_ind_depress.csv')