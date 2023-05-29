# Network Science (INFONPM) Final Assignment
# Timo Damm & Pascalle Doorn

# general preparations
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import ruptures as rpt

# Data import and preparations
# import
header=np.array(["t","i","j"])
data=pd.read_csv("/home/dot/PAOS_ResMa/complex_systems_track/CSRP/data/tij_InVS15.dat", sep = ' ', header = None) #main data with interactions
metadata=pd.read_csv("/home/dot/PAOS_ResMa/complex_systems_track/CSRP/data/metadata_InVS15.txt")
#data #optional, for inspection. 

def count_distinct(arr):
    distinct = set(arr)
    return len(distinct)
hist=data[0]
#print(count_distinct(hist))
#interactions = plt.hist(x=hist, bins= 18488, color='#0504aa',
                            #alpha=0.7, rwidth=0.85)
#plt.xlabel('Day')
#plt.ylabel('Number of Interactions')
#plt.title('Interaction Frequencies')

#splitting into days
#print(len(hist))
n=np.linspace(1,78249,num=78249)
seconds=plt.scatter(x=hist, y=n)


diff = [hist[i+1]-hist[i] for i in range(len(hist)-1)]

#find 9 largest jumps and their location
sorted(range(len(diff)), key=lambda x: diff[x])
cuts= sorted(range(len(diff)), key=lambda x: diff[x])[-9:]

#cut in these 9 places to have 10 days data separately
#print(cuts)
day1 = data.iloc[0:12161]
day2 = data.iloc[12162:20415]
day3 = data.iloc[20416:27319]
day4 = data.iloc[27320:35774]
day5 = data.iloc[35775:43512]
day6 = data.iloc[43513:50764]
day7 = data.iloc[50765:59893]
day8 = data.iloc[59894:64570]
day9 = data.iloc[64571:72511]
day10 = data.iloc[72511:]
day1 = day1.to_csv()
#print(day1)

# making graphs from the separate days
G1 = nx.read_edgelist("day1", data = ['1','2'])

# basic model (Louvain for all days, unweighted, no historical continutity)

# -----basic model (Louvain for all days, unweighted, no historical continutity)------
# Execute Louvain algorithm 100 times
num_runs = 100
best_partition = None
best_modularity = float('-inf')

for _ in range(num_runs):
    partition = community.best_partition(G1) #replace G1 with day subgraph of your day here to compute modularity for other days.
    modularity = community.modularity(partition, G1)

    if modularity > best_modularity:
        best_partition = partition
        best_modularity = modularity

# Print the best partition and modularity
print("Best Modularity:", best_modularity)
#print("Best Partition:", best_partition) # best partition optimal

# weighted edges model (Louvain with weighted edges for all days, no historical continutity)

# Best Combination of Local Communities (BCLC) (has historical continutity)
# ----Best Combination of Local Communities (BCLC) (has historical continutity)-----
# Run Louvain algorithm n times on G1 and G2
n = 100
results_1 = []
results_2 = []
for i in range(n):
    # Run Louvain algorithm on G1
    partition_1 = community.best_partition(G1)
    results_1.append(partition_1)

    # Run Louvain algorithm on G2
    partition_2 = community.best_partition(G2)
    results_2.append(partition_2)

# Compute cumulative modularity for all n*n partitions
max_modularity = float('-inf')
best_partition_pair = None

for partition_1 in results_1:
    for partition_2 in results_2:
        # Create a mapping between partitions
        mapping = {}
        for node in partition_1:
            mapping[node] = next((key for key, value in partition_2.items() if value == partition_1[node]), None)

        # Apply the mapping to partition_1
        mapped_partition_1 = {node: partition_2[mapping[node]] for node in partition_1 if mapping[node] is not None}

        # Compute modularity using mapped_partition_1 and partition_2
        modularity = community.modularity(mapped_partition_1, G1) + community.modularity(partition_2, G2)

        # Check if modularity is higher than current maximum
        if modularity > max_modularity:
            max_modularity = modularity
            best_partition_pair = (mapped_partition_1, partition_2)

# Access the best partition pairs
best_partition_1, best_partition_2 = best_partition_pair

# Print the results
print("Best partition pair:")
print("Partition 1:", best_partition_1)
print("Partition 2:", best_partition_2)
print("Overall modularity:", max_modularity)

# Day 3
# Run Louvain algorithm n times on G2 and G3
n = 100
results_2 = []
results_3 = []
for i in range(n):
    # Run Louvain algorithm on G2
    partition_2 = community.best_partition(G2)
    results_2.append(partition_2)

    # Run Louvain algorithm on G3
    partition_3 = community.best_partition(G3)
    results_3.append(partition_3)

# Compute cumulative modularity for all n*n partitions
max_modularity = float('-inf')
best_partition_pair = None

for partition_2 in results_2:
    for partition_3 in results_3:
        # Create a mapping between partitions
        mapping = {}
        for node in partition_2:
            if node in partition_3:
                mapped_node = next((key for key, value in partition_3.items() if value == partition_2[node]), None)
                if mapped_node is not None:
                    mapping[node] = mapped_node

        # Apply the mapping to partition_2
        mapped_partition_2 = {node: partition_3[mapping[node]] for node in mapping if partition_3[mapping[node]] in G3.nodes}

        # Check if all nodes in mapped_partition_2 are present in G2
        if all(node in G2.nodes for node in best_partition_2):
            # Compute modularity using best_partition_2 and mapped_partition_2
            modularity = community.modularity(best_partition_2, G2) + community.modularity(mapped_partition_2, G3)

            # Check if modularity is higher than current maximum
            if modularity > max_modularity:
                max_modularity = modularity
                best_partition_pair = (best_partition_2, mapped_partition_2)

# Access the best partition pairs
best_partition_2, best_partition_3 = best_partition_pair

# Print the results
print("Best partition pair:")
print("Partition 2:", best_partition_2)
print("Partition 3:", best_partition_3)
print("Overall modularity:", max_modularity)

#Day 4-10
# Create a list of graphs G4 to G10
graphs = [G4, G5, G6, G7, G8, G9, G10]

# Initialize the best partition pairs with best_partition_2 and None
best_partition_pairs = [(best_partition_2, None)] * len(graphs)

# Compute cumulative modularity for all n*n partitions
max_modularity = float('-inf')

# Loop over all graphs
for idx, graph in enumerate(graphs):
    results = []

    # Run Louvain algorithm n times on the current graph
    for _ in range(n):
        partition = community.best_partition(graph)
        results.append(partition)

    # Loop over all partitions of the current graph
    for partition in results:
        # Create a mapping between partitions
        mapping = {}
        for node in best_partition_3:
            mapped_node = next((key for key, value in partition.items() if value == best_partition_3[node]), None)
            if mapped_node is not None:
                mapping[node] = mapped_node

        # Apply the mapping to best_partition_3
        mapped_partition = {node: partition.get(mapping[node], -1) for node in best_partition_3}

        # Compute modularity using best_partition_3 and the mapped partition
        modularity = community.modularity(best_partition_3, G3) + community.modularity(mapped_partition, graph)

        # Check if modularity is higher than the current maximum
        if modularity > max_modularity:
            max_modularity = modularity
            best_partition_pairs[idx] = (best_partition_3, mapped_partition)

# Print the results
for idx, graph in enumerate(graphs):
    best_partition, mapped_partition = best_partition_pairs[idx]
    print(f"Best partition pair for G{idx+4}:")
    print("Partition 3:", best_partition)
    print(f"Partition {idx+4}:", mapped_partition)
    print("Overall modularity:", max_modularity)
    print()