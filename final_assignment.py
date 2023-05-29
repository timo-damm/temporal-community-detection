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