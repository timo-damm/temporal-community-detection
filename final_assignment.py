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
# ---- modified Best Combination of Local Communities (BCLC) (has historical continutity)-----
# Day 1 and 2
# Initialize variables
n = 100
results_1 = []
results_2 = []
modularity_1 = []
modularity_2 = []

# Run Louvain algorithm on G1
for i in range(n):
    partition_1 = community.best_partition(G1)
    modularity_1.append(community.modularity(partition_1, G1))
    results_1.append(partition_1)

# Run Louvain algorithm on G2
for i in range(n):
    partition_2 = community.best_partition(G2)
    modularity_2.append(community.modularity(partition_2, G2))
    results_2.append(partition_2)

# Compute cumulative modularity
cumulative_modularity = np.sum(modularity_1) + np.sum(modularity_2)

# Compute unstandardised Jaccard index
max_jaccard_index = 0
max_jaccard_partition_1 = None
max_jaccard_partition_2 = None

for partition_1 in results_1:
    for partition_2 in results_2:
        intersection = sum(1 for i, j in zip(partition_1.values(), partition_2.values()) if i == j)
        jaccard_index = intersection / len(partition_1)
        
        if jaccard_index > max_jaccard_index:
            max_jaccard_index = jaccard_index
            max_jaccard_partition_1 = partition_1
            max_jaccard_partition_2 = partition_2

# Print the results
print("Cumulative Modularity:", cumulative_modularity)
print("Modularity of Partition 1:", max(modularity_1))
print("Modularity of Partition 2:", max(modularity_2))
print("Max Jaccard Index:", max_jaccard_index)
print("Partition 1:", max_jaccard_partition_1) 
print("Partition 2:", max_jaccard_partition_2) #these two partitions maximise both Jaccard and modularity, therefore also historical continuity

#Day 3-10
# Create a list of graphs
graphs = [G3, G4, G5, G6, G7, G8, G9, G10]

# Initialize variables
n = 100
results = []
modularity = []
partition = []

# Run Louvain algorithm on G3
for i, G in enumerate(graphs):
    for j in range(n):
        part = community.best_partition(G)
        modularity.append(community.modularity(part, G))
        results.append(part)
        partition.append((i + 3, j))

# Compute cumulative modularity
cumulative_modularity = np.sum(modularity)

# Compute unstandardised Jaccard index
max_jaccard_index = 0
max_jaccard_partition = None

for part_1 in results:
    for part_2 in results:
        intersection = sum(1 for i, j in zip(part_1.values(), part_2.values()) if i == j)
        jaccard_index = intersection / len(part_1)
        
        if jaccard_index > max_jaccard_index:
            max_jaccard_index = jaccard_index
            max_jaccard_partition = part_1

# Find the partition that maximizes both modularity and Jaccard index for each graph
for i in range(3, 10):
    G = graphs[i - 3]
    max_modularity = -1
    max_jaccard_index = 0
    max_partition = None
    
    for j in range(n):
        part = community.best_partition(G)
        mod = community.modularity(part, G)
        intersection = sum(1 for a, b in zip(part.values(), max_jaccard_partition.values()) if a == b)
        jaccard_index = intersection / len(part)
        
        if mod > max_modularity and jaccard_index > max_jaccard_index:
            max_modularity = mod
            max_jaccard_index = jaccard_index
            max_partition = part
    
    max_jaccard_partition = max_partition

    print(f"Graph G{i}")
    print("Max Modularity:", max_modularity)
    print("Max Jaccard Index:", max_jaccard_index)
    print("Partition:", max_partition)
    print()
    
    # Set up figure and axes
    fig, ax = plt.subplots(figsize=(6, 6))

    # Draw the graph
    pos = nx.fruchterman_reingold_layout(G)
    nx.draw_networkx(G, pos, ax=ax, with_labels=False, node_color='darkblue', edge_color='gray', node_size=15)

    # Create the custom legend
    legend_text = "Modularity: {:.4f}\nJaccard Index: {:.4f}".format(max(modularity), max_jaccard_index)
    plt.text(0.05, 0.95, legend_text, transform=ax.transAxes, fontsize=12, verticalalignment='top')

    # Add a title indicating the day
    plt.title("Day {} Contact Network".format(i))

    # Adjust spacing and remove axis ticks
    plt.tight_layout()
    ax.tick_params(left=False, bottom=False)

    # Save the graph as a PNG file
    #plt.savefig("day{}_network.png".format(i), dpi=300)
    
# Visualisations (modify if needed for other days than 1)
#Set up figure and axes
fig, ax = plt.subplots(figsize=(6, 6))

# Draw the graph
pos = nx.fruchterman_reingold_layout(G1)
nx.draw_networkx(G1, pos, ax=ax, with_labels=False, node_color='darkblue', edge_color='gray', node_size = 15)

# Set figure title and axis labels
plt.title("Day 1 Contact Network")

# Create the legend
legend_text = "Modularity: {:.4f}\nJaccard Index: {:.4f}".format(max(modularity), max_jaccard_index)
plt.text(0.05, 0.95, legend_text, transform=ax.transAxes, fontsize=12, verticalalignment='top')

# Adjust spacing and remove axis ticks
plt.tight_layout()
ax.tick_params(left=False, bottom=False)

#plt.savefig("day1_network.png", dpi=300)
plt.show()