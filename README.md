# temporal-community-detection
This is a python script using temporal data for different community detection algorithms (Louvain, Louvain with weitghted edges, modified BCLC)
It is written for temporal networks and in this case uses a network of ten days worth of face-to-face interactions in a French public organisation. 
The data was collected and made available by the Sociopatterns collaboration and can be found (here)[http://www.sociopatterns.org/datasets/test/]. 
Essentially, we compare a naive implementation of the Louvain algorithm for community detection with one that uses weighted edges to encode repeated interactions and a modified version of the "Best Combination of Local Communities" algorithm.

Requires Python > 3.8 with pandas, numpy, networkx, matplotlib, ruptures and community.

The code is structured as follows: 
1. preparations and data import (including splitting the data into 10 temporal snapshots and making graphs for the separate days)
2. basic descriptive statistics
3. naive Louvain implementation
4. Louvain with weighted edges
5. modified BCLC algorithm
6. visualisations
