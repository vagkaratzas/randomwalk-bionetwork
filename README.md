# randomwalk-bionetwork
# UNDER CONSTRUCTION! TODO: add more stub networks
Starting from a specific disease, the user is allowed to command random walkers to explore a multi level network of their own input choosing and are returned an edgelist with the walked edges and their aggregate step count values.

Make sure Rscript is in the PATH environment variable and you have your input files in the current directory.

Please pre-install required libraries in order to run this script.
Required libraries: igraph

Input arguments:
1. starting disease node, string, must be contained in the Disease-Gene-DisGeNET.txt.
2. number of walkers, integer, default=1
3. number of steps per walker, integer, default=1000
4. restart alarm, integer, the walker restarts at the starting disease node after the designated number of steps is reached in the restart alarm, default=0
5. memory capacity, integer, remembers self and memory capcity-1 previous visited nodes, default=1
6. leviFlight, 0-1, if 1 each walker uses a cauchy distribution to decide the number of steps per round, default=0

Example: Rscript randomwalk-bionetwork.R "Idiopathic Pulmonary Fibrosis" 1 1000 50 2 0
