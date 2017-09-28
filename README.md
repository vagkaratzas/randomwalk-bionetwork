# randomwalk-bionetwork
Starting from a specific disease, the user is allowed to command random walkers to explore a multi level network of their own input choosing and are returned an edgelist with the walked edges and their aggregate step count values.

Please pre-install required libraries in order to run this script.
Required libraries: igraph

Input arguments:
1. starting disease node, string, must be contained in the Disease-Gene-DisGeNET.txt.
2. Character string containing the choices for stub-networks. Starting network contains Diseases and Genes and the user has the ability to add pathways and drugs. Character string is represented like: pathway-drug
3. number of walkers, integer, default=1
4. number of steps per walker, integer, default=1000
5. restart alarm, integer, the walker restarts at the starting disease node after the designated number of steps is reached in the restart alarm, default=0
6. memory capacity, integer, remembers self and memory capcity-1 previous visited nodes, default=1
7. leviFlight, 0-1, if 1 each walker uses a cauchy distribution to decide the number of steps per round, default=0

Example: Rscript randomwalk-bionetwork.R "Idiopathic Pulmonary Fibrosis" "pathway-drug" 1 1000 50 2 0

*Make sure Rscript is in the PATH environment variable and you have your input files in the current directory.
