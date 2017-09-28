# randomwalk-bionetwork
# randomwalk-bionetwork-multithread is the multithread version of the script, runnining parallel walkers depending on the number of cores of your processor and greatly reducing execution times. Its debug output (execution times of walkers) is appended on a file named parellel_log.txt.

Starting from a specific disease, the user is allowed to command random walkers to explore a multi level network of their own input choosing and are returned an edgelist with the walked edges and their aggregate step count values.

* Please pre-install required libraries in order to run this script.
* Required libraries: igraph
* Required libraries for multithread version: foreach, doParallel

Input arguments:
1. Starting Disease node (character string). Must be contained in the Disease-Gene-DisGeNET.txt. Selecting a starting disease is mandatory - no default value.
2. Character string containing the choices for stub-networks. Starting network contains Diseases and Genes and the user has the ability to add pathways and drugs (which will connect to the main network via Gene nodes). Character string is represented like: pathway-drug, pathwaydrug, etc.
3. Number of walkers (integer), default=1.
4. Number of steps per walker (integer), default=1000.
5. Restart alarm (integer). The walker restarts at the starting disease node after the designated number of steps is reached in the restart timer, default=0.
6. Memory capacity (integer). Remembers self + (memory capcity-1) previous visited nodes, default=1.
7. Levi Flight (0-1). If 1 each walker uses a cauchy distribution to decide the number of steps per round, default=0.

Example: Rscript randomwalk-bionetwork.R "Idiopathic Pulmonary Fibrosis" pathwaydrug 2 100 50 2 0

* Make sure Rscript is in the PATH environment variable and you have your input files in the current directory.
* All input files can be downloaded from this page.
