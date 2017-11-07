#!/usr/bin/env Rscript
# randomwalk-bionetwork
#Starting from a specific disease, the user is allowed to command random walkers to explore a multi level network of their own input choosing and are returned an edgelist with the walked edges and their aggregate step count values.

#*Please pre-install required libraries in order to run this script.
#*Required libraries: igraph

#Input arguments:
#1. Starting Disease node (character string). Must be contained in the Disease-Gene-DisGeNET.txt. Selecting a starting disease is mandatory - no default value.
#2. Character string containing the choices for stub-networks. Starting network contains Diseases and Genes and the user has the ability to add pathways and drugs (which will connect to the main network via Gene nodes). Character string is represented like: pathway-drug, pathwaydrug, etc.
#3. Number of walkers (integer), default=1.
#4. Number of steps per walker (integer), default=1000.
#5. Restart alarm (integer). The walker restarts at the starting disease node after the designated number of steps is reached in the restart timer, default=0.
#6. Memory capacity (integer). Remembers self + (memory capcity-1) previous visited nodes, default=1.
#7. Levi Flight (0-1). If 1 each walker uses a cauchy distribution to decide the number of steps per round, default=0.

#Example: Rscript randomwalk-bionetwork.R "Idiopathic Pulmonary Fibrosis" pathwaydrug 2 100 50 2 0

#*Make sure Rscript is in the PATH environment variable and you have your input files in the current directory.
#*All input files can be downloaded from this page.


### Functions
makeNetwork <- function(vectorOfFiles) { #function that creates the main graph out of the user's inputs for stub subnetworks
	edgeList <- read.delim(vectorOfFiles[1], header = FALSE, sep = "\t") #read first default file
	if (length(vectorOfFiles) > 1){ #if user choses more stub files than just the default (disease-genes)
		for (i in 2:length(vectorOfFiles)){ #for the rest of the network files chosen by user (except first which is already read)
			tempEdgelist <- read.delim(vectorOfFiles[i], header = FALSE, sep = "\t") #read next input edgelist
			edgeList <- rbind(edgeList, tempEdgelist) #combine to main edgelist
		}
	}
	netNodes <- cbind(as.character(edgeList[,1]),as.character(edgeList[,2])) #make sure all nodes will be interpreted as strings
	netNodes[is.na(netNodes)] <- "NA" #if value NA found, replace it with the string NA
	graph <- graph_from_edgelist(netNodes, directed = FALSE) #make graph from massive edgelist
	graph <- simplify(graph) #removes possible reoccuring edges, case sensitive
	return(graph)
}

randGenerator <- function() { #generates a random alpharethemtic string
	a <- do.call(paste0, replicate(5, sample(LETTERS, 1, TRUE), FALSE)) #1 string with 5 characters
	paste0(a, sprintf("%04d", sample(9999, 1, TRUE)), sample(LETTERS, 1, TRUE)) #paste with 4digits and 1 last character
}

memoryCheck <- function(memoryVector, memoryVectorIterator, nodeIndex) { #checks and iterates memory for given node
	if (memoryVectorIterator == length(memoryVector)){ #if memory full, start from beggining/overwrite first element
		memoryVector[1] <- nodeIndex
		memoryVectorIterator <- 1
	}
	else {
		memoryVectorIterator <- memoryVectorIterator + 1 #increment memory vector iterator by 1
		memoryVector[memoryVectorIterator] <- nodeIndex #put next node id in next memory spot
	}
	return(list(memV = memoryVector, memVI = memoryVectorIterator))
}

rwStep <- function(currentNode,memoryVector,memoryVectorIterator,fileConn) { #function to decide next node transition
	x <- 1 #default number of steps of round
	if (leviFlight==1) x <- ceiling(abs(rcauchy(1)))#if levi walk, sample number of steps for this round. This may create edges, previously unexisting in the final Edgelists
	spamRestarts <- 0 #will need to count memory restart attempts in order to restart from starting disease
	while (TRUE){ #repeat until iteratingNode not in RW's memory
		iteratingNode <- currentNode #need to remember starting node of this round in case of memory restart
		for (i in 1:x){	
			nodeIndex <- match(iteratingNode,V(graph)$name) #use current name to figure out the current index
			neighbors <- as.matrix(neighbors(graph, nodeIndex)) #use current index to find out neighbors
			neighborNodePosition <- sample(1:length(neighbors), 1) #pick the position of a uniformally-random neighbor
			iteratingNode <- V(graph)[neighbors[neighborNodePosition]]$name #receive name of the next node in case of Levi Flight iterations
		}
		nodeIndex <- match(iteratingNode,V(graph)$name) #use next name to figure out the next index
		if (!(nodeIndex %in% memoryVector) | (memoryCapacity<1)){ #continue only if node index not in memory and walker has memory capacity
			break #from while TRUE
		}
		else {
			if (spamRestarts<5) { #check if walker trapped, 5 times neighbor in memory
				spamRestarts <- spamRestarts + 1
			}
			else {
				#print("###Walker Trapped! Restarting at input Disease.###")
				#print(currentNode)
				cat(sprintf("%s\t%s\t1\n", currentNode, diseaseRestart), file = fileConn)
				memResult <- memoryCheck(memoryVector, memoryVectorIterator, nodeIndex) #put starting disease in memory again, if trapped and decided to restart there
				memoryVector <- memResult$memV
				memoryVectorIterator <- memResult$memVI
				return(list(itn = diseaseRestart, mv = memoryVector, mvi = memoryVectorIterator)) #the walker will jump to the restart point by brute force, even though it might be on the memory and create an edge
			}
		}
	}#end while TRUE
	if (memoryCapacity > 0){ #if memory exists, place in memory
		memResult <- memoryCheck(memoryVector, memoryVectorIterator, nodeIndex)
		memoryVector <- memResult$memV
		memoryVectorIterator <- memResult$memVI
	}#endif
	cat(sprintf("%s\t%s\t1\n", currentNode, iteratingNode), file = fileConn)
	return(list(itn = iteratingNode, mv = memoryVector, mvi = memoryVectorIterator))
}

randomWalk <- function() {
	library(igraph) #for reasons, this has to defined here as well
	randString <- randGenerator() #each walker will create his own steps and final edgelist files
	edgesFilename <- paste("temp_results/Steps_", randString, "_total", totalSteps, "_restart", restartAlarm, "_mem", memoryCapacity, "_levi", leviFlight, ".txt", sep="") #create unique edges filename with randstring and metadata
	currentNode <- diseaseRestart #start at disease chosen by user
	memoryVector <- vector(mode = "integer", length = 1) #initialization, if length remains 1, memory not used
	if (memoryCapacity > 0) { #create memory vector if selected and place first element
		memoryVector <- vector(mode="integer", length = memoryCapacity) #memory vector creation of unique node positions
		memoryVector[1] <- match(currentNode, V(graph)$name) #first element -> unique number of restart point
		memoryVectorIterator <- 1 #will be carried to rwStep together with the vector in order to correctly overwrite oldest nodes
	}
	fileConn<-file(edgesFilename, "a") #append
	restartTimer <- 0 #initialize restart counter, will restart after restartAlarm steps
	for(i in 1:totalSteps){
		if (restartTimer==restartAlarm){ #if time to restart
			currentNode <- diseaseRestart #hops to input Disease node without creating an edge and without counting this as a step
			restartTimer <- 0 #reinitialize restart timer
		}
		result <- rwStep(currentNode,memoryVector,memoryVectorIterator,fileConn) #list of variables as result
		currentNode <- result$itn
		memoryVector <- result$mv
		memoryVectorIterator <- result$mvi
		restartTimer <- restartTimer + 1
	}
	close(fileConn) #closing file connection
	return(list(efn = edgesFilename, rndStr = randString))
}

makeFinalEdgeList <- function(stepsFile, randStr) { #function that returns edgeist out of the total steps taken
	expfile <- stepsFile #paste("temp_results/", stepsFile, sep = "") #filename
	steps <- read.delim(expfile, FALSE, strip.white = TRUE, sep="\t") #read steps File
	file.remove(stepsFile) #deleting steps' file
	stepNodes <- cbind(as.character(steps[,1]), as.character(steps[,2])) #make sure edge names are characters
	stepNodes[is.na(stepNodes)] <- "NA" #convert to string NA to avoid errors
	g = graph.edgelist(stepNodes, directed = FALSE) #make graph from edgelist
	E(g)$weight = as.numeric(steps[,3]) #all weights are currently 1 but will be incremented with simplify, this step is needed for simplify
	g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "sum")) #increment same edges
	graphEdgeList <- get.edgelist(g) #final edgelist from graph, 2 columns, edges
	graphEdgeList <- cbind(graphEdgeList, rep(0, length(graphEdgeList[,1]))) #add extra 3rd column for weights
	graphEdgeList[,3] <- E(g)$weight #calculate and put weights in 3rd column
	graphEdgeList <- graphEdgeList[order(graphEdgeList[,3], decreasing = TRUE),] #order edges by decreasing weight
	resultEdgelistFilename <- paste("SuperBioNetworkEdgeList_", randStr, "_total", totalSteps, "_restart", restartAlarm, "_mem", memoryCapacity, "_levi", leviFlight, ".txt", sep="") #create file name for Result EdgeList
	expfile <- paste("temp_results/", resultEdgelistFilename, sep = "")
	write.table(graphEdgeList, expfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
	return(graphEdgeList)
}

parallelWalk <- function() { #main function used for parallel walkers
	results <- randomWalk() #no arguements, global variables (user's input) will be accessed
	makeFinalEdgeList(results$efn,results$rndStr)
	rm(results) #memory flush
}

aggregateWalkers <- function() {
	edgelists <- list.files("temp_results") #make a vector with all edgelists' filenames
	aggregatedEdgeList <- read.delim(paste("temp_results/", edgelists[1], sep = ""), header = FALSE, sep = "\t") #read first edgelist to initialize
	if (length(edgelists) > 1) {
		for (i in 2:length(edgelists)){ #for the rest of the network files chosen by user (except first which is already read)
			tempEdgelist <- read.delim(paste("temp_results/", edgelists[i], sep = ""), header = FALSE, sep = "\t") #read next input edgelist
			aggregatedEdgeList <- rbind(aggregatedEdgeList, tempEdgelist) #combine to main edgelist
		}
	}
	filesToDelete<-paste(getwd(), "/temp_results/", list.files("temp_results"), sep="") #variable containing filenames of each walker's steps, to be deleted
	unlink(filesToDelete) #deleting all step files
	netNodes <- cbind(as.character(aggregatedEdgeList[,1]),as.character(aggregatedEdgeList[,2])) #make sure all nodes will be interpreted as strings
	netNodes[is.na(netNodes)] <- "NA" #if value NA found, replace it with the string NA
	graph <- graph_from_edgelist(netNodes, directed = FALSE) #make graph from massive edgelist
	E(graph)$weight = as.numeric(aggregatedEdgeList[,3]) #all weights are currently 1 but will be incremented with simplify, this step is needed for simplify
	graph <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "sum")) #increment same edges
	graphEdgeList <- get.edgelist(graph) #final edgelist from graph, 2 columns, edges
	graphEdgeList <- cbind(graphEdgeList, rep(0, length(graphEdgeList[,1]))) #add extra 3rd column for weights
	graphEdgeList[,3] <- E(graph)$weight #calculate and put weights in 3rd column
	### order
	orderedGraphEdgeList <- graphEdgeList[order(as.numeric(graphEdgeList[,3]), decreasing = TRUE),] #order edges be decreasing weight
	### topPercentage
	#topRankedLines <- ceiling(topRankedPercentage / 100 * length(orderedGraphEdgeList[,1])) #calculate top Lines according to user input percentage
	#topRankedEdges <- head(orderedGraphEdgeList, topRankedLines)
	### writeToOutputFile
	randString <- randGenerator() #for unique edgelist name
	outfileName <- paste(randString,"aggregatedEdgelist.txt", sep="_")
	write.table(orderedGraphEdgeList, outfileName, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t") #write topRanked edges
}

### Main
### Libraries
library(igraph) #library for the main graph creation - edgelists (read and write)
library(foreach) #libary for multi threading
library(doParallel) #libary for multi threading
### User Inputs
args = commandArgs(trailingOnly = TRUE) #allow use of args
if (length(args) == 0) { #check if correct number of args
	stop("Need at least a starting disease node string as argument, as written in Disease-Gene-DisGeNET.txt", call.=FALSE) #if not disease name, throw error and stop
} else if (length(args) == 1) {
	# defaults
	args[2] <- ""
	args[3] <- 1
	args[4] <- 1000
	args[5] <- 0
	args[6] <- 1
	args[7] <- 0
}else if (length(args) == 2) {
	args[3] <- 1
	args[4] <- 1000
	args[5] <- 0
	args[6] <- 1
	args[7] <- 0
} else if (length(args) == 3) {
	args[4] <- 1000
	args[5] <- 0
	args[6] <- 1
	args[7] <- 0
} else if (length(args) == 4) {
	args[5] <- 0
	args[6] <- 1
	args[7] <- 0
} else if (length(args) == 5) {
	args[6] <- 1
	args[7] <- 0
} else if (length(args) == 6) {
	args[7] <- 0
}

#Random Walker Variable Values
diseaseRestart <- args[1] #"Idiopathic Pulmonary Fibrosis" #Disease Restart point
#Stub NetWork
networkFiles <- vector(mode="character") #initialization of vector with file names (edgelists) that will be iterated to create the main body of the SuperNetWork (graph)
networkFiles <- append(networkFiles, "Disease-Gene-DisGeNET.txt") #Default Network. Disease-Gene Edgelist from DisGeNET Curated Database
if(!file.exists("Disease-Gene-DisGeNET.txt")){ #checking if folder contains the subnetworks in the form of edgelists
		stop("Make sure you have downloaded the edgelist text files in the current directory", call.=FALSE) #need edgelists text files in the same directory
}
if (grepl("pathway", args[2])){
	networkFiles <- append(networkFiles, "Pathway_Gene_Reactome.txt") #Pathway-Gene (and additional molecules at the moment) from Reactome Pathways Gene Set
	print("Pathways will be added to the bionetwork.")
}
if (grepl("drug", args[2])){
	networkFiles <- append(networkFiles, "Drug-Gene-DGIdb.txt") #Drug-Gene (multiple interaction types) from DrugGeneInteraction database (DGIdb)
	print("Drugs will be added to the bionetwork.")
}
walkers <- as.integer(args[3]) #1 #number of walkers
if (walkers < 0) { #check if negative value for walkers
	print("Number of walkers must be a positive number. Choosing default value = 1")
	walkers <- 1
}
totalSteps <- as.integer(args[4]) #1000 #Total Steps
if (totalSteps < 0) { #check if negative value for walkers
	print("Number of total steps per walker must be a positive number. Choosing default value = 1000")
	totalSteps <- 1000
}
restartAlarm <- as.integer(args[5]) #0 #After how many steps to restart
if (restartAlarm < 0) { #check if negative value for walkers
	print("Restart Alarm counter must be a positive number. Choosing default value = 0")
	restartAlarm <- 0
}
memoryCapacity <- as.integer(args[6]) #1 #remembers current node + memoryCapacity-1 previously visited nodes
if (memoryCapacity < 0) { #check if negative value for walkers
	print("Memory capacity must be a positive number. Choosing default value = 1")
	memoryCapacity <- 1
}
leviFlight <- as.integer(args[7]) #0 #0 or 1, simple walk or Levi Flight (Cauchy distribution for number of steps per round)
if (leviFlight < 0) { #check if negative value for walkers
	print("Levi Flight choice must be either 0 or 1. Choosing default value = 0", call.=FALSE)
	leviFlight <- 0
}

#Execution
print("Creating SuperBioNetwork...")
graph <- makeNetwork(networkFiles) #construct main Graph
print("Graph ready!")
dir.create(file.path(getwd(), "temp_results")) #for each walkers results
filesToDelete<-paste(getwd(), "/temp_results/", list.files("temp_results"), sep="") #variable containing filenames of possible carbage in temp_results folder
unlink(filesToDelete) #deleting leftovers in temp_results folder
print('Starting Random Walks...')
no_cores <- detectCores() - 1 #Calculate the number of cores, -1 because the laptop will freeze if uses all cores for R, probably
cl <- makeCluster(no_cores, outfile = "parallel_log.txt") #Initiate cluster and open connection with outfile log
registerDoParallel(cl) #registering the backend
system.time(foreach(i=1:walkers, .combine=cbind) %dopar% { #system.time to time the whole program
	ptm <- proc.time() #start execution timer
	cat(sprintf("Starting Walker %i.\n", i)) #cat sprintf in order to write in output file of clusters, defined in makeCluster
	parallelWalk()
	elapsedtime <- proc.time() - ptm #execution timer
	cat(sprintf("Walker %i ran for %f seconds.\n",i,elapsedtime[3])) #print execution time of walker, elapsed time[3] is the elapsed time variable of the time vector
	rm(ptm) #memory flush
	gc() #garbage collector
})
stopCluster(cl)
aggregateWalkers() #create final output edgelist
