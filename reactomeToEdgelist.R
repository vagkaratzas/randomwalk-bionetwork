react <-  read.csv("ReactomePathways.gmt", header=FALSE, sep="\t") #length(react) -> max genes in a pathway
react_char <- sapply(react, as.character)
react_char[is.na(react_char)] <- "NA"
fileConn <- file("Pathway_Gene_Reactome.txt", "a") #a for append
for (i in 1:(length(react[,1]))){ #all rows
	for(j in 2:length(react)){ #all columns
		if (react_char[i,j] == "") break
		cat(sprintf("%s\t%s\n", react_char[i,1], react_char[i,j]), file = fileConn)
	}
}
close(fileConn)
