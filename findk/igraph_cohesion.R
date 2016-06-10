require(igraph)

#input_graph  <- "largest_4connected.net"
input_graph  <- "largest_5connected.net"
#input_graph  <- "erdos_renyi_1000_0.02.net"

mydate <- date(); cat(mydate, "\n")
g <- read.graph(file=input_graph, format="pajek")
k <- cohesion(g)
cat ("final k = ", k, "\n")
mydate <- date(); cat(mydate, "\n")
