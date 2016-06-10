require(igraph)
require(parallel)
cores=12 # Hard coding for now

#input_graph  <- "largest_0connected.net"
#input_graph  <- "largest_1connected.net"
#input_graph  <- "largest_2connected.net"
input_graph  <- "largest_3connected.net"
#input_graph  <- "largest_4connected.net"
#input_graph  <- "largest_5connected.net"
#input_graph  <- "erdos_renyi_1000_0.10.net"

g <- read.graph(file=input_graph, format="pajek")
k <- 10000 # Initialize to something ridiculously high

###### Functions #####
nvdp <- function(index, g)
{
  nodes	<- unlist(index)
  result <- vertex.disjoint.paths(g, nodes[1], nodes[2])
  return(result)
}

getpairs <- function(v, g)
{
  ns <- unlist(neighbors(g, v))
  nslen <- length(ns)
  pairs <- list()
  pc <- 0

  # Get number of vertex disjoint paths between
  # the neighbors of the candidate vertex
  for (i in ns) {
    for (j in ns) {
      if(i < j && are.connected(gp,i,j) == 0) {
        pc <- pc + 1
        pairs[[pc]] <- c(i,j)
      }
    }
  }
  return(pairs)
} 

##### End Functions #####


mydate <- date(); cat(mydate, "\n")

# Check to see if graph is disconnected
# If so, then graph is k=0
if (!is.connected(g)) {
  k <- 0
  cat ("final k = ", k, "\n")
  mydate <- date(); cat(mydate, "\n")
  q(save='no')
}

# Check to see if it contains articulation points
# If so, then graph is k=1
if (length(articulation.points(g)) > 0) {
  k <- 1
  cat ("final k = ", k, "\n")
  mydate <- date(); cat(mydate, "\n")
  q(save='no')
}

# At this point, we have at least bicomponent
# Do a quick test to get minimal node degree
kmax1 <- min(degree(g))
if (kmax1 == 2) {
  k <- 2
  cat ("final k = ", k, "\n")
  mydate <- date(); cat(mydate, "\n")
  q(save='no')
}

# Find isolated nodes using minimum clique size of three
excluded <- vector()
d <- degree(g)                                # Recalculate node degrees
cl <- maximal.cliques(g, min=3)               # Get all cliques of size >= 3
for (i in 1:length(cl)) {                     # Loop over cliques
  clique <- unlist(cl[[i]])                   # Get members of clique
  dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
  clsize <- length(clique)                    # Get size of clique
  worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes
  isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes

  # Maintain list of included vertices
  excluded <- append(excluded, isolated) #0.6.6
  #excluded <- append(excluded, as_ids(isolated)) #1.0.0

  nworldly  <- length(worldly)
  nisolated <- length(isolated)
  if (nworldly < k && nisolated > 0) {
    k <- length(worldly)
    # If clique has two worldly vertices, then we have bicomponent
    if (k == 2) {
      cat ("final k = ", k, "\n")
      mydate <- date(); cat(mydate, "\n")
      q(save='no')
    }
  }
}
kmax2 <- k

# Done with quick tests, report upper limit on k
k <- min(kmax1, kmax2)
cat("k-max from degrees = ", kmax1, "\n")
cat("k-max from clique  = ", kmax2, "\n")
cat("k-max              = ", k, "\n")

# Construct a new graph that excludes all of the isolated nodes
gp <- g - excluded
numV <- length(V(gp))
cat("#nodes excluded in simplified graph: ", length(excluded), "\n")

# Get list of all possible pairs to consider
pairs <- unlist(mclapply(V(gp), getpairs, gp, mc.cores=cores), recursive=FALSE)
cat("Start vertex disjoint path calculations\n")
cat("#pairs        = ", length(pairs), "\n")
pairs <- unique(pairs)
cat("#unique pairs = ", length(pairs), "\n")

results <- unlist(mclapply(pairs, nvdp, gp, mc.cores=cores))
kmin <- min(results)
cat("Minimum k based on vertex disjoint paths = ", k, "\n")
k <- min(k, kmin)

cat ("final k = ", k, "\n")
mydate <- date(); cat(mydate, "\n")
