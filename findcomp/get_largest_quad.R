require(igraph)

##### Parameters for k-components analysis
input_graph  <- "largest_triconnected.net"
output_graph <- "largest_quadconnected.net"
kcomponent   <- 4

##### Functions

find_all_2seps <-function(g) {
  # Find all 2-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  excluded <- vector()

  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=3)               # Get all cliques of size >= 3
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    excluded <- append(excluded, isolated)      # Maintain list of included vertices
    if (length(worldly) == 2) {                 #   Keep track of 2-separators
      pc <- pc + 1                              #   since they may not be separators 
      pl[[pc]] <- sort(worldly)                 #   after excluding isolated nodes
    }
  }

  nex <- length(excluded)
  gp <- g - excluded
  numV <- length(V(gp))

  cat("    Searching for all 2-separators: #excluded: ", nex, "\n")

  for (i in 1:numV) {
    h <- gp - i
    arts <- articulation.points(h)
    for (art in sort(arts)) {
      pc <- pc + 1
      pl[[pc]] <- sort(c(V(gp)$orig[i], V(h)$orig[art]))
    }
  }

  pl <- unique(pl)
  cat("    Found 2-separators: ", length(pl), "\n")
  return(pl)
}


find_all_3seps <- function(g) {
  # Find all 3-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  excluded <- vector()

  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=4)               # Get all cliques of size >= 4
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    excluded <- append(excluded, isolated)      # Maintain list of included vertices
    if (length(worldly) == 3) {                 #   Keep track of 3-separators
      pc <- pc + 1                              #   since they may not be separators 
      pl[[pc]] <- sort(worldly)                 #   after excluding isolated nodes
    }
  }

  nex <- length(excluded)
  gp <- g - excluded
  numV <- length(V(gp))

  cat("    Searching for all 3-separators: #excluded: ", nex, "\n")

  for (i in 1:(numV-1)) {
    for (j in (i+1):numV) {
      nodes <- c(i,j)
      h <- gp  - nodes
      arts <- articulation.points(h)
      for (art in sort(arts)) {
        pc <- pc + 1
        pl[[pc]] <- sort(c(V(gp)$orig[i], V(gp)$orig[j], V(h)$orig[art]))
      }
    }
  }

  pl <- unique(pl)
  cat("    Found 3-separators: ", length(pl), "\n")
  return(pl)
}

find_tough_3seps <- function(g) {
  # Find all 3-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  excluded <- vector()

  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=4)               # Get all cliques of size >= 4
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    excluded <- append(excluded, isolated)      # Maintain list of included vertices
    if (length(worldly) == 3) {                 #   Keep track of 3-separators
      pc <- pc + 1                              #   since they may not be separators 
      pl[[pc]] <- sort(worldly)                 #   after excluding isolated nodes
    }
  }

  nex <- length(excluded)
  gp <- g - excluded
  numV <- length(V(gp))

  cat("    Searching for tough 3-separators: #excluded: ", nex, "\n")

  psbl <- vector()
  for (v in V(gp)) {
    ns <- unlist(neighbors(gp, v))
    for (i in ns) {
      for (j in ns) {
        if(i > j && are.connected(gp,i,j) == 0) {
          vdp <- vertex.disjoint.paths(gp,i,j)
          if (vdp == 3) {
            psbl <- append(psbl, v)
          }
        }
      }
    }
  }
 
  psbl <- unique(psbl)
  cat("psbl: ", sort(psbl), "\n")
  np <- length(psbl)

  if (np > 0) {
    for (i in 1:(np-2)) {
      for (j in (i+1):(np-1)) {
        for (k in (j+1):np) {
          gpp <- gp - c(psbl[i], psbl[j], psbl[k])
          #result <- is.separator(gp, c(psbl[i], psbl[j], psbl[k]))
          if(!is.connected(gpp)) {
            pc <- pc + 1
            pl[[pc]] <- sort( c(V(gp)$orig[psbl[i]], V(gp)$orig[psbl[j]], V(gp)$orig[psbl[k]] ))
          }
        }
      }
    }
  }

  pl <- unique(pl)
  cat("    Found 3-separators: ", length(pl), "\n")
  return(pl)
}

find_easy_3seps <- function(g) {
  # Find all 3-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  numV <- length(V(g))

  exclude <- rep.int(0, numV)                   # Initialize filter
  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=4)               # Get all cliques of size >= 4
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    exclude[isolated] <- 1
  }

  nex <- length(exclude[which(exclude == 1)])

  cat("    Searching for easy 3-separators: #excluded: ", nex, "\n")

  for (e in E(g)) {
    #RSS nodes <- ends(g,e)
    nodes <- get.edge(g,e)
    if (exclude[nodes[1]] == 1 || exclude[nodes[2]] == 1) {
      next
    }
    h <- g - c(nodes[1],nodes[2])
    arts <- articulation.points(h)
    for (art in sort(V(h)$orig[arts])) {
      pc <- pc + 1
      pl[[pc]] <- sort(c(nodes[1],nodes[2],art))
    }
  }

  pl <- unique(pl)
  cat("    Found 3-separators: ", length(pl), "\n")
  return(pl)
}

kcore_delete_nodes <- function(g, kcomponent) {
  # Reduce a graph to k-cores by deleting nodes
  for (i in 1:100) {
    d <- degree(g)
    lowd_nodes <- V(g)[which(d < kcomponent)]
    g <- g - lowd_nodes
    if (length(lowd_nodes) == 0) {
      break
    }
  }
  return(g)
}

bicomponent_processing <- function(g, kcomponent, candidate) {
  # Find bicomponents
  bc <- biconnected.components(g)

  # Get size of largest bicomponent
  cutoff <- 0
  for(i in 1:bc$no) { 
    if(length(bc$components[[i]]) > cutoff) {
      cutoff <- length(bc$components[[i]])
    }
  }
  
  # Loop over bicomponents
  for (j in 1:bc$no) {
    gp <- induced.subgraph(g, bc$components[[j]])

    # Get size of bicomponent and skip if largest
    nverts_before <- length(V(gp))
    nedges_before <- length(E(gp))
    if(nverts_before == cutoff) {
      next
    }    

    ### Reduce to k-cores before writing out
    gp <- kcore_delete_nodes(gp, kcomponent)
    nverts <- length(V(gp))
    nedges <- length(E(gp))

    if(nverts <= kcomponent) {
    }
    if(nverts > kcomponent) {
      nedges_clique <- (nverts*(nverts-1))/2

      # Find size of minimum separator
      msep <- minimum.size.separators(gp)
      if(length(msep) > 0) {
        msep_min <- length(msep[[1]])
      } else {
        msep_min <- 0
      }

      if(nedges == nedges_clique) {
        candidate <- candidate + 1
	ksize <- length(V(gp))
        outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
        write.graph(gp, outfile, format="pajek")
      }
      else if(msep_min >= kcomponent) {
        candidate <- candidate + 1
        ksize <- length(V(gp))
        outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
        write.graph(gp, outfile, format="pajek")
      }
      else if(nedges != nedges_clique) {

        mwBlocks  <- cohesive.blocks(gp)
        blocks    <- blocks(mwBlocks)
        cohesion  <- cohesion(mwBlocks)
        parent    <- parent(mwBlocks)
        nblocks   <- length(blocks)

        for (ib in 1:nblocks) {
          if(cohesion[ib] >= kcomponent && cohesion[parent[ib]] < kcomponent) {
            gpp <- induced.subgraph(gp, blocks[[ib]])
            candidate <- candidate + 1
            ksize <- length(V(gpp))
            outfile <- paste("p", kcomponent, "_c", cohesion[ib], "_s", ksize, "_id", candidate, ".net", sep="")
            write.graph(gpp, outfile, format="pajek")
          }
        }

      }
    }
  }
  return(candidate)
}

reduce_to_largest_bicomponent <- function(g) {
  # Reduce graph to largest bicomponent
  bc <- biconnected.components(g)
  lbc <- vector(length=bc$no)
  for(i in 1:bc$no) { lbc[i] <- length(bc$components[[i]]) }
  maxbc <- which.max(lbc)
  g <- induced.subgraph(g, bc$components[[maxbc]])
  g <- set.vertex.attribute(g, "orig", V(g), V(g))
  return(g)
}

apply_separators <- function(g, candidate, sepsize, kcomponent, pl) {
  ### Use 2-separators to reduce graph
  w <- g   # Make a copy of graph to work with
  for (i in 1:length(pl)) {
    nodes <- unlist(pl[[i]])
    h <- w - nodes
    clust <- clusters(h)
    ncl <- clust$no
    cutoff <- max(clust$csize) + sepsize

    for(j in 1:ncl) {
      setm <- (V(h)$orig[which(clust$membership == j)])
      setp <- union(nodes, setm)

      if(length(setp) == cutoff) {
        next
      }

      g[setm, V(g)] <- FALSE
      gp <- induced.subgraph(w, unlist(setp))

      ### Reduce to k-cores before writing out
      nverts_before <- length(V(gp))
      nedges_before <- length(E(gp))
      gp <- kcore_delete_nodes(gp, kcomponent)
      nverts <- length(V(gp))
      nedges <- length(E(gp))

      if(nverts <= kcomponent) {
      }
      if(nverts > kcomponent) {
        nedges_clique <- (nverts*(nverts-1))/2

        # Find size of minimum separator
        msep <- minimum.size.separators(gp)
        if(length(msep) > 0) {
          msep_min <- length(msep[[1]])
        } else {
          msep_min <- 0
        }

        if(nedges == nedges_clique) {
          candidate <- candidate + 1
          ksize <- length(V(gp))
          outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
          write.graph(gp, outfile, format="pajek")
        }
        else if(msep_min >= kcomponent) {
          candidate <- candidate + 1
          ksize <- length(V(gp))
          outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
          write.graph(gp, outfile, format="pajek")
        }
        else if(nedges != nedges_clique) {
          mwBlocks  <- cohesive.blocks(gp)
          blocks    <- blocks(mwBlocks)
          cohesion  <- cohesion(mwBlocks)
          parent    <- parent(mwBlocks)
          nblocks   <- length(blocks)

          for (ib in 1:nblocks) {
            if(cohesion[ib] >= kcomponent && cohesion[parent[ib]] < kcomponent) {
              gpp <- induced.subgraph(gp, blocks[[ib]])
              candidate <- candidate + 1
              ksize <- length(V(gpp))
              outfile <- paste("p", kcomponent, "_c", cohesion[ib], "_s", ksize, "_id", candidate, ".net", sep="")
              write.graph(gpp, outfile, format="pajek")
            }
          }

        }
      }

    }
  }
  return(list(g, candidate))
}


##### Start main program

mydate <- date(); cat(mydate, "\n")

g             <- read.graph(file=input_graph, format="pajek")
g             <- set.vertex.attribute(g, "orig", V(g), V(g))
candidate     <- 0 # Initialize k-candidate count
lbc_size_old  <- 0

cat("\nSTEP 1 - iteratively trim isolated nodes, reduce to k-cores, reduce to bicomponents\n")
for (imaster in 1:100) {
  cat("---- STEP 1 iteration ", imaster, "\n")

  ###### Process cliques with isolated nodes
  for (iii in 1:100) {
    count <- 0                                    # Keep track of number of processed cliques
    d <- degree(g)                                # Recalculate node degrees
    cl <- maximal.cliques(g, min=kcomponent)      # Get all cliques of size >= kcomponent
    for (i in 1:length(cl)) {                     # Loop over cliques
      clique <- unlist(cl[[i]])                   # Get members of clique
      dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
      clsize <- length(clique)                    # Get size of clique
      isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
      worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes

      # If the clique contains exactly (k-1) worldly nodes, trim from main graph
      if(length(worldly) == (kcomponent-1) && length(isolated) > 0) {
        count <- count + 1

        # Write out the cliques that are not part of larger k-component
        if(clsize > kcomponent) { 
          #RSS gp <- induced.subgraph(g, as_ids(clique))
          gp <- induced.subgraph(g, clique)
	  msep_min <- length(V(gp)) - 1
          candidate <- candidate + 1
          ksize <- length(V(gp))
          outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
          write.graph(gp, outfile, format="pajek")
        }
        #RSS g[as_ids(isolated), as_ids(clique)] <- FALSE
        g[isolated, clique] <- FALSE
      }

    }

    d <- degree(g)
    size <- length(V(g)[which(d > 0)])
    cat("  Clique processing: Vertex count = ", size, "\n")

    if(count == 0) break
  }

  ### Reduce to k-cores by deleting edges
  xold <- 0
  alld_nodes <- V(g)
  for (jjj in 1:100) {
    d <- degree(g)
    lowd_nodes <- V(g)[which(d < kcomponent)]
    if (length(lowd_nodes) == xold) {
      break
    }
    xold <- length(lowd_nodes)
    #RSS g[as_ids(lowd_nodes), as_ids(alld_nodes)] <- FALSE
    g[lowd_nodes, alld_nodes] <- FALSE
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  k-coring: Vertex count = ", size, "\n")

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
  if(lbc_size == lbc_size_old) break
  lbc_size_old <- lbc_size
}

cat("\nSTEP 2 - iteratively find/use 2-separators, reduce to bicomponents\n")
for (imaster in 1:100) {
  cat("---- STEP 2 iteration ", imaster, "\n")

  ### Find 2-separators
  pl <- find_all_2seps(g)
  if(length(pl) == 0) break

  ### Apply 2-separators
  sepsize <- 2
  tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
  g         <- tmp[[1]]
  candidate <- tmp[[2]]

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  Apply 2-separators: Vertex count = ", size, "\n")

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
}

cat("\nSTEP 3 - Find/apply easy 3-separators\n")
for (imaster in 1:100) {
  cat("---- STEP 3 iteration ", imaster, "\n")

  ### Find and apply easy 3-separators
  pl <- find_easy_3seps(g)
  if(length(pl) == 0) break
  sepsize <- 3
  tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
  g         <- tmp[[1]]
  candidate <- tmp[[2]]

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  Apply easy 3-separators: Vertex count = ", size, "\n")

  ### Reduce to k-cores by deleting edges
  xold <- 0
  alld_nodes <- V(g)
  for (jjj in 1:100) {
    d <- degree(g)
    lowd_nodes <- V(g)[which(d < kcomponent)]
    if (length(lowd_nodes) == xold) {
      break
    }
    xold <- length(lowd_nodes)
    #RSS g[as_ids(lowd_nodes), as_ids(alld_nodes)] <- FALSE
    g[lowd_nodes, alld_nodes] <- FALSE
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  k-coring: Vertex count = ", size, "\n")

  # --- After k-coring, may have new lower degree separators ---

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")

  ### Find and apply 2-separators
  for (ijk in 1:100) {
    pl <- find_all_2seps(g)
    if(length(pl) == 0) {break}
    sepsize <- 2
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]

    d <- degree(g)
    size <- length(V(g)[which(d > 0)])
    cat("  Apply 2-separators: Vertex count = ", size, "\n")
    
    ### Process bicomponents and save largest for further processing
    candidate <- bicomponent_processing(g, kcomponent, candidate)
    g <- reduce_to_largest_bicomponent(g)
    lbc_size <- length(V(g))
    cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
  }
}


# --------------------------
cat("\nSTEP 4 - Find/apply all 3-separators\n")
for (imaster in 1:100) {
  cat("---- STEP 4 iteration ", imaster, "\n")

  n3sep <- 0

  ### Find and apply tough 3-separators
  pl <- find_tough_3seps(g)
  n3sep <- n3sep + length(pl)
  if(length(pl) > 0) {
    sepsize <- 3
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  Apply tough 3-separators: Vertex count = ", size, "\n")

  ### Reduce to k-cores by deleting edges
  xold <- 0
  alld_nodes <- V(g)
  for (jjj in 1:100) {
    d <- degree(g)
    lowd_nodes <- V(g)[which(d < kcomponent)]
    if (length(lowd_nodes) == xold) {
      break
    }
    xold <- length(lowd_nodes)
    #RSS g[as_ids(lowd_nodes), as_ids(alld_nodes)] <- FALSE
    g[lowd_nodes, alld_nodes] <- FALSE
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  k-coring: Vertex count = ", size, "\n")

  # --- After k-coring, may have new lower degree separators ---

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")

  ### Find and apply 2-separators
  for (ijk in 1:100) {
    pl <- find_all_2seps(g)
    if(length(pl) == 0) {break}
    sepsize <- 2
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]

    d <- degree(g)
    size <- length(V(g)[which(d > 0)])
    cat("  Apply 2-separators: Vertex count = ", size, "\n")
    
    ### Process bicomponents and save largest for further processing
    candidate <- bicomponent_processing(g, kcomponent, candidate)
    g <- reduce_to_largest_bicomponent(g)
    lbc_size <- length(V(g))
    cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
  }

  ### Find and apply easy 3-separators
  pl <- find_easy_3seps(g)
  n3sep <- n3sep + length(pl)
  if(length(pl) > 0) {
    sepsize <- 3
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  Apply easy 3-separators: Vertex count = ", size, "\n")

  ### Reduce to k-cores by deleting edges
  xold <- 0
  alld_nodes <- V(g)
  for (jjj in 1:100) {
    d <- degree(g)
    lowd_nodes <- V(g)[which(d < kcomponent)]
    if (length(lowd_nodes) == xold) {
      break
    }
    xold <- length(lowd_nodes)
    #RSS g[as_ids(lowd_nodes), as_ids(alld_nodes)] <- FALSE
    g[lowd_nodes, alld_nodes] <- FALSE
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  k-coring: Vertex count = ", size, "\n")

  # --- After k-coring, may have new lower degree separators ---

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")

  ### Find and apply 2-separators
  for (ijk in 1:100) {
    pl <- find_all_2seps(g)
    if(length(pl) == 0) {break}
    sepsize <- 2
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]

    d <- degree(g)
    size <- length(V(g)[which(d > 0)])
    cat("  Apply 2-separators: Vertex count = ", size, "\n")
    
    ### Process bicomponents and save largest for further processing
    candidate <- bicomponent_processing(g, kcomponent, candidate)
    g <- reduce_to_largest_bicomponent(g)
    lbc_size <- length(V(g))
    cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
  }

  if(n3sep == 0) break

}




write.graph(g, output_graph, format="pajek")
mydate <- date(); cat(mydate, "\n")
