# Install and load required packages
library(ape)
library(igraph)
library(openxlsx)

# Note that before starting, you need to organize the USalign output into a CSV file

# Read the USalign output CSV file
USalign_output_data <- read.csv("example.csv")

# Extract unique PDB chain names
pdb_chains <- unique(c(as.character(USalign_output_data$PDBchain1), as.character(USalign_output_data$PDBchain2)))

# Create an empty matrix
TM_score_matrix_data <- matrix(0, nrow = length(pdb_chains), ncol = length(pdb_chains), dimnames = list(pdb_chains, pdb_chains))

# Fill in the matrix data
for (i in 1:nrow(USalign_output_data)) {
  pdb_chain1 <- USalign_output_data$PDBchain1[i]
  pdb_chain2 <- USalign_output_data$PDBchain2[i]
  tm1 <- USalign_output_data$TM1[i]
  tm2 <- USalign_output_data$TM2[i]
  TM_score_matrix_data[pdb_chain1, pdb_chain2] <- tm1
  TM_score_matrix_data[pdb_chain2, pdb_chain1] <- tm2
}

# Modify row and column names, assuming they are predicted monomer structures, removing chain information ".pdb:A"
rownames(TM_score_matrix_data) <- gsub("\\.pdb:A", "", rownames(TM_score_matrix_data))
colnames(TM_score_matrix_data) <- gsub("\\.pdb:A", "", colnames(TM_score_matrix_data))

# Change the TM_score of protein against itself in the matrix to 1
TM_score_matrix_data[TM_score_matrix_data == 0] <- 1

# Output the TM_score matrix as a CSV file
write.csv(TM_score_matrix_data, file = "TM_score_matrix_data.csv")

# Convert the TM_score matrix to a distance matrix
dist_matrix <- dist(TM_score_matrix_data)

# Hierarchical clustering
hclust_tree <- hclust(dist_matrix, method = "ward.D2")

# Convert to phylo object
phylo_tree <- as.phylo(hclust_tree)

# Export the structure-based clustering tree in Newick format
write.tree(phylo_tree, file = "structure_based_hclust_tree.newick")

# Set the threshold for the structure similarity network
threshold <- 0.75

# Create an edge list based on the similarity threshold and filter out self-loops
overthresholdedges <- which(TM_score_matrix_data >= threshold, arr.ind = TRUE)
overthresholdedges <- overthresholdedges[overthresholdedges[, 1] != overthresholdedges[, 2], ]

# Create an empty graph object
graph <- graph.empty()

# Add nodes
nodes <- rownames(TM_score_matrix_data)
graph <- add.vertices(graph, nv = length(nodes), name = nodes)

# Add edges
for (i in 1:nrow(overthresholdedges)) {
  graph <- add.edges(graph, c(overthresholdedges[i, 1], overthresholdedges[i, 2]))
}

# Since TM_scores are directional, we need to convert the directed graph to an undirected graph
graph <- as.undirected(graph, mode = "collapse")

# clustering
clusters <- fastgreedy.community(graph)

# Get the size of each cluster
cluster_sizes <- sizes(clusters)

# Sort by cluster size in descending order
sorted_clusters <- clusters[order(cluster_sizes, decreasing = TRUE)]

# Get members of each cluster
cluster_members <- membership(clusters)

# Identify singleton nodes (clusters with only a single member)
singleton_nodes <- names(cluster_members[cluster_members %in% which(sizes(clusters) == 1)])

# Create Cytoscape export file
cytoscape_export <- createWorkbook()

# Create edge sheet
addWorksheet(cytoscape_export, sheetName = "Edges")
writeData(cytoscape_export, sheet = "Edges", x = "Source", startCol = 1, startRow = 1)
writeData(cytoscape_export, sheet = "Edges", x = "Target", startCol = 2, startRow = 1)

# Get edge list
edges <- get.edgelist(graph)

# Fill in the edge sheet data
if (nrow(edges) > 0) {
  writeData(cytoscape_export, sheet = "Edges", x = V(graph)[edges[, 1]]$name, startCol = 1, startRow = 2)
  writeData(cytoscape_export, sheet = "Edges", x = V(graph)[edges[, 2]]$name, startCol = 2, startRow = 2)
}

# Find the starting and ending rows of the current edge sheet
last_edge_row <- nrow(edges) + 1  # Since the first row is the header, data starts from the second row

# Add singleton nodes to the last of the first column
if (length(singleton_nodes) > 0) {
  writeData(cytoscape_export, sheet = "Edges", x = singleton_nodes, startCol = 1, startRow = last_edge_row + 1)
}

# Save the Cytoscape-readable input file in Excel format
saveWorkbook(cytoscape_export, "structure_based_similarity_network_cytoscape_export.xlsx", overwrite = TRUE)

# Create an empty data frame to store node and cluster information
export_sorted_clusters <- data.frame(protein = character(), cluster_name = character(), stringsAsFactors = FALSE)

# Iterate through sorted_clusters and add data to the data frame
for (cluster_name in names(sorted_clusters)) {
  proteins <- sorted_clusters[[cluster_name]]
  # Add each protein and corresponding cluster_name to the data frame
  for (protein in proteins) {
    # Check if the protein is in singleton_nodes
    if (protein %in% singleton_nodes) {
      cluster_name <- "singleton"  # Change cluster_name to "singleton"
    }
    export_sorted_clusters <- rbind(export_sorted_clusters, data.frame(protein = protein, cluster_name = cluster_name))
  }
}

# Export node and cluster information as a CSV file
write.csv(export_sorted_clusters, "structure_based_similarity_network_export_sorted_clusters.csv", row.names = FALSE, quote = TRUE)


# Set styles for nodes and edges, and perform simple visualization of the structure similarity network in R
vertex_attrs <- list(
  size = 5,
  color = "lightblue",
  frame.color = "gray",
  font.size = 8,
  label.cex = 0.7
)
edge_attrs <- list(
  width = E(graph)$weight * 1,
  color = "gray"
)
plot(graph, vertex.label = NA, vertex.color = vertex_attrs$color, 
     vertex.frame.color = vertex_attrs$frame.color, vertex.size = vertex_attrs$size,
     edge.width = edge_attrs$width, edge.color = edge_attrs$color,
     edge.arrow.size = 0.5, edge.width = 0.5)

