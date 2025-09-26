###################################################
# simulate temporal and snp data
###################################################
setwd("./R/network/transnetwork/")
source("sim_outbreak.R")
set.seed(42)

# simulate outbreak
outbreak <- simulate_outbreak(infection_rate = 3.0, removal_rate = 3.0, target_size = 10)
numcase <- nrow(outbreak)

# a tree plot for true transmissions
ttree_true <- ttree_from_transmission(outbreak, dateLastSample = 2008)
plot(ttree_true)
plot(ttree_true, type = "detailed", w.shape = 10, w.scale = 0.1)

# a time phylogenetic tree where split times are the infection times and tips are removal times
timetree <- build_nj_tree(outbreak, plot = TRUE)

# simulate a coalescent tree in mutation units based on the time tree
phytree <- sim_coaltree(timetree, theta = 1e-5, murate = 1e-4)
plot(phytree)

# simulate sequences based on the coalescent tree
write.tree(phytree, file = "tree.tre")
system("seq-gen -mHKY -l100000 tree.tre > seq")

# calculate snp distances
data <- read.dna("seq")
distance <- dist.dna(data, model = "N", as.matrix = TRUE)
index <- order(as.numeric(row.names(distance)))
newdistance <- distance[index, index]

# save temporal and snp data
write.csv(outbreak[, c(1, 5, 7)], file = "/Users/lliu/Library/CloudStorage/OneDrive-UniversityofGeorgia/Dropbox/Github/Julia/bestNet/data/sim_temporal.csv", row.names = FALSE)
write.csv(newdistance, file = "/Users/lliu/Library/CloudStorage/OneDrive-UniversityofGeorgia/Dropbox/Github/Julia/bestNet/data/sim_snp.csv", row.names = FALSE)

#######################################
# transphylo analysis
#######################################
plot(timetree)
axisPhylo(backward = FALSE)
ptree <- ptreeFromPhylo(timetree, dateLastSample = 2007.964)
plot(ptree)

w_shape <- 10
w_scale <- 0.1
datetime <- 2008

res <- inferTTree(ptree, mcmcIterations = 10000, w.shape = w_shape, w.scale = w_scale, dateT = datetime)

# transmission edges
matrixwiw <- computeMatWIW(res, burnin = 0.5)
netedges <- transmission_edges_matrixwiw(matrixwiw)
p_match <- length(which(as.numeric(netedges[, 2]) - outbreak[, 2] == 0)) / numcase
p_match

# Plot the results
med <- medTTree(res)
plot(med)
ttree <- extractTTree(med)
plot(ttree, type = "detailed", w_shape, w_scale)

##############################################
# bestnet analysis
##############################################
outbreak_est <- summary_bestnet("parameter_sim_50.csv", onsite_time = outbreak$onsite_time, removal_time = outbreak$removal_time)$transmission
diff2 <- outbreak_est - outbreak
accurracy <- length(which(diff2$parent_id == 0)) / (numcase - 1)

ttree_est <- ttree_from_transmission(outbreak_est, dateLastSample = 2008)
plot(ttree_est)
plot(ttree_est, type = "detailed", w_shape, w_scale)

transplot(outbreak_est, showlabel = F, style = 1, vertex_sizes = rep(5, 50), vertex_label_cex = rep(0.7, 50))
transplot(outbreak_est, showlabel = F, style = 2, vertex_sizes = rep(7, 50), vertex_label_cex = rep(0.7, 50))
transplot(outbreak_est, style = 3, vertex_sizes = rep(7, 50), vertex_label_cex = rep(0.7, 50))
transplot(outbreak_est, style = 4, vertex_sizes = rep(7, 50), vertex_label_cex = rep(0.7, 50))

##############################################
# network plots
###############################################

# Create a simple static network with 5 nodes
net <- network.initialize(5)
add.edge(net, 1, 2)
add.edge(net, 2, 3)
add.edge(net, 4, 5)

# Plot a timeline of the network's activity
el <- cbind(
    c(16, 13, 13, 10, 13, 16, 10, 13, 1, 10, 8, 1, 4, 4, 2, 2),
    1:16
)
# a vector of infection times
infectionTimes <- c(
    583, 494, 634, 40, 712, 701, 224, 719,
    674, 0, 749, 621, 453, 665, 709, 575
)
# make a network object, include the infection time
infTree <- network(el,
    vertex.attr = list(infectionTimes),
    vertex.attrnames = list("infectionTimes")
)

transmissionTimeline(infTree, time.attr = "infectionTimes")

# Create a simple graph
g <- graph_from_edgelist(
    matrix(c(1, 2, 2, 3, 3, 4, 4, 1), byrow = TRUE, ncol = 2)
)

# Add an edge attribute (e.g., weight or label)
E(g)$weight <- c(2.5, 1.2, 3.8, 2.0)
E(g)$label <- paste("w=", E(g)$weight)
V(g)$names <- as.character(1:4)

# Plot with edge labels
plot(
    g,
    edge.label = E(g)$label, # show edge attribute
    edge.width = E(g)$weight, # edge width proportional to weight
    edge.color = "green",
    edge.label.color = "black",
    edge.label.cex = 1.7,
    vertex.size = 25,
    vertex.color = "skyblue",
    vertex.label.color = "orange",
    vertex.label.cex = 2.7
)
