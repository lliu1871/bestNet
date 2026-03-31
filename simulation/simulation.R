###################################################
# simulation for transphylo
###################################################
setwd("./R/network/transnetwork/")
source("sim_outbreak.R")
library(TransPhylo)
set.seed(42)

numsim <- 10
transphylo_accurracy <- 1:numsim
for (i in 1:numsim) {
    # simulate outbreak
    outbreak <- simulate_outbreak(infection_rate = 3.0, removal_rate = 3.0, target_size = 10)
    numcase <- nrow(outbreak)

    # a transmission tree of the first 50 cases generated from outbreak
    timetree <- build_nj_tree(outbreak, add_coalescence = FALSE, theta = 10^-4, mu = 1, plot = FALSE)
    ptree <- ptreeFromPhylo(timetree, dateLastSample = 2007.964)
    res <- inferTTree(ptree, mcmcIterations = 10000, w.shape = 10, w.scale = 0.1, dateT = 2008)

    # transmission edges and accuracy
    matrixwiw <- computeMatWIW(res, burnin = 0.5)
    netedges <- transmission_edges_matrixwiw(matrixwiw)
    p_match <- length(which(as.numeric(netedges[, 2]) - outbreak[, 2] == 0)) / numcase
    transphylo_accurracy[i] <- p_match
}
mean(transphylo_accurracy)

###################################################
# simulation for bestnet
###################################################
numsim <- 10
set.seed(42)
bestnet_accurracy <- 1:numsim
for (i in 1:numsim) {
    # simulate outbreak
    outbreak <- simulate_outbreak(infection_rate = 3.0, removal_rate = 3.0, target_size = 10)
    write.csv(outbreak, paste("outbreak", i, ".csv", sep = ""))
    numcase <- nrow(outbreak)

    # save temporal data
    tempfile <- paste("sim_temporal", i, ".csv", sep = "")
    write.csv(outbreak[, c(1, 5, 7)], file = paste("/Users/lliu/Library/CloudStorage/OneDrive-UniversityofGeorgia/Dropbox/Github/Julia/bestNet/data/", tempfile, sep = ""), row.names = FALSE)

    # simulate DNA
    phytree <- build_nj_tree(outbreak, add_coalescence = FALSE, theta = 10^-4, mu = 10^-4, plot = FALSE)
    phytreestring <- write.tree(phytree)
    write(phytreestring, "tree.tre")
    seqfile <- paste("seq", i, sep = "")
    system(paste("seq-gen -mHKY -l1000000 tree.tre >", seqfile))

    # calculate snp distances
    data <- read.dna(seqfile)
    distance <- dist.dna(data, model = "N", as.matrix = TRUE)
    index <- order(as.numeric(row.names(distance)))
    newdistance <- distance[index, index]
    snpfile <- paste("sim_snp", i, ".csv", sep = "")
    write.csv(newdistance, file = paste("/Users/lliu/Library/CloudStorage/OneDrive-UniversityofGeorgia/Dropbox/Github/Julia/bestNet/data/", snpfile, sep = ""), row.names = FALSE)
}

for (i in 1:numsim) {
    inputfile <- paste("parameter", i, ".csv", sep = "")
    outbreak_est <- summary_bestnet(inputfile, plot = FALSE)$transmission
    outbreak <- read.csv(paste("outbreak", i, ".csv", sep = ""))[, -1]
    diff2 <- outbreak_est - outbreak
    bestnet_accurracy[i] <- length(which(diff2$parent_id == 0)) / (numcase - 1)
}
mean(bestnet_accurracy)

########################################
# accurracy plot
########################################
plot.new()
par(mfrow = c(1, 1))
inf_rate <- 1:3
tran <- c(0.22, 0.25, 0.35)
bestnet <- c(0.54, 0.56, 0.62)
data <- cbind(tran, bestnet)
row.names(data) <- inf_rate
barplot(t(data), col = c("#a1adbe", "#5ea3a9"), beside = T, ylab = "Accuraccy", xlab = "infection rate", ylim = c(0, 1), cex.axis = 1.2, cex.names = 1.2, , cex.lab = 1.5)
legend("topright", legend = c("transPhylo", "Bayesian Model"), cex = 1.7, col = c("#a1adbe", "#5ea3a9"), pch = 15)
