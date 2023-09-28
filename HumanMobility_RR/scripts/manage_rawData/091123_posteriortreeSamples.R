
library(ape)
library(lubridate)
library(dplyr)
s=68 ### random GPSC choice
resbd <- readRDS(paste0("./data/phylogenies/bd_GPSC",s))
################################################################
## Get tips and node heights from final BD tree
################################################################
Ntips = length(resbd$tree$tip.label)

## Tip dates
tip_dates_from_CI_bd = resbd$CI[1:Ntips,1]

## Node dates
nroot = length(resbd$tree$tip.label)+1 ## Checked and it's the root 
genetic_distance_mat = dist.nodes(resbd$tree) ## Pairwise time of differences inferred with BD
rownames(genetic_distance_mat) = c(resbd$tree$tip.label, Ntips+(1:(Ntips-1)))
colnames(genetic_distance_mat) = c(resbd$tree$tip.label, Ntips+(1:(Ntips-1)))
distance_to_root = genetic_distance_mat[nroot,]
root_height = tip_dates_from_CI_bd[which(resbd$tree$tip.label == names(distance_to_root[1]))] - distance_to_root[1]  ## Take one tip, doesn't matter which tip is used
nodes_dates_mean_from_CI_bd = root_height + distance_to_root[Ntips+(1:(Ntips-1))]
################################################################

################################################################
## Get mean tip dates and node heights from BD trace
################################################################
tip_dates_from_trace_bd = apply(resbd$record[,1:Ntips], MARGIN = 2, mean)
nodes_dates_mean_from_trace = apply(resbd$record[,(Ntips+1):(2*Ntips-1)], MARGIN = 2, mean)
################################################################

################################################################
## Compare the two
################################################################
par(mfrow=c(1,2))
plot(tip_dates_from_trace_bd, tip_dates_from_CI_bd)
abline(a = 0, b=1, lty=2)
plot(nodes_dates_mean_from_CI_bd, nodes_dates_mean_from_trace)
abline(a = 0, b=1, lty=2)
################################################################

################################################################
## Reconstruct one tree from the trace
################################################################
tree_ID = 155 ## I just chose a random number
tip_dates = resbd$record[tree_ID,1:Ntips]
node_dates = resbd$record[tree_ID,(Ntips+1):(2*Ntips-1)]
tipandnodes_dates = c(tip_dates, node_dates)

reconstructed_tree = resbd$tree ## Take BD tree for template
## Replace edge lengths, based on the tip and node dates that were extracted above
reconstructed_tree$edge.length = tipandnodes_dates[reconstructed_tree$edge[,2]] - tipandnodes_dates[reconstructed_tree$edge[,1]]

## Plot new tree
plot(reconstructed_tree)
################################################################

################################################################
## Compare distance matrices
################################################################
dist_reconstructed_tree = dist.nodes(reconstructed_tree)
dist_final_BD_tree = dist.nodes(resbd$tree)

plot(dist_final_BD_tree, dist_reconstructed_tree)
################################################################
