rm(list = ls())

library(factoextra)
library(cluster)
library(dendextend)
library(dplyr)

setwd("/home/aline/teletravail/tests_SC/")

# the two following files will be sorted in the same order
# Distance matrix data
distance <-read.csv("input_files/average_distance_table.tsv",header=T,row.names=1,sep="\t")

# sorting of the distance matrix 
distance_sorted <-distance[order(rownames(distance),decreasing=F),]
distance_sorted_transposed <- t(distance_sorted)
distance <- distance_sorted_transposed[order(rownames(distance),decreasing=F),]


# Binning data (the bin should be given only by a number)
binning_ref <- read.csv("input_files/contig_bin_correspondance_with_length.tsv",header=T,row.names=1,sep="\t")
binning_ref <- binning_ref[order(rownames(binning_ref),decreasing = F),]

# number of contigs
nb_tot_contigs = nrow(binning_ref)

# Number of selected bins (add one if the count begins by 0)
nb_bins = max(binning_ref[,1])+1

# the names of columns of the output are c("cluster", "neighbor", "sil_width"))
SC_REF<-as.data.frame(silhouette(as.integer(binning_ref[,1]), as.dist(distance))[,c(1,3)])

# Creation of a result dataframe with Bin names, bins number, sc ref and contig length as columns
results <-as.data.frame(cbind(binning_ref[,1],SC_REF[,2],binning_ref[,2]))
row.names(results) <-row.names(binning_ref)
names(results)= c("bins_ref","SC","contig_length")
binning_ref_SC_sort <- results[order(results[,2],decreasing=T),]

write.csv(binning_ref_SC_sort,"output_files/binning_ref_SC_sort.csv")

# sorting of the distance matrix by decreasing SC
distance_sort <- as.data.frame(cbind(results[,2],distance))
distance_sort <- distance_sort[order(distance_sort[,1],decreasing=T),]
distance_sort <- distance_sort[,c(2:ncol(distance_sort))]
distance_sort <- t(distance_sort)
distance_sort <- as.data.frame(cbind(results[,2],distance_sort))
distance_sort <- distance_sort[order(distance_sort[,1],decreasing=T),]
distance_sort <- distance_sort[,c(2:ncol(distance_sort))]
# write.csv(distance_sort,"output_files/average_distance_table_ref_SC_sort.csv")

# Distance matrix data
distance <-distance_sort

# Binning data
binning_ref <- binning_ref_SC_sort


# tours 1 to N

# modify here the number of 'tours'
N <- 5

for (i in 1:N) {

tour <- as.character(i)
  
# Threshold contig
nb_tot_contigs = nrow(binning_ref)
bad_SC = which(binning_ref$SC <0)
nb_contig_bad_SC <- length(bad_SC)
nb_contig_good_SC <- nb_tot_contigs-nb_contig_bad_SC


# number of bins
nb_bins = max(binning_ref[,1])

# Silhouette score associated computing for reference binning
SC_REF<-as.data.frame(silhouette(as.integer(binning_ref[,1]), as.dist(distance))[,c(1,3)])
sum_SC_REF <- sum(SC_REF[,2])
sum_SC_average_per_bin_REF <- sum(aggregate(SC_REF$sil_width, by=list(cluster=SC_REF$cluster), FUN=mean)[,2])
sum_SC_contig_length_REF <-sum(SC_REF[,2]*binning_ref[,3])


#Initialisation of the best som of average SC/bin corresponding to the reference
best_SC<-as.data.frame(silhouette(as.integer(binning_ref[,1]), as.dist(distance))[,c(1,3)])
best_sum_weigthed_contig_length <- sum(best_SC[,2]*binning_ref[,3])

#Initialisation of a new binning
binning_new <- binning_ref

for (contig in nb_contig_good_SC: nb_tot_contigs)
{
  best_bin = binning_ref[contig,1]

  for (bin_number in 0:nb_bins)
  {
    binning_new[contig,1] <- bin_number
    new_SC <- as.data.frame(silhouette(as.integer(binning_new[,1]), as.dist(distance))[,c(1,3)])
    new_sum_weighted_contig_length <- sum(new_SC[,2]*binning_ref[,3])
    if (new_sum_weighted_contig_length > best_sum_weigthed_contig_length)
    {
      best_sum_weigthed_contig_length = new_sum_weighted_contig_length
      best_SC = new_SC
      best_bin = bin_number
     }
    
    binning_new[contig,1] = best_bin
  }
}

  
#Computing of SC for the new binning
SC_new<-as.data.frame(silhouette(as.integer(binning_new[,1]), as.dist(distance))[,c(1,3)])
  
# creation of a new binning dataframe
results <-cbind(binning_ref[,1],SC_REF[,2],SC_new,binning_ref[,3])
row.names(results) <-row.names(binning_ref)
names(results)= c("Bins","SC_ref",paste0("bin_round", tour),"SC","contig_length")
binning_round_n <-results[,c(3:5)]
binning_round_n_SC_sort <- binning_round_n[order(binning_round_n[,2],decreasing=T),]


#write.csv(results,paste0("output_files/binning_round",tour,".csv"))
write.csv(binning_round_n_SC_sort,paste0("output_files/binning_round",tour,"_SC_sort.csv"))

# sorting of the distance matrix by decreasing SC (results[,4])
distance_sort <- as.data.frame(cbind(results[,4],distance))
distance_sort <- distance_sort[order(distance_sort[,1],decreasing=T),]
distance_sort <- distance_sort[,c(2:ncol(distance_sort))]
distance_sort <- t(distance_sort)
distance_sort <- as.data.frame(cbind(results[,4],distance_sort))
distance_sort <- distance_sort[order(distance_sort[,1],decreasing=T),]
distance_sort <- distance_sort[,c(2:ncol(distance_sort))]

#write.csv(distance_sort,paste0("output_files/average_distance_table_round",tour,"_SC_sort.tsv"))

# basic information for next tour
distance <- distance_sort
binning_ref <- binning_round_n_SC_sort  
}


