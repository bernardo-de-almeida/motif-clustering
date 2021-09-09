
#####
# Step 1: Prepare motif databases
#####

# TF motif models were downloaded from https://resources.aertslab.org/papers/iregulon/motifColl-10k-all-public.tar.gz in cluster-buster format
# (see here for details on the format: https://aertslab.org/#data-resources-all, under *SUPPLEMENTARY MATERIAL TO THE IREGULON PAPER*)

# create Markov Background Model - order 3
fasta-get-markov -n m /groups/stark/genomes/dm3/dm3.fa dm3.3-order.markov

# convert PWM models to MEME format
chen2meme singletons/bergman*cb -bg dm3.3-order.markov > all.dbs.meme
chen2meme singletons/cisbp*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/flyfactorsurvey*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/homer*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/jaspar*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/stark*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/idmmpmm*cb -bg dm3.3-order.markov >> all.dbs.meme

#####
# Step 2: Compute pair-wise motif similarity
#####

tomtom \
	-dist kullback \
	-motif-pseudo 0.1 \
	-text \
	-min-overlap 1 \
	all.dbs.meme all.dbs.meme \
> tomtom.all.txt

# remove last lines that have details
head -n -4 tomtom.all.txt > tomtom.all.treated.txt

#####
# Step 3: Hierarchically cluster motifs by similarity in R (Motif_clustering.R)
#####
