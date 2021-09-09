# Reference compendium of non-redundant TF motifs

Here is the general workflow for systematically clustering redundant motifs by similarity to remove redundancy. This approach was adapted from the awesome [Vierstra et al., Nature 2020](https://www.nature.com/articles/s41586-020-2528-x) and applied to 6,502 TF motif models from multiple species (mostly focused on *Drosophila* and human TFs).

## Included motif databases

6,502 TF motif models were obtained from [iRegulon](http://iregulon.aertslab.org/collections.html) covering the following databases:
- [Bergman](http://bergmanlab.genetics.uga.edu/?page_id=274) (version 1.1; [Down et al., 2007](https://www.ncbi.nlm.nih.gov/pubmed/17238282))
- [CIS-BP](http://cisbp.ccbr.utoronto.ca/) (version 1.02; [Weirauch et al., 2014](https://www.ncbi.nlm.nih.gov/pubmed/25215497))
- [FlyFactorSurvey](http://pgfe.umassmed.edu/ffs/) (2010; [Zhu et al., 2011](https://www.ncbi.nlm.nih.gov/pubmed/21097781))
- [HOMER](http://homer.salk.edu/homer/) (2010; [Heinz et al., 2010](https://www.ncbi.nlm.nih.gov/pubmed/20513432))
- [JASPAR](http://jaspar.genereg.net/) (version 5.0_ALPHA; [Mathelier et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/26531826))
- Stark (2007; [Stark et al., 2007](https://www.ncbi.nlm.nih.gov/pubmed/17994088))
- [iDMMPMM](http://autosome.ru/iDMMPMM/) (2009; [Kulakovskiy and Makeev, 2009](https://link.springer.com/article/10.1134/S0006350909060013))

TF motif models were downloaded from [here](https://resources.aertslab.org/papers/iregulon/motifColl-10k-all-public.tar.gz) in cluster-buster format (see here for details on the format: https://aertslab.org/#data-resources-all, under *SUPPLEMENTARY MATERIAL TO THE IREGULON PAPER*).

## Requirements

- R 3.5.1
  - data.table
  - motifStack
  - TFBSTools
- Meme (https://meme-suite.org/meme/)
- Tomtom (http://meme-suite.org/doc/download.html)

## Scripts used to create compendium of non-redundant TF motifs

**Motif_clustering_Drosophila.sh**
- Step 1: Prepare motif databases
- Step 2: Compute pair-wise motif similarity (using TOMTOM)
<br/><br/>

**Motif_clustering.R**
- Step 3: Hierarchically cluster motifs by similarity (distance: correlation, complete linkage)
Below is a heatmap representation of motifs clustered by simililarity and clusters identified cutting the dendrogram at height 0.8.
<img src="https://data.starklab.org/almeida/Motif_clustering/Clusters_heatmaps/All_motifs_hierarchically_clustered_heatmap_pairwise_similarity_scores.png" alt="drawing" width="500"/>

You can check the position of all motif clusters in the heatmap using the heatmaps at https://data.starklab.org/almeida/Motif_clustering/Clusters_heatmaps/ (cluster highlighted on x- and y-axis on red).
Example of [cluster 30 highlighted](https://data.starklab.org/almeida/Motif_clustering/Clusters_heatmaps/Highlight_cluster_30.png).
<br/><br/>

- Step 4: Annotation of motif clusters: clusters were manually curated and annotated with the respective motif types
- Step 5: Plot motif logos aligned for each motif cluster ([example of motifs from cluster 30 - GATA/1](https://data.starklab.org/almeida/Motif_clustering/Clusters_logos/Cluster30_GATA.1_44motifs.pdf)). All motif logos per cluster at at https://data.starklab.org/almeida/Motif_clustering/Clusters_logos/.
<br/><br/>

**Create_consensus_TF_motif_database.Rmd**
- Step 6: Curate metadata information with cluster information and save PWM models into single R object [TF_clusters_PWMs.RData](https://data.starklab.org/almeida/Motif_clustering/TF_clusters_PWMs.RData)
<br/><br/>

*TF_clusters_PWMs.RData* is a list containing:
- metadata table with information for each TF motif model
- list with all PWMs (log-odds, position weight matrix)
- list with all PWMs (position probability matrix) 

More information on PWM (Position weight matrix): https://en.wikipedia.org/wiki/Position_weight_matrix
<br/><br/>

**Scan genome using all motif models**  
These motif models can be used to scan any DNA sequence of interest in R as follows:
```
library(motifmatchr) (https://bioconductor.org/packages/release/bioc/html/motifmatchr.html)
# using PWM log-odds
motif_ix <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds,
                        Sequences,
                        genome = "BSgenome.Dmelanogaster.UCSC.dm3", p.cutoff = 1e-4, bg="genome", out = "scores")
```

## Enrichment of TF motifs in developmental and housekeeping Drosophila S2 enhancers
Here is an example of motif enrichment analysis in developmental and housekeeping Drosophila S2 enhancers over negative regions ([volcano plots](https://data.starklab.org/almeida/Drosophila_enhancers_motif_enrichment/Motif_enrichment_volcano_plots.pdf)) using this TF motif database. To remove motif redundancy, only the most significant TF motif per motif cluster was shown.

## Questions
If you have any questions/requests/comments please contact me at [bernardo.almeida94@gmail.com](mailto:bernardo.almeida94@gmail.com)

## Acnowledgments
We thank Gert Hulselmans and [Stein Aerts](https://aertslab.org/) for sharing the TF motif PWM collection. This approach was inspired on [Jeff Vierstra](https://www.vierstra.org/)'s work (see here for more details on his [Non-redundant TF motif matches genome-wide for mouse and human TFs](https://www.vierstra.org/resources/motif_clustering)).