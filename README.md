# Rshiny-RNAseq
Final Project for BF591-- interactive Rshiny application displaying RNAseq analysis data

exploratory data set:  Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration
  - This dataset profiled gene expression with RNASeq in murine cardiac tissue across different stages of development to identify genes associated with the loss of these cell’s capacity to regenerate.

Included components:
Sample Information Exploration:
  - Tab with a summary of the table that includes a summary of the type and values in each column
  - Tab with a data table displaying the sample information, with sortable columns
  - Tab containing plots of continuous variables

Counts Matrix Exploration:
  - Tab with text or a table summarizing the effect of the filtering, including:
      - number of samples
      - total number of genes
      - number and percentage of genes passing current filter
      - number and percentage of genes not passing current filter
- Tab with diagnostic scatter plots, where genes passing filters are marked in a darker color, and genes filtered out are lighter:
    - median count vs variance (consider log scale for plot)
    - median count vs number of zeros
- Tab with a clustered heatmap of counts remaining after filtering
    - consider enabling log-transforming counts for visualization
- Tab with a scatter plot of principal component analysis projections. You may either:
    - allowing the user to select which principal components to plot in a scatter plot (e.g. PC1 vs PC2)
    - allowing the user to plot the top N principal components as a beeswarm plot

Differential Expression:
  - Tab with sortable table displaying differential expression results
  - Volcano Plot displaying DEGs and log2FoldChange

Visualization of Individual Gene Expression:
  - Content displaying a plot of the selected type with the normalized gene counts for the selected gene split out by the categorical variable chosen
