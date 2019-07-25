# ABioTrans
A Biostatistical tool for Transcriptomics Analysis

If you find ABioTrans useful, please cite our paper:  
Zou Y, Bui TT and Selvarajoo K (2019) ABioTrans: A Biostatistical Tool for Transcriptomics Analysis. Front. Genet. 10:499. doi: 10.3389/fgene.2019.00499

Full text (free access) can be found at https://www.frontiersin.org/articles/10.3389/fgene.2019.00499/full 

## How to set up
1. Install R from https://cran.r-project.org/ 
2. Install Rstudio from https://www.rstudio.com/ 
3. Download ABioTrans-master.zip on GitHub and unzip it. Please do not modify www inside ABioTrans folder.
4. Open the ABioTrans.R file using RStudio and click `RunApp` button on the topright. When you run the code for the first time, installation of the required R packages can take up to 30 minutes. For subsequent runs, ABioTrans will only take 30s.

You can start your analysis now!

## How to do analysis
### Home
* Choose an RNA-Seq data file in comma-separated value (.csv) format. 
* If you input raw data (read counts), please make sure that the first column contains gene names, and the read counts of each genotype (conditions: wildtype, mutants, replicates, etc.) are in the following columns. Each genotype column should have a column name.
* Along with raw read counts, you can provide gene length (base pair) information in two-column .csv file, with the first column specifying gene names, which must match the gene names in raw data file, and the second column specifying gene length in base pair. Gene length file is required for normalization methods for sequencing depth and gene length: RPKM, FPKM, TPM
* List of negative control genes (spike-in or stably expressed genes accross all samples), if available, should be contained in one-column .csv file. Negative control genes are required for Remove Unwated Variation (RUV) normalziation.
* If you input a normalised data file, it should have gene names in rows and genotypes in columns, following the usual format of files deposited in the GEO database. 
* Finally, a metadata table matching column names of data file to experimental conditions should be given in two-column .csv formate. Metadata table is required for differential expression analysis
* Move on to the preprocessing and analysis methods once a datafile is loaded.

### Preprocessing
 Specify the cut-off expression values (same unit to your input data file - either raw read counts or normalized expression), and the minimum number of columns (samples) whose expression is above threshold value. Normalization methods are available upon your input of supporting data files (gene length and negative control genes). Relative Log Expression (RLE) plot of raw and normalized data are displayed to compare the effects of normalziation

### Scatter plot
 Choose the variable on x and y axis. Take log transformation on the values if needed. The colors refer to kernel density estimation (KDE). You can download a single plot as PDF. You can also download all pairs of samples scatter plot in one PDF file, which may take some time to run. 

### Distributions
 Choose the column you want to fit and compare with any combinations of the six statistical distributions. You can adjust the slider or input the range for x-axis to zoom in to see the fitted curve.

### Pearson and Spearman correlations
 Choose a correlation method and view it in either correlation heatmap, correlation plot in circle or correlation matrix.

### PCA and K-means clustering
 To take sample from all the genes, choose a gene sample size and a gene sample order. By defualt the full sample is taken. Then the variance percentage of all principal components,  2-D plot of any PC-axis combination and 3-D plot will be shown. K-means clustering is available for PCA 2-D and 3-D plots.

### Differential Expression Analysis
* ABioTrans provides 3 Differential Expression (DE) Analysis methods for multiple replicate dataset: edgeR, DESeq2 and NOISeq. For data with single replicate in all experiment condition, NOISeq method can simulate technical replicates to carry out DE test.
* For edgeR and DESeq2, raw read counts datafile must be provided. For NOISeq, gene expression should be normalized (by select normalization method in preprocessing tab if raw counts is inputted, or by directly providing normalized gene expression)  
* Metadata file is required for DE Analysis. Please make sure metadata contains all column names from input datafile and match them with experimental condition
* Specify DE methods, two conditions you wish to compare (condition 2 is compared against condition 1), and threshold criteria. By convention, DE genes are thresholded at 0.05 False Discovery Rate (FDR) and 2-fold change
* Volcano plot of DE result, dispersion plot of input data, and heatmap of DE genes are only available for edgeR and DESeq2 methods.

### Heatmap
* Heat map and hierachical clustering can be carried out on DE genes resulted from DE analysis tab by selecting `DE result`. Otherwise, you can specify the minimum fold change and minimum number of samples (columns in input data file) satisfying the fold change. Finally, you need to specify the number of clusters on rows (genes), then hit `Submit`.
* The gene names in each cluster are displayed in the `Gene clusters` panel, corresponding to the heatmap you just generated.

### Noise
* Select a desired noise plot according to your data. By default the name of the first replicate in each genotype is used as the name of the genotype. You can specify the names of the genotypes, only make sure that all names are different (due to the mapping mechanism of Plotly). Please be reminded that the consecutive columns will be regarded as of the same genotype.
* Noise here refers to the average transcriptome noise — the squared coefficient of variation — defined as the variance (σ2) of expression divided by the square mean expression (μ2). 
* If `replicates` is selected, the noise within one genotype will be computed. 
* If `genotypes (average of replicates)` is selected, the gene expression data within one genotype will first be taken average of, and the average value will be used as the expression value of that certain genotype. Then noise will be computed between the anchor genotype and all the other genotypes. 
* If `genotypes (no replicate)` is selected, each sample column is treated as a different genotype. Noise will directly be computed between the anchor genotype and the rest of the genotypes. 

### Entropy
* First specify whether your data is a time series data.
* Entropy here refers to Shannon entropy. It is calculated for each sample column.
* If the data is not in time series, entropy values of all samples will directly be displayed. 
* Otherwise, if it is a time series data, specify the number of time points. For example, if the number of time points is 6, then sample 1 to 6 will be regarded as genotype A time 1 - time 6, sample 7 to 12 will be genotype B time 1- time 6, and so on and so forth. For a line chart, the entropy of each genotype will be shown in one line. 

### Gene ontology
* Provide a list of genes for each cluster in NCBI gene ID format, (in one-column .csv format, or the .csv file downloaded from DE analysis tab), background gene set, select the over-representation test (`clusterProfiler`, `GOstats`, or `enrichR`), then select the proper organism and gene identifier.
* Click Plot button and wait for around 15 to 20 seconds for the result to be generated. Please avoid launching a new analysis when the system is still in a busy state. A list of enriched GO terms will be displayed in the `Table` tab.
* If clusterProfiler and GOStats are chosen, a pie chart showing relative size of all associated level-2 GO terms will be displayed in `Pie chart` tab. Please note that the GO terms in pie chart might not be significantly enriched since no over-representation test was carried out to generate the pie chart.
* For clusterProfiler method, graph visualization of user-specified GO terms are displayed in `Graph` tab

## Download instructions
1. Scatter plot, distribution fit, correlation plot and heatmap can be saved as PDF by clicking `Download as PDF`. You can name the file before saving it. Also, you can directly drag the plot from the GUI to a folder on the computer.

2. PCA, noise, entropy and gene ontology pie chart are supported by Plotly. If you are using the GUI in R window or Safari, click `Download as PNG` and it may take 5 seconds to save the graph. You cannot name the file yourself, so please be aware of overwriting issues. Even if the plot is saved successfully, there may be a 404 Not found message after downloading. Please go to your working directory and check. If you cannot save the graph successfully, please try the following code the Rstudio command line. You only need to run this code once. (You can also simply take a screenshot of the graph.)
```R
webshot::install_phantomjs()
```
If you are using R window or Safari, and you cannot save PCA-3D graph successfully, please create a plotly account (free). Please go to this website on the instructions of creating a new account and Find your authentication API keys in your online settings (https://plot.ly/r/getting-started/). After that please go to the first two lines in ABioTrans.R file, remove both of the "#" signs, and replace the username and api key with your own. Save the file before you run. However, please be reminded that this  `Download as PNG` button only downloads the PCA-3D plot from the default angle.

3. For Chrome users, the `Download as PNG` button may be disabled. Please put your mouse over the graph region and click the left most camera sign to save the plot. It will be saved in the download folder on your computer. You are suggested to use Chrome because it saves you trouble from creating an account and changing the code. More importantly, you can download the PCA-3D plot from any angle.
