# ABioTransPlus
A Biostatistical tool for Transcriptomics Analysis

[ABioTransPlus](http://combio-sifbi.org/ABioTrans/) is a web server for transcriptome data analysis and visualization. It supports the analysis of microarray and RNASeq data and performs ten different bio-statistical analyses that cover the common analytics for gene expression data. Furthermore, it gives the users access to several bioinformatics tools to perform 12 different bioinformatics analyses on gene/protein datasets.

ABioTransPlus is designed as a one-stop server that helps the users perform all tasks through an intuitive graphical user interface (GUI) that waves the hassle of coding, installing tools, packages or libraries and dealing with operating systems compatibility and versioning issues, some of the complications that make data analysis tasks more challenging for biologists. ABioTrans Plus is an open-source tool and the website is free and open to all users and there is no login requirement.

ABioTransPlus  is available at http://combio-sifbi.org/ABioTrans/

If you find ABioTrans useful, please cite our paper:  
Zou Y, Bui TT and Selvarajoo K (2019) ABioTrans: A Biostatistical Tool for Transcriptomics Analysis. Front. Genet. 10:499. doi: 10.3389/fgene.2019.00499

Full text (free access) can be found at https://www.frontiersin.org/articles/10.3389/fgene.2019.00499/full 

## Transcriptomic analysis
### Import data and pre-processing

ABioTransPlus supports two types of common data formats for gene expression analysis: RNA-Seq count matrix and Microarray CEL files.

#### RNA-Seq gene expression matrix format

##### Import data

Choose an RNA-Seq data file in comma-separated value (.csv) format. 

* If you input raw data (read counts), please make sure that the first column contains gene names, and the read counts of each genotype (conditions: wildtype, mutants, replicates, etc.) are in the following columns. Each genotype column should have a column name.

* Along with raw read counts, you can provide gene length (base pair) information in two-column .csv file, with the first column specifying gene names, which must match the gene names in raw data file, and the second column specifying gene length in base pair. Gene length file is required for normalization methods for sequencing depth and gene length: RPKM, FPKM, TPM

* Alternatively, if you input a normalized data file, it should have gene names in rows and genotypes in columns, following the usual format of files deposited in the GEO database. 

  ![Figure 1a: Required format for raw counts data file](https://github.com/buithuytien/ABioTrans/blob/master/Test%20data/Eg_raw.png)

  ![Figure 1b: Gene length file format](https://github.com/buithuytien/ABioTrans/blob/master/Test%20data/Eg_gene_length.png) 

List of negative control genes (spike-in or stably expressed genes accross all samples), if available, should be contained in one-column .csv file. Negative control genes are required for Remove Unwated Variation (RUV) normalziation.

![Figure 1c:Negative control gene file format ](https://github.com/buithuytien/ABioTrans/blob/master/Test%20data/Eg_negative_control_genes.png) 

Finally, a metadata table matching column names of data file to experimental conditions should be given in two-column .csv format. **Metadata table is required** for differential expression analysis

![Figure 1d: Metadata file format](https://github.com/buithuytien/ABioTrans/blob/master/Test%20data/Eg_metadata.png) 

##### Pre-processing

Preprocessing involves two steps: removing lowly expressed genes and normalizing the remaining gene expression. 

- For removing lowly expressed genes, the user need to specify threshold expression values (which must be in same unit to the input data file - either raw read counts or normalized expression), and the minimum number of data columns that must exceed the threshold value. 
- Normalization methods are available upon the availability of supporting data files: normalization for sequencing depth, including TPM and RPKM, requires gene length and normalization for sample variation, including RUV, requires negative control genes. User can download the filtered, normalized data in the `Data` tab.

Relative Log Expression (RLE) plots of raw and processed data are displayed to visualize the effects of normalization. Distribution of gene expression in each data column is visualized by violin plot.

![Figure 2: Preprocessing panel with RLE plots of raw data (upper figure) and filtered, RUV-normalized data (lower figure). Gene expression with minimum five counts in at least two columns are retained.](![1b_preprocess(2).PNG](https://github.com/buithuytien/ABioTrans/blob/online-version/Screenshots/1b_preprocess(2).PNG)

#### Microarray Preprocessing (new feature)

Microarrays can be used in many types of experiments including genotyping, epigenetics, translation profiling and gene expression profiling </br></br>
Raw Count data & Metadata can be extracted from the Microarray which can further be used to perform all the different analysis.<br>
The user needs to zip all the ```CEL``` files along with the ```SDRF```(Sample and Data Relationship Format) file ( e.g. E-MTAB-2967.sdrf.txt) & upload it.</br></br>
Once the processing is complete, the user can download the output in ```csv``` format.
![alt text](https://github.com/rahul799/ABioTrans/blob/update-readme/Screenshots/Microarray.png)

### Scatter plot

Scatter plot compares any 2 samples (or 2 replicates) by displaying the respective expression of all genes in 2D space. It is recommended to preform normalization for sequencing depth (TPM, RPKM, FPKM) for this step (and so does distribution fitting, correlation, hierarchical clustering, noise and entropy).

As gene expression data is naturally skewed towards very high expression level region, we recommend applying log-transformation to capture the whole data range. User can choose between log base 2, natural log, or log base 10. An option to add linear regression line is also available. Gene expression data is densely distributed in the lowly expressed region, making the dots usually indistinguishable in regular scatter plot. 

ABioTransPlus overlay a 2D kernel density estimation on the scatter plot to visualize the density of expression level. The user can choose to download each single scatter plot, or to download all pairs of samples scatter plot in one PDF file, which may take some time to run.

![Figure 3: Scatter plot](https://github.com/buithuytien/ABioTrans/blob/online-version/Screenshots/2a_scatter.PNG) 

### Distribution fitting
Choose the column you want to fit and compare with any combinations of the six statistical distributions. You can adjust the slider or input the range for x-axis to zoom in to see the fitted curve.

Distribution fitting compares the gene expression to a number of statistical continuous distributions, which can be used to validate the data. To visualize the comparison, ABioTransPlus displays the Cumulative Distribution Function of the preprocessed gene expression data with the user-selected theoretical distributions. Once it is confirmed that the gene set follow a distribution, it would be safe conclude the validity of the gene expression data. `AIC table` is also provided in `AIC table` tab to show the best fitted distribution in each sample

![Figure 3: Distribution fitting - Cumulative Distribution Function](https://github.com/buithuytien/ABioTrans/blob/online-version/Screenshots/2b_distfit.PNG) 

### Pearson and Spearman correlations
The Pearson correlation evaluates the linear relationship between two continuous variables. A relationship is linear when a change in one variable is associated with a proportional change in the other variable.

The Spearman correlation evaluates the monotonic relationship between two continuous or ordinal variables. In a monotonic relationship, the variables tend to change together, but not necessarily at a constant rate. The Spearman correlation coefficient is based on the ranked values for each variable rather than the raw data.

User shall choose a correlation method and view it in either correlation heatmap, correlation plot in circle or correlation matrix.

![Figure 4: Correlation heatmap](https://github.com/buithuytien/ABioTrans/blob/online-version/Screenshots/2c_corr.PNG) 

### PCA and K-means clustering
Principal Components Analysis (PCA) is an multivariate statistical technique for simplifying high-dimensional data sets ([Basilevsky 1994](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2669932/#R1)). Given *m* observations on *n* variables, the goal of PCA is to reduce the dimensionality of the data matrix by finding *r* new variables, where *r* is less than *n*. Termed principal components, these *r* new variables together account for as much of the variance in the original *n* variables as possible while remaining mutually uncorrelated and orthogonal. Each principal component is a linear combination of the original variables, and so it is often possible to ascribe meaning to what the components represent. A PCA analysis of transcriptomic data consider the genes as variables, creating a set of “principal gene components” that indicate the features of genes that best explain the experimental responses they produce.

![Figure 5: Principal Component plot](https://github.com/buithuytien/ABioTrans/blob/online-version/Screenshots/2d_pca(2).PNG) 

### Differential Expression Analysis

DE analysis identifies the genes that are statistically different in expression levels between the 2 selected conditions. Two important threshold are:

- The lower bound of expression **fold change** between the 2 selected conditions

- The upper bound of hypothesis test p-value

AbioTransPlus implements 3 popular methods to identify DE genes:

- [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

- [EdgeR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/)

- [NOISeq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4666377/)

For data with single replicate in all experiment condition, NOISeq method can simulate technical replicates to carry out DE test. Metadata file is required for DE Analysis. Please make sure metadata contains all column names from input data file and match them with experimental condition 

For edgeR and DESeq2, raw read counts data file must be provided. For NOISeq, gene expression should be normalized for sequencing depths (by select normalization method in Preprocessing tab if raw counts file is inputted, or by directly providing normalized gene expression) 

To carry out the analysis, first the user needs to specify DE methods, two conditions to compare (condition 2 is compared against condition 1), fold change cut-off value and False Discovery Rate (FDR or adjusted p-value) threshold, then hit the “Submit” button. By convention, DE genes are filtered at FDR = 0.05 and 2-fold change. 

When the computation finishes, table of DE genes, volcano plot of DE result and dispersion plot of input data are displayed in their respective tabs. Please note that volcano plot and dispersion plot are only available for edgeR and DESeq2 methods.

![Figure 6: Volcano plot summarizing differential expression analysis result](https://github.com/buithuytien/ABioTrans/blob/online-version/Screenshots/2e_volcano.PNG)

### Heatmap of gene expression and Hierarchical clustering
Hierarchical clustering is used to find the groups of co-expressed genes. The clustering is performed on normalized expressions of differentially expressed genes using Ward clustering method.

Heat map of gene expression and hierarchical clustering can be carried out on DE genes resulted from DE analysis tab by selecting `DE result`. Otherwise, you can specify the minimum fold change and minimum number of samples (columns in input data file) satisfying the fold change. Finally, you need to specify the number of clusters on rows (genes), then hit `Submit`.

The gene names in each cluster are displayed in the `Gene clusters` panel, corresponding to the heatmap you just generated.

![Figure 7: Heatmap of gene expression and hierarchical clustering for co-expressed genes](https://github.com/buithuytien/ABioTrans/blob/online-version/Screenshots/2f_headmap_cluster(2).PNG)

### Transcriptome-wide average noise

Transcriptome-wide average noise is used to quantify between gene expressions scatter of all replicates in one experimental condition. Formula is adapted from [Kumar et. al](https://www.nature.com/articles/srep07137) 's paper

* User to select a desired noise plot according to your data. By default the name of the first replicate in each genotype is used as the name of the genotype. You can specify the names of the genotypes, only make sure that all names are different (due to the mapping mechanism of Plotly). Please be reminded that the consecutive columns will be regarded as of the same genotype.
* Noise here refers to the average transcriptome noise — the squared coefficient of variation — defined as the variance (σ2) of expression divided by the square mean expression (μ2). 
* If `replicates` is selected, the noise within one genotype will be computed. 
* If `genotypes (average of replicates)` is selected, the gene expression data within one genotype will first be taken average of, and the average value will be used as the expression value of that certain genotype. Then noise will be computed between the anchor genotype and all the other genotypes. 
* If `genotypes (no replicate)` is selected, each sample column is treated as a different genotype. Noise will directly be computed between the anchor genotype and the rest of the genotypes. 

### Entropy

Shannon entropy ([Shannon, 1948](https://onlinelibrary.wiley.com/doi/10.1002/j.1538-7305.1948.tb01338.x)) measures the disorder of a high-dimensional system, where higher values indicate increasing disorder. ABioTransPlus adapts the Shannon entropy formula to quantify the variability of one sample of gene expression.

* First specify whether your data is a time series data.
* If the data is not in time series, entropy values of all samples will directly be displayed. 
* Otherwise, if it is a time series data, specify the number of time points. For example, if the number of time points is 6, then sample 1 to 6 will be regarded as genotype A time 1 - time 6, sample 7 to 12 will be genotype B time 1- time 6, and so on and so forth. For a line chart, the entropy of each genotype will be shown in one line. 

![Figure 8: Entropy](https://github.com/buithuytien/ABioTrans/blob/online-version/Screenshots/2g_entropy.PNG)

### t-distributed stochastic neighbor embedding (t-SNE)  (new feature)
t-SNE is a dimensionality-reduction approach that reduces the complexity of highly complex data such as the transcriptomic data. It visualizes the sample interrelations in a 2- or 3-dimensional visualization. This allows the identification of the close similarities between samples through the relative location of mapped points. Since t-SNE is nonlinear and able to control the trade-off between local and global relationships among points, its visualization of the clusters is usually more compelling when compared with the other methods [[Cieslak 2020](https://pubmed.ncbi.nlm.nih.gov/31784353/)]. 

ABioTrans Plus introduces an intuitive interface that allows performing t-SNE analysis on the processed untransformed transcriptomic data through entering three inputs perplexity value, the number of principal components (PC) and the number of clusters. The user can also choose to log transform the data before submission. 

### Clustering with random forest (new feature)
Random forest clustering is an unsupervised learning approach, where each sample is clustered into different classes, based on their similarity (usually based on Euclidean distance). 

ABioTrans Plus uses the random forest algorithm to generate a proximity matrix, a rough estimate of the distance between samples based on the proportion of times the samples end up in the same leaf node of the decision tree. The proximity matrix is converted to a dist matrix which is then inputted to the hierarchical clustering algorithm [[Chen 2012](https://pubmed.ncbi.nlm.nih.gov/22546560/)].

The ```diff_genes.csv``` file will be generated after each run, containing names of differential genes.
![alt text](https://github.com/SnowMelody/ABioTrans/blob/master/ABT_updated/screenshots/Random%20forest.png)

### Self-Organizing Map (SOM) (new feature)
SOM is a dimensionality reduction technique that produces a two-dimensional, discretized representation of the high-dimensional gene expression matrix. It uses a neighbourhood function to preserve the topological properties of the input gene expression matrix. Each data point (e.g. 1 sample) in the input gene expression matrix recognizes itself by competing for representation. SOM mapping steps start from initializing the weight vectors. From there, a sample vector is selected randomly and the map of weight vectors is searched to find which weight best represents that sample. Each weight vector has neighbouring weights that are close to it. The weight that is chosen is rewarded by being able to become more like that randomly selected sample vector. The neighbours of that weight are also rewarded by being able to become more like the chosen sample vector. This allows the map to grow and form different shapes. Most generally, they form square/rectangular/hexagonal/L shapes in 2D feature space [[Yin 2020](https://link.springer.com/chapter/10.1007%2F978-3-540-78293-3_17)].

ABioTrans Plus provides 5 types of SOM plots: </br>
```Property```: Uses values of codebook vectors (weight of gene vectors) and output as coloured nodes </br>
```Count```: Shows how many genes are mapped to each node </br>
```Codes```: Displays codebook vectors </br>
```Distance```: Shows how close genes are from each other when they are mapped </br>
```Cluster```: Uses hierarchical clustering to cluster the SOM </br>

There are 4 parameters to control the SOM plots. ```Samples used``` determines the sample chosen for the plots, either all samples or individual ones. ```No. of horizontal grids``` and ```No. of vertical grids``` changes the number of nodes used in the SOM. ```No. of clusters```  classifies the SOM nodes into the specified cluster size (for cluster plot).
![alt text](https://github.com/SnowMelody/ABioTrans/blob/master/ABT_updated/screenshots/SOM.png)

## Gene/Protein Set Analysis

Differentially expressed gene (DGE) analysis usually outputs a list of genes that are statistically determined as differentially expressed. Then, the list of DEGs is analyzed, interpreting and annotated to learn more about the functions, pathways and cellular processes that these genes are involved in. Most of the gene differential expression analysis tools do not include bioinformatics features for gene set analysis or include few basic analyses such as GO and pathways enrichment. 

In ABioTrans Plus, we designed the GO feature to be dynamic by reading the GO terms associated with the genes/proteins directly from UniProt Knowledgebase [[UniProt Consortium 2019](https://pubmed.ncbi.nlm.nih.gov/30395287/)] then visualize each of the three GO domains (cellular component, molecular function and biological process) in an independent tab in a bar chart and on a downloadable tabular format (Figure YY). Furthermore, we introduced 11 new bioinformatics analyses that can be performed on a given gene/protein dataset. 

### Pathways Enrichment Analysis

For a given gene or protein set, ABioTrans Plus uses g:Profiler [[Raudvere 2019](https://doi.org/10.1093/nar/gkz369)] to perform a pathway enrichment analysis and displays the results as a network where the nodes are the pathways and the edges are the overlap between the pathways (Figure X). We use Cytoscape.JS for the network visualization [[Franz 2016](https://academic.oup.com/bioinformatics/article/32/2/309/1744007)] and, therefore, the network properties such as the colour and layout can be changed from the left panel. 

The users can choose from nine different network layouts and five different network colour schemes. The overlap between the pathways can be also changed from the provided controllers so that the user can choose an overlap cutoff or overlap range. The enrichment results can also be downloaded as a CSV file.      

### Protein-Protein Interaction

Investigating PPIs is one of the essential steps in systems biology studies. PPI databases are growing in size with improved accuracy of PPI curations. Furthermore, several recent large-scale projects provide experimentally-confirmed interaction such as the reference map of the human binary protein interactome (HuRI) project [[Luck 2020](https://www.nature.com/articles/s41586-020-2188-x)]. The PPI feature in ABioTrans Plus provides the users with an interface where they can upload a set of proteins (UniProt accessions) and get all the interactions associated with them. The interactions are visualized as a network where the nodes are the proteins and the edges are the interactions and the node size corresponds to the number of interactors of the protein. We use Cytoscape.JS for PPI visualization [[Franz 2016](https://academic.oup.com/bioinformatics/article/32/2/309/1744007)]. We provide users with five network styles and nine network layouts to customize their results. The results are also displayed as an interactions table and can be downloaded as a network or an interaction table. 

### Complex Enrichment 

The identification of the subunits of the protein complexes is important to understand the functions and the formation of these macromolecular machines. We provide the user with a complex enrichment feature that allows the identification of proteins in the provided dataset that are part of a known protein complex. This feature uses CORUM databases, which contain curated complex information for mammalian proteins [[Giurgiu 2019](https://pubmed.ncbi.nlm.nih.gov/30357367/)]. This feature provides the user with complex-forming proteins and the complex information in the submitted dataset.   

### Protein Function

UniProt Knowledgebase provides a detailed function for thousands of proteins. The protein function feature retrieves the protein function information from UniProt of a given protein set (UniProt accessions). The retrieved protein functions are displayed in a downloadable tabular format. 

### Protein Subcellular Localization

The protein localization critically affects the protein function. The change of the protein localization is a cellular approach to regulate biological processes. UniProt Knowledgebase provides curated subcellular localization information for the proteins. The protein subcellular localization feature provides the user with an interface to get the subcellular localization information for a given list of proteins (UniProt accessions) and display the results in a downloadable tabular format.

### Protein Domains

The protein domains are the functional subunits of the proteins that contribute to their overall function. Each domain is responsible for a certain function or interaction. ABioTrans Plus provides the users with a protein domain feature that connects to UniProt Knowledgebase and retrieves the domain information associated with each protein in a given list of UniProt accession. 

### Tissue Expression

The distinct expression profile of genes and proteins per tissue is what gives the different tissues the suitability for their functions. The tissue expression feature in ABioTrans Plus provides the user with the tissue expression for each protein in a given protein list (UniProt accessions) through retrieving this information from UniProt Knowledgebase. The result is displayed in a downloadable tabular format. 8) Gene Co-expression: The co-expression analysis is a common analysis that assesses the expression level of different genes to identify simultaneously expressed genes. The resultant co-expression networks are used to identify functionally related genes or genes being controlled by the same transcriptional mechanism [[Vella 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5359264/)]. ABioTrans Plus provides the users with an interface where they can submit a co-expression query to GeneMANIA [[Franz 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6030815/)] then shows the results at GeneMANIA’s website in a new tab. Currently, we support queries for nine model organisms including human, yeast, E. coli, C. elegans, Arabidopsis, Drosophila, zebrafish, mouse and rat.9) Protein Physicochemical Properties: For a given set of proteins (UniProt accessions), this feature provides the user with the complete sequences of them in a single FASTA file and allows the user to investigate their physicochemical properties. The physicochemical analysis includes sequence charge, GRAVY index [[Kyte 1982](https://pubmed.ncbi.nlm.nih.gov/7108955/)] and hydrophobicity (Figure YA).

### Protein Evolutionary Analysis

For a given set of proteins (UniProt accessions), this feature provides the user with a phylogenetic and evolutionary analysis that include multiple sequence alignment (MSA) of the protein sequences, clustering based on the amino acid sequences, chromosomal location or gene tree (Figure YB).     

### Protein Pathological Analysis

Several diseases are associated with the malfunction of certain genes or proteins. The disease-protein association is collected in different online resources such as OMIM databases [[Amberger 2019](https://pubmed.ncbi.nlm.nih.gov/30445645/)], DisProt [[Hatos 2020](https://academic.oup.com/nar/article/48/D1/D269/5622715)] and DisGeNET [[Pinero 2019](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1021/5611674)]. ABioTrans Plus provides the users with an interface that retrieves the disease-protein association from online databases for a given list of proteins (UniProt accessions). The disease-protein association is visualized as bubble chats that shows the distribution of the proteins among the disease (each bubble is protein and the bubble size is the number of associated disease) or the distribution of proteins among the diseases (each bubble is a disease and the bubble size is the number of associated proteins) (Figure YC).    

The features that communicate with UniProt Knowledgebase use UniProtR, an R package for data retrieval and visualization from UniProt [[Soudy 2020](https://pubmed.ncbi.nlm.nih.gov/31843688/)]. Since all the bioinformatics features only accept gene names (gene symbol) or UniProt accessions, we provide the users on each page with links to two ID converters UniProt ID mapping [[UniProt Consortium 2019](https://pubmed.ncbi.nlm.nih.gov/30395287/)] and g:Convert [[Raudvere 2019](https://doi.org/10.1093/nar/gkz369)] to convert their identifiers to gene names or UniProt accessions. 



## Download instructions
1. Scatter plot, distribution fit, correlation plot and heatmap can be saved as PDF by clicking `Download as PDF`. You can name the file before saving it. Also, you can directly drag the plot from the GUI to a folder on the computer.

2. PCA, noise, entropy and gene ontology pie chart are supported by Plotly. If you are using the GUI in R window or Safari, click `Download as PNG` and it may take 5 seconds to save the graph. You cannot name the file yourself, so please be aware of overwriting issues. Even if the plot is saved successfully, there may be a 404 Not found message after downloading. Please go to your working directory and check. If you cannot save the graph successfully, please try the following code the Rstudio command line. You only need to run this code once. (You can also simply take a screenshot of the graph.)
```R
webshot::install_phantomjs()
```
If you are using R window or Safari, and you cannot save PCA-3D graph successfully, please create a plotly account (free). Please go to this website on the instructions of creating a new account and Find your authentication API keys in your online settings (https://plot.ly/r/getting-started/). After that please go to the first two lines in ABioTrans.R file, remove both of the "#" signs, and replace the username and api key with your own. Save the file before you run. However, please be reminded that this  `Download as PNG` button only downloads the PCA-3D plot from the default angle.

3. For Chrome users, the `Download as PNG` button may be disabled. Please put your mouse over the graph region and click the left most camera sign to save the plot. It will be saved in the download folder on your computer. You are suggested to use Chrome because it saves you trouble from creating an account and changing the code. More importantly, you can download the PCA-3D plot from any angle.

## Docker Support
Now we also have a docker image containing all the dependencies needed to run ABioTrans Plus. It also support shiny-server opensource which can be exposed in port 3838.

It can be pulled using 
```docker pull jaktab/genecloudomics-webserver:latest```
