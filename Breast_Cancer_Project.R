
----------######################################################----------
          ######         Analysis of RNA-Seq Data         ######         
          ######        Building a Simple Classifier      ######  
          ######    Data Source: The Cancer Genome Atlas  ######
          ######         Author: Anusha Bailoor           ######
          ######         Date: 21st February 2024         ######
----------######################################################----------

                ##### Part 1: Introduction #####
# Cancer is a disease that plagues millions of people around the world. It is a widely studied topic and yet devoid of any cures whatsoever. For the year 2020 alone, Cancer accounted for more than 10 million deaths worldwide (Sharma et al., 2010). Breast Cancer is one of the leading causes of Cancer in women, especially in the US, with a predicted estimate that almost 12% of the female population in the US will be diagnosed with it over the course of their lifetimes (Waks & Winer, 2019).    
# RNASeq Data is a data mine of information which could prove to be beneficial when looking at bioinformatics analysis, specifically or any biological processes for that matter. It enables us to detect mutations in single nucleotides, does not require a strong hold of genomic sequence knowledge and the most important benefit that we are going to exploit further in this script is that it gives information about quantitative expression levels (Castillo et al., 2017; Costa-Silva et al., 2017).  
# This research topic interests me, not only because it is a disease that hits close to home, with more than a few of my family members either being plagued by or having had lost their battle against it, but also because from a scientific viewpoint, even after tremendous research efforts having been made in the field, the complete eradication of this disease is still a distant dream. The prospect that working on this project will help me understand the mechanisms that govern the progression of cancer, while strengthening my base of conceptual understanding of bioinformatics is the reason why I want to work on a project that involves data wrangling of RNASeq Data collected specifically from cancer patients. Here, my main aim is to use raw RNASeq data, normalize it, use this clean data to build a simple classifier based on differentially expressed genes in Breast Cancer and make hierarchical clusters of the top differentially-expressed genes.
# The data that will be used in this script is sourced from The Cancer Genome Atlas (TCGA) (Tomczak et al., 2015) (URL: https://www.cancer.gov/ccg/research/genome-sequencing/tcga). We will be looking at the Breast Cancer dataset that is identified by TCGA as BRCA. The TCGA-BRCA dataset provides information for about 1095 fully accessible patients and includes a variety of information ranging from clinical, expression, DNA Methylation, proteomic and genotypic data. For the scope of this script, we will be focusing on the Transcriptome Profiling data that contains the RNASeq data (Berger et al., 2018). 
# In order to do complete justice to my own personal motivations and the quality of my work, the script below is an adaptation of the work carried out by the Costa Lab, specifically, the PhD students, Fabio Ticoni and Tiago Maie (URL: https://www.costalab.org/wp-content/uploads/2021/11/handout_day41.html) and is inspired by the work by Dr Shraddha Pai in the Cancer Analysis Workshop by OICR in 2021 (URL: https://pailab.oicr.on.ca/CBW_CAN_DataIntegration_2021/module-11-lab-1-three-way-classifier-from-four-genomic-layers.html#get-and-prepare-data) 

                 ##### Part 2: Package Loading #####
# Let us load all the necessary packages that will be required to undertake an analysis of this sort

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("edgeR")
# BiocManager::install("limma")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(limma)
# install.packages("tidyverse")
library(tidyverse)
# install.packages("glmnet")
library(glmnet)
# install.packages("gplots")
library(gplots)
# install.packages("caret")
library(caret)
# install.packages("RColorBrewer")
library(RColorBrewer)
 install.packages("gprofiler2")
library(gprofiler2)

                ##### Part 3: Data Review #####
# In this section, we will look at the data that is available from TGCA and  will review it

GDCprojects = getGDCprojects()

head(GDCprojects[c("project_id", "name")])

TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

query_TCGA = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts"
)

df_BRCA = getResults(query_TCGA)
View(df_BRCA)
colnames(df_BRCA)

table(df_BRCA$sample_type)

#Now to keep the global environment clean let us get rid of query_TCGA and df_BRCA.
rm(query_TCGA)
rm(df_BRCA)

                ##### Part 4: Data Retrieval #####
# Now that we have a preliminary understanding of the data, let us download it. As per the tabulated details of sample types, only 14 samples are metastatic, in order to simplify the analysis, let us avoid involving those samples into our study. Also, let us download these files on to our device so that the analysis is carried out smoothly. 

query_TCGA_brca = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query_TCGA_brca)

BRCA_data <- GDCprepare(query_TCGA_brca)

dim(BRCA_data)

# The SummarizedExperiment package from Bioconductor helps wrangle some data, particularly when looking at Experimental data. Let's use the colData function from the package to look at all the columns within this Large Ranged Summarized Experiment matrix that has just been created called BRCA_data.

colnames(colData(BRCA_data))


          ##### Part 5: Preliminary Data Exploration ####
# Now that we have a lot of information from TCGA regarding Breast Cancer data, let us explore some of the data from the Large Ranged Summarized Experiment data.

table(BRCA_data@colData$tumor_descriptor)
table(BRCA_data@colData$tissue_or_organ_of_origin)
table(BRCA_data@colData$age_at_index)
table(BRCA_data@colData$sample_type)
table(BRCA_data@colData$gender)
table(BRCA_data@colData$race)
table(BRCA_data@colData$ethnicity)

# From the above data exploration step, I would like to visualize the data as it would give us a better understanding. 

race_data <- table(BRCA_data@colData$race)

race_data_df <- as.data.frame(race_data)

ggplot(race_data_df, aes(x = "Var1", y = "Freq", fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Distribution of Breast Cancer Cases by Race",
       fill = "Race") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Now, let us visualize the distribution of the cancer cases by their gender

gender_data <- table(BRCA_data@colData$gender)
gender_data_df <- as.data.frame(gender_data)

ggplot(data = gender_data_df, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(breaks = seq(0, 30, by = 5)) +
  geom_text(aes(label = Freq),
            stat = "identity",
            vjust = -0.5,
            size = 3) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(title = "Distribution of Breast Cancer Cases by Gender", x = "Gender", y = "Case Count", fill= "Gender")

     ##### Part 6: Pan-optic Data Exploration and Filtering #####
# Now with some initial data exploration and visualization out of the way, let us get into the data filtering and normalization of the RNASeq pipeline that will be used to build the classifier. The normalization step is crucial as it helps us to perform differential expression (DE) analysis without the batch effects that is often associated with raw data. The pipeline will consist of a function that encapsulates all the necessary functions into one workflow. The R package limma will be used for this.Also, the voom method or the mean-variance modelling at the observational level is used, this method transforms RNASeq Data to log2 counts/million (logCPM), estimates the mean variance relationship and, uses this information to compute appropriate observational-level weights. 
# Here, the function limma_RNA_workflow, uses the package limma to clean our data using feature selection. In this function three arguments have to inputted by the user, first, the data they want to normalize, secondly, the feature they want to select and lastly the reference group that they want to look at. Next the design matrix must be created, which will take in the conditions that the user inputs, the ~ symbol on line 147 indicates that a formula is being constructed. Now, to remove the genes that have low levels of differential expression, a DGEList is created, here the default of removal of genes with less than 10 reads is maintained. Code on line 157 normalizes the data to minimize batch effects and technical variation using the Trimmed-mean of M-values (TMM) normalization method. The plot=TRUE code on line 158 ensures the generation of a Mean-variance trend everytime this function is run. But, a plot without any kind of fitting will be incomplete, so to fit the plot, lmfit is used, along with eBayes to produce an object that holds a myriad of statistics that we can leverage to rank the differentially expressed genes. The topTable function helps us get a table of the top 100 genes, as we have set the parameter to 100. Finally, the function returns three variables, the voomObj that generates the plot, fit that contains a myriad of statistics based on Bayesian concepts and a table, topGene containing the top 100 differentially expressed genes.

limma_RNA_workflow = function(
    BRCA_data,
    condition_variable,
    reference_group=NULL){
  
  design_factor = colData(BRCA_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(BRCA_data),
                samples=colData(BRCA_data),
                genes=as.data.frame(rowData(BRCA_data)))
  
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  topGenes = topTable(fit, coef=ncol(design), number=100, sort.by="p")
  dim(design)
  return(
    list(
      voomObj=v, 
      fit=fit, 
      topGenes=topGenes
    )
  )
}

# Now, let us call the function and look the top 100 differentially expressed genes comparing two different types of cancer definitions, either Primary Tumor or Solid Tissue Normal and since Solid Tissue Normal are the healthy tissues, we will use that as our reference group.  

limma_BRCA_res = limma_RNA_workflow(
  BRCA_data=BRCA_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)

# Now to look at the Principal Component Analysis Plot for the voomObj. The first step to do this is to extract the data and the condition variable, followed by the most important step of perfroming the PCA. Let us also add in colors in order to make the PCA plot pop, to do that we can make a custom palette with the desired colors. The next step is to plot a graph with the PCA1 and PCA2 with our desired color palette and ensuring we have a legend that tells us which color stands for which cancer definition type. Finally we save the resultant PCA. 

group <- factor(limma_BRCA_res$voomObj$targets[, "definition"])
expression_data <- t(limma_BRCA_res$voomObj$E)

pca <- prcomp(expression_data)

plot_colors <- c("blue", "purple")
palette(plot_colors)

plot(pca$x[, 1:2], col = as.numeric (group), pch = 19)

legend("bottomright", inset = 0.01, levels(group), pch = 19, col = 1:length(levels(group)))

res_pca <- pca

# Additionally we can look at making a hierarchical cluster or a heatmap of the differentially expressed genes that we collected using the RNASeq Normalization function limma_RNA_workflow. To do this, we must first make a normalized expression matrix from our limma_BRCA_res object, followed by getting a list of gene names and renaming the genes within our expression matrix. Finally, we want to look at the topGenes that are already sorted by adjusted p-value in our list limma_BRCA_res, we could use head function to get the top 20 of them. Using only the top 20 will help us understand the clustering better also, let's use sample type within the TCGA BRCA to group the data within the heatmap

BRCA_em = as.matrix(t(limma_BRCA_res$voomObj$E))

gene_names = limma_BRCA_res$voomObj$genes[,"gene_name"]

colnames(BRCA_em) = gene_names

head(limma_BRCA_res$topGenes,20)

top_genes = limma_BRCA_res$topGenes$gene_name[1:20]

sample_type = factor(limma_BRCA_res$voomObj$targets$sample_type)

# Now, to assign a few parameters, the color is set using the brewer.pal function from the RColorBrewer package. The hclust function is a part of the stats package and the function, clust_func performs complete linkage clustering whereas the function, dist_func uses the inverse of correlation as distance.

hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)

clust_func = function(x) hclust(x, method="complete")

dist_func = function(x) as.dist((1-cor(t(x)))/2)

# Next, let us construct the a heatmap using the heatmap.2 function from gplots. The visualization of a heatmap is based on the tinkering a variety of parameters like using the hmcol character to define the color map, setting the dendogram to show both the axis, calling on the functions clust_func and dist_func to define the hierarchical clustering method and using the correlation coefficient for distance respectively, to touch upon a few important parameters. The rest of the parameters pertain to the aesthetic of the heatmap.

gene_heatmap = heatmap.2(
  t(BRCA_em[,top_genes]),
  scale="row",          
  density.info="none",  
  trace="none",         
  col=hmcol,            
  labCol=FALSE,         
  ColSideColors=as.character(as.numeric(sample_type)), 
  dendrogram="both",  
  hclust = clust_func,  
  distfun = dist_func,  
  cexRow=1.0,           
  keysize = 1.25,       
  margins=c(1,8)        
)

            ##### Part 7: Building the Classifier #####
# Now that we have visualized the TCGA BRCA data a bit, let us get on with building a classifier from the same. This classifier when functioning correctly must be able to classify a sample as a tumor or not. We can do this by first building a simple linear model and then an Elastic Net Model. The absolute first step to any machine learning model is to generate a training and testing set that will help us train and validate our data. Let's start by extracting the data we need, since we want to classify on the basis of sample type, that is, whether it is a Primary Tumor or a Normal Tissue, let us look at the normalized expression values from the limma_BRCA_res$voomObj$E and need to specify 'definition'column as the factor.

class_mat <- as.matrix(t(limma_BRCA_res$voomObj$E))
class_fact <- as.factor(limma_BRCA_res$voomObj$targets$definition)

# Now to build the training and testing model. Here I will be using the function createDataPartition which is a part of the caret package to do this. It typically partitions the data into 75-25, wherein, 75% is for the training model and 25% is for the testing model.

set.seed(80)
train_mod = createDataPartition(class_fact, p=0.75, list=FALSE)

train_1 = class_mat[train_mod, ]
test_1  = class_mat[-train_mod, ]

train_2 = class_fact [train_mod]
test_2  = class_fact [-train_mod]

#Now to train an Elastic Net Model. Typically, the Elastic Net Model combines the two other models, namely LASSO and Ridge Regression. Here, the parameter alpha can determine whether the model skews more towards a LASSO regression (by setting alpha to 1) or towards Ridge (by setting alpha to 0). For now, let us set it to 0.5. 

ENM_res = cv.glmnet(
  x = train_1,
  y = train_2,
  alpha = 0.5,
  family = "binomial"
)

       ##### Part 8: Machine Learning Model Validation #####
# The models have been generated but we must validate the same so that we can check its efficacy. In order to do that, we can use our test data. 

pred = predict(ENM_res, newx=test_1, type="class", s="lambda.min")

# One method of evaluating our model is to use a confusion matrix that returns a simple table with all the predictions v/s the true values. This allows us to look at the true positives, true negatives, false positives and false negatives. 

confusion_mat <- table(pred, test_2)
print(confusion_mat)

# Now let us check the some of the other parameters that are important when looking at a machine learning model, like the sensitivity, specificity and precision. All these parameters are part of the caret package. 

print(paste0("Sensitivity: ",sensitivity(confusion_mat)))
print(paste0("Specificity: ",specificity(confusion_mat)))
print(paste0("Precision: ",precision(confusion_mat)))

# This means that the model seems to be working very well, but we are also interested in looking at the differentially expressed genes that help the model make these accurate classifications. This type of information can also be looked at by taking a look at the coefficients (in this case, genes) that the Elastic Net Model selected to make the predictions.

ENM_pred_gene <- coef(ENM_res, s = "lambda.min")
dim(ENM_pred_gene)
head(ENM_pred_gene)

# We can see a lot of genes pop on in this, but we are most interested in the non-zero ones, so let's get some information regarding the non-zero ones.

ENM_pred_gene = ENM_pred_gene[ENM_pred_gene[,1] != 0,]
head(ENM_pred_gene)

# Now, let's try to get a count of all th relevant genes that have contributed in this model and are non-zero. Also, we should get rid of the first gene or the Intercept, as that is a variable of the model itself
ENM_pred_gene = ENM_pred_gene[-1]

relevant_genes = names(ENM_pred_gene)
length(relevant_genes)

# Okay, now that we have the counts, let us look at their names as well

head(relevant_genes)

# The Ensembl annotation on this is not that understandable, perhaps looking at the common gene name would be more helpful. We can do that by using common gene names from limma_BRCA_res$voomObj$genes. Lets assign it to a new data frame altogether.

common_gene_names <- limma_BRCA_res$voomObj$genes[relevant_genes,"gene_name"]
head(common_gene_names)

# A nice way to understand the changes would be to generate a hierarchical cluster of the relevant genes from the Elastic Net Model. But we have also looked at a similar cluster for our limma model, so we could also try to see whether the same genes have been relevant for both models as well. The parameters for the heat map generation remain the same, just the color palette is slightly adjusted. Also the purple shows the genes that are common between the two heatmaps while the white shows the genes that are unique to the Elastic Net Model.

hmcol = colorRampPalette(rev(brewer.pal(4, "RdYlBu")))(256)


colorLimmaGenes = ifelse(
  (relevant_genes %in% limma_BRCA_res$topGenes$gene_id),
  "purple",
  "white"
)

gene_heatmap_2 = heatmap.2(
  t(class_mat[,relevant_genes]),
  scale="row",          
  density.info="none",  
  trace="none",         
  col=hmcol,            
  labRow=common_gene_names, 
  RowSideColors=colorLimmaGenes,
  labCol=FALSE,         
  ColSideColors=as.character(as.numeric(class_fact)), 
  dendrogram="both",   
  hclust = clust_func,       
  distfun = dist_func, 
  keysize = 1.25,
  cexRow=.5,            
  margins=c(1,5)        
)

hc = as.hclust(gene_heatmap$rowDendrogram)

# Let us run some preliminary functional enrichment analysis using gprofiler2 package on R. To do this first, cut the tree into 2 groups, up-regulated in tumor and up-regulated in control
clusters <- cutree(hc, k = 2)
table(clusters)

gprofiler_cols <- c("significant", "p_value","term_size", "intersection_size", "term_id", "term_name")

#For Group 1, up regulated in tumor tissue
query_genes <- names(clusters[clusters %in% 1])
result <- gost(query_genes, organism = "hsapiens")$result

# Check the column names in the result
colnames(result)

# Select the desired columns and to account for any discrepancies between the column names and our desired columns
if (all(gprofiler_cols %in% colnames(result))) {
  selected_result <- result[, gprofiler_cols]
  selected_result
} else {
  missing_cols <- gprofiler_cols[!gprofiler_cols %in% colnames(result)]
  cat("The following columns are missing in the result:", missing_cols, "\n")
}

# Now, redo this for Group 2, to find the up regulated genes in normal tissue
query_genes_normal <- names(clusters[clusters %in% 2])
result_normal <- gost(query_genes_normal, organism = "hsapiens")$result

if (all(gprofiler_cols %in% colnames(result_normal))) {
  selected_result_normal <- result_normal[, gprofiler_cols]
  selected_result_normal
} else {
  missing_cols <- gprofiler_cols[!gprofiler_cols %in% colnames(result_normal)]
  cat("The following columns are missing in the result:", missing_cols, "\n")
}

# Here, again it is clearly visible that there is a clear demarcation between the normal tissue samples and the breast cancer tissue samples. Also, it is noteworthy that within the cancer tissue fewer up regulated genes were identified when compared to the normal tissue. 

                   ##### Part 9: Discussion #####
# The primary goal of this project was to use open-source Cancer RNA-Seq data, find differentially expressed genes in it and use this data as a factor to classify possibly unseen data as Tumor tissue type or normal based on the training models. In this regard, if we look at the PCA plot, it shows a clear demarcation in terms of RNA expression profile between the two sample groups, that is, the Primary Tumor Tissue and the Solid Normal Tissue. The first heatmap generated by the limma workflow shows the top 20 genes within the samples and the ones in red are the Normal tissues and the ones in black are the Primary Tumor Tissue. Now to move on to the machine learning model itself, the   accuracy, that is the specificity and sensitivity of the model show that the model works well and will be able to classify unseen data effectively as well. 
# Now, I also realize the various biases that this analysis might have. In the context of sample size, I believe that the sample size was decent to undertake a study of this scale but it must be noted that some of the data within the open-source data portal is hidden and unavailable without prior authority. Another noteworthy point of TCGA portal is that not all types of cancers data available on the portal have the same treatment, that is, not all types of cancers have a wide variety of information available and/or have a small sample size perhaps due to the aggressive nature of the cancer or the low occurrence of that particular type of cancer(5-6). Comparing both the heatmaps, I can infer that the algorithms used by both limma model and the Elastic Net Model are completely different,the inference is drawn by the fact that top genes within the two models are different, in fact the top genes of LIMMA seem to be somewhere towards the end in the Elastic Net Model heatmap which begs the question which one is better? In my opinion, the heatmap generated by the LIMMA is better as the Elastic Net Model contains certain biases which leads to the low or almost no overlap between the two heatmaps and leads to a more streamlined demarcation between the Primary Tumor Tissue (highlighted by black at the top) and the Solid Normal Tissue (highlighted by red at the top) in the heatmap generated using the Elastic Net Model (7-9). 
# Now to touch upon the next steps in this research, I feel the scope is immense as there is a lot of information that is provided by The Cancer Genome Atlas for Breast Cancer. If I had more time to work on this, I would be interested in a few major types of analysis, the first being Survival Analysis which as the name suggests, enables us to understand patient outcomes and survival chances using Kaplan-Meier Plots and R packages like survival and survminer. Secondly, I would like to work on Proteomics-based Breast Cancer research on RStudio using BioConductor vignettes like RforProteomics. Finally, I would also like to work with the gprofiler2 package to generate even more visualizations like Manhattan Plots and perform Gene Ontology Analysis on the differentially expressed genes found in this analysis. (All of these packages are referenced in the Reference section under Package references) 


                 ##### Part 10: References #####
# Paper References:
# Berger, A. C., Korkut, A., Kanchi, R. S., Hegde, A. M., Lenoir, W., Liu, W., Liu, Y., Fan, H., Shen, H., Ravikumar, V., Rao, A., Schultz, A., Li, X., Sumazin, P., Williams, C., Mestdagh, P., Gunaratne, P. H., Yau, C., Bowlby, R., … Mariamidze, A. (2018). A Comprehensive Pan-Cancer Molecular Study of Gynecologic and Breast Cancers. Cancer Cell, 33(4), 690-705.e9. https://doi.org/10.1016/J.CCELL.2018.03.014/ATTACHMENT/7388D231-E5C1-484B-925B-2D8D1F899297/MMC5.XLSX
#Castillo, D., Gálvez, J. M., Herrera, L. J., Román, B. S., Rojas, F., & Rojas, I. (2017). Integration of RNA-Seq data with heterogeneous microarray data for breast cancer profiling. BMC Bioinformatics, 18(1), 1–15. https://doi.org/10.1186/S12859-017-1925-0/FIGURES/11
# Costa-Silva, J., Domingues, D., & Lopes, F. M. (2017). RNA-Seq differential expression analysis: An extended review and a software tool. PLoS ONE, 12(12). https://doi.org/10.1371/JOURNAL.PONE.0190152
# Das, J., Gayvert, K. M., Bunea, F., Wegkamp, M. H., & Yu, H. (2015). ENCAPP: Elastic-net-based prognosis prediction and biomarker discovery for human cancers. BMC Genomics, 16(1), 1–13. https://doi.org/10.1186/S12864-015-1465-9/FIGURES/6
# Han, H. (2015). Diagnostic biases in translational bioinformatics. BMC Medical Genomics, 8(1), 1–17. https://doi.org/10.1186/S12920-015-0116-Y/TABLES/6
# Lee, J. S., Cho, Y., Lee, J., Ko, G. W., & Yu, D. (2022). A study on bias effect of elastic net penalized regression for model selection†. The Korean Data & Information Science Society, 33(1), 67–86. https://doi.org/10.7465/JKDI.2022.33.1.67
# Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. https://doi.org/10.1093/NAR/GKV007
# Sharma, G. N., Dave, R., Sanadya, J., Sharma, P., & Sharma, K. K. (2010). VARIOUS TYPES AND MANAGEMENT OF BREAST CANCER: AN OVERVIEW. Journal of Advanced Pharmaceutical Technology & Research, 1(2), 109. /pmc/articles/PMC3255438/
# Tomczak, K., Czerwińska, P., & Wiznerowicz, M. (2015). The Cancer Genome Atlas (TCGA): an immeasurable source of knowledge. Contemporary Oncology (Poznan, Poland), 19(1A), A68–A77. https://doi.org/10.5114/WO.2014.47136
#Waks, A. G., & Winer, E. P. (2019). Breast Cancer Treatment: A Review. JAMA, 321(3), 288–300. https://doi.org/10.1001/JAMA.2018.19323
#Wartmann, H., Heins, S., Kloiber, K., & Bonn, S. (2021). Bias-invariant RNA-sequencing metadata annotation. GigaScience, 10(9), 1–13. https://doi.org/10.1093/GIGASCIENCE/GIAB064

# Package References:
citation("TCGAbiolinks")
citation("SummarizedExperiment")
citation ("edgeR")
citation("limma")
citation("tidyverse")
citation("glmnet")
citation("gplots")
citation("caret")
citation("survival")
citation("gprofiler2")

