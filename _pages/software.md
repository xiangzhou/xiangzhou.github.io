---
layout: archive
title: "Software"
permalink: /software/
author_profile: true
---


## Bayesian Analysis for Spatial Segmentation (BASS)

BASS is a method for multi-scale and multi-sample analysis in spatial transcriptomics. BASS performs multi-scale transcriptomic analyses in the form of joint cell type clustering and spatial domain detection, with the two analytic tasks carried out simultaneously within a Bayesian hierarchical modeling framework. For both analyses, BASS properly accounts for the spatial correlation structure and seamlessly integrates gene expression information with spatial localization information to improve their performance. In addition, BASS is capable of multi-sample analysis that jointly models multiple tissue sections/samples, facilitating the integration of spatial transcriptomic data across tissue samples.

* The package is currently available on <a href="https://github.com/zhengli09/BASS">github</a>.
* Citation: Zheng Li, and Xiang Zhou (2022). Multi-scale and multi-sample analysis enables accurate cell type clustering and spatial domain detection in spatial transcriptomic studies. Genome Biology. 23: 168.
* Contact <a href="mailto:zlisph@umich.edu">Zheng Li</a> with any questions, comments, or bugs reports.


## Conditional AutoRegressive model based Deconvolution (CARD)

CARD is the software that leverages cell type specific expression information from single cell RNA sequencing (scRNA-seq) for the deconvolution of spatial transcriptomics. A unique feature of CARD is its ability to model the spatial correlation in cell type composition across tissue locations, thus enabling spatially informed cell type deconvolution. Modeling spatial correlation allows us to borrow the cell type composition information across locations on the entire tissue to accurately infer the cell type composition on each individual location, achieve robust deconvolution performance in the presence of mismatched scRNA-seq reference, impute cell type compositions and gene expression levels on unmeasured tissue locations, and facilitate the construction of a refined spatial tissue map with a resolution much higher than that measured in the original study.

* The software is currently available on <a href="https://github.com/YingMa0107/CARD/">github</a>, with a tutorial available <a href="https://yingma0107.github.io/CARD/">here</a>.
* All scripts for reproducing the results presented in the paper is available on <a href="https://github.com/YingMa0107/CARD-Analysis">github</a>.
* Citation: Ying Ma, and Xiang Zhou (2022). Spatially informed cell type deconvolution for spatial transcriptomics. Nature Biotechnology. 40: 1349–1359.
* Contact <a href="mailto:ying_ma@brown.edu">Ying Ma</a> with any questions, comments, or bugs reports.


## COmposite likelihood-based COvariance regression NETwork model (CoCoNet)

CoCoNet is a composite likelihood-based covariance regression network model for identifying trait-relevant tissues or cell types. CoCoNet integrates tissue-specific gene co-expression networks constructed from either bulk or single cell RNA sequencing studies with association summary statistics from genome-wide association studies. CoCoNet relies on a covariance regression network model to express gene-level effect sizes for the given GWAS trait as a function of the tissue-specific co-expression adjacency matrix. With a composite likelihood-based inference algorithm, CoCoNet is scalable to tens of thousands of genes.

* The package is currently available <a href="http://lulushang.org/docs/Projects/CoCoNet">here</a>.
* All analysis scripts used in the paper is available <a href="http://lulushang.org/docs/Projects/CoCoNet/Reproduce">here</a>.
* Citation: Lulu Shang, Jennifer A. Smith and Xiang Zhou (2020). Leveraging gene co-expression patterns to infer trait-relevant tissues in genome-wide association studies. PLOS Genetics. 16: e1008734.
* Contact <a href="mailto:shanglu@umich.edu">Lulu Shang</a> with any questions, comments, or bugs reports.


## Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)

DBSLMM is an accurate and scalable method for constructing polygenic scores in large biobank scale data sets. DBSLMM relies on a flexible modeling assumption on the effect size distribution to achieve robust and accurate prediction performance across a range of genetic architectures. DBSLMM also relies on a simple deterministic search algorithm to yield an approximate analytic estimation solution using summary statistics only, which, when paired with further algebraic innovations, resulting in substantial computational savings.

* The package is currently available on <a href="https://github.com/biostat0903/DBSLMM">github</a>.
* Scripts to fit the other PGS methods are available at <a href="https://biostat0903.github.io/DBSLMM/Scripts.html">here</a>.
* Citation: Sheng Yang and Xiang Zhou (2020). Accurate and scalable construction of polygenic scores in large biobank data sets. American Journal of Human Genetics. 106: 679-693.
* Contact <a href="mailto:yangsheng@njmu.edu.cn">Sheng Yang</a> with any questions, comments, or bugs reports.


## latent Dirichlet Process Regression (DPR)

DPR is a software package implementing the latent Dirichlet process regression method for genetic prediction of complex traits. DPR relies on the Dirichlet process to assign a prior on the effect size distribution itself and is thus capable of inferring an effect size distribution from the data at hand. Effectively, DPR uses infinitely many parameters a priori to character the effect size distribution, and with such a flexible modeling assumption, DPR is capable of adapting to a broad spectrum of genetic architectures and achieves robust predictive performance across a wide range of complex traits.

* The package is currently available on <a href="https://github.com/biostatpzeng/DPR">github</a>.
* Citation: Ping Zeng and Xiang Zhou (2017). Non-parametric genetic prediction of complex traits with latent Dirichlet process regression models. Nature Communications. 8: 456.
* Contact <a href="mailto:zpstat@xzhmu.edu.cn">Ping Zeng</a> with any questions, comments, or bugs reports.


## Effect size Correlation for COnfounding determination (ECCO)

ECCO is computationally efficient approach for determining the optimal number of PEER factors for eQTL mapping analysis. ECCO requires the availability of an outcome phenotype in addition to the usual genotype and expression data required for eQTL mapping studies. With the outcome phenotype, ECCO estimates the gene expression effect on the phenotype for one gene at a time through two different analyses: a differential expression regression analysis and a Mendelian randomization (MR) analysis. By computing and examining the correlation between the estimated effect sizes from the two different analyses, ECCO can subsequently determine the optimal number of PEER factors for eQTL mapping analysis.

* The package, along with its GTEx eQTL mapping results, is currently available on <a href="https://github.com/fanyue322/ECCO">github</a>.
* Scripts to reproduce all results in the manuscript are available at <a href="https://github.com/fanyue322/ECCOreproduce">here</a>.
* Citation: Fan Yue, Huanhuan Zhu, Yanyi Song, Qinke Peng and Xiang Zhou (2021). Efficient and effective control of confounding in eQTL mapping studies through joint differential expression and Mendelian randomization analyses. Bioinformatics. 37: 296–302.
* Contact <a href="mailto:xafanyue@stu.xjtu.edu.cn">Yue Fan</a> with any questions, comments, or bugs reports.


## Genetic and Environmental Covariance estimation by composite-liKelihood Optimization (GECKO)

GECKO is a computational method for estimating both genetic and environmental covariances using GWAS summary statistics. GECKO improves estimation accuracy of method of moments algorithms while keeping computation in check. GECKO relies on composite likelihood, is scalable computationally, uses only on summary statistics, provides accurate genetic and environmental covariance estimates across a range of scenarios, and accommodates SNP annotation stratified covariance estimation.

* The software is currently available on <a href="https://github.com/borangao/GECKO">github</a>.
* Citation: Boran Gao, Can Yang, Jin Liu, and Xiang Zhou (2021). Accurate genetic and environmental covariance estimation with composite likelihood in genome-wide association studies. PLOS Genetics. e1009293.
* Contact <a href="mailto:borang@umich.edu">Boran Gao</a> with any questions, comments, or bugs reports.


## Genome-wide Efficient Mixed Model Association (GEMMA: LMM, mvLMM, BSLMM, and MQS)

GEMMA is the software implementing the Genome-wide Efficient Mixed Model Association algorithm for a standard linear mixed model and some of its close relatives for genome-wide association studies (GWAS):

* It fits a univariate linear mixed model (LMM) for marker association tests with a single phenotype to account for population stratification and sample structure, and for estimating the proportion of variance in phenotypes explained (PVE) by typed genotypes (i.e. "chip heritability").
* It fits a multivariate linear mixed model (mvLMM) for testing marker associations with multiple phenotypes simultaneously while controlling for population stratification, and for estimating genetic correlations among complex phenotypes.
* It fits a Bayesian sparse linear mixed model (BSLMM) using Markov chain Monte Carlo (MCMC) for estimating PVE by typed genotypes, predicting phenotypes, and identifying associated markers by jointly modeling all markers while controlling for population structure.
* It estimates variance component/SNP heritability, and partitions it by different SNP functional categories. In particular, it uses HE regression or REML AI algorithm to estimate variance components when individual-level data are available. It uses MQS to estimate variance components when only summary statisics are available. 

It is computationally efficient for large scale GWAS and uses freely available open-source numerical libraries.

* The software is currently available on <a href="https://github.com/genetics-statistics/GEMMA">github</a>, with a <a href="GEMMAmanual.pdf">User Manual</a> (last edited on 05/18/2016) and a draft updated <a href="demo.txt">demo.txt</a> for variance component estimation.
* An old version, <a href="gemma-0.94.tar.gz">version 0.94</a>, is compiled on 01/12/2014.
* Citations:
  * Software tool and univariate linear mixed model: Xiang Zhou and Matthew Stephens (2012). Genome-wide efficient mixed-model analysis for association studies. Nature Genetics. 44: 821–824.
  * Multivariate linear mixed model: Xiang Zhou and Matthew Stephens (2014). Efficient multivariate linear mixed model algorithms for genome-wide association studies. Nature Methods. 11(4): 407–409.
  * Bayesian sparse linear mixed model: Xiang Zhou, Peter Carbonetto and Matthew Stephens (2013). Polygenic modeling with Bayesian sparse linear mixed
models. PLoS Genetics. 9(2): e1003264.
  * Variance component estimation: Xiang Zhou (2017). A unified framework for
variance component estimation with summary statistics in genome-wide
association studies. Annals of Applied Statistics. 11(4): 2027-2051.
* Click <a href="https://github.com/genetics-statistics/GEMMA/issues">here</a> if you have any questions, comments, or bugs reports.


## Gene-based Integrative Fine-mapping through conditional TWAS (GIFT)

GIFT is a Gene-based Integrative Fine-mapping for performing conditional TWAS analysis. GIFT examines one genomic region at a time, jointly models the GReX of all genes residing in the focal region, and carries out TWAS conditional analysis in a maximum likelihood framework. In the process, GIFT explicitly models the gene expression correlation and cis-SNP LD across different genes in the region and accounts for the uncertainty in the constructed GReX. As a result, GIFT provides effective type I error control, refines marginal TWAS associations into a much smaller set of putatively causal associations, and yields high statistical power with reduced false discoveries.

* The software is currently available on <a href="https://yuanzhongshang.github.io/GIFT/">github</a>.
* Citation: Lu Liu, Ran Yan, Ping Guo, Jiadong Ji, Weiming Gong, Fuzhong Xue, Zhongshang Yuan, and Xiang Zhou (2023). Conditional transcriptome-wide association study for fine-mapping causal genes. Nature Genetics. in press.
* Contact <a href="mailto:luliuu@umich.edu">Lu Liu</a> with any questions, comments, or bugs reports.


## integrative Differential expression and gene set Enrichment Analysis (iDEA)

iDEA is a method for performing joint differential expression (DE) and gene set enrichment (GSE) analysis. iDEA builds upon a hierarchical Bayesian model for joint modeling of DE and GSE analyses. It uses only summary statistics as input, allowing for effective data modeling through complementing and pairing with various existing DE methods. It relies on an efficient expectation-maximization algorithm with internal Markov Chain Monte Carlo steps for scalable inference. By integrating DE and GSE analyses, iDEA can improve the power and consistency of DE analysis and the accuracy of GSE analysis over common existing approaches. 

* The package is currently available on <a href="https://xzhoulab.github.io/iDEA/">github</a>.
* All analysis scripts used in the paper is available <a href="https://github.com/xzhoulab/iDEA-Analysis">here</a>.
* Citation: Ying Ma, Shiquan Sun, Xuequn Shang, Evan T. Keller, Mengjie Chen and Xiang Zhou (2020). Integrative differential expression and gene set enrichment analysis using summary statistics for single cell RNAseq studies. Nature Communications. 11: 1585.
* Contact <a href="mailto:sqsunsph@xjtu.edu.cn">Shiquan Sun</a> or <a href="mailto:ying_ma@brown.edu">Ying Ma</a> with any questions, comments, or bugs reports.


## Integrative Methylation Association with GEnotypes (IMAGE)

IMAGE is a method that performs methylation quantitative trait locus (mQTL) mapping in bisulfite sequencing studies. IMAGE jointly accounts for both allele-specific methylation information from heterozygous individuals and non-allele-specific methylation information across all individuals, enabling powerful ASM-assisted mQTL mapping. In addition, IMAGE relies on an over-dispersed binomial mixed model to directly model count data, which naturally accounts for sample non-independence resulting from individual relatedness, population stratification, or batch effects that are commonly observed in sequencing studies. IMAGE relies on a penalized quasi-likelihood (PQL) approximation-based algorithm for scalable model inference.

* The package is currently available on <a href="https://github.com/fanyue322/IMAGE">github</a>.
* All analysis scripts used in the paper are available <a href="https://github.com/fanyue322/IMAGEreproduce">here</a>.
* Citation: Yue Fan, Tauras P. Vilgalys, Shiquan Sun, Qinke Peng, Jenny Tung and Xiang Zhou (2019). High-powered detection of genetic effects on DNA methylation using integrated methylation QTL mapping and allele-specific analysis. Genome Biology. 20: 220.
* Contact <a href="mailto:xafanyue@stu.xjtu.edu.cn">Yue Fan</a> with any questions, comments, or bugs reports.