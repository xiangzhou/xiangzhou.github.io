---
layout: archive
title: "Software"
permalink: /software/
author_profile: true
---


## Bayesian Analysis for Spatial Segmentation (BASS)

BASS is a method for multi-scale and multi-sample analysis in spatial transcriptomics. BASS performs multi-scale transcriptomic analyses in the form of joint cell type clustering and spatial domain detection, with the two analytic tasks carried out simultaneously within a Bayesian hierarchical modeling framework. For both analyses, BASS properly accounts for the spatial correlation structure and seamlessly integrates gene expression information with spatial localization information to improve their performance. In addition, BASS is capable of multi-sample analysis that jointly models multiple tissue sections/samples, facilitating the integration of spatial transcriptomic data across tissue samples.

* The package is currently available on <a href="https://github.com/zhengli09/BASS">github</a>
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
* Scripts to fit the other PGS methods are also available at <a href="https://biostat0903.github.io/DBSLMM/Scripts.html">github</a>.
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
* Scripts to reproduce all results in the manuscript are also available at <a href="https://github.com/fanyue322/ECCOreproduce">github</a>.
* Citation: Fan Yue, Huanhuan Zhu, Yanyi Song, Qinke Peng and Xiang Zhou (2021). Efficient and effective control of confounding in eQTL mapping studies through joint differential expression and Mendelian randomization analyses. Bioinformatics. 37: 296–302.
* Contact <a href="mailto:xafanyue@stu.xjtu.edu.cn">Yue Fan</a> with any questions, comments, or bugs reports.


## Genetic and Environmental Covariance estimation by composite-liKelihood Optimization (GECKO)

GECKO is a computational method for estimating both genetic and environmental covariances using GWAS summary statistics. GECKO improves estimation accuracy of method of moments algorithms while keeping computation in check. GECKO relies on composite likelihood, is scalable computationally, uses only on summary statistics, provides accurate genetic and environmental covariance estimates across a range of scenarios, and accommodates SNP annotation stratified covariance estimation.

* The software is currently available on <a href="https://github.com/borangao/GECKO">github</a>.
* Citation: Boran Gao, Can Yang, Jin Liu, and Xiang Zhou (2021). Accurate genetic and environmental covariance estimation with composite likelihood in genome-wide association studies. PLOS Genetics. e1009293.
* Contact <a href="mailto:borang@umich.edu">Boran Gao</a> with any questions, comments, or bugs reports.

