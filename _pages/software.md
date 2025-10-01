---
layout: archive
title: "Software"
permalink: /software/
author_profile: true
---


## Bayesian Analysis for Spatial Segmentation (BASS)

BASS is a method for multi-scale and multi-sample analysis in spatial transcriptomics. BASS performs multi-scale transcriptomic analyses in the form of joint cell type clustering and spatial domain detection, with the two analytic tasks carried out simultaneously within a Bayesian hierarchical modeling framework. For both analyses, BASS properly accounts for the spatial correlation structure and seamlessly integrates gene expression information with spatial localization information to improve their performance. In addition, BASS is capable of multi-sample analysis that jointly models multiple tissue sections/samples, facilitating the integration of spatial transcriptomic data across tissue samples.

* The software package is currently available on <a href="https://github.com/zhengli09/BASS">github</a>.
* Citation: Zheng Li, and Xiang Zhou (2022). Multi-scale and multi-sample analysis enables accurate cell type clustering and spatial domain detection in spatial transcriptomic studies. Genome Biology. 23: 168.
* Contact <a href="mailto:zlisph@umich.edu">Zheng Li</a> with any questions, comments, or bugs reports.


## CAsual relationship identificatioN using ONe sample instrumental variable model (Canon)

Canon (CAsual relationship identificatioN using ONe sample instrumental variable model) is a one-sample IV analysis method specifically designed to systematically identify genes potentially causally influenced by perturbed target genes across different sc-CRISPR platforms. Canon examines one gene pair at a time, employs a probabilistic model to automatically select valid IVs among gRNAs, relies on a scalable sampling-based algorithm for calibrated p-values computation, while explicitly accounting for the correlation between exposure and outcome genes inherent to the one-sample design and directly modeling the off-target effects of gRNAs on the outcome genes. As a result, Canon provides effective type I error control and high power across different sc-CRISPR platforms. 

* The software is currently available on <a href="https://github.com/pekjoonwu/Canon">github</a>.
* Citation: Peijun Wu, and Xiang Zhou (2025). Canon: Causal inference of downstream genes in single cell CRISPR studies via instrumental variable analysis. 
* Contact <a href="mailto:pekjoonw@umich.edu">Peijun Wu</a> with any questions, comments, or bugs reports.


## Conditional AutoRegressive model based Deconvolution (CARD)

CARD is the software that leverages cell type specific expression information from single cell RNA sequencing (scRNA-seq) for the deconvolution of spatial transcriptomics. A unique feature of CARD is its ability to model the spatial correlation in cell type composition across tissue locations, thus enabling spatially informed cell type deconvolution. Modeling spatial correlation allows us to borrow the cell type composition information across locations on the entire tissue to accurately infer the cell type composition on each individual location, achieve robust deconvolution performance in the presence of mismatched scRNA-seq reference, impute cell type compositions and gene expression levels on unmeasured tissue locations, and facilitate the construction of a refined spatial tissue map with a resolution much higher than that measured in the original study.

* The software is currently available on <a href="https://github.com/YingMa0107/CARD/">github</a>, with a tutorial available <a href="https://yingma0107.github.io/CARD/">here</a>.
* All scripts for reproducing the results presented in the paper is available on <a href="https://github.com/YingMa0107/CARD-Analysis">github</a>.
* Citation: Ying Ma, and Xiang Zhou (2022). Spatially informed cell type deconvolution for spatial transcriptomics. Nature Biotechnology. 40: 1349–1359.
* Contact <a href="mailto:ying_ma@brown.edu">Ying Ma</a> with any questions, comments, or bugs reports.



## CELl type-specific spatially variable gene IdentificatioN Analysis (CELINA)

CELINA is the software that can be used to systematically identify cell type specific spatially variable genes (ct-SVGs) across a variety of spatial transcriptomics platforms. Celina examines one gene at a time and uses a spatially varying coefficient model to explicitly and accurately model gene’s spatial expression pattern in relation to the cell type distribution across tissue locations. As a result, Celina provides effective type I error control and high statistical power in both single cell and spot resolution spatial transcriptomics.

* The software is currently available on the <a href="https://lulushang.org/Celina_Tutorial/index.html">website</a> with source code available on <a href="https://github.com/pekjoonwu/CELINA">github</a>.
* Citation: Lulu Shang, Peijun Wu, and Xiang Zhou (2025). Statistical identification of cell type-specific spatially variable genes in spatial transcriptomics. Nature Communications. in press.
* Contact <a href="mailto:lshang@mdanderson.org">Lulu Shang</a> or <a href="mailto:pekjoonw@umich.edu">Peijun Wu</a> with any questions, comments, or bugs reports.


## COmposite likelihood-based COvariance regression NETwork model (CoCoNet)

CoCoNet is a composite likelihood-based covariance regression network model for identifying trait-relevant tissues or cell types. CoCoNet integrates tissue-specific gene co-expression networks constructed from either bulk or single cell RNA sequencing studies with association summary statistics from genome-wide association studies. CoCoNet relies on a covariance regression network model to express gene-level effect sizes for the given GWAS trait as a function of the tissue-specific co-expression adjacency matrix. With a composite likelihood-based inference algorithm, CoCoNet is scalable to tens of thousands of genes.

* The package is currently available <a href="http://lulushang.org/docs/Projects/CoCoNet">here</a>.
* All analysis scripts used in the paper is available <a href="http://lulushang.org/docs/Projects/CoCoNet/Reproduce">here</a>.
* Citation: Lulu Shang, Jennifer A. Smith and Xiang Zhou (2020). Leveraging gene co-expression patterns to infer trait-relevant tissues in genome-wide association studies. PLOS Genetics. 16: e1008734.
* Contact <a href="mailto:lshang@mdanderson.org">Lulu Shang</a> with any questions, comments, or bugs reports.


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


## subcelluar Expression LocaLization Analysis (ELLA)

ELLA is a statistical method for modeling the subcellular localization of mRNAs and detecting genes that display spatial variation within cells in high-resolution spatial transcriptomics. ELLA utilizes a nonhomogeneous Poisson process to model the spatial count data within cells, creates a unified cellular coordinate system to anchor diverse shapes and morphologies across cells, and relies on an expression intensity function to capture the subcellular spatial distribution of mRNAs. ELLA can be applied to an arbitrary number of cells and detect a wide variety of subcellular localization patterns across diverse spatial transcriptomic techniques, while producing effective control of type I error and yielding high statistical power. With a computationally efficient algorithm, ELLA is scalable to tens of thousands of genes across tens of thousands of cells.

* The package is currently available on <a href="https://github.com/jadexq/ELLA/tree/main">github</a> with a tutorail available <a href="https://jadexq.github.io/ELLA/">here</a>.
* Citation: Jade Xiaoqing Wang, and Xiang Zhou (2024). ELLA: Modeling subcellular spatial variation of gene expression within cells in high-resolution spatial transcriptomics.
* Contact <a href="mailto:jadewang@umich.edu">Jade Wang</a> with any questions, comments, or bugs reports.


## Fine-mApping of causal genes for BInary Outcomes (FABIO)

FABIO is a transcriptome-wide association study (TWAS) fine-mapping method specifically designed for binary traits that is capable of modeling all genes jointly on an entire chromosome. FABIO relies on a probit model to directly relate multiple GReX to binary outcome. Additionally, it jointly models all genes located on a chromosome to account for the correlation among GReX arising from cis-SNP LD and expression correlation across genomic regions. As a result, FABIO effectively controls false discoveries while offering substantial power gains over existing TWAS fine-mapping approaches.

* The software is currently available on <a href="https://github.com/superggbond/FABIO/">github</a>.
* Citation: Haihan Zhang, Kevin He, Lam C. Tsoi, and Xiang Zhou (2024). FABIO: A TWAS fine-mapping method for prioritizing causal genes in binary traits. PLOS Genetics. e1011503.
* Contact <a href="mailto:hhzhang@umich.edu">Haihan Zhang</a> with any questions, comments, or bugs reports.


## Fast Cell-Cell Communication analysis (FastCCC)

FastCCC is a highly scalable, permutation-free statistical toolkit tailored to identify critical cell-cell communications (CCCs) in the form of ligand-receptor interactions (LRIs) and uncover novel biological insights in single-cell transcriptomics studies. FastCCC presents an analytic solution for computing p-values in CCC analysis, enabling scalable analysis without the need for computationally intensive permutations. It introduces a modular communication score computation framework that calculates various communication scores through a range of algebraic operations between ligand and receptor expression levels, capturing a broad spectrum of CCC patterns and ensuring robust analysis. Additionally, FastCCC not only enables the analysis of large-scale datasets containing millions of cells, but also introduces reference-based CCC analysis, where large-scale datasets are treated as reference panels to substantially improve CCC analysis on user-collected datasets. 

* The software is currently available on the <a href="https://svvord.github.io/FastCCC/">website</a>, with source code available on <a href="https://github.com/Svvord/FastCCC">github</a>.
* Citation: Siyu Hou, Wenjing Ma, and Xiang Zhou (2025). FastCCC: A permutation-free framework for scalable, robust, and reference-based cell-cell communication analysis in single cell transcriptomics studies. Nature Communications.
* Contact <a href="mailto:siyuh@umich.edu">Siyu Hou</a> with any questions, comments, or bugs reports.



## Fast Genotype-Environment interaction analysis (fastGxE)

fastGxE is a scalable and effective genome-wide multi-environment genotype-environment interaction (GxE) method designed to identify genetic variants that interact with one or multiple environmental factors. fastGxE controls for both polygenic effects and polygenic interaction effects, is robust to the number of environmental factors contributing to GxE interactions, and ensures scalability for genome-wide analysis in large biobank studies. fastGxE is accompanied with the algorithm mmSuSiE, which is an extension of SuSiE specifically designed for mixed-model analysis and aims to identify the environmental factors driving the detected GxE interactions by fastGxE. 

* The software is currently available on <a href="https://chaoning.github.io/fastGxE">github</a>.
* Citation: Chao Ning, and Xiang Zhou (2025). fastGxE: Powering genome-wide detection of genotype-environment interactions in biobank studies.
* Contact <a href="mailto:chaon@umich.edu">Chao Ning</a> with any questions, comments, or bugs reports.



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

* The software is currently available on <a href="https://github.com/genetics-statistics/GEMMA">github</a>, with a <a href="GEMMAmanual.pdf">User Manual</a> (last edited on 05/18/2016) and a sample code <a href="demo.txt">demo.txt</a> for variance component estimation.
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
* Citation: Lu Liu, Ran Yan, Ping Guo, Jiadong Ji, Weiming Gong, Fuzhong Xue, Zhongshang Yuan, and Xiang Zhou (2023). Conditional transcriptome-wide association study for fine-mapping causal genes. Nature Genetics. 56: 348–356.
* Contact <a href="mailto:luliuu@umich.edu">Lu Liu</a> with any questions, comments, or bugs reports.


## integrative Differential expression and gene set Enrichment Analysis (iDEA)

iDEA is a method for performing joint differential expression (DE) and gene set enrichment (GSE) analysis. iDEA builds upon a hierarchical Bayesian model for joint modeling of DE and GSE analyses. It uses only summary statistics as input, allowing for effective data modeling through complementing and pairing with various existing DE methods. It relies on an efficient expectation-maximization algorithm with internal Markov Chain Monte Carlo steps for scalable inference. By integrating DE and GSE analyses, iDEA can improve the power and consistency of DE analysis and the accuracy of GSE analysis over common existing approaches. 

* The package is currently available on <a href="https://xzhoulab.github.io/iDEA/">github</a>.
* All analysis scripts used in the paper is available <a href="https://github.com/xzhoulab/iDEA-Analysis">here</a>.
* Citation: Ying Ma\*, Shiquan Sun\*, Xuequn Shang, Evan T. Keller, Mengjie Chen and Xiang Zhou (2020). Integrative differential expression and gene set enrichment analysis using summary statistics for single cell RNAseq studies. Nature Communications. 11: 1585.
* Contact <a href="mailto:sqsunsph@xjtu.edu.cn">Shiquan Sun</a> or <a href="mailto:ying_ma@brown.edu">Ying Ma</a> with any questions, comments, or bugs reports.


## Integrative Methylation Association with GEnotypes (IMAGE)

IMAGE is a method that performs methylation quantitative trait locus (mQTL) mapping in bisulfite sequencing studies. IMAGE jointly accounts for both allele-specific methylation information from heterozygous individuals and non-allele-specific methylation information across all individuals, enabling powerful ASM-assisted mQTL mapping. In addition, IMAGE relies on an over-dispersed binomial mixed model to directly model count data, which naturally accounts for sample non-independence resulting from individual relatedness, population stratification, or batch effects that are commonly observed in sequencing studies. IMAGE relies on a penalized quasi-likelihood (PQL) approximation-based algorithm for scalable model inference.

* The software is currently available on <a href="https://github.com/fanyue322/IMAGE">github</a>.
* All analysis scripts used in the paper are available <a href="https://github.com/fanyue322/IMAGEreproduce">here</a>.
* Citation: Yue Fan, Tauras P. Vilgalys, Shiquan Sun, Qinke Peng, Jenny Tung and Xiang Zhou (2019). High-powered detection of genetic effects on DNA methylation using integrated methylation QTL mapping and allele-specific analysis. Genome Biology. 20: 220.
* Contact <a href="mailto:xafanyue@stu.xjtu.edu.cn">Yue Fan</a> with any questions, comments, or bugs reports.


## integrative MApping of Pleiotropic association (iMAP)

iMAP performs integrative mapping of pleiotropic association and functional annotations using penalized Gaussian mixture models. iMAP relies on a multinomial logistic regression model to incorporate a large number of binary and continuous SNP annotations, and, with a sparsity-inducing penalty term, is capable of selecting a small, informative set of annotations. In addition, iMAP directly models summary statistics from GWASs and uses a multivariate Gaussian distribution to account for phenotypic correlation between traits. As a result, iMAP is capable of integrating both binary and continuous SNP annotations, selecting informative annotations from a large set of potentially non-informative ones and using GWAS summary statistics while simultaneously accounting for phenotypic correlation between traits. 

* The software is currently available on <a href="https://github.com/biostatpzeng/iMAP">github</a>.
* An early version, <a href="iMAP-master.zip">version 1.00alpha</a>, was last modified on 10/04/2017.
* Citation: Ping Zeng, Xingjie Hao and Xiang Zhou (2018). Pleiotropic mapping and annotation selection in genome-wide association studies with penalized Gaussian mixture models. Bioinformatics. 34: 2797-2807.
* Contact <a href="mailto:zpstat@xzhmu.edu.cn">Ping Zeng</a> with any questions, comments, or bugs reports.


## Integrative and Reference-Informed tissue Segmentation (IRIS)

IRIS (Integrative and Reference-Informed tissue Segmentation) is a method for spatial domain detection in spatially resolved transcriptomics (SRT). IRIS models multiple tissue slices jointly and segments each tissue slice into multiple spatial domains. IRIS also accounts for the spatial correlation structure commonly observed across locations on each tissue slice and explicitly models the similarity in cell type composition underling similar spatial domains across tissue slices. A unique feature of IRIS is its ability to incorporate a scRNA-seq data to serve as the reference for domain detection, which allows IRIS to seamlessly integrate the cell type specific transcriptomic profiles from the scRNA-seq reference to the SRT dataset to substantially improve the accuracy in spatial domain detection. As a result, IRIS is accurate, scalable, and robust for spatial domain detection across a range of SRT technologies with distinct spatial resolutions.

* The software is currently available on <a href="https://github.com/YingMa0107/IRIS/">github</a>, with a tutorial available <a href="https://yingma0107.github.io/IRIS/">here</a>.
* Citation: Ying Ma, and Xiang Zhou (2024). Integrative and reference-informed spatial domain detection for spatial transcriptomics. Nature Methods. 21: 1231-1244.
* Contact <a href="mailto:ying_ma@brown.edu">Ying Ma</a> with any questions, comments, or bugs reports.


## Local genetic correlation estimation using summary statistics (Logica)

Logica (LOcal GenetIc Correlation across Ancestries) is a method specifically designed to estimate local genetic correlations across ancestries. Logica employs a bivariate linear mixed model that explicitly accounts for diverse LD patterns across ancestries, operates on GWAS summary statistics, and utilizes a maximum likelihood framework for robust inference. Logica enhances the accuracy of MoM estimates, produces well-controlled false discovery rates, and demonstrates greater power in detecting local genetically correlated regions. An important byproduct of Logica is its reformulation of a joint heritability test across ancestries, which explicitly accounts for the composite natural of the null hypothesis, resulting in well-calibrated P-values. 

* The software is currently available <a href="https://github.com/borangao/Logica">github</a>.
* Citation: Boran Gao, and Xiang Zhou (2025). Logica: A likelihood framework for cross-ancestry local genetic correlation estimation using summary statistics. 
* Please contact <a href="mailto:gao824@purdue.edu">Boran Gao</a> with any questions, comments, or bug reports.



## Mixed model Association for Count data via data AUgmentation (MACAU)

MACAU is the software implementing the Mixed model Association for Count data via data AUgmentation algorithm. It fits a binomial mixed model to perform differential methylation analysis for bisulfite sequencing studies. It fits a Poisson mixed model to perform differential expression analysis for RNA sequencing studies. It is computationally efficient for large scale sequencing studies and uses freely available open-source numerical libraries. 

* Downloads:
  * <a href="macau">Executable</a>, Version 1.40, for x86 64-bit Linux, compiled on 08/09/2019.
  * <a href="macau-v1.40.rar">Source Code</a>, Version 1.40, last modified on 08/09/2019.
  * <a href="macau-example.tar.gz">Example Data Sets</a>: a bisulfite sequencing data with 438,865 methylation sites for 50 baboons; and	an RNA squencing data with 12,018 genes for 63 baboons.
  * <a href="MACAUmanual.pdf">User Manual</a>, last edited on 03/15/2017.
  * <a href="bb.r">R script</a> to fit a beta-binomial model.
* Citations: 
  * Binomial mixed models for bisulfite sequencing analysis: Amanda J. Lea, Jenny Tung, and Xiang Zhou (2015). A flexible, effcient binomial mixed model for identifying differential DNA methylation in bisulfite sequencing data. PLoS Genetics. 11: e1005650.
  * Poisson mixed models for RNA sequencing analysis: Shiquan Sun, Michelle Hood, Laura Scott, Qinke Peng, Sayan Mukherjee, Jenny Tung, and Xiang Zhou (2017). Differential expression analysis for RNAseq using Poisson mixed models. Nucleic Acids Research. 45(11): e106.
* Contact <a href="mailto:sqsunsph@xjtu.edu.cn">Shiquan Sun</a> with any questions, comments, or bugs reports.


## MArginal ePIstasis Test (MAPIT)

MAPIT is the software implementing the new strategy for mapping epistasis: instead of directly identifying individual pairwise or higher-order interactions, MAPIT focuses on mapping variants that have non-zero <i>marginal epistatic effects</i> --- the combined pairwise interaction effects between a given variant and all other variants. By testing marginal epistatic effects, MAPIT can identify candidate variants that are involved in epistasis without the need to identify the exact partners with which the variants interact, thus potentially alleviating much of the statistical and computational burden associated with standard epistatic mapping procedures. MAPIT is based on a variance component model, and relies on a recently developed variance component estimation method for efficient parameter inference and p-value computation.

* The software is currently available on <a href="https://github.com/lorinanthony/MAPIT">github</a>.
* Citation: Lorin Crawford, Ping Zeng, Sayan Mukherjee, and Xiang Zhou (2017). Detecting epistasis with the marginal epistasis test in genetic mapping studies of quantitative traits. PLOS Genetics. e1006869.
* Contact <a href="mailto:lorin_crawford@brown.edu">Lorin Crawford</a> with any questions, comments, or bugs reports.


## Multi-ancestry Sum of the Single Effects Model (MESuSiE)

MESuSiE is a method for multi-ancestry fine-mapping analysis in genome-wide association studies. MESuSiE explicitly models both shared and ancestry-specific causal variants across ancestries, properly accounts for the diverse LD pattern observed in different ancestries, relies on GWAS marginal summary statistics as input, and extends the recent scalable variational inference algorithm SuSiE, which was developed for ancestry-specific fine-mapping, towards scalable multi-ancestry fine-mapping.

* The software is currently available on <a href="https://github.com/borangao/MESuSiE">github</a>.
* Citation: Boran Gao, and Xiang Zhou (2022). MESuSiE: Multi-ancestry fine-mapping for scalable and powerful discovery of shared and ancestry-specific putative causal variants in genome-wide association studies. Nature Genetics. 56: 170-179.
* Contact <a href="mailto:borang@umich.edu">Boran Gao</a> with any questions, comments, or bugs reports.


## Multi-ancEstry TRanscriptOme-wide analysis (METRO)

METRO is a method that leverages expression data collected from multiple genetic ancestries to enhance the power of TWAS. METRO incorporates expression prediction models constructed from multiple ancestries through a joint likelihood-based inference framework, allowing us to account for the uncertainty in the prediction models constructed in each expression study. In addition, METRO is capable of inferring the contribution of expression prediction models in different genetic ancestries towards explaining and informing the gene-trait association, allowing us to interrogate the ancestry-dependent transcriptomic mechanisms underlying gene-trait association.

* The package is currently available at <a href="https://github.com/zhengli09/METRO">github</a>.
* Citation: Zheng Li, Wei Zhao, Lulu Shang, Thomas H Mosley, Sharon LR Kardia, Jennifer A. Smith, and Xiang Zhou (2021). METRO: Multi-ancestry transcriptome-wide association studies for powerful gene-trait association detection. American Journal of Human Genetics. 109: 783-801.
* Contact <a href="mailto:zlisph@umich.edu">Zheng Li</a> with any questions, comments, or bugs reports.


## Mendelian Randomization with Automated Instrument Determination (MRAID)

MRAID is a software for carrying out Mendelian randomization analysis. MRAID borrows ideas from fine-mapping approaches to model an initial set of candidate SNP instruments that are in potentially high LD with each other and automatically selects among them the suitable instruments for MR analysis. MRAID also explicitly models two types of horizontal pleiotropic effects that are either uncorrelated or correlated with the instrumental effects on the exposure to ensure effective control of horizontal pleiotropy. MRAID achieves both analytic tasks through a joint likelihood inference framework and relies on a scalable sampling-based algorithm to compute calibrated p-values for causal inference. As a result, MRAID provides calibrated type I error control for causal effect testing in the presence of horizontal pleiotropy, reduces false positives and, as a by-product, estimates the proportion of SNPs exhibiting uncorrelated or correlated horizontal pleiotropy.

* The software is currently available on <a href="https://github.com/yuanzhongshang/MRAID">github</a>.
* Citation: Zhongshang Yuan, Lu Liu, Ping Guo, Ran Yan, Fuzhong Xue, and Xiang Zhou (2022). Likelihood based Mendelian randomization analysis with automated instrument selection and horizontal pleiotropic modeling. Science Advances. 8: eabl5744.
* Contact <a href="mailto:yuanzhongshang@sdu.edu.cn">Zhongshang Yuan</a> with any questions, comments, or bugs reports.


## Multi-trait assisted Polygenic Scores (mtPGS)

mtPGS is a method that constructs accurate PGS for a target trait of interest through leveraging multiple traits relevant to the target trait. Specifically, mtPGS borrows SNP effect size similarity information between the target trait and its relevant traits to improve the effect size estimation on the target trait, thus achieving accurate PGS. In the process, mtPGS flexibly models the shared genetic architecture between the target and the relevant traits to achieve robust performance, while explicitly accounting for the environmental covariance among them to accommodate different study designs with various sample overlap patterns. In addition, mtPGS uses only summary statistics as input and relies on a deterministic algorithm with several algebraic techniques for scalable computation.

* The software is currently available on <a href="https://github.com/xuchang0201/mtPGS">github</a>.
* Citation: Chang Xu, Santhi Ganesh, and Xiang Zhou (2023). mtPGS: Leverage multiple correlated traits for accurate polygenic score construction.
* Contact <a href="mailto:xuchang@umich.edu">Chang Xu</a> with any questions, comments, or bugs reports.


## Omnigenic Mendelian Randomization (OMR)

OMR is the software that explores the benefits of the omnigenic architecture for MR analysis. OMR builds upon the omnigenic modeling assumption on SNP effect sizes and use all genome-wide SNPs to serve as instrumental variables without any instrumental variable pre-selection. OMR imposes a general modeling assumption on the horizontal pleiotropic effects and relies on a scalable composite likelihood framework for causal effect inference.

* The software is currently available on <a href="https://github.com/wanglu205/OMR">github</a>.
* All scripts for reproducing the results presented in the paper is available <a href="https://github.com/wanglu205/OMRreproduce">here</a>.
* Citation: Lu Wang, Boran Gao, Yue Fan, Fuzhong Xue, and Xiang Zhou (2021). Mendelian randomization under the omnigenic architecture. Briefings in Bioinformatics. 22: bbab322.  
* Contact <a href="mailto:willa0205@yeah.net">Lu Wang</a> with any questions, comments, or bugs reports.


## Probabilistic Mendelian Randomization for TWAS (PMR-Egger and moPMR-Egger)

The PMR software package implements two methods:

PMR-Egger, which is a method that fits probabilistic Mendelian randomization with an Egger regression assumption on horizontal pleiotropy for transcriptome-wide association studies (TWASs). PMR-Egger relies on a new MR likelihood framework that unifies many existing TWAS and MR methods, accommodates multiple correlated instruments, tests the causal effect of gene on trait in the presence of horizontal pleiotropy, directly performs genome-wide test of horizontal pleiotropy, and, with a newly developed parameter expansion version of the expectation maximization algorithm, is scalable to hundreds of thousands of individuals.

moPMR-Egger, which extends PMR-Egger towards analyzing multiple outcome traits in TWAS applications. moPMR-Egger examines one gene at a time, relies on its cis-SNPs that are in potential linkage disequilibrium with each other to serve as instrumental variables, and tests its causal effects on multiple traits jointly. A key feature of moPMR-Egger is its ability to test and control for potential horizontal pleiotropic effects from instruments, thus maximizing power while minimizing false associations for TWASs. moPMR-Egger provides calibrated type I error control for both causal effects testing and horizontal pleiotropic effects testing and is more powerful than existing univariate TWAS approaches in detecting causal associations.

* The package is currently available on <a href="https://github.com/yuanzhongshang/PMR">github</a>.
* Citations:
  * Zhongshang Yuan, Huanhuan Zhu, Ping Zeng, Sheng Yang, Shiquan Sun, Can Yang, Jin Liu, and Xiang Zhou (2020). Testing and controlling for horizontal pleiotropy with the probabilistic Mendelian randomization in transcriptome-wide association studies. Nature Communications. 11: 3861.
  * Lu Liu, Ping Zeng, Fuzhong Xue, ZhongshangYuan, and Xiang Zhou (2021). Multi-trait transcriptome-wide association studies with probabilistic Mendelian randomization. American Journal of Human Genetics. 108: 240-256.
* Contact <a href="mailto:yuanzhongshang@sdu.edu.cn">Zhongshang Yuan</a> or <a href="mailto:luliuu@umich.edu">Lu Liu</a> with any questions, comments, or bugs reports.


## Penalized QuasiLikelihood for genomic Sequencing count data (PQLseq)

PQLseq is a method that fits generalized linear mixed models for analyzing RNA sequencing and bisulfite sequencing data. It estimates gene expression or methylation heritability for count data. It performs differential expression analysis in the presence of individual relatedness or population stratificaiton.

* The software is available on <a href="https://github.com/xzhoulab/PQLseq">github</a> and <a href="https://cran.r-project.org/web/packages/PQLseq/index.html">CRAN</a>, with a <a href="PQLseqManual.pdf">User Manual</a> (last edited on 03/20/2018) and <a href="pqlseq_SA_update.R">an updated R script</a> for accessing the residuals, tau1 and tau2 (modified on 04/05/2023).
* Citation: Shiquan Sun\*, Jiaqiang Zhu\*, Sahar Mozaffari, Carole Ober, Mengjie Chen and Xiang Zhou (2018). Heritability estimation and differential analysis with generalized linear mixed models in genomic sequencing studies. Bioinformatics. 35: 487-496.
* Contact <a href="mailto:sqsunsph@xjtu.edu.cn">Shiquan Sun</a> or <a href="mailto:jiaqiang@umich.edu">Jiaqiang Zhu</a> with any questions, comments, or bugs reports.


## PGS-based phenotype prediction interval (PredInterval)

PredInterval is designed to quantify phenotype prediction uncertainty through the construction of well-calibrated prediction intervals. PredInterval is non-parametric in natural and extracts information based on quantiles of phenotypic residuals through cross-validations, thus achieving well-calibrated coverage of true phenotypic values across a range of settings and traits with distinct genetic architecture. In addition, the PredInterval framework is general and can be paired with any PGS method. 

* The software is available on <a href="https://github.com/xuchang0201/PredInterval">github</a>.
* Citation: Chang Xu, Santhi Ganesh, and Xiang Zhou (2024). Statistical construction of calibrated prediction intervals for polygenic score based phenotype prediction.
* Contact <a href="mailto:xuchang@umich.edu">Chang Xu</a> with any questions, comments, or bugs reports.


## Sex-dimorphic mapping with the sum of the single effects model (sdSuSiE)

sdSuSiE is a software designed for powerful and scalable sex-dimorphic fine-mapping. sdSuSiE explicitly models sex differences in genetic effects while properly accounting for correlations between different types of genetic effects induced by LD and sample imbalance between sexes. It uses summary statistics from sex-stratified GWASs as inputs and extends the recent scalable variational inference algorithm SuSiE to fine-map sex-dimorphic effects. As a result, sdSuSiE achieves well-calibrated credible set coverage and posterior inclusion probabilities (PIPs), demonstrating high statistical power with controlled false discovery rates, all while maintaining computational efficiency. 

* The software is available on <a href="https://yuanzhongshang.github.io/sdSuSiE/">github</a>.
* Citation: Lu Liu, Zhongshang Yuan, and Xiang Zhou (2025). Sex-dimorphic fine-mapping with summary statistics in genome-wide association studies. 
* Contact <a href="mailto:luliuu@umich.edu">Lu Liu</a> with any questions, comments, or bugs reports.


## Scalable Multiple Annotation integration for trait-Relevant Tissue identification (SMART)

SMART is a software implementing the Scalable Multiple Annotation integration for trait-Relevant Tissue identification and usage. It extends the commonly used linear mixed model to relate variant effect sizes to variant annotations by introducing variant specific variance components that are functions of multiple annotations. It quantifies and evaluates the joint contribution of multiple annotations to genetic effect sizes by performing parameter inference using the widely used generalized estimation equation (GEE). The GEE-based algorithm in SMART allows for the use of summary statistics and naturally accounts for the correlation among summary statistics due to linkage disequilibrium. With GEE statistics, SMART applies mixture models to classify tissues into two categories—those that are relevant to the trait and those that are not—thus formulating the task of identifying trait-relevant tissues into a classification problem. 

* Software: <a href="SMART_0.1.0.tar.gz">version 0.10</a>, released on 10/10/2017, with a <a href="SMART_example.rar">Data Example</a>.
* Manual: Use install.packages() command to install the package. Use ?SMART or ?LDSM or ?EMmix to see the description of these functions and some simple scripts to run the data example. A detailed manual will appear soon.
* Citation: Xingjie Hao, Ping Zeng, Shujun Zhang and Xiang Zhou (2018). Identifying and exploiting trait-relevant tissues with multiple functional annotations in genome-wide association studies. PLOS Genetics. e1007186.
* Contact <a href="mailto:xingjiegenetics@qq.com">Xingjie Hao</a> with any questions, comments, or bugs reports.


## Spatial Domain Transition detection (SpaDOT)

SpaDOT a deep-learning-based computational method designed to accurately identify spatial domains and infer their transitions and relationships across time points in spatiotemporal transcriptomics data. SpaDOT employs a variational autoencoder (VAE) framework integrated with a Gaussian Process prior and graph neighbor information to capture the low-dimensional manifold underlying the data, incorporating hidden clustering variables to characterize distinct spatial domains. In the process, SpaDOT explicitly accounts for both global and local spatial correlation structures among spatial locations within the tissue at each time point and leverages optimal transport (OT) analysis to infer the relationships between spatial domains across time points, thereby enhancing accuracy and interpretation.

* The software is currently available on <a href="https://github.com/marvinquiet/SpaDOT">github</a>, with a tutorial available <a href="https://marvinquiet.github.io/SpaDOT/">here</a>.
* Citation: Wenjing Ma, Siyu Hou, Lulu Shang, Jiaying Lu, and Xiang Zhou (2025). Optimal transport modeling uncovers spatial domain dynamics in spatiotemporal transcriptomics studies. 
* Contact <a href="mailto:wenjinma@umich.edu">Wenjing Ma</a> with any questions, comments, or bugs reports.


## Spatiall aware probabilistic Principal Component Analysis (SpatialPCA)

SpatialPCA is the software that perform spatially aware dimension reduction for spatial transcriptomics. SpatialPCA explicitly models the spatial correlation structure across tissue locations to preserve the neighboring similarity of the original data in the low dimensional manifold. The low dimensional components obtained from SpatialPCA thus contain valuable spatial correlation information and can be directly paired with existing computational tools developed in scRNA-seq for effective and novel downstream analysis in spatial transcriptomics. In particular, the low-dimensional components from SpatialPCA can be paired with scRNA-seq clustering methods to enable effective de novo tissue domain detection and can be paired with scRNA-seq trajectory inference methods to enable effective developmental trajectory inference on the tissue. Because of the data generative nature of SpatialPCA and its explicit modeling of spatial correlation, it can also be used to impute the low dimensional components on new and unmeasured spatial locations, facilitating the construction of a refined spatial map with a resolution much higher than that measured in the original study.

* The software is currently available on <a href="https://github.com/shangll123/SpatialPCA">github</a>, with a tutorial available <a href="http://lulushang.org/SpatialPCA_Tutorial">here</a>.
* All scripts for reproducing the results presented in the paper is available
  <a href="http://lulushang.org/docs/Projects/SpatialPCA">here</a>.
* Citation: Lulu Shang, and Xiang Zhou (2022). Spatially aware dimension reduction for spatial transcriptomics. Nature Communications. 13: 7203.
* Contact <a href="mailto:lshang@mdanderson.org">Lulu Shang</a> with any questions, comments, or bugs reports.


## Spatial PAttern Recognition via Kernels (SPARK and SPARK-X)

SPARK and SPARK-X are methods for detecting genes with spatial expression patterns in spatially resolved transcriptomic studies. SPARK directly models count data generated from various spatial resolved transcriptomic techniques through generalized spatial linear models. While SPARK-X relies on a scalable non-parametric testing framework to model a wide variety of spatial transcriptomics data collected through different technologies. Both SPARK and SPARK-X rely on algebraic innovations for scalable computatation as well as newly developed statistical formulas for hypothesis testing, producing well-calibrated p-values and yielding high statistical power. Both SPARK and SPARK-X are implemented in the SPARK software.

* The software is currently available on <a href="https://xzhoulab.github.io/SPARK/">github</a>.
* All analysis scripts used in the two papers are available <a href="https://xzhoulab.github.io/SPARK/03_experiments/">here</a>.
* Citations:
  * Shiquan Sun\*, Jiaqiang Zhu\*, and Xiang Zhou (2020). Statistical analysis of spatial expression pattern for spatially resolved transcriptomic studies. Nature Methods. 17: 193-200. 
  * Jiaqiang Zhu, Shiquan Sun, and Xiang Zhou (2021). Non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies. Genome Biology. 22: 184.
* Contact <a href="mailto:jiaqiang@umich.edu">Jiaqiang Zhu</a> or <a href="mailto:sqsunsph@xjtu.edu.cn">Shiquan Sun</a> with any questions, comments, or bugs reports.


## Spatially Resolved Transcriptomics simulator (SRTsim)

SRTsim is a software simulator for generating synthetic spatially resolved transcriptomic (SRT) data based on a wide variety of SRT techniques. SRTsim incorporates spatial localization information to simulate SRT expression count data in a reproducible and scalable fashion, thus facilitating SRT experimental design and methodology development. A key benefit of SRTsim is its ability to not only maintain various location-wise and gene-wise SRT count properties but also preserve the spatial expression patterns of the SRT data on the tissue, thus making it feasible to evaluate SRT method performance for various SRT-specific analytic tasks using the synthetic data.

* The software is currently available on <a href="https://xzhoulab.github.io/SRTsim/">github</a>.
* Citation: Jiaqiang Zhu*, Lulu Shang*, and Xiang Zhou (2022). SRTsim: spatial pattern preserving simulations for spatially resolved transcriptomics.
* Contact <a href="mailto:jiaqiang@umich.edu">Jiaqiang Zhu</a> or <a href="mailto:lshang@mdanderson.org">Lulu Shang</a> with any questions, comments, or bugs reports.


## Variational Inference based Probabilistic Canonical Correlation Analysis (VIPCCA)

VIPCCA is a software that relies based on a non-linear probabilistic canonical correlation analysis, for effective and scalable single cell data alignment. VIPCCA leverages cutting-edge techniques from deep neural network for non-linear modeling of single cell data, thus allowing users to capture the complex biological structures from integration of multiple single-cell datasets across technologies, data types, conditions, and modalities. In addition, VIPCCA relies on variational inference for scalable computation, enabling efficient integration of large-scale single cell datasets with millions of cells. Importantly, VIPCCA can transform multi-modalities into lower dimensional space without any post-hoc data processing, a unique and desirable feature that is in direct contrast to existing alignment methods.

* The software is currently available on <a href="https://github.com/jhu99/vipcca">github</a>.
* All analysis scripts used in the paper are available <a href="https://github.com/jhu99/vipcca_paper">here</a>. </p>
* Citation: Jialu Hu, Mengjie Chen, and Xiang Zhou (2022). Effective and scalable single-cell data alignment with non-linear canonical correlation analysis. Nucleic Acids Research. 50: e21. 
* Contact <a href="mailto:hujialu.xd@gmail.com">Jialu Hu</a> with any questions, comments, or bugs reports.


## Variability preserving ImPutation for Expression Recovery (VIPER)

VIPER is a method that performs Variability Preserving ImPutation for Expression Recovery in single cell RNA sequencing studies. VIPER is based on nonnegative sparse regression models and is capable of progressively inferring a sparse set of local neighborhood cells that are most predictive of the expression levels of the cell of interest for imputation. A key feature of VIPER is its ability to preserve gene expression variability across cells after imputation.

* The package is currently available on <a href="https://github.com/ChenMengjie/VIPER">github</a>.
* All analysis scripts used in the paper is available <a href="https://github.com/ChenMengjie/Vpaper2018">here</a>.
* Citation: Mengjie Chen and Xiang Zhou (2018). VIPER: variability-preserving imputation foraccurate gene expression recovery insingle-cell RNA sequencing studies. Genome Biology. 19:196.
* Contact <a href="mailto:mengjiechen@uchicago.edu"> Mengjie Chen </a> with any questions, comments, or bugs reports.


## Variant-set test INtegrative TWAS for GEne-based analysis (VINTAGE)

VINTAGE is a unified statistical framework for integrative analysis of GWAS and eQTL mapping studies to identify and decipher gene-trait associations. VINTAGE explicitly quantifies and tests the proportion of genetic effects on a trait potentially mediated through gene expression using a local genetic correlation test, and further leverages such information to guide the integration of gene expression mapping study towards gene association mapping in GWAS through a genetic variance test. The explicit quantification of local genetic correlation in VINTAGE allows its gene association test to unify two seemingly unrelated methods, SKAT and TWAS, into the same analytic framework and include both as special cases, thus achieving robust performance across a range of scenarios. 

* The package is currently available on <a href="https://github.com/zhengli09/VINTAGE">github</a>, with both a tutorial and code for reproducing all real data analysis results available <a href="https://zhengli09.github.io/VINTAGE-analysis/">here</a>.
* Citation: Zheng Li, Boran Gao, and Xiang Zhou (2024). VINTAGE: A unified framework integrating gene expression mapping studies with genome-wide association studies for detecting and deciphering gene-trait associations.
* Contact <a href="mailto:zlisph@umich.edu"> Zheng Li </a> with any questions, comments, or bugs reports.


## Paternity Inference from Low-Coverage Sequencing Data (WHODAD)

WHODAD is a software package implementing the WHODAD method for paternity inference from low-coverage sequencing data.

* Software <a href="whodad.tar">version 1.00alpha</a> (compiled on 10/16/2015), with a <a href="WHODADmanual.pdf">User Manual</a> (last edited on 10/20/2015).
* Citation: Noah Snyder-Mackler, William H Majoros, Michael L Yuan, Amanda O Shaver, Jacob B Gordon, Gisela H Kopp, Stephen A Schlebusch, Jeffrey D Wall, Susan C Alberts, Sayan Mukherjee, Xiang Zhou, and Jenny Tung (2016). Efficient genome-wide sequencing and low-coverage pedigree analysis from non-invasively collected
samples. Genetics. 203: 699-714.
