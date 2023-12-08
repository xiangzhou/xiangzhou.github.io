---
layout: archive
title: "Software"
permalink: /software/
author_profile: true
---


## Bayesian Analysis for Spatial Segmentation (BASS)

BASS is a method for multi-scale and multi-sample analysis in spatial transcriptomics. BASS performs multi-scale transcriptomic analyses in the form of joint cell type clustering and spatial domain detection, with the two analytic tasks carried out simultaneously within a Bayesian hierarchical modeling framework. For both analyses, BASS properly accounts for the spatial correlation structure and seamlessly integrates gene expression information with spatial localization information to improve their performance. In addition, BASS is capable of multi-sample analysis that jointly models multiple tissue sections/samples, facilitating the integration of spatial transcriptomic data across tissue samples.

* The package is currently available at <a href="https://github.com/zhengli09/BASS">github</a>
* Citation: Zheng Li, and Xiang Zhou (2022). Multi-scale and multi-sample analysis enables accurate cell type clustering and spatial domain detection in spatial transcriptomic studies. Genome Biology. 23: 168.
* Contact <a href="mailto:zlisph@umich.edu"> Zheng Li </a> with any questions, comments, or bugs reports.


## Conditional AutoRegressive model based Deconvolution (CARD)

CARD is the software that leverages cell type specific expression information from single cell RNA sequencing (scRNA-seq) for the deconvolution of spatial transcriptomics. A unique feature of CARD is its ability to model the spatial correlation in cell type composition across tissue locations, thus enabling spatially informed cell type deconvolution. Modeling spatial correlation allows us to borrow the cell type composition information across locations on the entire tissue to accurately infer the cell type composition on each individual location, achieve robust deconvolution performance in the presence of mismatched scRNA-seq reference, impute cell type compositions and gene expression levels on unmeasured tissue locations, and facilitate the construction of a refined spatial tissue map with a resolution much higher than that measured in the original study.

* The software is currently available on <a href="https://github.com/YingMa0107/CARD/">github</a>, with a tutorial available <a href="https://yingma0107.github.io/CARD/">here</a>.
* All scripts for reproducing the results presented in the paper is available on <a href="https://github.com/YingMa0107/CARD-Analysis">github</a>.
* Citation: Ying Ma, and Xiang Zhou (2022). Spatially informed cell type deconvolution for spatial transcriptomics. Nature Biotechnology. 40: 1349–1359.
* Contact <a href="mailto:ying_ma@brown.edu"> Ying Ma </a> with any questions, comments, or bugs reports.

## COmposite likelihood-based COvariance regression NETwork model (CoCoNet)

CoCoNet is a composite likelihood-based covariance regression network model for identifying trait-relevant tissues or cell types. CoCoNet integrates tissue-specific gene co-expression networks constructed from either bulk or single cell RNA sequencing studies with association summary statistics from genome-wide association studies. CoCoNet relies on a covariance regression network model to express gene-level effect sizes for the given GWAS trait as a function of the tissue-specific co-expression adjacency matrix. With a composite likelihood-based inference algorithm, CoCoNet is scalable to tens of thousands of genes.

* The package is currently available <a href="http://lulushang.org/docs/Projects/CoCoNet">here</a>.
* All analysis scripts used in the paper is available <a href="http://lulushang.org/docs/Projects/CoCoNet/Reproduce">here</a>.
* Citation: Lulu Shang, Jennifer A. Smith and Xiang Zhou (2020). Leveraging gene co-expression patterns to infer trait-relevant tissues in genome-wide association studies. PLOS Genetics. 16: e1008734.
* Contact <a href="mailto:shanglu@umich.edu"> Lulu Shang </a> with any questions, comments, or bugs reports.