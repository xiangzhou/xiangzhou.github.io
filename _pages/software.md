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