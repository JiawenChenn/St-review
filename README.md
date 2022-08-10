# St-review
This repo includes codes and processed data used to produce results in 

* [A comprehensive comparison on cell type composition inference for spatial transcriptomics data](https://academic.oup.com/bib/article-abstract/23/4/bbac245/6618233?redirectedFrom=fulltext)

## Repository structure
* ``processed_data``: This folder contains processed data for the three datasets we used in our analysis: mouse olfactory bulb (MOB), human developing heart, and mouse cortex. Both spatial transcriptomics data that are subject to cell type deconvolution inference and reference single cell RNA-seq data are released. Details are described in Table 2 of our [review paper](https://www.biorxiv.org/content/10.1101/2022.02.20.481171v1).
* ``scripts``: 
  * ``evaluation_metric.R``: codes used to calculate evaluation metrics between inferred and true cell type proportions. We used the following three metrics: RMSE, distance correlation, and difference.
  * other folders named after the name of each method include the code used to perform inference in the MOB data using the internal reference. Please refer to the following github/websites of each method for the most up-to-date pipeline.

### Related link
* Adroit https://github.com/TaoYang-dev/AdRoit
* cell2location https://cell2location.readthedocs.io/en/latest/index.html
* destVI https://docs.scvi-tools.org/en/stable/user_guide/models/destvi.html
* RCTD https://github.com/dmcable/spacexr
* STdeconvolve https://jef.works/STdeconvolve/
* stereoscope https://github.com/almaan/stereoscope
* SpatialDWLS http://giottosuite.com/articles/analyses_deconvolution_Oct2021.html
* SPOTlight https://github.com/MarcElosua/SPOTlight
* DSTG https://github.com/Su-informatics-lab/DSTG
* Tangram https://github.com/broadinstitute/Tangram  https://tangram-sc.readthedocs.io/en/latest/index.html

## Reference

* Yang, T., et al., AdRoit is an accurate and robust method to infer complex transcriptome composition. Communications Biology, 2021. 4(1): p. 1218.
* Kleshchevnikov, V., et al., Cell2location maps fine-grained cell types in spatial transcriptomics. Nature Biotechnology, 2022.
* Lopez, R., et al., DestVI identifies continuums of cell types in spatial transcriptomics data. Nature Biotechnology, 2022.
* Cable, D.M., et al., Robust decomposition of cell type mixtures in spatial transcriptomics. Nature Biotechnology, 2021.
* Miller, B.F., et al., Reference-free cell-type deconvolution of multi-cellular pixel-resolution spatially resolved transcriptomics data. bioRxiv, 2021: p. 2021.06.15.448381.
* Andersson, A., et al., Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography. Communications Biology, 2020. 3(1): p. 565.
* Dong, R. and G.-C. Yuan, SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biology, 2021. 22(1): p. 145.
* Elosua-Bayes, M., et al., SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Research, 2021. 49(9): p. e50-e50.
* Song, Q. and J. Su, DSTG: deconvoluting spatial transcriptomics data through graph-based artificial intelligence. Briefings in Bioinformatics, 2021. 22(5): p. bbaa414.
* Dries, R., et al., Giotto: a toolbox for integrative analysis and visualization of spatial expression data. Genome Biology, 2021. 22(1): p. 78.
* Elosua-Bayes, M., et al., SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Research, 2021. 49(9): p. e50-e50.
* Song, Q. and J. Su, DSTG: deconvoluting spatial transcriptomics data through graph-based artificial intelligence. Briefings in Bioinformatics, 2021. 22(5): p. bbaa414.
* Biancalani, T., et al., Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. Nature Methods, 2021. 18(11): p. 1352-1362.
* Asp, M., et al., A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart. Cell, 2019. 179(7): p. 1647-1660.e19.
* Tepe, B., et al., Single-Cell RNA-Seq of Mouse Olfactory Bulb Reveals Cellular Heterogeneity and Activity-Dependent Molecular Census of Adult-Born Neurons. Cell Reports, 2018. 25(10): p. 2689-2703.e3.
* Yao, Z., et al., A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. Cell, 2021. 184(12): p. 3222-3241.e26.
* Stickels, R.R., et al., Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2. Nature Biotechnology, 2021. 39(3): p. 313-319.
* Eng, C.-H.L., et al., Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH+. Nature, 2019. 568(7751): p. 235-239.


## Contact
For questions/comments, please contact Jiawen Chen (jiawenn@email.unc.edu). The processed MOB scRNA-seq external reference with all genes, mouse SSp slide-seqV2, mouse SSp external reference with all genes exceed the file size limit, which can be accessed upon request.





