# inTRINSiC

Glioblastoma multiforme (GBM) is a highly malignant form of cancer that lacks effective treatment options or well‐defined strategies for personalized cancer therapy. The disease has been stratified into distinct molecular subtypes; however, the underlying regulatory circuitry that gives rise to such heterogeneity and its implications for therapy remain unclear. We developed a modular computational pipeline, Integrative Modeling of Transcription Regulatory Interactions for Systematic Inference of Susceptibility in Cancer (inTRINSiC), to dissect subtype‐specific regulatory programs and predict genetic dependencies in individual patient tumors. Using a multilayer network consisting of 518 transcription factors (TFs), 10,733 target genes, and a signaling layer of 3,132 proteins, we were able to accurately identify differential regulatory activity of TFs that shape subtype‐specific expression landscapes. Our models also allowed inference of mechanisms for altered TF behavior in different GBM subtypes. Most importantly, we were able to use the multilayer models to perform an in silico perturbation analysis to infer differential genetic vulnerabilities across GBM subtypes and pinpoint the MYB family member MYBL2 as a drug target specific for the Proneural subtype.

This repository holds source code for the regression modeling in the paper published in Molecular Systems Biology:
Integrated regulatory models for inference of subtype‐specific susceptibilities in glioblastoma
https://www.embopress.org/doi/full/10.15252/msb.20209506

Please refer to regressionPipeline.R for a step-by-step guide through the regression pipeline that assigns magnitude and directionality to edges in the input backbone transcription regulatory network.

The input data you need for this pipeline are:

1. Gene expression matrix (CSV/TSV format with row and column names)
2. Binary transcription regulatory network edge list: first column contains target genes and second column contains the corresponding TFs

The pipeline outputs weighted adjacency matrices for the models built from each subtype/cell type. The weights correspond to the F values for the regulatory edges.

We will update this repository soon with the source code for the next steps in the inTRINSiC pipeline: protein activity inference and gene essentiality prediction.
