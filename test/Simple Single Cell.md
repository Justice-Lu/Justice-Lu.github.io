---
layout: default
title: Simple Single Cell 
nav_order: 2
---

## Simple Single Cell 
{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>


## Single-cell RNA sequencing Data (scRNA-seq)

Single-cell sequencing is a powerful and efficient technology that profiles the whole transcriptome of a large pool of cells. And as a result, such data is rich and may be manipulated to ask a variety of questions via different computational methods. 

A single-cell data post-processed is usually read in with a metadata that describes the samples and individual cells.
	
<!--
orig.ident  nCount_RNA  nFeature_RNA specie seurat_clusters cluster condition
GGTATTAATCTC_0105	Moe_0105	2859.618	1725	mouse	1	mOSNs	Mouse
GTCCTAGGTAGA_0105	Moe_0105	3261.495	2171	mouse	4	mOSNs/SUS	Mouse
AACAGCGTAAGC_0105	Moe_0105	3126.335	1960	mouse	24	iOSNs	Mouse
TCATCAGAGATC_0105	Moe_0105	2193.834	960	mouse	1	mOSNs	Mouse
CAGCTGATGTCC_0105	Moe_0105	1508.292	538	mouse	1	mOSNs	Mouse
AGGCTCATGACA_0105	Moe_0105	1707.622	650	mouse	4	mOSNs/SUS	Mouse
GGTTATTCGTGG_0105	Moe_0105	1784.866	701	mouse	0	mOSNs	Mouse
ACCTGAAGGAGT_0105	Moe_0105	2607.098	1298	mouse	12	iOSNs	Mouse
CAAGAGATAGCC_0105	Moe_0105	2219.591	1066	mouse	1	mOSNs	Mouse
TGCTGCTCCCAG_0105	Moe_0105	2481.906	1261	mouse	2	iOSNs	Mouse
-->

And the actual single cell data is simply a data table of a all the labeled cell barcodes on y axis, with the corresponding gene expression counts on the x axis. 

One of the most straight forward way to visualize single cell data is through UMAP. Where each cells can be represented on a series of principal component space projected in 2D. Where the closer the points are to each other, the more similar their gene expression profile are to each other and vice versa. 

And through known gene markers and clustering, we can predict certain clusters to represent a certain cell types. 

![](/resource/1_Single_Cell/SIF_Umap_cluster.png)


## Cluster composition 

For this example the scRNA-seq data is conducted on mouse, human, and an covid-patient. And one of the initial comparison we can look at is how does the cells from each of those samples compare to each other in terms of cluster composition and expression level. To do that we can simply look at the total counts of a cluster, and the percentage of the cells the cluster made up of using a dot plot. 

![](/resource/1_Single_Cell/SIF_percent_table.png)

Dotplot is extremely a very useful visualization tool especially at comparing multiple dimensions between categories. However, keep in mind that it may get messy and less intuitive when not applied approriately. 

























