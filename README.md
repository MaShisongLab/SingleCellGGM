# SingleCellGGM

**An algorithm for single-cell gene co-expression network analysis**

<b>SingleCellGGM</b> (single-cell graphical Gaussian model) is a MATLAB algorithm for single-cell gene co-expression network analysis based on the graphical Gaussian model (GGM). The algorithm is modified from a previously published method ([Ma *et al*, 2007](#References)) and a previous MATLAB algorithm [rgsGGM](https://github.com/MaShisongLab/rgsGGM) ([Wang *et al*, 2023](#References)), both for conducting GGM co-expression network analysis on bulk transcriptome dataset.

SingleCellGGM takes a single-cell gene expression matrix as input and uses a process consisting of ~20,000 iterations to calculate partial correlation coefficients (<i>pcors</i>) between gene pairs. For each gene pair, SingleCellGGM also calculates how many cells the gene pair are co-expressed in. SingleCellGGM then takes the gene pairs with <i>pcor</i> >= a selected cutoff value and co-expressed in >= a selected number of cells to construct a GGM gene co-expression network. For more details, please refer to [Xu et al, 2023](#References). 

## Table of Contents
- [Install](#Install)
- [Usage](#Usage)
- [References](#References)

## Install
This algorithm requires [MATLAB](https://www.mathworks.com/products/matlab.html). Copy the file `SingleCellGGM.m` to your working folder or the MATLAB scripts folder and start using it.

## Usage

The algorithm takes a log-normalized gene expression matrix, the number of iterations, the names of the genes, and the name of dataset as inputs. The expression matrix should have samples in rows and genes in columns. The sample numbers should be large and the low-expression genes should be filtered out first. 

<B>ggm = SingleCellGGM(`expression_matrix`,`number_of_iterations`,`gene_names`,`dataset_name`)</B>

`expression_matrix` - the gene expression matrix, samples in rows, genes in columns <br/>
`number_of_iterations` - the number of iterations used for *pcor* calculation, usually 20000<br/>
`gene_names` - the names for the genes in the matrix <br/>
`dataset_name` - the name of the dataset


Below, we use a mouse single-cell gene expression matrix obtained from the MCA project ([Han *et al*, 2018](#References)) as an example to demonstrate how to conduct single-cell GGM gene co-expression network analysis via SingleCellGGM. The matrix file "MCA_Figure2-batch-removed.txt.tar.gz" can be downloaded from [Figureshare](https://figshare.com/ndownloader/files/10351110?private_link=865e694ad06d5857db4b) as provided by MCA. Unzip and place the file "Figure2-batch-removed.txt" into the MATLAB working folder. We also obtained the Ensembl gene IDs for the genes within the matrix and saved it in a file "data/MCA.ensembl.gene.ids.txt".  

```matlab
% MATLAB code
% Read in the single-cell gene expression matrix
expression_matrix = readtable('Figure2-batch-removed.txt','ReadRowNames',true);
expression_matrix = table2array(expression_matrix);

% Log normalization
expression_matrix = log2( expression_matrix ./ sum( expression_matrix ) * 10000 + 1 );
expression_matrix = expression_matrix';

% Read in gene names
gene = readcell('data/MCA.ensembl.gene.ids.txt');

% Out of 25133 genes in the matrix, 24802 have Ensembl gene IDs. Only the 
% genes with Ensembl gene IDs will be used for network construction.
idx = contains(gene,'ENSMUSG');

expression_matrix = expression_matrix(:,idx);
gene = gene(idx);

% Conduct single-cell gene co-expression analysis via SingleCellGGM
ggm = SingleCellGGM( expression_matrix, 20000, gene, 'mca')

% Examine the results
ggm
ggm.SigEdges(1:5,:)

% Save all gene pairs to a file for gene co-expression construction
writetable(ggm.SigEdges(:,1:3),'mca.ggm.network.txt','Delimiter','tab','WriteVariableNames',FALSE)
```

Here is a glimpse into the results:
```shell
ggm

ggm =

  SingleCellGGM with properties:

                gene_num: 24802
               gene_name: {24802x1 cell}
                pcor_all: [24802x24802 single]
       pcor_sampling_num: [24802x24802 int16]
    coexpressed_cell_num: [24802x24802 int32]
             samples_num: 61637
             RoundNumber: 20000
                SigEdges: [127229x9 table]


ggm.SigEdges(1:5,:)

ans =

  10x9 table

           GeneA                   GeneB              Pcor      SamplingTime       r       Cell_num_A    Cell_num_B    Cell_num_coexpressed    Dataset
    ____________________    ____________________    ________    ____________    _______    __________    __________    ____________________    _______

    'ENSMUSG00000109644'    'ENSMUSG00000057228'    0.033203        132         0.32947       660            406               191              'mca'
    'ENSMUSG00000109644'    'ENSMUSG00000021228'    0.042609        143         0.20348       660            149                71              'mca'
    'ENSMUSG00000109644'    'ENSMUSG00000021751'    0.030988        128         0.17124       660            162                61              'mca'
    'ENSMUSG00000109644'    'ENSMUSG00000089678'    0.038486        141         0.20272       660            120                65              'mca'
    'ENSMUSG00000109644'    'ENSMUSG00000109311'    0.042568        136         0.45847       660            692               321              'mca'
```

The network can then be clustered via network clustering algorithm to obtain gene co-expression module and used for down-stream analysis.

## References

Han X, Wang R, Zhou Y, Fei L, Sun H, Lai S, Saadatpour A, Zhou Z, Chen H, Ye F, Huang D, Xu Y, Huang W, Jiang M, Jiang X, Mao J, Chen Y, Lu C, Xie J, Fang Q, Wang Y, Yue R, Li T, Huang H, Orkin SH, Yuan GC, Chen M, and Guo G. 2018. Mapping the Mouse Cell Atlas by Microwell-Seq. *Cell* 172: 1091-1107.

Ma S, Gong Q, and Bohnert HJ. 2007. An Arabidopsis gene network based on the graphical Gaussian model. *Genome Research* 17:1614-1625.

Xu Y, Wang Y, and Ma S. 2023. SingleCellGGM enables gene expression program identification from single-cell transcriptomes and facilitates universal cell label transfer. *submitted*

Wang Y, Zhang Y, Yu N, Li B, Gong J, Mei Y, Bao J, and Ma S. 2023. Decoding transcriptional regulation via a human gene expression predictor. *Journal of Genetics and Genomics* https://doi.org/10.1016/j.jgg.2023.01.006
 

