# SingleCellGGM

**An algorithm for single-cell gene co-expression network analysis**

<b>SingleCellGGM</b> (single-cell graphical Gaussian model) is a MATLAB algorithm for single-cell gene co-expression network analysis based on the graphical Gaussian model (GGM). The algorithm is modified from a previously published method ([Ma *et al*, 2007](#References)) and a previous MATLAB algorithm [rgsGGM](https://github.com/MaShisongLab/rgsGGM) ([Wang *et al*, 2023](#References)), both for conducting GGM co-expression network analysis on bulk transcriptome datasets.

SingleCellGGM takes a single-cell gene expression matrix as input and uses a process consisting of ~20,000 iterations to calculate partial correlation coefficients (<i>pcors</i>) between gene pairs. For each gene pair, SingleCellGGM also calculates in how many cells the gene pair are co-expressed. SingleCellGGM then takes the gene pairs with <i>pcor</i> >= a selected cutoff value and co-expressed in >= a selected number of cells to construct a GGM gene co-expression network. For more details, please refer to [Xu et al, 2024](#References). 

## Table of Contents
- [Install](#Install)
- [Usage](#Usage)
- [Tutorial](#Tutorial)
- [References](#References)

## Install
This algorithm requires [MATLAB](https://www.mathworks.com/products/matlab.html). Copy the files `SingleCellGGM.m`, `fdr_control.m`, and `adjust_cutoff.m` to your working folder or the MATLAB scripts folder to begin utilization.

## Usage

The algorithm takes a log-normalized gene expression matrix, the number of iterations, the names of the genes, and optionally, the name of dataset as inputs. The expression matrix should have cells in rows and genes in columns. The cell numbers should be large (>5000 recommended) and low-expression genes should be filtered out first. 

<B>ggm = SingleCellGGM(`expression_matrix`,`number_of_iterations`,`gene_names`,`dataset_name`)</B>

`expression_matrix` - gene expression matrix (double format), cells in rows, genes in columns<br/>
`number_of_iterations` - number of iterations used for *pcor* calculation, typically set to 20000<br/>
`gene_names` - names of genes in the matrix <br/>
`dataset_name` (optional) - name of the dataset (default is 'na')<br><br>
An alternative for <i>number_of_iterations</i> is calculated as <i>round(gene_number * (gene_number - 1) / 39980)</i>, ensuring that each gene pair is sampled, on average, in 100 iterations.

<B>fdr = fdr_control(`expression_matrix`,`ggm`,`permutation_fraction`)</B><br>

`expression_matrix` - gene expression matrix
<br>`ggm` - result from SingleCellGGM function
<br>`permutation_fraction` - fraction of genes to be permutated (default is 1)

<B>ggm_adjusted = adjust_cutoff(`ggm`,`pcor_cutoff`,`coexpressed_cell_cutoff`)</B>

`ggm` - result from SignleCellGGM function
<br>`pcor_cutoff` - cutoff for pcor to be adjusted (default is 0.03)
<br>`coexpressed_cell_cutoff` - cutoff for number of coexpressed cells to be adjusted (default is 10)

Below, we use a mouse single-cell gene expression matrix obtained from the MCA project ([Han *et al*, 2018](#References)) as an example to demonstrate how to conduct single-cell GGM gene co-expression network analysis via SingleCellGGM. The matrix file "MCA_Figure2-batch-removed.txt.tar.gz" can be downloaded from [Figureshare](https://figshare.com/ndownloader/files/10351110?private_link=865e694ad06d5857db4b) as provided by MCA. Unzip and place the file "Figure2-batch-removed.txt" into the MATLAB working folder. We also obtained the Ensembl gene IDs for the genes within the matrix and saved it in a file "data/MCA.ensembl.gene.ids.txt".  

```matlab
% MATLAB code
% Read in the single-cell gene expression matrix
expression_matrix = readtable('Figure2-batch-removed.txt','ReadRowNames',true);
expression_matrix = table2array(expression_matrix);

% Ensure the matrix is in double format
expression_matrix = double(expression_matrix);

% Log normalization
expression_matrix = log2( expression_matrix ./ sum( expression_matrix ) * 10000 + 1 );
expression_matrix = expression_matrix';

% Read in gene names
gene = readcell('data/MCA.ensembl.gene.ids.txt');

% Out of 25133 genes in the matrix, 24802 have Ensembl gene IDs. Only the 
% genes with Ensembl gene IDs will be used for network construction
idx = contains(gene,'ENSMUSG');

expression_matrix = expression_matrix(:,idx);
gene = gene(idx);

% Filter out genes expressed in less than 10 cells
idx_expression_filter = sum(expression_matrix > 0) >= 10;
expression_matrix = expression_matrix(:,idx_expression_filter);
gene = gene(idx_expression_filter);

% Convert to sparse matrix format to save memory
expression_matrix = sparse( expression_matrix);

% Conduct single-cell gene co-expression analysis via SingleCellGGM
ggm = SingleCellGGM( expression_matrix, 20000, gene, 'mca')

% Examine the results
ggm
ggm.SigEdges(1:5,:)

% Save all gene pairs to a file for gene co-expression network construction
writetable(ggm.SigEdges(:,1:3),'mca.ggm.network.txt','Delimiter','tab','WriteVariableNames',false)

% FDR Control - control the false discovery rate with different pcor cutoff value
fdr = fdr_control(expression_matrix, ggm)

% Inspect FDR rate
fdr
fdr.fdr

% Adjust pcor cutoff to 0.02
ggm_adjusted = adjust_cutoff(ggm, 0.02)

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
             DatasetName: 'mca'


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
```shell
fdr

fdr = 

  fdr_control with properties:

       gene_num: 24802
      gene_name: {24802x1 cell}
    samples_num: 61637
    RoundNumber: 20000
       SigEdges: [1185x9 table]
            fdr: [91x5 table]
    DatasetName: 'mca'

fdr.fdr

ans =

  91x5 table

    Pcor     SigEdgeNum        FDR         SigPermutatedEdgeNum    SigPermutatedEdgeProportion
    _____    __________    ____________    ____________________    ___________________________

     0.01     1645160      '     0.283'           464795                     0.00151          
    0.011     1279237      '     0.203'           259572                    0.000844          
    0.012     1019825      '     0.141'           143319                    0.000466          
    0.013      832700      '    0.0939'            78154                    0.000254          
    0.014      694468      '    0.0613'            42545                    0.000138          
    0.015      589634      '    0.0391'            23077                     7.5e-05          
    0.016      508505      '    0.0246'            12532                    4.07e-05          
    0.017      444704      '    0.0153'             6783                    2.21e-05          
    0.018      392243      '   0.00962'             3775                    1.23e-05          
    0.019      348797      '   0.00604'             2108                    6.85e-06          
     0.02      312383      '   0.00379'             1185                    3.85e-06          
    0.021      281495      '   0.00244'              687                    2.23e-06          
    0.022      254676      '   0.00147'              374                    1.22e-06          
    0.023      231776      '  0.000962'              223                    7.25e-07          
    0.024      211169      '  0.000578'              122                    3.97e-07          
    0.025      193134      '  0.000373'               72                    2.34e-07          
    0.026      176888      '  0.000232'               41                    1.33e-07          
    0.027      162583      '  0.000148'               24                     7.8e-08          
    0.028      149481      '    0.0001'               15                    4.88e-08          
    0.029      137858      '  7.98e-05'               11                    3.58e-08          
     0.03      127229      '  4.72e-05'                6                    1.95e-08          
    0.031      117748      '   3.4e-05'                4                     1.3e-08          
    0.032      109148      '  2.75e-05'                3                    9.75e-09          
    0.033      101235      '  2.96e-05'                3                    9.75e-09          
    0.034       94054      '  1.06e-05'                1                    3.25e-09          
    0.035       87357      '  1.14e-05'                1                    3.25e-09          
    0.036       81282      '  1.23e-05'                1                    3.25e-09          
    0.037       75650      '  1.32e-05'                1                    3.25e-09          
    0.038       70563      '  1.42e-05'                1                    3.25e-09          
    0.039       65973      '  1.52e-05'                1                    3.25e-09          
     0.04       61606      '< 1.52e-05'                0                           0          
    0.041       57559      '< 1.52e-05'                0                           0          
    0.042       53880      '< 1.52e-05'                0                           0          
    0.043       50509      '< 1.52e-05'                0                           0          
    0.044       47420      '< 1.52e-05'                0                           0          
    0.045       44399      '< 1.52e-05'                0                           0                   

```
```shell
ggm_adjusted

ggm_adjusted = 

  adjust_cutoff with properties:

                gene_num: 24802
               gene_name: {24802x1 cell}
                pcor_all: [24802x24802 single]
       pcor_sampling_num: [24802x24802 int16]
                     PCC: [24802x24802 single]
    coexpressed_cell_num: [24802x24802 int32]
             samples_num: 61637
             RoundNumber: 20000
                SigEdges: [312383x9 table]
             DatasetName: 'mca'

```
The network can then be clustered via network clustering algorithm to obtain gene co-expression module and used for down-stream analysis.

## Tutorial

Check out [here](https://github.com/MaShisongLab/SingleCellGGM_Network_Analysis_Tutorial) for a tutorial on SingleCellGGM network analysis and downstream co-expression module analysis.

## References

Han X, Wang R, Zhou Y, Fei L, Sun H, Lai S, Saadatpour A, Zhou Z, Chen H, Ye F, Huang D, Xu Y, Huang W, Jiang M, Jiang X, Mao J, Chen Y, Lu C, Xie J, Fang Q, Wang Y, Yue R, Li T, Huang H, Orkin SH, Yuan GC, Chen M, and Guo G. 2018. Mapping the Mouse Cell Atlas by Microwell-Seq. *Cell* 172: 1091-1107.

Ma S, Gong Q, and Bohnert HJ. 2007. An Arabidopsis gene network based on the graphical Gaussian model. *Genome Research* 17:1614-1625.

Xu Y, Wang Y, and Ma S. 2024. SingleCellGGM enables gene expression program identification from single-cell transcriptomes and facilitates universal cell label transfer. *submitted*

Wang Y, Zhang Y, Yu N, Li B, Gong J, Mei Y, Bao J, and Ma S. 2023. Decoding transcriptional regulation via a human gene expression predictor. *Journal of Genetics and Genomics* 50:305-317.

