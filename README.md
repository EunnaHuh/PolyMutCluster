# PolyMutCluster
PolyMutCluster provides an analysis pipeline for quantifying the similarity between multi-dimensional protein signaling profiles. It is designed to be widely applicable to any mutagenesis and drug respone data with improved speed using parallel computing.<br/>
<hr style="border:2px solid gray"> </hr>

## Publication
Huh E, Gallion J, Lichtarge O. PolyMutCluster: Unsupervised Clustering of Multi-Dimensional Data in Signaling Profiles <br/>
(submitted)<br/><br/>
If you have any questions or comments, please feel free to contact Eunna Huh (heo2399@gmail.com)
<hr style="border:2px solid gray"> </hr>

## Content
* [Download code](#Download-Code)
* [Installation](#Installation)
* [Run PolyMutCluster](#Run-PolyMutCluster)
* [File Formats](#File-Formats)
* [Acknowledgement](#Acknowledgement)
<hr style="border:2px solid gray"> </hr>

### Installation
#### Requirements
* [Anaconda](https://docs.anaconda.com/anaconda/install/) or [MiniConda](https://docs.conda.io/en/latest/miniconda.html)
* python >= 3.7

#### Setup
```bash
conda create -n PolyMutCluster python=3.7.1
source activate PolyMutCluster
pip install -r env.txt
```
<hr style="border:2px solid gray"> </hr>

### Run PolyMutCluster
```bash
source activate PolyMutCluster
python PolyMutCluster.py 
```
#### Arguments <br/>
`-In` Directory of assay readouts file in csv format <br/>
`-Err`  Error type. SEM:standard error of the mean or SD: standard deviation<br/>
`-Rep`  Directory of experiment replication file in csv format when error type is SEM. Default is 4 <br/>
`-Norm` If you have wild type (WT) values, we recommand to type ND (Normalized Difference against WT). If you don't have WT, type MinMax which rescales each experiment (column) based on minimum and maximum values<br/>
<br/>
Example 1:
```bash
python PolyMutCluster.py -In <assay readouts file directory> -Err 'SEM' -Rep <replicate file directory> -Norm 'MinMax'
```
Example 2:
```bash
python PolyMutCluster.py -In <assay readouts file directory> -Err 'SD'  -Norm 'MinMax'
```

#### Options <br/>
`-p`  The number of CPU to be used. Default is to assign all available CPUs <br/>
`-itr`  The nummber of error-propagated matrix generated. Default is 500 <br/>
`-k`  The maximum number of K. Default is 3 <br/>
`-l`  Select linkage method [scipy.cluster.hierarchy.linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html). Default is average <br/>
`-pdist`  A measurement for a pairwise correlation of cluster frequency matrix. Options: braycurtis, canberra, chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard, jensenshannon, kulsinski, mahalanobis, matching, minkowski, rogerstanimoto, russellrao, seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule, pearson, kendall, spearman. Default is euclidean <br/>
`-log`  If you want to change scale of measurement from log to linear in the process of normalization, type the names of experiments <br/>
`-S`  If you want to save every normalized matrices, type 'ON' <br/> 
<br/>
Example:
```bash
python PolyMutCluster.py -In <assay readouts file directory> -Err 'SEM' -Rep <replicate file directory> -Norm 'ND' -p 3 -itr 1000 -k 7 -l 'ward' -pdist 'euclidean' -log 'EC50' 'tka'  -S 'ON'
```

<hr style="border:2px solid gray"> </hr>

### File Formats
* The required information on Inputdata.csv: 
  - Rows: List of mutations
  - Columns: List of experiments. Column displaying average is written as the name of assay and parameter, seperated by underscore "_". For instance, "cAMP_Emax" or "cAMP_EC50". Column diaplaying either standard deviation (SD) or standard error of the mean (SEM) is written after the names of assay and parameter. For instance "cAMP_Emax_SEM" or "cAMP_EC50_SD"
  - Values: Average and standard deviation(SE)/standard error of the mean(SEM) from multiple assay readouts.

For reference, see `tutorial/data/inputdata.csv` <br/> <br/>

* The required information on Replicate.csv :
  - When standard error of the mean (SEM) is inputted as an experimental error, PolyMutCluster need information on the number of repetition performed for each assay to convert standard error of the mean (SEM) to standard deviation (SD).
  - when standard deviation (SD) is inputted as an experimental error, you don't need to create Replicate.csv 
  - Rows: List of mutations
  - Columns: List of experiments. Column diaplaying standard error of the mean (SEM) is written after the names of assay and parameter. For instance "cAMP_Emax_SEM" or "cAMP_EC50_SEM"
  - Values: The number of repetition conducted per experiment

For reference, see `tutorial/data/Replicate.csv` <br/> <br/>
<hr style="border:2px solid gray"> </hr>

### Acknowledgement
[Benredjem-Gallion](https://github.com/JonathanGallion/Benredjem-Gallion) Original dvelopment for PolyMutCluster


