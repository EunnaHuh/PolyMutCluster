# PolyMutCluster
PolyMutCluster provides an analysis pipeline for quantifying the similarity between multi-dimensional protein signaling profiles. It is designed to be widely applicable to any mutagenesis and drug respone data with improved speed using parallel computing.<br/>
This code is based on .<br/>
<hr style="border:2px solid gray"> </hr>

## Publication
Huh E, Gallion J, Lichtarge O. PolyMutCluster: Unsupervised Clustering of Multi-Dimensional Data in Signaling Profiles <br/>
(submitted)<br/><br/>
If you have any questions or comments, please feel free to reach out Eunna Huh (heo2399@gmail.com)
<hr style="border:2px solid gray"> </hr>

## Content
* [Download code](#Download-Code)
* [Installation](#Installation)
* [Run PolyMutCluster](#Run-PolyMutCluster)
* [File Formats](#Input-File-Formats)
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
`-In` Directory of input file in csv format <br/>
`-Err`  Type of error "SEM":standard error of the mean or "SD":standard deviation<br/>
`-Rep`  Directory of experiment replication file in csv format when type of error is SEM. Default is 4<br/>
`-Norm` If you have wild type (WT) values, we recommand to type ND (Normalized Difference against WT). If you don't have WT, type MinMax which rescales each experiment (column) based on minimum and maximum values<br/>
<br/>
Example:
```bash
python PolyMutCluster.py -In '../exp/Inputdata.csv' -Err 'SEM' -Rep '../exp/Replicate.csv' -Norm 'MinMax'
```
#### Options <br/>
`-p`  The number of CPU assigned. Default is to use all available CPU <br/>
`-itr`  The nummber of error propagated matrix generated. Default is 500 <br/>
`-k`  The maximum number of K. Default is 3 <br/>
`-l`  Select linkage method [scipy.cluster.hierarchy.linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html) <br/>
`-pdist`  Select distance function [scipy.spatial.distance.pdist](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html) <br/>
`-log`  If the measurements are in log scale and you want to change them into linear scale, type the names of experiments <br/>
`-Save_suppl`  If you want to save all error-propagated matrices, type 'ON' <br/> 
<br/>
Example:
```bash
python PolyMutCluster.py -In <assay readouts file directory> -Err 'SEM' -Rep < replicate file directory> -Norm 'MinMax' -p 4 -iter 1000 -k 7 -l 'ward' -pdist 'euclidean' -log 'EC50' 'tka'  -Save_suppl 'ON'
```

<hr style="border:2px solid gray"> </hr>

### File Formats
* The required information in Inputdata.csv: 
  - Rows: List of mutations
  - Columns: List of experiments. Column displaying average is written as the name of assay and parameter, seperated by underscore "_". For instance "cAMP_Emax" or "cAMP_EC50". Column diaplaying either standard deviation (SD) or standard error of the mean (SEM) is written after the names of assay and parameter. For instance "cAMP_Emax_SEM" or "cAMP_EC50_SD"
  - Values: Average and standard deviation(SE)/standard error of the mean(SEM) from multiple assay readouts.

For reference, see `tutorial/data/inputdata.csv` <br/> <br/>

* The required information in Replicate.csv :
  - when standard error of the mean (SEM) is inputted as an experimental error, PolyMutCluster need information on the number of repetition performed for each assay to convert standard error of the mean (SEM) tostandard deviation (SD).
  - when standard deviation (SD) is inputted as an experimental error, you don't need to create Replicate.csv 
  - Rows: List of mutations
  - Columns: List of experiments. Column diaplaying standard error of the mean (SEM) is written after the names of assay and parameter. For instance "cAMP_Emax_SEM" or "cAMP_EC50_SEM"
  - Values: The number of repetition conducted for a experiment

For reference, see `tutorial/data/Replicate.csv` <br/> <br/>
<hr style="border:2px solid gray"> </hr>

### Acknowledgement
[Benredjem-Gallion](https://github.com/JonathanGallion/Benredjem-Gallion) Original dvelopment for PolyMutCluster


