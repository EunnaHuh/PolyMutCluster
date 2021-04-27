# PolyMutCluster
PolyMutCluster provides an analysis pipeline for quantifying the similarity between multi-dimensional protein signaling profiles. <br/>
This code is re-wrote from [Benredjem-Gallion](https://github.com/JonathanGallion/Benredjem-Gallion) to make it widely applicable to any mutagenesis and drug respone data and improve running speed by parallel computing. 
<hr style="border:2px solid gray"> </hr>

## Publication
Huh E, Gallion J, Lichtarge O. PolyMutCluster: Unsupervised Clustering of Multi-Dimensional Data in Signaling Profiles <br/>
~~(submitted)~~<br/><br/>
If you have any questions or comments, please feel free to reach out Eunna Huh (heo2399@gmail.com)
<hr style="border:2px solid gray"> </hr>

## Content
* [Download code](#Download-Code)
* [Installation](#Installation)
* [Run PolyMutCluster](#Run-PolyMutCluster)
* [Input File Formats](#Input-File-Formats)
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
`-In`  Directory of input file in csv format <br/>
`-Err`  Type of error "SEM":standard error of the mean or "SD":standard deviation <br/>
`-Rep`  Directory of experiment replication file in csv format when type of error is SEM. Default is 4 <br/>
`-Norm`  If you have wild type (WT) values, we recommand to type ND (Normalized Difference against WT). If you don't have WT, type MinMax which rescales each experiment (column) based on minimum and maximum values. <br/><br/>
Example:
```bash
python PolyMutCluster.py -In '../exp/Inputdata.csv' -Err 'SEM' -Rep '../exp/Replicate.csv' -Norm 'MinMax'
```
#### Options <br/>
`-p`  The number of CPU assigned. Default is all available CPU 
`-itr`  The nummber of error propagated matrix generated. Default is 500 <br/>
`-k`  The maximum number of K. Default is 3 <br/>
`-l`  Select linkage method [scipy.cluster.hierarchy.linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html) <br/>
`-pdist`  Select distance function [scipy.spatial.distance.pdist](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html) <br/>
`-log`  If the measurements are in log scale and you want to change them to linear scale, type the names of experiments <br/>
`-Save_suppl`  If you want to save all error propagated matrix please type 'ON' <br/> <br/>

<hr style="border:2px solid gray"> </hr>

### Input File Formats


