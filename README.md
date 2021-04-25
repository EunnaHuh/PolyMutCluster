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
* Arguments
`-In`  Please type input file directory <br/>
`-Err` Please select type of error either "SEM":standard error of the mean or "SD":standard deviatio
<hr style="border:2px solid gray"> </hr>

### Input File Formats
