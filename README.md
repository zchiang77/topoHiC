TopoHiC is a computational tool for Hi-C analysis. It enables simultaneous identification of chromatin loops and topologically associating domains (TADs).

## Download

First, you need to install `pydory` package, which enables high efficient computation of H1 loops in very large datasets. Click [here](https://github.com/nihcompmed/Dory) to move to the github page for installation help. Note that it requires gcc-10.2 with openMP support. 

We recommend that you create a virtual environment prior to installation to avoid potential version conflicts.

After that, download all the python scripts in our repository to the same directory. 

## Data preparation

TopoHiC takes Hi-C files in `.cool` format as input.


, and then modify the parameter settings in `main.py`.

`parant_dir` is the directory that contains your Hi-C 
