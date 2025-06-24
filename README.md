TopoHiC is a computational tool for Hi-C analysis. It enables simultaneous identification of chromatin loops and topologically associating domains (TADs).

## Download

First, you need to install `pydory` package, which enables high efficient computation of H1 loops in very large datasets. Use `pip install pydory` to install. If you need further assistance with the installation of `pydory`, please refer to [this page](https://github.com/nihcompmed/Dory). Note that it requires gcc-10.2 with openMP support. We recommend that you create a virtual environment prior to installation to avoid potential version conflicts.

After that, download all the python scripts in our repository to the same directory. 

## Data preparation

TopoHiC takes Hi-C files in `.cool` format as input. We recommend use Hi-C files at 5kb or even higher resolution. If you want to use Hi-C at higher resolution, such as 1kb, you might need to check the scripts and change the file names in `main.py` and `domain.py`. In case of 5kb resolution (default), we recommend that you rename your Hi-C files as `NAME-5kb.cool` format, such as `replicate1-5kb.cool` and `replicate2-5kb.cool`. You may put your Hi-C files under a PARENT_DIR, such as `hicfiles/`.

## Parameter settings

After the preparation of Hi-C files, you need to modify the parameter settings in `main.py`.

- `parant_dir` is the directory that contains your Hi-C, for example, `hicfiles/`
- `exp_list` is a list of NAME of your Hi-C files, for example, `['replicate1', 'replicate2']`
- `weight_col` is the name of the column that stores the weighted or raw contact frequencies which you want to use. If you use `cooler balance` to perform normalization with default parameters, then you need to set this to `'weight'`. In other cases, `'KR'`, `'VC'`, etc are all available.
- `chrom_list` is a list of all the names of chromosomes you want to calculate, for example, `['chr1', 'chr2']` or `['1', '2']`. Note that it must be consistent with the chromosome names in your Hi-C files.

