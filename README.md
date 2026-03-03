TopoHiC is a computational tool for Hi-C analysis. It enables simultaneous identification of chromatin loops and topologically associating domains (TADs).

## Download

First, you need to install `pydory` package, which enables high efficient computation of H1 loops in very large datasets. Use `pip install pydory` to install. If you need further assistance with the installation of `pydory`, please refer to [Dory](https://github.com/nihcompmed/Dory). Note that it requires gcc-10.2 with openMP support. We recommend that you create a virtual environment prior to installation to avoid potential version conflicts.

Then, the following software packages are required:

- python = 3.12.2
- cooler = 0.10.3
- numpy = 1.26.4
- pandas = 2.2.1
- h5py = 3.13.0
- matplotlib = 3.9.1
- networkx = 3.4.2
- pydot = 1.4.2
- sklearn = 1.3.2
- scipy = 1.14.0
- numba = 0.59.1

This pipeline also requires the import of the following packages: `pickle`, `itertools`, `copy`, `gc`, `collections`, `math`, `subprocess`. These packages are part of the standard Python library and require no additional installation.

To get started, download all the python scripts in our repository to the same directory. 

## Data preparation

TopoHiC takes Hi-C files in `.cool` format as input. We recommend use Hi-C files at 5kb or even higher resolution. If you want to use Hi-C at higher resolution, such as 1kb, you might need to check the scripts and change the file names in `main.py` and `domain.py`. In case of 5kb resolution (default), we recommend that you rename your Hi-C files as `NAME-5kb.cool` format, such as `replicate1-5kb.cool` and `replicate2-5kb.cool`. You may put your Hi-C files under a PARENT_DIR, such as `hicfiles/`.

## Parameter settings

After the preparation of Hi-C files, you need to modify the parameter settings in `main.py`.

- `parent_dir` is the directory that contains your Hi-C, for example, `hicfiles/`
- `exp_list` is a list of NAME of your Hi-C files, for example, `['replicate1', 'replicate2']`
- `weight_col` is the name of the column that stores the weighted or raw contact frequencies which you want to use. If you use `cooler balance` to perform normalization with default parameters, then you need to set this to `'weight'`. In other cases, `'KR'`, `'VC'`, etc are all available.
- `chrom_list` is a list of all the names of chromosomes you want to calculate, for example, `['chr1', 'chr2']` or `['1', '2']`. Note that it must be consistent with the chromosome names in your Hi-C files.
- `resolution` is the resolution of your Hi-C files, such as `5000`
- `edge_limit_bp` is the maximum distance in the 3D point cloud. Points that are separated by a distance exceeding this limit are regarded as being infinitely distant from each other. You can set this value to `-1` in order to consider all distances.

We advise against modifying other parameters unless you possess a comprehensive understanding of persistent homology.

## Calculate loop scores

Use command `python main.py` in order to calculate loop scores for each Hi-C file. The results will be stored by default under `parent_dir/topohic` directory.

## Find TADs

First, modify the parameter settings, including `parent_dir`, `exp_list`, `chrom_list` and `resolution`, in `domain.py` similar to `main.py`. The parameter `domain_limit` is set to be 50000 by default which means that candidate TADs whose sizes are below 50kb will be filtered out.

By default, `pair_merge` is set to be `True` which means topoHiC will find TADs using the merged approach, in which case the `score_threshold` is set to be 0.5. However, if you want to implement topoHiC using the individual approach, you need to change `pair_merge` to `False` and `score_threshold` to 0.

In merged approach, the resulting TADs will be stored under `parent_dir/topohic/merge` directory. Otherwise, they will be stored under separate directories.

## Example

In the following steps, we will walk you through downloading and processing the example dataset from the GEO database.

- **Step1:** Download the folder named `example`​ to local. Then copy all the scripts you have just downloaded into the `example/topohic/​` directory, placing them together with `run1.py`​ and `domain.py`.
- **Step2:** Click this link [GSM1551569](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1551569&format=file&file=GSM1551569%5FHIC020%2Ehic) to download the file with the .hic suffix to your local device. This is the biological replicate 1 in Dataset1 in our paper. You may rename it to `rep1.hic` for convenience. Then, install [HiCExplorer](https://github.com/deeptools/HiCExplorer). Run `hicConvertFormat` command to convert the `.hic` format to `.cool`. For example,
```
hicConvertFormat -m rep1.hic -o rep1.cool --inputFormat hic --outputFormat cool -r 5000
```
The output file should be named as `rep1_5000.cool`. You can create a folder, for instance, "Rao2014", and put the cooler file in the new folder. After that, run `cooler balance` to perform ICE normalization. Remember to use the -p parameter to specify the core number used. For example,
```
cooler balance --ignore-diags 0 -p 64 --max-iters 500 -f rep1_5000.cool
```

- **Step3:** We assume that at this step, your data is stored in `example/Rao2014/`, and your Python scripts are located in `example/topohic/`. Now, navigate back to the script directory and open `run1.py`. You might want to modify the following values in the script:
```
parent_dir = "../Rao2014/"
exp_list = ['rep1']
chrom_list = ['22']
threads = 64
```

The `parent_dir` is the folder in which you store your data. The `exp_list` contains the labels of your datasets. The `chrom_list` contains the chromosomes on which you wish to apply the pipeline. The `threads` parameter specifies the number of CPU cores the pipeline will use. After that, run `python3 run1.py`. The pipeline will compute H1 loops and compute loop scores. By default, it will create a new folder in `example/Rao2014/topohic/` to store all the results. 

You can find the computed loops in `example/Rao2014/topohic/rep1/topohic-22/rep1-sc4/scratch.txt`. To process the result file, we provide a Python script named `topohicMultiple.py`. Through straightforward text processing, it converts the results in `scratch.txt`​ into a structured format for easier reading and analysis. This script will convert the genomic bins to chromosome coordinates; therefore it requires a reference file named as `binchr5.txt` which can be obtained through 
```
cooler dump --header -t bins rep1-5kb.cool > binchr5.txt
```
The processed results will be output to a folder named `topohicMultiple​` in the same directory as the script itself by default. Please ensure this folder is created. The code is relatively simple, users can easily customize the input and output settings as needed. (Please refer to [Dory](https://github.com/nihcompmed/Dory) and its original paper for detailed explanation of the results of persistent homology)

- **Step4:** If the chromatin loop is all you want, then you may stop here. Since we here only have one replicate, we set `pair_wise=False` in the `domain.py` to extract TADs through individual approach. You may want to modify the following values in `domain.py` as we did before in the previous step:
```
parent_dir = "../Rao2014/"
exp_list = ['rep1']
chrom_list = ['22']
```

After that, run `python3 domain.py`. The computed TADs are in `example/Rao2014/topohic/rep1/topohic-22/rep1-sc4/domain.csv`. The output CSV file contains five columns. **Note that the last two columns correspond to features that are currently under development and not yet complete; please ignore and do not use them.**

## Help

If you need any help, feel free to contact zhenchiang@sjtu.edu.cn.

