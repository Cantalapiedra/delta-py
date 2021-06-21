# delta-py
Python port of the R delta_statistic

Based on:
https://academic.oup.com/bioinformatics/article/35/11/1862/5144670
https://github.com/mrborges23/delta_statistic

It uses PastML MPPA for ACR:
https://academic.oup.com/mbe/article/36/9/2069/5498561
https://github.com/evolbioinfo/pastml

## Requirements

- python 3
- numpy
- pandas
- pastml (https://github.com/evolbioinfo/pastml)

## Usage

`python delta.py tree.nw data.csv "colname"`

tree.nw is a phylogenetic tree in newick format
data.csv is a comma-separated file, with a header, with first column the tree leaves IDs and next columns the data to be used for phylogenetic signal.
"colname" is the name of the column (or a comma-separated list of columns) to be used from data.csv to compute phylogenetic signal.

Other parameters:
CPUs: by default (0) uses all available CPUs.
p-value repetitions: number of repetitions to compute p-value. By default (0) p-value is not computed.
data separator: the field separator to parse data.csv. By default is ",".
ID index: the column of data.csv with the tree leaves IDs. By default (0), IDs are in the first column. It is 0-based, so position 3, for example, is ID index = 2.

Example using all parameters:

`python delta.py tree.nw data.csv "cluster" 10 100 "\t" 1`

## API usage

```
from delta import compute_delta_pvalue                                                                                                                  
                                                                                                                                                        
tree="../tree_trait_examples_claudia/cluster_138_tree.nw"                                                                                               
data="../tree_trait_examples_claudia/cluster_138_trait.csv"                                                                                             
columns=["cluster"]                                                                                                                                     
threads=10                                                                                                                                              
pval_reps=10                                                                                                                                            
data_sep=","                                                                                                                                            
id_index=0                                                                                                                                              
                                                                                                                                                        
delta, p_value = compute_delta_pvalue(tree, data, columns, threads, pval_reps, data_sep, id_index, verbose = True)                                      
                                                                                                                                                        
print(delta)                                                                                                                                            
print(p_value)
```

see [test_api.py](test_api.py)



