# ImmunoPepper
ImmunopPepper is a script which predicts the pepetide sequency from the RNA sequence splicegraph. You can specify the sample and output using the command line below.

```
--output_dir
./tests/test1
--ann_path
./tests/test1/data/test1pos.gtf
--ref_path
./tests/test1/data/test1pos.fa
--splice_path
./tests/test1/data/posgraph/spladder/genes_graph_conf3.merge_graphs.pickle
--count_path
./tests/test1/data/posgraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5
--maf_path
./tests/test1/data/test1pos.maf
--vcf_path
./tests/test1/data/test1pos.vcf
--samples
test1pos
--mutation_mode
somatic_and_germline
--debug
--filter_redundant
```
### Immunopepper Installation

It is recommended to setup a separate virtual or conda environment.

The basic immunopepper package can be installed from source using:
```
pip install -r requirements.txt -r requirements_dev.txt
make install
```

After that, consult the output of the following command for usage options
```
immunopepper -h
```
