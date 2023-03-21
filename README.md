# Phylogenomics

## merge_assembly_pieces.py
Genome sequences of some genomes, specially the Refseq genomes (starts with GCF_), are stored in muliple multiple-fasta files as shown below:

```
./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr1.fna    
./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr2.fna    
./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr3.fna    
./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr4.fna    
./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/unplaced.fna
```
Fragmented storage brings inconvenience to subsequent analysis. Therefore, I wrote this script to merge these fragmented genome sequences into one FASTA file.


**Input** is a tab-delimited file as shown below:
```
GCA_021398005.1    ./GCA_021398005.1/ncbi_dataset/data/GCA_021398005.1/GCA_021398005.fna
GCA_021436885.1    ./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr1.fna
```

**Output** is a FASTA file.   
``output_directory/GCA_021398005.1.fna``

**Recommended usage**    
``
find . -name '*.fna' | grep -v -e 'cds' -e 'rna' -e 'protein' | awk -F '/' '{print $2"\t"$0}' | merge_assembly_pieces.py - -o output_directory
``

