# ![Title](images/Title.png)

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

## summary_BUSCO_results.py
When we are carrying out a big project including thousands of genomes, we should summarize all BUSCO results and further to remove bad assemblies. 
For this purpose, I wrote ``summary_BUSCO_results.py``.

**Recommended usage**    
``
find . -name 'short_summary.specific*txt' | summary_BUSCO_results.py - | sed 's/genome_label/Assembly/;s/.fna//' > busco_statistics.tsv
``

**Results** as shown below:
```
Assembly	C	S	D	F	M	Complete BUSCOS (C)	Complete and single-copy BUSCOs (S)	Complete and duplicated BUSCOs (D)	Fragemented BUSCOs (F)	Missing BUSCOs (M)   Total BUSCO groups searched
GCA_013390195.1	97.8%	97.5%	0.3%	0.2%	2.0%	1669	1664	5	3	34	1706
GCA_003568745.1	97.0%	96.3%	0.7%	0.8%	2.2%	1655	1643	12	13	38	1706
GCA_014705165.1	97.6%	97.1%	0.5%	0.4%	2.0%	1665	1656	9	6	35	1706
```

# single_copy_BUSCO_dataset.py

