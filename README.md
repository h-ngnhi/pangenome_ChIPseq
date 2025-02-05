# Analysis of ZNF146/OZF and ZNF507 targeting LINE-1 sequences, applied on Human Pangenome

## Resources

The analysis will be done following [Creamer et al., 2022](https://pubmed.ncbi.nlm.nih.gov/35100360/)'s analysis. We will be using the data from [ENCODE](https://www.encodeproject.org/), more specifically raw ChIP sequencing data of [ZNF146](https://www.encodeproject.org/experiments/ENCSR689YFA/) and [ZNF507](https://www.encodeproject.org/experiments/ENCSR598TIR/).  

[GenPipes](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_chipseq.html) tool is used to analyze the data on **GRCh37/hg19** reference genome.

## Project

**Ultimate goal**: New studies of transposable elements and the improvement of the pangenome.

**Project goal**: The same analysis will be redone by mapping the data on the new Human Pangenome and comparing to the results generated from the old genome reference.

**Proposed workflow**:  

1. Download raw ChIP-seq data from ENCODE.
2. Run Genpipes to map the data on **hg19** (paper) and **GRCh38** and later on **Pangenome** and view the peaks.
3. Interselect peaks with [RepeatMasker](https://www.repeatmasker.org/) ([Tutorial](https://www.clementgoubert.com/post/a-simple-pipeline-for-te-annotation-in-an-assembled-genome#:~:text=A%20simple%20pipeline%20to%20annotate%20transposable%20elements%20%28TEs%29,has%20homology%20to%20different%20consensus%20of%20the%20library.)) annotated repeats for different TE classes (to see if it agrees with the paper that ZNF146 and ZNF507 only associate with LINE-1).
4.


