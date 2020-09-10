# Tempus Bioinformatics Challenge
### Overview
The challenge is to construct a prototype a variant annotation tool that annotates vcf files. The tool accepts vcf files as input and outputs a tsv file with annotations. A jupyter notebook is included as an extended guide that explains the code with test data included. 

### Requirements
The tool is a python script, so Python 3.7 or greater must be installed.

The script utilizes the following libraries: pysam, json, and requests. 
Json and requests are native to python, so there is no need to install them. Requests will be used to access the Exome Aggregation Consortium (ExAC) API (http://exac.hms.harvard.edu/) in order to get data about the variants. The json package will parse data from the ExAC API.

Pysam can be installed with:
```{r eval=FALSE,echo=TRUE}
pip install pysam
```
or if you have Anaconda (https://anaconda.org/bioconda/pysam) installed:
```{r eval=False, echo=True}
conda install -c bioconda pysam
```
The pysam library is used to parse vcf files. More info here: https://pysam.readthedocs.io/en/latest/.
### Usage
Either of these commands will work on the command line:
```{r eval=FALSE,echo=TRUE}
./tempus_bioi_challenge.py -i /path/to/in.vcf -o /path/to/out.tsv
```
```{r eval=FALSE,echo=TRUE}
python3 tempus_bioi_challenge.py -i /path/to/in.vcf -o /path/to/out.tsv
```
### Input VCF File
The input file must be a VCF file. The following fields must also be included in the VCF file in addition to mandatory header:

**INFO FIELD**
- INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
- INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">
- INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">

**FORMAT FIELD**
- FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
- FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">

### Output TSV file
The TSV file contains the following annotations for each variant: 
- reference allele i.e. AT
- alternate allele(s) i.e. ATA
- variant type i.e. deletion
- effect i.e. missense 
- total read depth of all samples
- number of reads supporting variants
- percentage of reads supporting the variant
- allele frequency
- read depth, number of reads supporting variant, and percentage of reads supporting the variant for each sample
