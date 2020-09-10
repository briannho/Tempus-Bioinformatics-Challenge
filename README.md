# Tempus Bioinformatics Challenge
The challenge is to construct a prototype a variant annotation tool that annotates vcf files. The tool accepts vcf files as input. The output is a tsv file containing the following information for each variant: reference allele, alternate allele(s), variant type, effect, total read depth of all samples, number of reads supporting variant, percentage of reads supporting the variant, allele frequency, and then read depth, number of reads supporting variant, and percentage of reads supporting the variant for each sample. 

# Requirements
The tool is a python script, so Python 3.7 or greater must be installed.

The script utilizes the following libraries: pysam, json, and requests. Json and requests are native to python, so there is no need to install them. Pysam can be installed with:
```{r eval=FALSE,echo=TRUE}
pip install pysam
```
or if you have Anaconda (https://anaconda.org/bioconda/pysam) installed:
```{r eval=False, echo=True}
conda install -c bioconda pysam
```
The pysam library is used to parse vcf files. More info here: https://pysam.readthedocs.io/en/latest/.

# Usage
Either of these commands will work on the command line:
```{r eval=FALSE,echo=TRUE}
./tempus_bioi_challenge.py -i /path/to/in.vcf -o /path/to/out.tsv
```
```{r eval=FALSE,echo=TRUE}
python3 tempus_bioi_challenge.py -i /path/to/in.vcf -o /path/to/out.tsv
```


