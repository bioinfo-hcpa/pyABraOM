# pyABraOM

![](https://img.shields.io/badge/python-3.x-purple)

## Summary

- [Introduction](#introduction)
- [Installation](#installation)
- [Search Types](#search-types)
    - [Search by gene](#search-by-gene)
    - [Search by region](#search-by-genomic-region)
    - [Search by variant](#search-by-variant)
- [Batch search](#batch-search)

## Introduction
   PyABraOM is an API developed to facilitate the access to human variant data from ABraOM database to access data from GRCh37/hg19 and GRCh38/hg38 genome version. The package offers support to three kinds of different searches and the possibility to search in batches.For plotting the data, please take a look at the [BIOVARS package](https://github.com/bioinfo-hcpa/biovars).
   
## Installation

   Currently there is not a PyPI version for the PyABraOM API, so the installation needs that you clone this repository and install it as local package.
```
$ git clone https://github.com/bioinfo-hcpa/pynoma.git
$ pip install -e pynoma #outside the directory
```

## Search Types

   There are 3 kinds of search by using pyABraOM and all returned in the same format: Gene, Genomic Region and Variant search.
   
   Since the ABraOM database has the two latest human genome versions, the user must specify the version in GRCh37/hg19 or GRCh38/hg38 to retrieve data from the database for the specific genome version. Besides, when searching by  genes, genomic regions or variant for chromosomes X or Y, a new column "Number of Hemizygous" will be added to the output dataframe, so the user should have caution when performing pandas concatenation operations, resulting in table cells with NaN values. Futhermore, the GATK and CEGH filter provided from ABraOM are available, but latter is not in the default output dataframe.
   
   
### Search by gene

The search method expects the gene name to be considered in the function Search_gene(version:str,gene:str,CEGH_Filter= False,Variant_ID= False, Process= True).

* version (str):the reference genome version (either "GRCh37/hg19" or "GRCh38/hg38")
* gene (str): the gene name to search in database
* CEGH_Filter (bool): whether to provide in-house quality control filter CEGH-USP
* Variant_ID (bool): whether to provide variant identifier
* Process (bool): whether to provide data in JSON format or pandas dataframe

```python
import pyabraom
from pyabraom import Seach_gene
df = Search_gene('hg38', 'ACE2')
```
### Search by genomic region

The genomic region method expects the genome region to be considered in the function Search_region(version:str,chromosome,start,end,CEGH_Filter= False,Variant_ID= False, Process= True). The X and Y chromosome must be provided as character in the chromosome parameter.

* version (str):the reference genome version (either "GRCh37/hg19" or "GRCh38/hg38")
* chromosome(obj): chromosome 
* start (int): where the region of interest starts 
* end (int):where the region of interest ends
* CEGH_Filter (bool): whether to provide in-house quality control filter CEGH-USP
* Variant_ID (bool): whether to provide variant identifier
* Process (bool): whether to provide data in JSON format or pandas dataframe

```python
import pyabraom
from pyabraom import Search_region
df = Search_region('hg38',4,980883,984868)
```

### Search by variant

The variant search method expects the variant identifier rsID i the function Variant_ID(version:str,variant:str,CEGH_Filter= False,Process= True):

* version (str):the reference genome version (either "GRCh37/hg19" or "GRCh38/hg38")
* variant (str): the variant identifier rsID to search in database
* CEGH_Filter (bool): whether to provide in-house quality control filter CEGH-USP
* Variant_ID (bool): whether to provide variant identifier
* Process (bool): whether to provide data in JSON format or pandas dataframe

```python
import pyabraom
from pyabraom import Variant_ID
df = Variant_ID('hg38','rs115134980')
```

### Batch Search 

The batch search method expects as parameter a list of one provided type of search in the function Searches(searches_list:list)

```python
import pyabraom
df = Searches([Search_gene('hg38','ACE2'),Search_gene('hg38','TMPRSS2')])
```
## Features

Not available yet.
* Variant Filter

## Additional Information

## References

Naslavsky, Scliar, Yamamoto, et al., 2022. Whole-genome sequencing of 1,171 elderly admixed individuals from São Paulo, Brazil.
