# pyABraOM

![](https://img.shields.io/badge/python-v1.x-purple)

## Introduction


## Search Types

   There are 3 kinds of search by using pyABraOM and all returned in the same format: Gene, Genomic Region and Variant search.
   
   Since the ABraOM database has the two latest human genome versions, the user must specify the version in GRCh37/hg19 or GRCh38/hg38 to retrieve data from the database for the specific genome version. Besides, when searching by  genes, genomic regions or variant for chromosomes X or Y, a new column "Number of Hemizygous" will be added to the output dataframe, so the user should have caution when performing pandas concatenation operations, resulting in table cells with NaN values. Futhermore, the GATK and CEGH filter provided from ABraOM are available, but latter is not in the default output dataframe.
   
   
### Search by gene

Search_gene(version:str,gene:str,CEGH_Filter= False,Variant_ID= False, Process= True):

```python
import pyabraom
from pyabraom import Seach_gene
df = Search_gene('hg38', 'ACE2')
```
### Search by  genomic region

Search_region(version:str,chromosome,start,end,CEGH_Filter= False,Variant_ID= False, Process= True)

```python
import pyabraom
from pyabraom import Search_region
df = Search_region('hg38',4,980883,984868)
```

### Search by  variant

Variant_ID(version:str,variant:str,CEGH_Filter= False,Process= True):

```python
import pyabraom
from pyabraom import Variant_ID
df = Variant_ID('hg38','rs115134980')
```

### Batch Search 

Searches(searches_list:list)

```python
import pyabraom
df = Searches([Search_gene('hg38','ACE2'),Search_gene('hg38','TMPRSS2')])
```
## Features

Not available yet.
* Variant Filter

## Additional Information

## References

```
