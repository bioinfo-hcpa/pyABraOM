import requests
import pandas as pd
import sys
import urllib3
from tqdm import tqdm

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def Genome_version(version:str):
 
          if version=="hg38" or "GRCh38":
              version = "abraomdbgenoma"
 
          elif version=="hg19" or "GRCh37":
              version= "abraomdb"
          else:
             raise Exception("Error:This version is not present in AbraOM")
 
          return version                                 
 
 
def Request(version:str,query:str,GATK_PASS=False):
   url ="https://abraom.ib.usp.br/script.php" 

   version = Genome_version(version)

   if not GATK_PASS:
      response= requests.post(url, data={"table": version, "str": query},
                              verify=False,
                              timeout=None) 

   else:
      response= requests.post(url, data={"table": version, 
                                         "str": query,
                                         "gatk":'PASS'},
                              verify=False,
                              timeout=None) 
   return response.json()
 

def Process_data(df,CEGH_Filter,Variant_ID):

    Chr = []
    start =[]
    ref = []
    alt =[]
    gene_list = []
    Filter =[]
    predfunc=[]
    allele_number=[]
    allele_alt_count=[]
    homozygous_alt_allele=[]
    frequency =[]
    CEGH=[]
    Variant_ID_list=[]
    hemizygous =[]
    data= pd.DataFrame()
    for i in range(0,len(df["data"])):
        Chr.append(df["data"][i]["Chr"])

        start.append(df["data"][i]["Start"])

        ref.append(df["data"][i]["Ref"])

        alt.append(df["data"][i]["Alt"])

        gene_name = (df["data"][i]["Gene_refGene"]).split(";")
        gene_list.append(gene_name[0])

        Filter.append(df["data"][i]["FILTER"])

        predfunc.append(df["data"][i]["PredictedFunc_refGene"].replace(';','_'))
        
        allele_number.append(df["data"][i]["Allele_number"])

        allele_alt_count.append(df["data"][i]["Allele_ALT_count"])

        homozygous_alt_allele.append(df["data"][i]["HomozygousALT_count"])

        frequency.append(df["data"][i]["Frequencies"])

        if CEGH_Filter:
           CEGH.append(df["data"][i]["CEGH_Filter"])

        if df["data"][i]["Chr"] in ['X','Y']:
           hemizygous.append(df["data"][i]["Hemizygous_count"])

        ID= df["data"][i]["avsnp147"].split()
        Variant_ID_list.append(ID[4].replace('</span><div',''))

    data["Chromosome"] = Chr

    data["Position"]=start

    data["Reference"]= ref

    data["Alternative"]=alt

    data["Annotation"]=predfunc

    data["Gene"]=gene_list

    if Variant_ID:
       data['rsID'] =  [i if i != 'NA' else '-' for i in Variant_ID_list]

    data["GATK Filter"]=Filter

    if CEGH_Filter:
       data["CEGH Filter"]=CEGH

    data["Number of Homozygotes"]= homozygous_alt_allele

    if len(hemizygous)!= 0:
        data["Number of Hemizygotes"]= hemizygous

    data["Allele Number"]= allele_number

    data["Allele Count"]= allele_alt_count

    data["Allele Frequency"]= frequency

    return data



def Dataframe_adjust(dataframe):
    pd.set_option('display.max_column',None)
    pd.set_option('display.max_rows',None)
    pd.set_option('display.max_seq_items',None)
    pd.set_option('display.max_colwidth', 20)
    pd.set_option('expand_frame_repr', True)
    pd.set_option('display.colheader_justify','left')

    return dataframe

def Searches(lista:list,verbose=True):
     appended_data = []
     searches_data = tqdm(lista, desc="Batch searching") if verbose else lista
     for data in searches_data:
         appended_data.append(data)

     final_data = pd.concat(appended_data, ignore_index=True)
     if 'Number of Hemizygotes' in final_data.columns:
         final_data['Number of Hemizygotes'] = final_data['Number of Hemizygotes'].fillna(0)
     if 'Variant_ID' in final_data.columns:
         final_data['Variant_ID'] = final_data['Variant_ID'].fillna('-')

     return final_data


def Search_gene(version:str,gene:str,CEGH_Filter= False,GATK_PASS= False,Variant_ID=False, Process= True):
        
        gene= gene.upper()
        response = Request(version,gene,GATK_PASS)
        
        if Process:
            result = Process_data(response,CEGH_Filter,Variant_ID)
            data= Dataframe_adjust(result)
            data= data.loc[data['Gene'] == gene]
            data.reset_index(inplace=True,drop=True)
        else:
           data= response

        return data
       

def Search_region(version:str,chromosome,start,end,CEGH_Filter= False,Variant_ID=False, Process= True):

      assert start < end, "End needs to be greater than start position."
      if isinstance(chromosome, int):
         region = "%d"":""%d""-""%d" %(chromosome, start, end)
      else:
         region = "%s"":""%d""-""%d" %(chromosome, start, end)

      response = Request(version,region)
      if Process:
         result = Process_data(response, CEGH_Filter, Variant_ID)
         data = Dataframe_adjust(result)
      else:
         data = response

      return data



def Variant_ID(version:str,variant:str,CEGH_Filter= False,GATK_PASS= False,Process= True):

    response = Request(version,variant,GATK_PASS)
    Variant_ID=False
    if Process:
       result = Process_data(response,CEGH_Filter,Variant_ID)
       data= Dataframe_adjust(result)

    else:
       data= response

    return data


def genome_ref_info(chromosome,version,position):
       
    if version in ["hg38","GRCh38"]:
       server = "https://rest.ensembl.org/sequence/region/human/%s:%d..%d:1?" %(chromosome,position,position)
    if version in ["hg19", "GRCh37"]:
       server=  "https://grch37.rest.ensembl.org/sequence/region/human/%s:%d..%d:1?" %(chromosome,position,position)
       
    request = requests.get(server, headers={ "Content-Type" : "text/plain"})
    if not request.ok:
       request.raise_for_status()
       sys.exit()
    return(request.text)




def Variant_ID_biovars(df,version):
    list1= []

    for i in range(0,len(df)): 
       if df['Reference'][i] != '-' :
          variant_ID_biovar= df['Chromosome'][i]+'-'+str(df['Position'][i])+'-'+df['Reference'][i]+'-'+df['Alternative'][i]
          list1.append([variant_ID_biovar,df['Allele Frequency'][i]])
          
       elif df['Reference'][i] == '-':
          prev_nt=genome_ref_info(df['Chromosome'][i], version, int(df['Position'][i]))
          variant_ID_biovar = df['Chromosome'][i]+'-'+str(df['Position'][i])+'-'+prev_nt+'-'+prev_nt+df['Alternative'][i]
          list1.append([variant_ID_biovar,df['Allele Frequency'][i]])


    return list1
