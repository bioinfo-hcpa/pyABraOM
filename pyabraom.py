import requests
import pandas as pd
import time
from pandas import json_normalize
import pkg_resources

def genome_version(version:str):
 
          if version=="hg38" or "GRCh38":
              version = "abraomdbgenoma"
 
          elif version=="hg19" or "GRCh37":
              version= "abraomdb"
          else:
             raise Exception("Error:This version is not present in AbraOM")
 
          return version                                 
 
 
def Request(version:str,query:str,GATK_PASS=False):
 
        url ="https://abraom.ib.usp.br/script.php" 
 
        version = genome_version(version)
        certificate_path = pkg_resources.resource_filename('pyabraom', 'CertBundle.pem')

        if GATK_PASS==False:
           response= requests.post(url, data={"table":version,"str":query},verify=certificate_path,timeout =None) 
 
        else:
           response= requests.post(url, data={"table":version,"str":query,"gatk":'PASS'},verify=certificate_path,timeout =None) 
 
        return response.json()
 
def process_data(df,CEGH_Filter):
    chr = []
    start =[]
    ref = []
    alt =[]
    gene = []
    filter =[]
    prediction_func=[]
    allele_number=[]
    allele_alt_count=[]
    homozygous_alt_allele=[]
    frequencie =[]
    CEGH=[]
    Variant_ID=[]
    data= pd.DataFrame()
    for i in range(0,len(df["data"])):
        chr.append(df["data"][i]["Chr"])
 
        start.append(df["data"][i]["Start"])
        if len(df["data"][i]["Ref"])<=10:
           ref.append(df["data"][i]["Ref"])
        else:
           ref_object =df["data"][i]["Ref"]
           ref.append(ref_object[:3])
 
        if len(df["data"][i]["Alt"])<=3:
           alt.append(df["data"][i]["Alt"])
        else:
           alt_object =df["data"][i]["Alt"]
           alt.append(alt_object[:3])
 
        if len(df["data"][i]["Gene_refGene"]) ==1:
           gene.append(df["data"][i]["Gene_refGene"])
 
        else:
           gene = (df["data"][i]["Gene_refGene"]).split()
           gene = gene[0]
 
        filter.append(df["data"][i]["FILTER"])
 
        prediction_func.append(df["data"][i]["PredictedFunc_refGene"])
        
        allele_number.append(df["data"][i]["Allele_number"])
 
        allele_alt_count.append(df["data"][i]["Allele_ALT_count"])
 
        homozygous_alt_allele.append(df["data"][i]["HomozygousALT_count"])
 
        frequencie.append(df["data"][i]["Frequencies"])
        if CEGH_Filter==True:
           CEGH.append(df["data"][i]["CEGH_Filter"])
 
        ID= df["data"][i]["avsnp147"].split()
        Variant_ID.append(ID[4].replace('</span><div',''))
 
    data["Chr"] = chr
    data["Start"]=start
    data["Ref"]= ref
    data["Alt"]=alt
    if Variant_ID==True:
       data["ID"]=Variant_ID  
    data["Gene"]=gene
    data["GATK Filter"]=filter
    if CEGH_Filter==True:
       data["CEGH"]=CEGH
    data["Annotation"]=prediction_func
    data["Homozygotes"]= homozygous_alt_allele
    data["Allele Number"]= allele_number
    data["Allele Alt"]= allele_alt_count
    data["Frequencie"]= frequencie
 
    return data
'''
def style(dataframe):
    dataframe = dataframe.style.set_table_styles(
    [{'selector':'th', 
      'props':[('font_family',' Sans-serif'),('text-align', 'center'),('background',' #FFFFFF'),("font-size", "10pt"),('font-weight','bold')]},
     
   {'selector':'td','props':[('font_family',' Times New Rowman'),('text-align', 'center'),("font-size", "10pt"),('width', '350px')]},
 
   {'selector':'tr:nth-of-type(odd)','props':[('background',' #f3f2f2')]},
 
   {'selector':'td:nth-child(1)','props':[('font-weight','bold')]},
 
   {'selector':'tr:hover','props':[('background','#bfbfbf')]}]).hide_index()
 
    return dataframe
'''
def dataframe_adjust(dataframe):
    pd.set_option('display.max_column',None)
    pd.set_option('display.max_rows',None)
    pd.set_option('display.max_seq_items',None)
    pd.set_option('display.max_colwidth', 600)
    pd.set_option('expand_frame_repr', True)
    pd.set_option('display.colheader_justify','right')
 
 
    return dataframe
 
 
def Search_gene(version:str,gene:str,CEGH_Filter=False,GATK_PASS=False,Process_data=True):
         
        response = Request(version,gene,GATK_PASS)
        if Process_data==True:
            result = process_data(response,CEGH_Filter)
            data= dataframe_adjust(result)
        else:
           data= response
 
        return data
       
 
def Search_region(version:str,chromosome,start,end,CEGH_Filter=False,Process_data=True):
        region = "%d"":""%d""-""%d" %(chromosome,start,end)
        response = Request(version,region)
        if Process_data==True:
           result = process_data(response,CEGH_Filter)
           dataframe= dataframe_adjust(result)
        else:
           data= response
 
        return dataframe
 
 
 
 
def Variant_ID(version:str,variant:str,CEGH_Filter=False,GATK_PASS=False,Process_data=True):
    response = Request(version,variant,GATK_PASS)
    if Process_data==True:
       result = process_data(response,CEGH_Filter)
       data= dataframe_adjust(result)
    else:
       data= response
    
    return data
