import requests
import pandas as pd
import time


def Genome_version(version:str):

          if version in ["hg38","GRCh38"]:
              version = "abraomdbgenoma"

          elif version in ["hg19", "GRCh37"]:
              version= "abraomdb"
          else:
             raise Exception("Error: This genome version is not part of ABraOM project.Try search for hg19 or hg38 version.")

          return  version                                

 
def Request(version:str,query:str,GATK_PASS=False):

        url = "http://abraom.ib.usp.br/script.php"
 
        version = genome_version(version)
        if GATK_PASS==False:
           response= requests.post(url, data={"table":version,"str":query},timeout =None) 

        else:
           response= requests.post(url, data={"table":version,"str":query,"gatk":'PASS'},timeout =None) 

        return response.json()

def Process_data(df,CEGH_Filter,Variant_ID,gene):

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
      if gene.upper() in df["data"][i]["Gene_refGene"]:
        Chr.append(df["data"][i]["Chr"])

        start.append(df["data"][i]["Start"])

        ref.append(df["data"][i]["Ref"])

        alt.append(df["data"][i]["Alt"])

        gene_name = (df["data"][i]["Gene_refGene"]).split(";")
        gene_list.append(gene_name[0])

        Filter.append(df["data"][i]["FILTER"])

        predfunc.append(df["data"][i]["PredictedFunc_refGene"])
        
        allele_number.append(df["data"][i]["Allele_number"])

        allele_alt_count.append(df["data"][i]["Allele_ALT_count"])

        homozygous_alt_allele.append(df["data"][i]["HomozygousALT_count"])

        frequency.append(df["data"][i]["Frequencies"])

        if CEGH_Filter==True:
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

    if Variant_ID==True:
       data['ID'] =  [i if i != 'NA' else '-' for i in Variant_ID_list]

    data["GATK Filter"]=Filter

    if CEGH_Filter==True:
       data["CEGH"]=CEGH

    data["Number of Homozygotes"]= homozygous_alt_allele

    if len(hemizygous)!= 0:
        data["Number of Hemizygous"]= hemizygous

    data["Allele Number"]= allele_number

    data["Allele Count"]= allele_alt_count

    data["Allele Frequency"]= frequency

    return data

def Variant_ID_biovar(dataframe):
    list1= []
    for i in range(0,len(dataframe)):
      if len(dataframe['Alternative'][i])==1: 
        variant_ID_biovar = dataframe['Chromosome'][i]+'-'+dataframe['Position'][i]+'-'+dataframe['Reference'][i]+'-'+dataframe['Alternative'][i]
        frequency = dataframe['Allele Frequency'][i]
        list1.append([variant_ID_biovar,frequency])
    return list1


def Dataframe_adjust(dataframe):
    pd.set_option('display.max_column',None)
    pd.set_option('display.max_rows',None)
    pd.set_option('display.max_seq_items',None)
    pd.set_option('display.max_colwidth', 20)
    pd.set_option('expand_frame_repr', True)
    pd.set_option('display.colheader_justify','left')

    return dataframe

def Searches(lista:list):
     appended_data = []

     for i in lista:
         data = i
         appended_data.append(data)

     final_data = pd.concat(appended_data, ignore_index=True)

     final_data['Hemizygous'] = final_data['Hemizygous'].fillna(0)

     final_data['ID'] = final_data['ID'].fillna('-')

     return final_data


def Search_gene(version:str,gene:str,CEGH_Filter= False,GATK_PASS= False,Variant_ID= False, Process= True):
         
        response = Request(version,gene,GATK_PASS)

        if Process==True:
            result = Process_data(response,CEGH_Filter,Variant_ID,gene)
            data= Dataframe_adjust(result)

        else:
           data= response

        return data
       

def Search_region(version:str,chromosome,start,end,CEGH_Filter= False,Variant_ID= False, Process= True):

        if start > end:
            raise Exception('End needs to be greater than start position.')

        else:
            region = "%s"":""%d""-""%d" %(chromosome,start,end)
            response = Request(version,region)

            if Process==True:
              result = Process_data(response,CEGH_Filter,Variant_ID)
              dataframe= Dataframe_adjust(result)

            else:
              data= response

        return dataframe


def Variant_ID(version:str,variant:str,CEGH_Filter= False,GATK_PASS= False,Process= True):

    response = Request(version,variant,GATK_PASS)

    if Process==True:
       result = Process_data(response,CEGH_Filter)
       data= Dataframe_adjust(result)

    else:
       data= response

    return data
    
