import requests
import pandas as pd
from pandas import json_normalize
import time
from Bio import Entrez, SeqIO
import pkg_resources

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
        certificate_path = pkg_resources.resource_filename('pyabraom', 'CertBundle.pem')

        if GATK_PASS==False:
           response= requests.post(url, data={"table":version,"str":query},verify=certificate_path,timeout =None) 
 
        else:
           response= requests.post(url, data={"table":version,"str":query,"gatk":'PASS'},verify=certificate_path,timeout =None) 
 
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
       data['rsID'] =  [i if i != 'NA' else '-' for i in Variant_ID_list]

    data["GATK Filter"]=Filter

    if CEGH_Filter==True:
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

def Searches(lista:list):
     appended_data = []

     for i in lista:
         data = i
         appended_data.append(data)

     final_data = pd.concat(appended_data, ignore_index=True)
     if 'Number of Hemizygotes' in final_data.columns:
         final_data['Number of Hemizygotes'] = final_data['Number of Hemizygotes'].fillna(0)
     if 'Variant_ID' in final_data.columns:
         final_data['Variant_ID'] = final_data['Variant_ID'].fillna('-')

     return final_data


def Search_gene(version:str,gene:str,CEGH_Filter= False,GATK_PASS= False,Variant_ID=False, Process= True):
         
        response = Request(version,gene,GATK_PASS)

        if Process==True:
            result = Process_data(response,CEGH_Filter,Variant_ID)
            data= Dataframe_adjust(result)

        else:
           data= response

        return data
       

def Search_region(version:str,chromosome,start,end,CEGH_Filter= False,Variant_ID=False, Process= True):

        if start > end:
            raise Exception('End needs to be greater than start position.')

        else:
            region = "%d"":""%d""-""%d" %(chromosome,start,end)
            response = Request(version,region)

            if Process==True:
              result = Process_data(response,CEGH_Filter,Variant_ID)
              dataframe= Dataframe_adjust(result)

            else:
              data= response

        return dataframe


def Variant_ID(version:str,variant:str,CEGH_Filter= False,GATK_PASS= False,Process= True):

    response = Request(version,variant,GATK_PASS)
    Variant_ID=False
    if Process==True:
       result = Process_data(response,CEGH_Filter,Variant_ID)
       data= Dataframe_adjust(result)

    else:
       data= response

    return data

def genome_ref_info(chromosome,version,position):
      chromosomes_hg38 = {1:568815597,
                    2:568815596,
                    3:568815595,
                    4:568815594,
                    5:568815593,
                    6:568815592,
                    7:568815591,
                    8:568815590,
                    9:568815589,
                    10:568815588,
                    11:568815587,
                    12:568815586,
                    13:568815585,
                    14:568815584,
                    15:568815583,
                    16:568815582,
                    17:568815581,
                    18:568815580,
                    19:568815579,
                    20:568815578,
                    21:568815577,
                    22:568815576,
                    'X':568815575,
                    'Y':568815574}


      chromosomes_hg37 = {1:224589800,
                    2:224589811,
                    3:224589815,
                    4:224589816,
                    5:224589817,
                    6:224589818,
                    7:224589819,
                    8:224589820,
                    9:224589821,
                    10:224589801,
                    11:224589802,
                    12:224589803,
                    13:224589804,
                    14:224589805,
                    15:224589806,
                    16:224589807,
                    17:224589808,
                    18:224589809,
                    19:224589810,
                    20:224589812,
                    21:224589813,
                    22:224589814,
                    'X':224589822,
                    'Y':224589823}

      if version in ["hg38","GRCh38"]:
        handle = Entrez.efetch(db="nucleotide", 
                              id= str(chromosomes_hg38[chromosome]) ,
                              rettype="fasta", 
                              strand=1, 
                              seq_start=position, 
                              seq_stop=position)
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)
      if version in ["hg19", "GRCh37"]:
        handle = Entrez.efetch(db="nucleotide", 
                              id= str(chromosomes_hg37[chromosome]) ,
                              rettype="fasta", 
                              strand=1, 
                              seq_start=position, 
                              seq_stop=position)
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)

def Variant_ID_biovar(dataframe,version):
    list1= []
    for i in range(0,len(dataframe)):
      if dataframe['Alternative'][i]=='-': 
        next_nt=genome_ref_info(int(dataframe['Chromosome'][i]),version, int(dataframe['Position'][i])+1)
        variant_ID_biovar = dataframe['Chromosome'][i]+'-'+dataframe['Position'][i]+'-'+dataframe['Reference'][i]+'-'+next_nt
        frequency = dataframe['Allele Frequency'][i]
        list1.append([variant_ID_biovar,frequency])

      elif dataframe['Reference'][i]=='-': 
        prev_nt=genome_ref_info(int(dataframe['Chromosome'][i]),version, int(dataframe['Position'][i]))
        variant_ID_biovar = dataframe['Chromosome'][i]+'-'+dataframe['Position'][i]+'-'+prev_nt+'-'+prev_nt+dataframe['Alternative'][i]
        frequency = dataframe['Allele Frequency'][i]
        list1.append([variant_ID_biovar,frequency])
      elif dataframe['Reference'][i]!='-' and dataframe['Alternative'][i]!='-':
        variant_ID_biovar = dataframe['Chromosome'][i]+'-'+dataframe['Position'][i]+'-'+dataframe['Reference'][i]+'-'+dataframe['Alternative'][i]
        frequency = dataframe['Allele Frequency'][i]
        list1.append([variant_ID_biovar,frequency])


    return list1
