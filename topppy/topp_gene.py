import csv
import json
import os
from pandas import DataFrame
import requests

def get_topp_cats()->list:
    """
    Get a list of ToppFun categories

    Returns: a list

    Examples:

        from topppy import get_topp_cats
        get_topp_cats()


    """
    toppCats = ["GeneOntologyMolecularFunction","GeneOntologyBiologicalProcess",
    "GeneOntologyCellularComponent","HumanPheno","MousePheno","Domain","Pathway",
    "Pubmed","Interaction","Cytoband","TFBS","GeneFamily","Coexpression","CoexpressionAtlas",
    "ToppCell","Computational","MicroRNA","Drug","Disease"]
    return toppCats


def topp_save(topp_data: DataFrame, filename:str = None, save_dir:str = None, split:bool = False, format:str = "xlsx")->None:
    """
    Save topp_data results (optionally) split by celltype/cluster

    Args:
        topp_data: Results from toppFun as a dataframe
        filename: File name prefix for each split file. Default: Current working directory
        save_dir: The directory to save files. Default: $HOME
        split: Boolean, whether to split the dataframe by celltype/cluster. Default: False
        format: Saved file format, one of ["xlsx", "csv", "tsv"]. Default: "xlsx"

    Returns: None

    Examples:

        from topppy import topp_save, topp_data
        topp_save(topp_data, filename="toppFun_results", split = True, format = "xlsx")

   """



    if save_dir is None:
        save_dir=os.getcwd()
    os.makedirs(save_dir, exist_ok=True)

    if not split:
        # Not grouped
        if format== 'xlsx':
            if filename is None:
                filename= 'toppData.xlxs'
            else:
                filename= f'{filename}.xlsx'
            path=os.path.join(save_dir, filename)
            topp_data.to_excel(path,header=True,index=False,sheet_name='toppData')
            print('Saving file:', filename, '\n')

        elif format== 'csv':
            if filename is None:
                filename= 'toppData.csv'
            else:
                filename= f'{filename}.csv'
            path = os.path.join(save_dir, filename)
            topp_data.to_csv(path,sep=',',quoting=csv.QUOTE_MINIMAL,header=True,index=False)
            print('Saving file:', filename, '\n')

        elif format== 'tsv':
            if filename is None:
                filename= 'toppData.tsv'
            else:
                filename= f'{filename}.tsv'
            path = os.path.join(save_dir, filename)
            topp_data.to_csv(path,sep='\t',quoting=csv.QUOTE_MINIMAL,header=True,index=False)
            print('Saving file:', filename, '\n')
    else:
        # Grouped
        # Traverse the cluster
        for gr in topp_data['Cluster'].unique():
            tmp_toppdata=topp_data.loc[topp_data.Cluster==gr]
            if format == 'xlsx':
                if filename is None:
                    current_filename = f'toppData_{gr}.xlsx'
                else:
                    current_filename = f'{filename}_{gr}.xlsx'
                path = os.path.join(save_dir, current_filename)
                tmp_toppdata.to_excel(path,header=True,index=False,sheet_name='toppData')
                print('Saving file:', current_filename, '\n')

            elif format == 'csv':
                if filename is None:
                    current_filename = f'toppData_{gr}.csv'
                else:
                    current_filename = f'{filename}_{gr}.csv'
                # Write the data into that format
                path = os.path.join(save_dir, current_filename)
                tmp_toppdata.to_csv(path, sep=',', quoting=csv.QUOTE_MINIMAL, header=True, index=False)
                print('Saving file:', current_filename, '\n')

            elif format == 'tsv':
                if filename is None:
                    current_filename = f'toppData_{gr}.tsv'
                else:
                    current_filename = f'{filename}_{gr}.tsv'
                path = os.path.join(save_dir, current_filename)
                tmp_toppdata.to_csv(path, sep='\t', quoting=csv.QUOTE_MINIMAL, header=True, index=False)
                print('Saving file:', current_filename, '\n')

def get_entrez(genes:list)->list:
    """
    Convert genes into Entrez format

    Args:
        genes:A list of genes

    Returns: a vector of genes in Entrez format

    Examples:

        from topppy import get_entrez
        get_entrez(genes)

    """
    url = "https://toppgene.cchmc.org/API/lookup"
    payload={'Symbols':genes}
    headers = {'Content-Type': 'application/json'}
    r=json.loads(requests.post(url,json=payload,headers=headers).text)
    new_gene_list=[]
    for g in r['Genes']:
        new_gene_list.append(g['Entrez'])
    return new_gene_list