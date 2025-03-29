import csv
import os


def get_topp_cat()->list:
    """
    Get a list of ToppFun categories

    Returns: a list

    Examples:

        get_topp_cat()


    """
    toppCats = ["GeneOntologyMolecularFunction","GeneOntologyBiologicalProcess",
    "GeneOntologyCellularComponent","HumanPheno","MousePheno","Domain","Pathway",
    "Pubmed","Interaction","Cytoband","TFBS","GeneFamily","Coexpression","CoexpressionAtlas",
    "ToppCell","Computational","MicroRNA","Drug","Disease"]
    return toppCats


def toppsave(topp_data,filename:str,save_dir:str,split:bool,format_1:str)->None:
    """
    Save toppData results (optionally) split by celltype/cluster

    Args:
        topp_data: Results from toppFun as a dataframe
        filename: filename prefix for each split file
        save_dir: the directory to save files
        split: Boolean, whether to split the dataframe by celltype/cluster
        format_1: Saved file format, one of ["xlsx", "csv", "tsv"]

    Returns: None

    Examples:

        toppsave(topp_data, filename="toppFun_results", split = TRUE, format_1 = "xlsx")

   """



    if save_dir is None:
        save_dir=os.getcwd()
        os.makedirs(save_dir, exist_ok=True)

    if not split:
        #不分组
        if format_1=='xlsx':
            if filename is None:
                filename='toppData.xlxs'
            else:
                filename=f'{filename}.xlsx'
            path=os.path.join(save_dir,filename)
            topp_data.to_excel(path,header=True,index=False,sheet_name='toppData')
            print('Saving file:',filename,'\n')

        elif format_1=='csv':
            if filename is None:
                filename='toppData.csv'
            else:
                filename=f'{filename}.csv'
            path = os.path.join(save_dir, filename)
            topp_data.to_csv(path,sep=',',quoting=csv.QUOTE_MINIMAL,header=True,index=False)
            print('Saving file:',filename,'\n')

        elif format_1=='tsv':
            if filename is None:
                filename='toppData.tsv'
            else:
                filename=f'{filename}.tsv'
            path = os.path.join(save_dir, filename)
            topp_data.to_csv(path,sep='\t',quoting=csv.QUOTE_MINIMAL,header=True,index=False)
            print('Saving file:',filename,'\n')
    else:
        #分组
        # 遍历cluster
        for gr in topp_data['Cluster'].unique():
            tmp_toppdata=topp_data.loc[topp_data.Cluster==gr]
            if format_1 == 'xlsx':
                if filename is None:
                    current_filename = f'toppData_{gr}.xlsx'
                else:
                    current_filename = f'{filename}_{gr}.xlsx'
                path = os.path.join(save_dir, current_filename)
                tmp_toppdata.to_excel(path,header=True,index=False,sheet_name='toppData')
                print('Saving file:', current_filename, '\n')

            elif format_1 == 'csv':
                if filename is None:
                    current_filename = f'toppData_{gr}.csv'
                else:
                    current_filename = f'{filename}_{gr}.csv'
                # 将数据写入该格式中
                path = os.path.join(save_dir, current_filename)
                tmp_toppdata.to_csv(path, sep=',', quoting=csv.QUOTE_MINIMAL, header=True, index=False)
                print('Saving file:', current_filename, '\n')

            elif format_1 == 'tsv':
                if filename is None:
                    current_filename = f'toppData_{gr}.tsv'
                else:
                    current_filename = f'{filename}_{gr}.tsv'
                path = os.path.join(save_dir, current_filename)
                tmp_toppdata.to_csv(path, sep='\t', quoting=csv.QUOTE_MINIMAL, header=True, index=False)
                print('Saving file:', current_filename, '\n')