import csv
import os
from data import topp_data

def get_topp_cat()->list:
    """
    Returns:Returns a list of all categories in the ToppGene database
    """
    topp_data_1=topp_data['category'].drop_duplicates()
    topp_cats=list(topp_data_1)
    return topp_cats

def toppsave(topp_data,filename:str,save_dir:str,split:bool,format_1:str)->None:
    '''

    Args:
        topp_data:Results from toppFun as a dataframe
        filename:filename prefix for each split file
        save_dir:the directory to save files
        split:Boolean, whether to split the dataframe by celltype/cluster
        format_1: Saved file format, one of c("xlsx", "csv", "tsv")

    Returns:Save toppData results (optionally) split by celltype/cluster

    '''

    if save_dir is None:
        save_dir=os.getcwd()

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
            topp_data.to_csv(path,sep=',',quoting=csv.QUOTE_NONE,header=True,index=False)
            print('Saving file:',filename,'\n')

        elif format_1=='tsv':
            if filename is None:
                filename='toppData.tsv'
            else:
                filename=f'{filename}.tsv'
            path = os.path.join(save_dir, filename)
            topp_data.to_csv(path,sep='\t',quoting=csv.QUOTE_NONE,header=True,index=False)
            print('Saving file:',filename,'\n')
    else:
        #分组
        # 遍历cluster
        for gr in topp_data['Cluster'].unique():
            tmp_toppdata=topp_data.loc[topp_data.Cluster==gr]
            if format_1 == 'xlsx':
                if filename is None:
                    filename = f'toppData_{gr}.xlsx'
                else:
                    filename = f'{filename}_{gr}.xlsx'
                path = os.path.join(save_dir, filename)
                tmp_toppdata.to_excel(path,header=True,index=False,sheet_name='toppData')
                print('Saving file:', filename, '\n')

            elif format_1 == 'csv':
                if filename is None:
                    filename = f'toppData_{gr}.csv'
                else:
                    filename = f'{filename}_{gr}.csv'
                # 将数据写入该格式中
                path = os.path.join(save_dir, filename)
                tmp_toppdata.to_csv(path, sep=',', quoting=csv.QUOTE_NONE, header=True, index=False)
                print('Saving file:', filename, '\n')

            elif format_1 == 'tsv':
                if filename is None:
                    filename = f'toppData_{gr}.tsv'
                else:
                    filename = f'{filename}_{gr}.tsv'
                path = os.path.join(save_dir, filename)
                tmp_toppdata.to_csv(path, sep='\t', quoting=csv.QUOTE_NONE, header=True, index=False)
                print('Saving file:', filename, '\n')