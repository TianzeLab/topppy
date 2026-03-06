import csv
import json
import os
import sys
import pandas as pd
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
        topp_data: Results from topp_fun as a dataframe
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


def process_markers(markers:DataFrame,cluster_col:str,gene_col:str,p_val_col:str,logFC_col:str,
                   num_genes:int = 1000,pval_cutoff:float = 0.5,fc_cutoff:float = 0,fc_filter:str ='ALL',genes_submit_cutoff=1000)->dict:

    markers_list=dict()
    for cl in markers[cluster_col].unique():
        all_cl_markers=markers[(markers[cluster_col]==cl)&(markers[p_val_col]<pval_cutoff)]
        if fc_filter=='ALL':
            all_cl_markers=all_cl_markers[abs(all_cl_markers[logFC_col])>fc_cutoff]
            all_cl_markers=all_cl_markers.sort_values(by=logFC_col,key=lambda x:abs(x),ascending=False)[[gene_col]]
        elif fc_filter=='UPREG':
            all_cl_markers=all_cl_markers[(all_cl_markers[logFC_col]>fc_cutoff)]
            all_cl_markers=all_cl_markers.sort_values(by=logFC_col,ascending=False)[[gene_col]]
        elif fc_filter=='DOWNREG':
            all_cl_markers=all_cl_markers[(all_cl_markers[logFC_col]<fc_cutoff)]
            all_cl_markers=all_cl_markers.sort_values(by=logFC_col)[[gene_col]]

        if num_genes:
            if len(all_cl_markers[gene_col])>num_genes:
                markers_list[cl]=all_cl_markers[gene_col].head(num_genes).astype(str).tolist()
            else:
                markers_list[cl]=all_cl_markers[gene_col].astype(str).tolist()
        else:
            markers_list[cl] = all_cl_markers[gene_col].astype(str).tolist()
    return markers_list


def get_topp(gene_list:list,key_type:str,topp_categories:list,
             max_results:int = 10,min_genes:int = 5,max_genes:int = 1500,pval_cutoff:float = 0.05,correction:str = 'FDR')->DataFrame:

    if key_type != 'ENTREZ':
        new_gene_list = get_entrez(gene_list)
    else:
        new_gene_list = gene_list
    payload = dict()
    payload['Genes'] = new_gene_list
    category_list = []
    for i in range(len(topp_categories)):
        cat_dict={
            'Type':topp_categories[i],
            'PValue':pval_cutoff,
            'MinGenes':min_genes,
            'MaxGenes':max_genes,
            'MaxResults':max_results,
            'Correction':correction
        }
        category_list.append(cat_dict)

    payload['Categories']=category_list
    url = 'https://toppgene.cchmc.org/API/enrich'
    headers = {'Content-Type': 'application/json'}
    r = json.loads(requests.post(url,json=payload,headers=headers).text)
    response_data = r['Annotations']
    keepers = ['Category','ID','Name','PValue','QValueFDRBH','QValueFDRBY','QValueBonferroni',
             'TotalGenes','GenesInTerm','GenesInQuery','GenesInTermInQuery','Source','URL']
    return_df = []
    for da in response_data:
        row={}
        for k in keepers:
            row[k]=da.get(k)
        return_df.append(row)
    return DataFrame(return_df)


def topp_fun(markers:DataFrame,topp_categories:list =None,cluster_col:str = 'cluster',gene_col:str = 'gene',p_val_col:str = 'adj_p_val_col',
             logFC_col:str = 'avg_logFC',num_genes:int = 1000,pval_cutoff:float = 0.5,fc_cutoff:float = 0,fc_filter:str = 'ALL',clusters:list = None,
             correction:str = 'FDR',key_type:str = 'SYMBOL',min_genes:int = 2,max_genes:int = 1500,max_results:int = 50)->DataFrame:
    """

    The topp_fun() function takes a DataFrame and selects genes to use in querying ToppGene.

    Args:
        markers: A vector of markers or DataFrame with columns as cluster labels
        topp_categories: A string or vector with specific toppfun categories for the query
        cluster_col: Column name for the groups of cells (e.g. cluster or celltype)
        gene_col: Column name for genes (e.g. gene or feature)
        p_val_col: Column name for the p-value or adjusted p-value (preferred)
        logFC_col: Column name for the avg log FC column
        num_genes: Number of genes per group to use for toppGene query
        pval_cutoff: (adjusted) P-value cutoff for filtering differentially expressed genes
        fc_cutoff: Avg log fold change cutoff for filtering differentially expressed genes
        fc_filter: Include "ALL" genes, or only "UPREG" or "DOWNREG" for each cluster
        clusters: Which clusters to include in toppGene query
        correction: P-value correction method ("FDR" is "BH")
        key_type: Gene name format
        min_genes: Minimum number of genes to match in a query
        max_genes: Maximum number of genes to match in a query
        max_results: Maximum number of results per cluster

    Returns: DataFrame

    Examples:

        from topppy import *
        toppdata=topp_fun(ifnb_de,topp_categories=None,cluster_col='celltype',gene_col='gene',p_val_col='p_val_adj',logFC_col='avg_log2FC')

    """


    msg='''
    This function returns data generated from ToppGene (https://toppgene.cchmc.org/)
    Any use of this data must be done so under the Terms of Use and citation guide established by ToppGene.
    Terms of Use: https://toppgene.cchmc.org/navigation/termsofuse.jsp
    Citations: https://toppgene.cchmc.org/help/publications.jsp'''
    print(msg)
    markers=DataFrame(markers)
    if cluster_col not in markers.columns:
        raise ValueError(f'Cluster column {cluster_col} not found in data.Please specify')
    if clusters is not None:
        markers=markers[markers[cluster_col].isin(clusters)]
    if fc_filter not in ['ALL','UPREG','DOWNREG']:
        raise ValueError("please select one of ['ALL','UPREG','DOWNREG'] for fc_filter")
    if isinstance(markers,DataFrame):
        marker_list=process_markers(
            markers=markers,
            cluster_col=cluster_col,
            gene_col=gene_col,
            p_val_col=p_val_col,
            logFC_col=logFC_col,
            num_genes=num_genes,
            pval_cutoff=pval_cutoff,
            fc_cutoff=fc_cutoff,
            genes_submit_cutoff=1000,
            fc_filter=fc_filter
        )
    else:
        raise ValueError("data format not recognized")
    if correction not in ['none','FDR','Bonferroni']:
        raise ValueError("invalid P-value correction method,please select either 'none','FDR',or 'Bonferroni'")
    if topp_categories is None:
        topp_categories=get_topp_cats()
    big_df=DataFrame()
    missing_clusters=[]
    for col in marker_list.keys():
        if col not in ['rank','X']:
            gene_list=marker_list[col]
            if len([g for g in gene_list if pd.notna(g)])>=min_genes:
                print('Working on cluster:',col)
                d=get_topp(
                    gene_list=gene_list,
                    topp_categories=topp_categories,
                    key_type='SYMBOL',
                    pval_cutoff=pval_cutoff,
                    min_genes=min_genes,
                    max_genes=max_genes,
                    max_results=max_results,
                    correction=correction
                )
                if d.shape[0]==0:
                    missing_clusters.append(col)
                elif len(markers.columns)==1:
                    big_df=pd.concat([big_df,d],ignore_index=True)
                else:
                    d['Cluster']=col
                    big_df=pd.concat([big_df,d],ignore_index=True)
            else:
                missing_clusters.append(col)
    if len(missing_clusters)>0:
        warn_msg=f'WARNING:no results found for cluster {missing_clusters}'
        print(warn_msg,file=sys.stderr)
    return big_df
