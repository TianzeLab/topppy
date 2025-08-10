import pandas as pd
from topppy import get_topp_cats
from topppy import topp_save
from topppy import get_entrez
from topppy import process_markers,get_topp,topp_fun
import shutil
import os

import requests
import json

def test_get_topp_cat():
    assert get_topp_cats() == ["GeneOntologyMolecularFunction", "GeneOntologyBiologicalProcess",
    "GeneOntologyCellularComponent","HumanPheno","MousePheno","Domain","Pathway",
    "Pubmed","Interaction","Cytoband","TFBS","GeneFamily","Coexpression","CoexpressionAtlas",
    "ToppCell","Computational","MicroRNA","Drug","Disease"]

def test_toppsave():
    data = {
        "Gene": ["GeneA", "GeneB", "GeneC"],
        "Score": [0.1, 0.5, 0.3],
        "Cluster": ["CellType1", "CellType1", "CellType2"]
    }
    data = pd.DataFrame(data)

    shutil.rmtree(path=os.path.expanduser('~') + os.sep + 'topppy', ignore_errors=True)
    os.mkdir(os.path.expanduser('~') + os.sep + 'topppy')
    topp_save(topp_data=data, filename='test_toppsave', save_dir=os.path.join(os.path.expanduser('~'), 'topppy'), split=False, format='xlsx')
    topp_save(topp_data=data, filename='test_toppsave', save_dir=os.path.join(os.path.expanduser('~'), 'topppy'), split=False, format='csv')
    topp_save(topp_data=data, filename='test_toppsave', save_dir=os.path.join(os.path.expanduser('~'), 'topppy'), split=False, format='tsv')

    assert os.path.exists(os.path.expanduser('~') + os.sep + 'topppy' + os.sep + "test_toppsave.xlsx") and os.path.getsize(os.path.expanduser('~') + os.sep + 'topppy' + os.sep + "test_toppsave.xlsx")
    assert os.path.exists(os.path.expanduser('~') + os.sep + 'topppy' + os.sep + "test_toppsave.csv") and os.path.getsize(os.path.expanduser('~') + os.sep + 'topppy' + os.sep + "test_toppsave.csv")
    assert os.path.exists(os.path.expanduser('~') + os.sep + 'topppy' + os.sep + "test_toppsave.tsv") and os.path.getsize(os.path.expanduser('~') + os.sep + 'topppy' + os.sep + "test_toppsave.tsv")

    topp_save(topp_data=data, filename='test_toppsave', save_dir=os.path.join(os.path.expanduser('~'), 'topppy'), split=True, format='xlsx')
    topp_save(topp_data=data, filename='test_toppsave', save_dir=os.path.join(os.path.expanduser('~'), 'topppy'), split=True, format='csv')
    topp_save(topp_data=data, filename='test_toppsave', save_dir=os.path.join(os.path.expanduser('~'), 'topppy'), split=True, format='tsv')
    for gr in data['Cluster'].unique():
        assert os.path.exists(
            os.path.expanduser('~') + os.sep + 'topppy' + os.sep + f"test_toppsave_{gr}.xlsx") and os.path.getsize(
            os.path.expanduser('~') + os.sep + 'topppy' + os.sep + f"test_toppsave_{gr}.xlsx")
        assert os.path.exists(
            os.path.expanduser('~') + os.sep + 'topppy' + os.sep + f"test_toppsave_{gr}.csv") and os.path.getsize(
            os.path.expanduser('~') + os.sep + 'topppy' + os.sep + f"test_toppsave_{gr}.csv")
        assert os.path.exists(
            os.path.expanduser('~') + os.sep + 'topppy' + os.sep + f"test_toppsave_{gr}.tsv") and os.path.getsize(
            os.path.expanduser('~') + os.sep + 'topppy' + os.sep + f"test_toppsave_{gr}.tsv")


def test_get_topp_gene():
    url = "https://toppgene.cchmc.org/API/lookup"
    headers = {'Content-Type': 'application/json'}
    data = {"Symbols":["FLDB","APOE","ENSG00000113196","ENSMUSG00000020287"]}
    response = json.loads(requests.post(url, json=data, headers=headers).text)


def test_get_entrez():
    genes=["IFNG", "FOXP3"]
    a=get_entrez(genes)
    assert isinstance(a,list)


def test_topp_fun():
    ds=pd.DataFrame({
        'celltype': ['A', 'A', 'B', 'B'],
        'gene': ['G1', 'G2', 'G3', 'G4'],
        'p_val_adj': [0.01, 0.02, 0.03, 0.04],
        'avg_log2FC': [1.5, -1.2, 0.8, 2.0]
    })
    a=topp_fun(ds,topp_categories=None,cluster_col="celltype",gene_col="gene",p_val_col="p_val_adj",logFC_col="avg_log2FC")
    assert isinstance(a,pd.DataFrame)
    assert 'Cluster' in a.columns