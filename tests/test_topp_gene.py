import pandas as pd
from topppy import get_topp_cats
from topppy import topp_save
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
    data=pd.DataFrame(data)
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