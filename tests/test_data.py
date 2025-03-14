from pandas import DataFrame

from topppy import *


def test_data():
    assert type(ifnb_de) is DataFrame
    assert type(pbmc_markers) is DataFrame
    assert type(topp_data) is DataFrame
