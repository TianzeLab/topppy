"""

ifnb_de:

pbmc_markers:

topp_data:
"""
from pandas import read_csv
import os


project_dir = os.path.split(os.path.realpath(__file__))[0]

ifnb_de = read_csv(project_dir + os.sep + 'ifnb_de.csv', index_col='index')
pbmc_markers = read_csv(project_dir + os.sep + 'pbmc_markers.csv', index_col='index')
topp_data = read_csv(project_dir + os.sep + 'topp_data.csv', index_col='index')