"""

ifnb_de:A dataframe of differentially expressed genes generated using the FindMarkers function for each cluster from the Kang 2018 IFNB dataset

pbmc_markers:A dataframe of marker genes generated using the

topp_data:A dataframe of of sample toppData results
"""
from pandas import read_csv
import os


project_dir = os.path.split(os.path.realpath(__file__))[0]

ifnb_de = read_csv(project_dir + os.sep + 'ifnb_de.csv', index_col='index')
pbmc_markers = read_csv(project_dir + os.sep + 'pbmc_markers.csv', index_col='index')
topp_data = read_csv(project_dir + os.sep + 'topp_data.csv', index_col='index')
