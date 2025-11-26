from topppy import *

def test_topp_plot():
    topp_plot(topp_data,category="GeneOntologyMolecularFunction",clusters=0,save=True,file_prefix="MF_cluster0")

def test_topp_balloon():
    topp_balloon(topp_data, balloons = 3, save = True, filename = "Balloon_plot")