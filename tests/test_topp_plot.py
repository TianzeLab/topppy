from topppy import *

def test_topp_plot():
    topp_plot(topp_data,category="GeneOntologyMolecularFunction",clusters=0,save=False,file_prefix="MF_cluster0")
    plot_list_combined = topp_plot(topp_data, category="GeneOntologyMolecularFunction", clusters=None, save=False, combine=True,
                          file_prefix="GO_molecular_function", ncols=3)
    plot_list = topp_plot(topp_data, category="GeneOntologyMolecularFunction", clusters=None, save=False, combine=False,
                          file_prefix="GO_molecular_function", ncols=3)
    plot_list[0]|plot_list[1]
    

def test_topp_balloon():
    topp_balloon(topp_data, balloons = 3, save = True, filename = "Balloon_plot")