import math
import numpy as np
import os
from pandas import DataFrame, Categorical
from plotnine import ggplot, geom_segment, aes, geom_point, scale_color_cmap, theme_bw, ylab, ggtitle, scale_y_discrete, \
    labs, theme, element_text, xlab, scale_x_discrete, element_rect
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def topp_plot(toppdata: DataFrame, category: str, clusters: list | int | str = None, cluster_col: str = "Cluster",
              num_terms: int = 10, p_val_adj: str = "BH", p_val_display: str = "log",
              save: bool = False, save_dir: str = None, width: int = 5, height: int = 6, file_prefix: str = None,
              y_axis_text_size: int = 8, combine: bool = False, ncols: int = None) -> ggplot | dict:
    '''
    Create a dotplot from toppdata results

    Args:
        toppdata: A toppdata results dataframe
        category: The topp categories to plot
        clusters: The cluster(s) to plot
        cluster_col: The column name for clusters (default: "Cluster")
        num_terms: The number of terms from the toppdata results to be plotted, per cluster
        p_val_adj: The P-value correction method: "BH", "Bonferroni", "BY", or "none"
        p_val_display: If "log", display the p-value in terms of -log10(p_value)
        save: Whether to save the file automatically
        save_dir: Directory to save file
        width: width of the saved file (inches)
        height: height of the saved file (inches)
        file_prefix: file prefix if saving the plot - the cluster name is also added automatically
        y_axis_text_size: Size of the Y axis text - for certain categories, it's helpful to decrease this
        combine: If TRUE and multiple clusters selected, return a patchwork object of all plots; if FALSE return list of plots
        ncols: If patchwork element returned, number of columns for subplots

    Returns: A ggplot object or a dict where the keys are the names of clusters and the values are ggplot objects

    '''
    test_cols=["Category","Name","PValue","GenesInTerm","GenesInQuery","GenesInTermInQuery"]
    for t in test_cols:
        if t not in toppdata.columns:
            raise ValueError("The column",t,"is missing from the toppData,please correct and retry." )
    GROUPBY_COL=cluster_col
    if GROUPBY_COL not in toppdata.columns:
        raise ValueError("Invalid cluster column:",GROUPBY_COL,"not found in toppdata.Select an existing column with cluster_col=")

    if category not in toppdata["Category"].unique():
        raise ValueError("Category",category,"not found in the data.Please select one of ",", ".join(toppdata["Category"].unique()))
    elif category is None:
        raise ValueError("Please select one of these categories:",", ".join(toppdata["Category"].unique()))

    if clusters is None:
        if "Cluster" in toppdata.columns:
            clusters=toppdata["Cluster"].unique()

    tmp_data=toppdata
    if p_val_adj not in ["BH","Bonferroni","BY","none","None","log"]:
        print("P value adjustment not found - using 'BH' by default. For no adjustment, use p_val_adj = 'none'.")

    if p_val_adj=="BH":
        p_val_col="QValueFDRBH"
    elif p_val_adj=="Bonferroni":
        p_val_col ="QValueBonferroni"
    elif p_val_adj =="BY":
        p_val_col ="QvalueFDRBY"
    elif p_val_adj =="none":
        p_val_col ="PValue"
    elif p_val_adj =="None":
        p_val_col ="PValue"
    else:
        p_val_col ="QValueFDRBH"

    if p_val_display=='log':
        color_label="-log10(FDR)"
        tmp_data["nlog10_fdr"]=-np.log10(tmp_data[p_val_col])
        p_val_display_column="nlog10_fdr"
    elif p_val_col in ["QValueFDRBH","QValueBonferroni","QvalueFDRBY"]:
        color_label="Adj. P-value"
        p_val_display_column=p_val_col
    else:
        color_label="P-value"
        p_val_display_column = p_val_col

    if save is True:
        if save_dir is None:
            output_dir=os.getcwd()
        else:
            output_dir=save_dir

    if isinstance(clusters, int) or isinstance(clusters, str):
        clusters = [clusters]

    if len(clusters)>1:
        if combine is True:
            print("Multiple clusters entered: function returns a list of ggplots\n")
        overall_plot_list=dict()
        for c in clusters:
            df_C=tmp_data[tmp_data[GROUPBY_COL]==c]
            df_C=df_C[df_C["Category"]==category]
            df_C["geneRatio"]=df_C["GenesInTermInQuery"]/df_C["GenesInTerm"]
            df_C=df_C.sort_values(by=p_val_display_column,ascending=False).head(num_terms)

            p=(ggplot(df_C,aes(x="geneRatio",y="reorder(Name,geneRatio)"))
               + geom_segment(aes(xend=0,yend="reorder(Name, geneRatio)"))
               + geom_point(mapping=aes(size="GenesInTermInQuery",color=p_val_display_column))
               + scale_color_cmap(name="-Log10(FDR)",cmap_name="plasma")
               + theme_bw()
               + ylab(category)
               + ggtitle(f"Cluster {c}")
               + theme(axis_text_y=element_text(size=y_axis_text_size))
               + scale_y_discrete(labels=lambda x:[textwrap.fill(i,20) for i in x])
               + labs(color=color_label,size="Genes from Query\n in Gene Set"))
            overall_plot_list[c] = p
            if save is True:
                if file_prefix is None:
                    save_filename=f"{category}_{c}_toppDotPlot.pdf"
                else:
                    save_filename=f"{file_prefix}_{c}.pdf"
                p.save(os.path.join(output_dir,save_filename),width=width,height=height)
        if combine is True:
            if ncols is None:
                ncols=min(3,len(overall_plot_list))
            nrows = math.ceil(len(overall_plot_list) / ncols)
            plot_values=list(overall_plot_list.values())
            combined_plots=None
            for i in range(0,len(plot_values),ncols):
                row_plots=plot_values[i:i+ncols]
                current_row=row_plots[0]
                for p_next in row_plots[1:]:
                    current_row=current_row|p_next
                if combined_plots is None:
                    combined_plots=current_row
                else:
                    combined_plots=combined_plots/current_row
            plt.ioff()
            fig=combined_plots.draw(show=False)
            fig.set_size_inches(width*ncols, height*nrows)
            fig.suptitle(category,fontsize=16)
            if save:
                save_filename = f"{file_prefix}_combined.pdf"
                combined_plots.save(os.path.join(output_dir, save_filename), width=width*ncols, height=height*nrows)
            return combined_plots
        else:
            return overall_plot_list
    elif len(clusters)==1:
        c=clusters[0]
        df_C1 = tmp_data[tmp_data[GROUPBY_COL] == c]
        df_C1 = df_C1[df_C1["Category"] == category]
        df_C1["geneRatio"] = df_C1["GenesInTermInQuery"] / df_C1["GenesInTerm"]
        df_C1 = df_C1.sort_values(by=p_val_display_column, ascending=False).head(num_terms)

        single_plot = (ggplot(df_C1, aes(x="geneRatio", y="reorder(Name, geneRatio)"))
             + geom_segment(aes(xend=0, yend="reorder(Name, geneRatio)"))
             + geom_point(mapping=aes(size="GenesInTermInQuery", color=p_val_display_column))
             + scale_color_cmap(name="-Log10(FDR)",cmap_name="plasma")
             + theme_bw()
             + ylab(category)
             + ggtitle(f"Cluster {c}")
             + theme(axis_text_y=element_text(size=y_axis_text_size))
             + scale_y_discrete(labels=lambda x: [textwrap.fill(i, 20) for i in x])
             + labs(color=color_label, size="Genes from Query\n in Gene Set"))

        if save is True:
            if file_prefix is None:
                save_filename = f"{category}_{c}_toppDotPlot.pdf"
            else:
                save_filename = f"{file_prefix}_{category}_{c}_toppDotPlot.pdf"
            single_plot.save(os.path.join(output_dir, save_filename), width=width, height=height)
        return single_plot
    else :
        df_C2 = tmp_data[tmp_data["Category"] == category]
        df_C2["geneRatio"] = df_C2["GenesInTermInQuery"] / df_C2["GenesInTerm"]
        df_C1 = df_C2.sort_values(by=p_val_display_column, ascending=False).head(num_terms)
        single_plot = (ggplot(df_C1, aes(x="geneRatio", y="reorder(Name, geneRatio)"))
                       + geom_segment(aes(xend=0, yend="reorder(Name, geneRatio)"))
                       + geom_point(mapping=aes(size="GenesInTermInQuery", color=p_val_display_column))
                       + scale_color_cmap(name="-Log10(FDR)",cmap_name="plasma")
                       + theme_bw()
                       + ylab(category)
                       + theme(axis_text_y=element_text(size=y_axis_text_size))
                       + scale_y_discrete(labels=lambda x: [textwrap.fill(i, 20) for i in x])
                       + labs(color=color_label, size="Genes from Query\n in Gene Set"))

        if save is True:
            if file_prefix is None:
                save_filename = f"{category}_toppDotPlot.pdf"
            else:
                save_filename = f"{file_prefix}_{category}_toppDotPlot.pdf"
            single_plot.save(os.path.join(output_dir, save_filename), width=width, height=height)
        return single_plot


def topp_balloon(toppdata: DataFrame, categories: list = None, balloons: int = 3, x_axis_text_size: int = 6,
                 cluster_col: str = "Cluster", filename: str = None, save: bool = False, height: int = 5,
                 width: int = 10) -> ggplot|dict:
    '''
    Create a balloon plot from toppdata results

    Args:
        toppdata: A toppdata results dataframe
        categories: The topp categories to plot
        balloons: Number of balloons per group to plot
        x_axis_text_size: Size of the text on the x axis
        cluster_col: The column name for clusters (default: "Cluster")
        filename: Filename of the saved balloon plot
        save: Save the balloon plot if TRUE
        height: Height of the saved balloon plot
        width: Width of the saved balloon plot

    Returns: A ggplot object or a dict where the keys are the names of clusters and the values are ggplot objects

    '''
    if categories is None:
        categories=toppdata["Category"].unique()
    GROUPBY_COL=cluster_col
    if GROUPBY_COL not in toppdata.columns:
        raise ValueError(f"Invalid cluster column:{GROUPBY_COL} not found in toppdata")
    balloon_list={}
    for cat in categories:
        print("Balloon Plot:",cat)
        df=toppdata[toppdata["Category"]==cat]
        df["nlog10_fdr"]=-np.log10(df["QValueFDRBH"])
        df["geneRatio"]=df["GenesInTermInQuery"]/df["GenesInTerm"]
        df=df.groupby(GROUPBY_COL, group_keys=False).apply(lambda x: x.nlargest(balloons, "nlog10_fdr")).reset_index(drop=True)
        df["geneRatio"] = df["GenesInTermInQuery"] / df["GenesInTerm"]
        order_Names = df.sort_values(GROUPBY_COL)["Name"].dropna().unique()
        df["Name"]=Categorical(df["Name"], categories=order_Names, ordered=True)
        p=(ggplot(df,aes(x="Name",y=GROUPBY_COL))
           + geom_point(aes(size="geneRatio",color="nlog10_fdr"))
           + scale_color_cmap(name="-Log10(FDR)",cmap_name="plasma")
           + labs(color="-Log10(FDR)",size="Gene Ratio")
           + xlab(cat)
           + theme_bw()
           + scale_x_discrete(labels=lambda x: [textwrap.fill(str(t), 30) for t in x])
           + theme(axis_text_x=element_text(size=x_axis_text_size,angle=60,ha="right"),panel_border=element_rect(color=None)))
        if filename is None:
            plot_filename=f"toppBalloon_{cat}.pdf"
        else:
            plot_filename=f"{filename}_{cat}.pdf"
        if save is True:
            p.save(plot_filename,height=height,width=width)
            print("save to:",plot_filename)
        print("\n")
        balloon_list[cat]=p
    if len(balloon_list)==1:
        return list(balloon_list.values())[0]
    else:
        return balloon_list



