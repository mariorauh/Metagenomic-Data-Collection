import argparse as ap
import pandas as pd
from pandas import read_csv
import matplotlib.pyplot as plt
import numpy as np
import ast

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def command_line():
    parser = ap.ArgumentParser("Metagenomic Data Collection (via MG-Rast) - Metadata Analysis")

    parser.add_argument("-i", "--input", help="Metadata csv file that shall be analyzed.",
                        required=True)
    #parser.add_argument("-o", "--output", help="Output File name", default="metadata_analyzer.png")
    #parser.add_argument("-f", "--format", help="Choose the output's format. Default: pdf", default="png")

    parser.add_argument("--keyword_pie", help="Use this to create a pie chart containing all keyword information.",
                        action='store_true')
    parser.add_argument("--keyword_pie_out", help="Optional: Set output file name for the keyword pie chart. Default: "
                                                  "piechart.pdf",
                        default="piechart", type=str)
    parser.add_argument("--keyword_bar", help="Use this to create a pie chart containing all keyword information.",
                        action='store_true')
    parser.add_argument("--keyword_bar_out", help="Optional: Set output file name for the keyword bar chart. Default: "
                                                  "barchart.pdf",
                        default="barchart", type=str)
    parser.add_argument("--alpha_diversity", help="Compute a Boxplot based on the alpha"
                                                  " diversity for each group of keywords",
                        action='store_true')
    parser.add_argument("--alpha_diversity_out", help="Optional: Set output file name for the alpha diversity boxplot."
                                                    " Default: alpha_div.pdf",
                        default="alpha_div", type=str)
    parser.add_argument("--rarefaction_curve", help="Compute a Boxplot based on the rarefaction curve "
                                                    "for each group of keywords.",
                        action='store_true')
    parser.add_argument("--rarefaction_curve_out", help="Optional: Set output file name for the rarefaction curve "
                                                        "boxplot. Default: rc_box.pdf",
                        default="rc_box", type=str)
    parser.add_argument("--sequence_count_raw", help="Compute a Boxplot based on the raw sequence count "
                                                    "for each group of keywords.",
                        action='store_true')
    parser.add_argument("--sequence_count_raw_out", help="Optional: Set output file name for the boxplot of "
                                                         "raw sequence count. "
                                                         "Default: seq_count.pdf",
                        default="seq_count", type=str)
    parser.add_argument("--species_count", help="Compute a Boxplot on the species count for each group of keywords.",
                        action='store_true')
    parser.add_argument("--species_count_out", help="Optional: Set output file name for the boxplot of the species"
                                                    " count. Default: species_count.pdf")

    return vars(parser.parse_args())


def keyword_graphs(df:pd.DataFrame, keyword_pie:bool, keyword_bar:bool, pie_out:str, bar_out:str):
    '''
    
    :param df: Dataframe containing all metadata
    :param keyword_pie: boolean that determines if the user would like to get a pie chart from the data
    :param keyword_bar: boolean that determines if the user would like to get a bar chart from the data 
    :return: file(s) containing the graph
    '''
    keywords = df['keyword'].tolist()
    counts = {}
    for key in keywords:
        if key in counts:
            counts[key] = counts[key] + 1
        else:
            counts[key] = 1


    #keys = [str(x) for x in counts.keys()]
    keys = []
    for x in list(counts.keys()):
        temp = ""
        x = ast.literal_eval(x)
        for i in x:
            temp+=f"{str(i[1])}, "
        keys.append(temp)

    values = np.fromiter(counts.values(), dtype=int)
    counter = 0
    # if keyword_pie true => create pie chart of keywords
    if keyword_pie:
        plt.pie(values, labels=keys)
        plt.tight_layout()
        plt.savefig(f"{pie_out}.pdf")
        print(f"{bcolors.OKGREEN}Pie Chart created successfully and saved to {pie_out}.pdf{bcolors.ENDC}")
        counter+=1
        #plt.show()

    if keyword_bar:
        plt.barh(keys, values)
        plt.xlabel("Counts")
        plt.tight_layout()
        plt.savefig(f"{bar_out}.pdf")
        print(f"{bcolors.OKGREEN}Bar Chart created successfully and saved to {bar_out}.pdf{bcolors.ENDC}")
        counter+=1
        #plt.show()

    return counter

## To DO: Create Boxplot for each keyword group for alpha diversity , species count & sequence count
def alpha_diversity(df:pd.DataFrame, alpha_out:str):
    '''

    :param df: Dataframe containing all metadata
    :param alpha_diversity: boolean that determines if the user would like to get an overview of alpha diversity
    :return: a file containing the graphs
    '''

    alpha_values = df["alpha_diversity_shannon"].tolist()
    keywords = df['keyword'].tolist()

    alpha_dict = {}

    for i, key in enumerate(keywords):
        #print(type(alpha_values[i]))
        if str(key) in alpha_dict:
            temp = alpha_dict[key]
            temp.append(alpha_values[i])
            alpha_dict[str(key)] = temp
        else:
            alpha_dict[str(key)] = [alpha_values[i]]

    #print(alpha_dict)
    values = [x for x in alpha_dict.values()]
    #print(values)

    plt.boxplot(values, vert=False)
    val_length = len(values)
    ticks = [i for i in range(val_length+1)]
    tick_names = [y for y in alpha_dict.keys()]

    keys = [""]
    for x in list(alpha_dict.keys()):
        temp = ""
        x = ast.literal_eval(x)
        for i in x:
            temp+=f"{str(i[1])}, "
        keys.append(temp)

    plt.yticks(ticks, keys)
    plt.xlabel("Alpha Diversity")
    plt.tight_layout()
    plt.savefig(f"{alpha_out}.pdf")
    print(f"{bcolors.OKGREEN}Alpha Diversity Boxplot created successfully and saved to {alpha_out}.pdf{bcolors.ENDC}")


def rarefaction_analyses(df:pd.DataFrame, rc_out:str):

    slope_values = df["RC_slope"].tolist()
    keywords = df['keyword'].tolist()

    slope_dict = {}

    for i, key in enumerate(keywords):
        #print(type(alpha_values[i]))
        if str(key) in slope_dict:
            temp = slope_dict[key]
            temp.append(slope_values[i])
            slope_dict[str(key)] = temp
        else:
            slope_dict[str(key)] = [slope_values[i]]

    #print(alpha_dict)
    values = [x for x in slope_dict.values()]
    #print(values)

    plt.boxplot(values, vert=False)
    val_length = len(values)
    ticks = [i for i in range(val_length+1)]
    tick_names = [y for y in slope_dict.keys()]

    keys = [""]
    for x in list(slope_dict.keys()):
        temp = ""
        x = ast.literal_eval(x)
        for i in x:
            temp+=f"{str(i[1])}, "
        keys.append(temp)

    plt.yticks(ticks, keys)
    plt.xlabel("RC Slopes")
    plt.tight_layout()
    plt.savefig(f"{rc_out}.pdf")
    print(f"{bcolors.OKGREEN}Rarefaction Curve Slopes boxplot created successfully and saved to {rc_out}.pdf"
          f"{bcolors.ENDC}")


def seq_count_raw(df:pd.DataFrame, seq_out:str):

    seq_values = df["sequence_count_raw"].tolist()
    keywords = df['keyword'].tolist()

    seq_dict = {}

    for i, key in enumerate(keywords):
        #print(type(alpha_values[i]))
        if str(key) in seq_dict:
            temp = seq_dict[key]
            temp.append(seq_values[i])
            seq_dict[str(key)] = temp
        else:
            seq_dict[str(key)] = [seq_values[i]]

    #print(alpha_dict)
    values = [x for x in seq_dict.values()]
    #print(values)

    plt.boxplot(values, vert=False)
    val_length = len(values)
    ticks = [i for i in range(val_length+1)]
    tick_names = [y for y in seq_dict.keys()]

    keys = [""]
    for x in list(seq_dict.keys()):
        temp = ""
        x = ast.literal_eval(x)
        for i in x:
            temp+=f"{str(i[1])}, "
        keys.append(temp)

    plt.yticks(ticks, keys)
    plt.xlabel("Sequence Count Raw")
    plt.tight_layout()
    plt.savefig(f"{seq_out}.pdf")
    print(f"{bcolors.OKGREEN}Sequence Count Raw boxplot created successfully and saved to {seq_out}.pdf{bcolors.ENDC}")


def species_count_boxplot(df:pd.DataFrame, species_out:str):

    species_values = df["species_count"].tolist()
    keywords = df['keyword'].tolist()

    species_dict = {}

    for i, key in enumerate(keywords):
        #print(type(alpha_values[i]))
        if str(key) in species_dict:
            temp = species_dict[key]
            temp.append(species_values[i])
            species_dict[str(key)] = temp
        else:
            species_dict[str(key)] = [species_values[i]]

    #print(alpha_dict)
    values = [x for x in species_dict.values()]
    #print(values)

    plt.boxplot(values, vert=False)
    val_length = len(values)
    ticks = [i for i in range(val_length+1)]
    tick_names = [y for y in species_dict.keys()]

    keys = [""]
    for x in list(species_dict.keys()):
        temp = ""
        x = ast.literal_eval(x)
        for i in x:
            temp+=f"{str(i[1])}, "
        keys.append(temp)

    plt.yticks(ticks, keys)
    plt.xlabel("Species Count")
    plt.tight_layout()
    plt.savefig(f"{species_out}.pdf")
    print(f"{bcolors.OKGREEN}Sequence Count Raw boxplot created successfully and saved to {species_out}.pdf{bcolors.ENDC}")

def main():

    print("Analyzing data ... ")

    args = command_line()
    df, keyword_pie, keyword_bar = read_csv(args["input"]), args["keyword_pie"], args["keyword_bar"]
    alpha_div = args["alpha_diversity"]
    rc = args["rarefaction_curve"]
    seq_count = args["sequence_count_raw"]
    pie_out, bar_out, alpha_out = args["keyword_pie_out"], args["keyword_bar_out"], args["alpha_diversity_out"]
    rc_out, seq_out = args["rarefaction_curve_out"], args["sequence_count_raw_out"]
    species_count, species_out = args["species_count"], args["species_count_out"]
    counter = keyword_graphs(df, keyword_pie, keyword_bar, pie_out, bar_out)
    if alpha_div:
        alpha_diversity(df, alpha_out)
        counter+=1

    if rc:
        rarefaction_analyses(df, rc_out)
        counter+=1

    if seq_count:
        seq_count_raw(df, seq_out)
        counter+=1

    if species_count:
        species_count_boxplot(df, species_out)
        counter+=1

    print(f"{bcolors.OKGREEN}{counter} File(s) was/were created!{bcolors.ENDC}")

if __name__ == '__main__':
    main()