import csv
import argparse as ap
from pandas import read_csv, concat, DataFrame


def command_line():
    parser = ap.ArgumentParser("Metagenomic Data Collection (via MG-Rast) - Check CSV files for duplicates")

    parser.add_argument("-i", "--input", help="List of input files containing metadata information in the format "
                                              "provided in the previous Data Collection Step.",
                        nargs='+',
                        required=True)
    parser.add_argument("-o", "--output", help="The output file's name", default="metadata.csv")

    return vars(parser.parse_args())

def export_metagenome_ids(files:list):

    all_metagenomes = []
    all_dfs = []
    for file in files:
        df = read_csv(file)
        #print(df)
        mgms = df['metagenome_id'].tolist()
        counter = 0
        print_duplicates(mgms)
        duplicate_rows = []
        for mgm in mgms:
            if mgm in all_metagenomes:
                #print(f"{mgm} is already present.")
                duplicate_rows.append(counter)
                #df = df[df.metagenome_id == mgm]
            else:

                all_metagenomes.append(mgm)

            counter+=1

        df = df.drop(duplicate_rows)
        #print(df)
        all_dfs.append(df)
    print_duplicates(all_metagenomes)
    final_df = concat(all_dfs, ignore_index=True)
    #print(final_df)
    return final_df

def print_duplicates(all_mgms:list):

    temp = []
    for mgm in all_mgms:
        if mgm in temp:
            print(f"{mgm} is a duplicate.")
        else:
            temp.append(mgm)

    return temp

def main():
    print("Checking for common samples")
    args = command_line()
    files, output = args["input"], args["output"]
    print(files)
    final_df = export_metagenome_ids(files)
    final_df.to_csv(output, index=False)
    print(f"File saved to {output}")

    
if __name__ == '__main__':
    main()