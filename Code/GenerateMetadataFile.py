import json as js
import csv
import os
import subprocess
from pathlib import Path
import argparse as ap
import time
import ast
import numpy as np
import requests


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
    parser = ap.ArgumentParser("Metagenomic Data Collection (via MG-Rast)")

    # adapted from stackoverflow.com
    def check_positive(value):
        ivalue = int(value)
        if ivalue < 0:
            raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
        return ivalue

    parser.add_argument("-o", "--output", help="Name of the Output csv file. Default: metadata", default="metadata")
    parser.add_argument("-l", "--limit", help="Maximum number of API Search results from MG-Rast. "
                                              " Multiple input values possible. See Readme for Example. Default: 40",
                        default=[40], type=list, action='append', nargs='+')
    parser.add_argument("-m", "--metadata", help="More information in Readme. Multiple Text words should be surrounded"
                                                 " by \"...\". Ensure, that the metadata field keywords match the key"
                                                 " words from the official MG-Rast API Search: "
                                                 "(https://www.mg-rast.org/mgmain.html?mgpage=searchapi). "
                                                 "Please enter the desired keywords in the following form "
                                                 "(Note: Add multiple requests by adding multiple \"-m\"): "
                                                 "-m <metadata field 1> <text 1> <metadata field 2> <text 2> ...",
                        type=list, required=True, action='append', nargs='+')
    parser.add_argument("-pd", "--public_data", help="Include Searching public data.",
                        action='store_true')
    parser.add_argument("--desc", help="Include to sort results in descending order. Default: Ascending",
                        action='store_true')
    parser.add_argument("--order_field", help="Select a metadata field that determines the ordering."
                                              "Default: \"created_on\"", default="created_on")
    parser.add_argument("-j", "--json", help="Possible Labeling of the json files. NOT Required. Note: This only works"
                                             ", if # of json arguments matches the number of \"-m\"."
                                             " Default: request_[1-n]", default=["request"], type=list, nargs='+')
    parser.add_argument("-r", "--rarefaction_threshold", help="Assign a threshold for the rarefaction curve to"
                                                              " increase the sensitivity of the files to be downloaded."
                                                              " Note: The smaller the threshold, the better the"
                                                              " sequencing depth.",
                        default=0.5, type=float
                        )
    parser.add_argument("--min_species_count", help="Set a minimum number of species count. Any dataset with less "
                                                    " species count will be discarded. ",
                        default=1000, type=check_positive)
    parser.add_argument("--set_min_readNumber", help="Set a minimum number of reads that each dataset shall consist of."
                        , default=1000000, type=check_positive)
    parser.add_argument("--phylogeny", help="If this option is selected, 16S rRNA will be considered as well, else not."
                        , action='store_true')
    parser.add_argument("--ignore_rc", help="If this option is selected, the slope of the rarefaction curve will not "
                                            "be considered anymore in the acceptance process for a dataset. "
                                            "Find out more in the Readme."
                        , action='store_true')
    parser.add_argument("--no_duplicate_proj", help="If this option is selected, the resulting metadata file will "
                                                    "only contain each project, once.",
                        action='store_true')

    return vars(parser.parse_args())


# create a dictionary that stores all the user defined requests
def metadata_converter(metadata:list):
    '''
    :param metadata: User input from the --metadata command
    :return: dict. keys = number of requests. values = metadata information for the request.
    '''
    requests = {}
    counter = 1
    for i in metadata:
        all = []
        for j in i:
            all.append(''.join(j))

        sorted_all = []
        for z,m in enumerate(all):
            if z % 2 == 1:
                sorted_all.append([all[z-1],all[z]])

        requests[counter] = sorted_all
        counter+=1

    return requests

def generate_curl_request(metadata_fields:list, limit:int, desc:bool, spd:bool, ordered_by:str):
    '''
    :param metadata_fields: list of metadata fields
    :param limit: int
    :param desc: boolean. True => Descending order
    :param spd: boolean. True => Search Public data
    :param ordered_by: string => metadata field that determines the ordering
    :return: a string containing the final curl request.
    Example: curl -F "limit=5" -F "order=created_on" -F "direction=asc" -F "public=yes" "https://api.mg-rast.org/search"
    '''
    curl = "curl  -F "
    curl+=f"\"limit=5\" -F \"order={ordered_by}\" "
    if desc:
        curl+=f"-F \"direction=desc\" "
    else:
        curl += f"-F \"direction=asc\" "

    if spd:
        curl+=f"-F \"public=yes\" "
    else:
        curl += f"-F \"public=no\" "

    for z,fields in enumerate(metadata_fields):
        curl += f"-F \"{fields[0]}={fields[1]}\" "
        #for j in fields:
        #    print(j)
        #    curl+=f"-F \"{j[0]}={j[1]}\" "

    curl+="\"https://api.mg-rast.org/search\""

    return curl

def generate_all_curls(metadata_fields:dict, limits:list, desc:bool, spd:bool, ordered_by:str):
    '''
    :param metadata_fields: dict containing all metadata fields for all requests
    :param limit: int
    :param desc: boolean. True => Descending order
    :param spd: boolean. True => Search Public data
    :param ordered_by: string => metadata field that determines the ordering
    :return: a dict containing all final curl request.
    '''
    all_curls = {}
    for i,key in enumerate(metadata_fields.keys()):
        all_curls[key] = generate_curl_request(metadata_fields[key], limits[i], desc, spd, ordered_by)

    return all_curls

def run_curls(all_curls:dict, json:list):
    '''
    Example: curl  -F "limit=5" -F "order=created_on" -F "direction=asc" -F "public=yes" -F "all=soil" "https://api.mg-rast.org/search"
    :param all_curls: dict containing all curls
    :return: creates json files based on the requests and names it: request_<n>.json
    '''
    file_list = []
    counter = 1
    for key in all_curls.keys():
        print(all_curls[key]+f" -o {json[key-1]}.json")
        os.system(f"{all_curls[key]} -o {json[key-1]}")
        print(f"{bcolors.OKGREEN}File saved to {json[key-1]}{bcolors.ENDC}")
        file_list.append(f"{json[key-1]}")
        counter+=1

    return file_list


def check_next(temp:dict):

    try:
        next_curl = temp['next']  # need this because curl is limited to 1000 datasets information
        n = True
        return next_curl, n
    except:
        n = False
        next_curl = "No more next"
        return next_curl, n

# import a given json file into the program
def import_json(file):


    if not file.startswith('.'):    # need this to ignore hidden files in MacOS

        all_results = {}
        with open(file) as json_file:

            temp = js.loads(json_file.read())

            next_curl, n = check_next(temp)

            for p in temp['data']:

                b = False

                for i in p:
                    if '16S' in str(p[i]) or '16s' in str(p[i]):
                        b = True

                if b:
                    continue

                try:
                    results = [p['metagenome_id'],p['project_name'],p['project_id'],p['biome'],p['country'],p['material'],p['feature'],
                                p['sequence_type'],p['seq_meth'],p['sequence_count_raw'],p['alpha_diversity_shannon']]
                except:
                    continue
                try:
                    results.append(p['env_package_name'])
                except:
                    results.append('None')

                all_results[p['metagenome_id']] = results

        return all_results, n, next_curl


def create_metadata(file_list:list, output:str, threshold:float, limits:list, min_species:int, metadata:dict, p:bool, ignore_slope:bool, min_reads:int, no_dup_proj:bool):
    '''
    Import all the previously created json files
    :param file_list: list of json files that were created previously
    :return: a file with metagenome information and a list with all unique metagenomic ids
    '''

    with open(f'{output}.csv', 'w', newline='') as csvfile:
        csvfile_writer = csv.writer(csvfile, delimiter=',')
        csvfile_writer.writerow(['metagenome_id','project_name','project_id','biome','country','material','feature',
                                 'sequence_type','seq_meth','sequence_count_raw','alpha_diversity_shannon',
                                 'env_package_name','species_count','RC_slope','keyword'])
        counter = 1
        metagenome_ids = []
        good_ids = []
        project_ids = []
        for k, file in enumerate(file_list):

            print(file)
            temp_ids = []
            d, n, next_curl = import_json(file)
            while d is not None:
                for key in d.keys():
                    b = False
                    for i in d[key]:
                        if not p and ("16S" in str(i) or "16s" in str(i)):
                            b = True
                            break
                        else:
                            mgm = d[key][0]
                            #print(mgm)
                            project_n = d[key][2]
                            # change this to mgm not in metagenome ids if projects can be the same
                            if mgm not in metagenome_ids:   # ensure that each metagenome only appears once
                                if no_dup_proj:
                                    if project_n not in project_ids:
                                        project_ids.append(project_n)

                                        metagenome_ids.append(mgm)
                                        curl = f"https://api-ui.mg-rast.org/metagenome/{mgm}" \
                                               f"?verbosity=stats&detail=rarefaction"
                                        rarefactions = os.popen(
                                            f"curl \"{curl}\"").read()  # creates String in list format
                                        rarefactions = ast.literal_eval(
                                            rarefactions)  # turn nested string list into actual list
                                        # print(type(rarefactions))
                                        # rarefactions = os.system(f'curl \"{curl}\"')

                                        try:
                                            r_co, grad, species_count = check_rarefaction(rarefactions, threshold, min_species, min_reads, ignore_slope)
                                        except:
                                            break
                                        # rarefaction curve coefficient

                                        if r_co:
                                            print(r_co)
                                            temp_ids.append(mgm)
                                            good_ids.append(mgm)
                                            d_key = d[key]
                                            d_key.append(species_count)
                                            d_key.append(grad)
                                            print(metadata)
                                            d_key.append(metadata[counter])

                                            csvfile_writer.writerow(d_key)

                                            if len(temp_ids) == limits[k]:
                                                break

                                else:
                                    metagenome_ids.append(mgm)
                                    curl = f"https://api-ui.mg-rast.org/metagenome/{mgm}" \
                                           f"?verbosity=stats&detail=rarefaction"
                                    rarefactions = os.popen(f"curl \"{curl}\"").read() # creates String in list format
                                    rarefactions = ast.literal_eval(rarefactions) # turn nested string list into actual list
                                    #print(type(rarefactions))
                                    #rarefactions = os.system(f'curl \"{curl}\"')

                                    r_co, grad, species_count = \
                                        check_rarefaction(rarefactions, threshold, min_species, min_reads, ignore_slope)
                                # rarefaction curve coefficient

                                    if r_co:

                                        #curl2 = f"curl \"https://api-ui.mg-rast.org/download/{mgm}?file=299.1\" > "


                                        print(r_co)
                                        temp_ids.append(mgm)
                                        good_ids.append(mgm)
                                        d_key = d[key]
                                        d_key.append(species_count)
                                        d_key.append(grad)
                                        print(metadata)
                                        d_key.append(metadata[counter])

                                        csvfile_writer.writerow(d_key)

                                        if len(temp_ids) == limits[k]:
                                            break

                    if len(temp_ids) == limits[k]:
                        return

                # if there exists a next, then do the following
                if n:
                    os.system(f"curl \"{next_curl}\" -o {file}")
                    d, n, next_curl = import_json(file)
                    #d = dict(os.popen(f"curl \"{next_curl}\""))
                    #print(f"This is d:{d}")

                else:
                    break

            counter+=1


    return good_ids


def check_rarefaction(r:list, threshold:float, min_species:int, min_reads:int, ignore_slope:bool):
    '''
    #print(f"This is the rarefactions number {rarefactions}")
    # we now use the numpy function gradient to compare all y values of r.
    r = np.array(r)
    #x = r[:,0]
    y = r[:,1]  # for each sublist in r, take the second value of the list (representing the vy value and assign it to y
    grad = np.gradient(y,3)
    last = grad[-1]
    if last<threshold: return True, last
    else: return False, last

    #print(number_rarefactions)
    '''
    try:
        r = np.array(r)

        x = r[:,0]
        y = r[:,1]

        y2 = y[-1]
        y1 = y[-3]

        max_numb_read = x[-1]

        if y2 < y1:
            y2 = y[-2]
            y1 = y[-4]
        print(min_species)
        slope = y2 - y1
        if y2 < min_species:
            return False, slope, y2
        elif max_numb_read < min_reads:
            return False, slope, y2
        else:
        # if ignore_slope is true => return True, slope, y2
            if ignore_slope:
                return True, slope, y2

            if slope < threshold:
                return True, slope, y2
            else:
                return False, slope, y2
    except:
        return False, 0,0



def json_name_converter(json:list, metadata:list):

    res = []
    for i in json:
        temp = ''.join(i)
        res.append(temp)

    all_res = []
    if res[0] == "request":
        for j in range(len(metadata)):
            all_res.append(f"request_{j}.json")

        return all_res
    else:
        return res

def limit_config(limits:list, metadata:dict):

    limits_number = len(limits)
    metadata_number = len(list(metadata.keys()))

    if metadata_number > 1:

        if limits_number == 1:  # if no limit was set
            all_limits = [limits[0] for i in range(0,metadata_number)]
            #print(all_limits)   # return all limits set to the same limit because no limit was chosen
            return all_limits

        elif limits_number < metadata_number:

            all_limits = [int(''.join(item)) for i, item in enumerate(limits[1])]
            for i in range(0,metadata_number-limits_number):
                all_limits.append(40)

            return all_limits

        else:
            all_limits = [int(''.join(item)) for i,item in enumerate(limits[1]) if i < metadata_number]
            return all_limits

    else:

        if limits_number == 1:
            return limits

        else:
            limit = [int(''.join(item)) for i,item in enumerate(limits[1]) if i < 1]
            return limit



def main():

    start = time.time()
    args = command_line()
    output, metadata, spd, limit, desc, ordered_by, min_species_count = args["output"], metadata_converter(args["metadata"]),\
                                                     args["public_data"], args["limit"], args["desc"],\
                                                     args["order_field"], args["min_species_count"]
    min_reads = args["set_min_readNumber"]
    phylogeny = args["phylogeny"] # will be False by default. If --phylogeny is included => True
    ignore_slope = args["ignore_rc"] # will be False by default. If --ignore_rc is include => True
    no_dup_proj = args["no_duplicate_proj"]
    json = json_name_converter(args["json"], list(metadata.keys()))
    threshold = args["rarefaction_threshold"]
    print(limit)
    print(metadata)
    print(min_species_count)
    limits = limit_config(limit, metadata)

    # generate all curls
    all_curls = generate_all_curls(metadata,limits,desc,spd,ordered_by)

    # run all curls and save the results to the desired json files
    all_files = run_curls(all_curls, json)
    print(limits)
    # create file with metagenmic information and returna list with all metagenomic ids
    metagenomic_ids = create_metadata(all_files, output, threshold, limits, min_species_count, metadata, phylogeny,
                                      ignore_slope, min_reads, no_dup_proj)



    stop = time.time()
    print(f"Overall Time: {round(stop-start,2)}s")

if __name__ == '__main__':
    main()
