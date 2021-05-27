# Metagenomic-Data-Collection

Here, we offer you different tools to download, process, and analyze metagenomic metadata. Please check below, to see what the different tools can be used for in more detail. For Metadata Collection from [MG-Rast](https://www.mg-rast.org/) it is possible to use the [GenerateMetadataFile.py](#Metagenomic-Data-Collection). For elimination of duplicates, we offer the [CSV_Check.py](#CSV-Checker). In order to get an overview over the metadata, the usage of [DataAnalysis.py]() is recommended.

## Metagenomic Data Collection

This command line program offers a pipeline of creating a metadata file based on user-defined input parameters. Also, we offer the possibility to download the desired metagenomic files. Furthermore, the user has a chance of selecting a threshold to only receive those metagenomic data sets with a proper rarefaction curve. 
All data sets will be downloaded from the official [MG-Rast Website](https://www.mg-rast.org/).

### Prerequisites

In order to run this program without problem *Python Version 3.5+* is required. This should include all necessary modules. If you encounter any problems running the programm, please contact [Mario Rauh](mailto:mario.rauh@student.uni-tuebingen.de?subject=[GitHub]%20MasterThesis-PGPT).

### Usage

From the command line, get an overview over all possible commands using the `-h` flag. In the following, all options will be explained in more detail.

#### Output

```
-o OUTPUT, --output OUTPUT
```

Name of the output file that will later contain all the metadata information. Per default, the name of the output file will be *metadata.csv*.
Note: It is **not** necessary to assign the *.csv* extension to the name.

#### Limit

```
-l LIMIT, --limit LIMIT
```
Adding this flag to the command offers the possibility for a user to self define a number of results that will be saved in the metadata file if they fullfil all other requirements (e.g. Rarefaction Curve, minimum number of reads, etc.). Per default, this value is set to 10.
*Multiple Metadata*: If multiple metadata arguments are used (see [Metadata](#Metadata) for more information), it is possible to define multiple limits. *Note*: The number of limits must **exactly** fit the number of metadata arguments (the number of *-m*), else the remaining metadata will have default limit(s).

*Examples*:
Single Metadata Argument:
```
-l 3
```
The resulting metadata file will contain a maximum of 3 dataset references.

2 Metadata Arguments and 2 limits:
```
-l 3 5
```
The resulting metadata file will contain two different metadata information. One will have a maximum of 3 dataset references and the other one a maximum of 5 dataset references.

2 Metadata Arguments and 1 limit:
```
-l 3
```
The resulting metadata file will contain two different metadata information. The first one will have a maximum of 3 dataset references and the other one a maximum of 10 dataset references, because 10 is the default limit.

**Caution**: This limit is **not** the same limit as on the API Search of MG-Rast

#### Metadata

```
-m METADATA [METADATA ...], --metadata METADATA [METADATA ...]
```

The most important flag of this program. By using it, it is possible to specify the API Search. Furthermore, this program offers an execution of multiple API Searches with just one command. Do so, by adding another *-m* flag for each API Search. The following rules need to be considered for each of the API Searches:

* Alternating, the metadata field and its keyword(s).
* Multiple keywords must be put into quotes.

The following is an example for 2 API search requests that can be executed in the same run:

```
-m project_name hmp -m country usa all "soil plant"
```

In this case, the first search will look for the metadata field *project_name* with the property *hmp*. The second search has the metadata field *country* with the property *usa* as well as a metadata field *all* with multiple properties *soil* and *plant*.

#### Public Data

```
-pd, --public_data
```

Include this flag to allow searching public data. If not included, public data will not be searched.

#### Sorting

```
--desc
```

If included, the results will be saved in *descending* order, *ascending* else.

#### Ordering

```
--order_field ORDER_FIELD
```

Determine a metadata field for the ordering. The default metadata field is *created_on*.

#### JSON

```
-j JSON [JSON ...], --json JSON [JSON ...]
```

One step includes the API Search results to be saved in *.json* files. This flag offers the possibility to keep a better overview over those files. For each *-m* flag, add a string to name the search results.
This only works, if the number of names is equal to number of API Search requests. 
If the number is not equal to each other or the JSON flag is not included in the command, the *.json* files are named *request_1.json*, *request_2.json* etc.

One example could be:

```
-j hmp soil
```

According to the above example from [Metadata](#metadata), two API Search requests were performed. Based on our input, the created *.json* files will be named *hmp.json* and *soil.json*.

#### Rarefaction Threshold

```
-r RAREFACTION_THRESHOLD, --rarefaction_threshold RAREFACTION_THRESHOLD
```

A value greater 0 that shall indicate the sequencing depth. The lower the rarefaction threshold, the deeper the sequencing, the better the results. Default value is set to 0.5.

#### Minimum Species Count

```
--min_species_count MIN_SPECIES_COUNT
```

Including this option, will allow you to set your own species count that a dataset at least has to have to be considered. Needs to be larger than 0. The higher this value is chosen, the less results can be provided. By default, this value is set to 1000.

Note: Minimum Species Count **≠** α-Diversity

#### Set Minimum Read Number

```
--set_min_readNumber SET_MIN_READNUMBER
```

Including this options, lets the you choose a minimum number of reads that a dataset needs to consist of to be included in the final metadata. By default, this value is set to 1,000,000.

#### Phylogeny

```
--phylogeny
```

If this flag is included in the command, 16S rRNA datasets will be considered in the results, as well. Because 16S rRNA offers the identification of species but does not provide the any more genomic information, this option is disabled by default.

#### Ignore Rarefaction Curve

```
--ignore_rc
```

If this option is included, the slope of the rarefaction curve won't have any influence on the resulting dataset collection. Therefore, the [--rarefaction_threshold](#Rarefaction-Threshold) can be ignored. 


## CSV Checker

For different input *csv* files, this program will check all files for duplicates. To do so, it focuses on the metagenomic IDs that can be found in all *csv* files.

**Note**: This program was designed to work on files that were generated using the [Metadata File Generator](#metagenomic-data-collection-via-mg-rast)

### Prerequisites for the Checker

In order to run this program without problem *Python Version 3.5+* is required. Additionally, Python's [pandas](https://pandas.pydata.org/) module needs to be installed. If you encounter any problems running the programm, please contact [Mario Rauh](mailto:mario.rauh@student.uni-tuebingen.de?subject=[GitHub]%20MasterThesis-PGPT).

### Usage for the Checker

#### Input of the Checker

```
-i INPUT [INPUT ...], --input INPUT [INPUT ...]
```

The input of this program is a list of file names which shall be checked for duplicates. 

Example:

```
-i metadata1.csv metadata2.csv
```

No limit is set for the number of input *csv* files.

#### Output of the Checker

```
-o OUTPUT, --output OUTPUT
```

The checker will save a *csv* file to the desired output name. It will contain all metagenomic metadata information from all the input files. Metagenomic IDs that were present in multiple files, will only be added once to the output file.
