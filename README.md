# Database_comparator
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A program for searching and analyzing databases using various algorithms.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Configuration file](#Configuration-file)
  - [Syntax of configuration file](#Syntax-of-config-file)
  - [Inserting configuration file to program](#Inserting-config-file-to-program)
- [Usage](#Usage)
  - [Exact match](#Exact-match)
  - [Aligner](#Aligner)
  - [BLAST](#BLAST)
  - [Hamming distances](#Hamming-distances)

## Overview
The program is for comparing and analyzing databases using various methods.

It utilizes the provided configuration to perform exact matching, sequence alignment, 
BLAST searches, and calculates Hamming distances between sequences. The class allows for exporting the results to 
different file formats, such as Excel, CSV, and Markdown.

Configuration of program is given by **config_file.txt**

## Installation
```shell
# Example installation command
pip install Database_comparator
```

## Configuration file
```text
# Databases
QUERY HEDIMED__230620_Hedimed_1_22_basic--table_EF_predelana.xlsx part3

DB Databases/Nakayama.csv CDR3b [Clone/SequenceID, Epitope]
DB Databases/McPAS-TCR-filtred.csv CDR3.beta.aa [PubMed.ID, Pathology, Additional.study.details]
DB Databases/vdjdb.csv cdr3 [antigen.gene, antigen.species, mhc.a, gene]
DB Databases/TCRdb_all_sequnces.csv AASeq [TCRDB_project_ID, RunId, cloneFraction]

# Smith–Waterman algorithm
SWA_tolerance 0.9
SWA_gap_score -1000
SWA_mismatch_score 0
SWA_match_score 1

# Blastp Algorithm
BLAST_e_value 0.05
BLAST_database_name clip_seq_db
BLAST_output_name blastp_output.txt

# Hamming distance
HD_max_distance 1

# Multiprocessing
number_of_processors 3
```
#### Syntax of config file:
```
# QUERY - query database 
QUERY >Name of query database< >Name of column with sequence<

# DB - Databases with the data we want to analyze
DB >Name of data database< >Name of column with sequence< >identifiers of sequence<

# SWA_tolerance - tolerance of Smith Waterman algorithm (score/max_score)
SWA_tolerance >float<

# Smith Waterman scoring
SWA_gap_score >int<
SWA_mismatch_score >int<
SWA_match_score >int<

BLAST_e_value >float<
BLAST_database_name >the name of the blast database that will be created if needed<
BLAST_output_name >name of output file<

HD_max_distance >Maximum Hamming distance(int)<

number_of_processors >number of processors for multprocessing(int)<
```
#### Inserting config file to program:
```python
from Database_comparator import db_compare

cfg_file = 'path_to_config_file.txt'
db = db_compare.DB_comparator(cfg_file)
```
## Usage
```python3
from Database_comparator import db_compare

cfg_file = 'path_to_config_file.txt'
db = db_compare.DB_comparator(cfg_file)

# Modules
db_exact_match = db.exact_match # Used fot exact match search
db_aligner = db.aligner # Used for Smith Waterman algorithm
db_blast = db.blast # Used for BLAST search
db_hamming = db.hamming_distances # Used for finding Hamming distances between sequences
```

```python3
# Exporting results
from Database_comparator import db_compare

cfg_file = 'path_to_config_file.txt'
db = db_compare.DB_comparator(cfg_file)

# Data computing....

db.export_data_frame(output_file="Results.xlsx", data_format="xlsx")
db.export_data_frame(output_file="Results.csv", data_format="csv")

```
### Exact match
```python3
from Database_comparator import db_compare

cfg_file = 'path_to_config_file.txt'
db = db_compare.DB_comparator(cfg_file)

# Program will search in single database for exact match with query database
db.exact_match.exact_match_search_in_single_database(database_index=0)
#Multiprocessing
db.exact_match.exact_match_search_in_single_database(database_index=0, parallel=True)
db.exact_match.exact_match_search_in_single_database_MULTIPROCESSING(database_index=0)
# Program will search all given databases for exact match with query database
db.exact_match.exact_match_search_in_all_databases()

# User can also use multiprocessing
db.exact_match.exact_match_search_in_all_databases(parallel=True)
# or
db.exact_match.exact_match_search_in_all_databases_MULTIPROCESSING()

```
### Aligner
```python
from Database_comparator import db_compare

cfg_file = 'path_to_config_file.txt'
db = db_compare.DB_comparator(cfg_file)

#Single core
db.aligner.smithWatermanAlgorithm_match_search_in_single_database(database_index=0)
db.aligner.smithWatermanAlgorithm_match_search_in_all_databases()

#Multiprocessing
db.aligner.smithWatermanAlgorithm_match_search_in_single_database(database_index=0, parallel=True)
db.aligner.smithWatermanAlgorithm_match_search_in_single_database_MULTIPROCESSING()
db.aligner.smithWatermanAlgorithm_match_search_in_all_databases(parallel=True)
db.aligner.smithWatermanAlgorithm_match_search_in_all_databases_MULTIPROCESSING()
```
### BLAST
```python
from Database_comparator import db_compare
cfg_file = 'path_to_config_file.txt'
db = db_compare.DB_comparator(cfg_file)
# Info about BLAST
db.blast.blast_database_info()

# Making BLAST database
db.blast.blast_make_database(name="BLAST_Database")

db.blast.blast_search_for_match_in_database() #Query is input database
db.blast.analyze_matches_in_database() #BLAST output will be analyzed with aligner

# User can also use this function.
db.blast.blast_search_and_analyze_matches_in_database()
```
### Hamming distances
```python
from Database_comparator import db_compare
cfg_file = 'path_to_config_file.txt'
db = db_compare.DB_comparator(cfg_file)

db.hamming_distances.find_hamming_distances_for_single_database(database_index=0)
db.hamming_distances.find_hamming_distances_for_all_databases()

# User can also do this
db.hamming_distances.find_hamming_distances_for_single_database(database_index=0, analyze=False)
db.hamming_distances.analyze_single_hamming_matrix(database_index=0)

db.hamming_distances.find_hamming_distances_for_all_databases(analyze=False)
db.hamming_distances.analyze_all_hamming_matrices()

# Hamming matrices are stored in >hamming_matrices_for_all_databases<
db_matrices = db.hamming_distances.hamming_matrices_for_all_databases
```

