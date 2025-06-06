<!--- jupyterlab:open-with:Markdown Preview --->
# Database comparator

[![PyPI version](https://badge.fury.io/py/Database-comparator.svg)](https://badge.fury.io/py/Database-comparator)
[![Downloads](https://pepy.tech/badge/database-comparator)](https://pepy.tech/project/database-comparator)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The Database Comparator is a versatile tool designed for searching, analyzing, and comparing biological sequence databases. It supports various algorithms, including exact matching, sequence alignment, BLAST searches, and Hamming distance calculations, facilitating comprehensive analysis of DNA and protein sequences. The program is highly customizable, allowing users to adjust parameters to suit their specific needs. It also supports multiprocessing, enabling faster processing of large datasets. The Database Comparator is a valuable resource for researchers, bioinformaticians, and anyone working with biological sequence data.

## Table of Contents

- [Installation](#installation)
- [Docker](#docker)
- [Configuration file](#configuration-file)
  - [Configuration file .txt format](#configuration-file-txt-format)
  - [Syntax of config file](#syntax-of-config-file)
  - [Configuration file .xlsx format](#configuration-file-xlsx-format)
  - [Inserting config file to program](#inserting-config-file-to-program)

- [Loging](#loging)
- [Exporting results](#exporting-results)
- [Usage](#usage)
  - [Exact match](#exact-match)
  - [Aligner](#aligner)
  - [BLAST](#blast)
  - [Hamming distances](#hamming-distances)


## Installation
Use the following command to install the program:
```shell
pip install Database-comparator
```

or clone the repository and install the program manually:
```shell
git clone https://github.com/preislet/Database_comparator.git
cd Database_comparator
```
### Blast
#### Windows:
BLAST needs to be installed [manually](https://www.ncbi.nlm.nih.gov/books/NBK569861/).

#### Linux:
```shell
sudo apt update
sudo apt install ncbi-blast+
```

## Docker
To download prepared docker image use:

```shell
docker pull ghcr.io/preislet/database_comparator:latest
```

Optionaly docker file is provided in repository. To build the image from the Dockerfile, follow these steps:

#### Step 1:
Run the following command to build the Docker image. 
```shell
docker build -t db-comparator .
```
#### Step 2:
After the image is successfully built, you can run a container from it:

```shell
docker run -p 8888:8888 db-comparator
```

#### Step 3:
Open a web browser and go to http://localhost:8888.

This will open the Jupyter Notebook interface where you can run the Database Comparator. Example notebooks are provided in the `notebooks` folder. You can also upload your own configuration file and data files to the container.

## Configuration file

The configuration file is used to adjust the program properly to the data that the user wants to analyze. The configuration folder contains all the information from the database query and the databases against which we want to compare the query. Optionally, internal parameters for the Smith Waterman algorithm, BLAST, etc. can be set. If these parameters are not specified, they will be set to the default value. Configuration file can be in **.txt** or .**xlsx** format. **We highly recommend using .xlsx format** because it is more user-friendly.

### Configuration file .txt format

The table below describes all available configuration options for the **Database Comparator**.

| Option Name            | Description                                                       | Type    | Default Value                | Example Values                   |
|------------------------|-------------------------------------------------------------------|---------|------------------------------|----------------------------------|
| `DB`                  | Defines a database path, sequence column, and result column.     | String  | None                         | `DB path/to/db.csv seq_col result_col [identifiers]` |
| `QUERY`               | Specifies the query file path and sequence column name.          | String  | None                         | `QUERY path/to/query.csv seq_col` |
| `SWA_tolerance`       | Tolerance for Smith-Waterman alignment.                          | Float   | `0.93`                        | `0.95`, `0.9`                   |
| `SWA_gap_score`       | Gap penalty for Smith-Waterman alignment.                        | Float   | None                          | `-2.0`, `-3.0`                  |
| `SWA_mismatch_score`  | Mismatch penalty for Smith-Waterman alignment.                   | Float   | None                          | `-1.0`, `-2.0`                  |
| `SWA_match_score`     | Match reward for Smith-Waterman alignment.                       | Float   | None                          | `2.0`, `3.0`                    |
| `SWA_matrix`          | Substitution matrix for alignment.                               | String  | None                          | `BLOSUM62`, `PAM250`            |
| `SWA_mode`            | Alignment mode (`global` or `local`).                            | String  | None                          | `local`, `global`               |
| `BLAST_e_value`       | E-value threshold for BLAST searches.                            | Float   | `0.05`                        | `1e-5`, `0.01`                  |
| `BLAST_database_name` | Name of the BLAST database.                                      | String  | `"clip_seq_db"`               | `"Any_name"`           |
| `BLAST_output_name`   | Name of the BLAST output file.                                   | String  | `"blastp_output.txt"`         | `"output.txt"`, `"results.tsv"` |
| `HD_max_distance`     | Maximum allowed Hamming distance.                               | Integer | `1`                           | `2`, `5`, `10`                  |
| `number_of_processors`| Number of CPU cores to use for multiprocessing.                 | Integer | `1`                           | `2`, `4`, `8`                   |
| `separator`           | Separator for results in the input DataFrame.                   | String  | `"\n"`                        | `";"`, `","`, `" "`             |

#### Notes:

- The `DB` and `QUERY` parameters are required in the configuration file.
- Some parameters (like `SWA_tolerance`, `SWA_match_score`, etc.) are specific to **Smith-Waterman alignment**.
- The `BLAST_*` parameters configure **BLAST sequence searches**.
- `HD_max_distance` is used for **Hamming distance calculations**.
- `separator` determines how multiple results are stored in the output file.


Example of configuration file:

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

```text
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
SWA_matrix >name of scoring matrix<
SWA_mode >local | global<

BLAST_e_value >float<
BLAST_database_name >the name of the blast database that will be created if needed<
BLAST_output_name >name of output file<

HD_max_distance >Maximum Hamming distance(int)<

number_of_processors >number of processors for multprocessing(int)<
```

#### Notes:
if you want to use the **default value** for some parameter, you can skip it in the configuration file. Default values are shown in the table above. 

If the configuration file is not correct, the program will raise an error and inform the user about the problem. 

**If the configuration file is not provided, the program will switch into testing only mode.**

### Configuration file .xlsx format

.xlsx format is more user-friendly and allows for easier configuration of the program. The Default .xlsx file is provided in github repository. The user can modify it according to their needs. The .xlsx file contains several sheets, each with a different purpose. All tables are predefined, and the user only needs to fill in the necessary data. Cells the yellow color are only cells that the user can modify. **If the user wants to use the default value for some parameter, they can leave the cell empty.**

The first sheet is the **Query** sheet, where the user can specify the query database and the databases against which they want to compare the query. It also contains the **Sepataor** parameter, which determines how multiple results are stored in the output file and **Number of processors** parameter, which determines the number of CPU cores to use for multiprocessing. The **Aligner** sheet is used to set parameters for the Smith-Waterman algorithm, such as tolerance, gap score, mismatch score, match score, scoring matrix, and alignment mode. The **BLAST** sheet is used to configure BLAST searches, including the E-value threshold, database name, and output file name. The **Hamming_distance** sheet is used to set the maximum allowed Hamming distance.

#### Notes:

- The **Query** sheet is required in the configuration file.
- The **Aligner**, **BLAST**, and **Hamming_distance** sheets are optional.

### Inserting config file to program:

```python
from Database_comparator import db_compare

if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file)
```


## Loging
Program is capable of logging all the actions that are performed. The Log folder (***DatabaseComparatorLogs***) is created in the directory where the program is run. User can specify log_roject and log_tag for better organization of the logs. Logs name also contains timestamp of the log creation.

```python
from Database_comparator import db_compare
if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file, log_project="MyProject", log_tag="MyTag")

```

Loging folder then can look like this:

```
DatabaseComparatorLogs/
├── Project1
│   ├── DB_comparator_run_2023-10-01_12-00-00_MyTag1.log
│   ├── DB_comparator_run_2025-10-01_12-30-00_MyTag2.log
│   └── DB_comparator_run_2025-10-01_13-00-00_MyTag3.log
└── Project2
    ├── DB_comparator_run_2025-10-01_14-00-00_MyTag1.log
    ├── DB_comparator_run_2025-10-01_14-30-00_MyTag2.log
    └── DB_comparator_run_2020-10-01_15-00-00_MyTag3.log
```

## Exporting results

The program is also capable of exporting the results to a .csv/.xlsx file. The user can specify the path to the output file and the separator used to separate the multiple results for same sequence. The separator is defined in the **configuration file**. If the file exceeds the maximum allowed size for .xlsx files, the program will automatically generate .csv backup file.

### Example of exporting results:

```python
from Database_comparator import db_compare

if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file, log_project="MyProject", log_tag="MyTag")

  # Data computing...

  db.export_data_frame(output_file="MyAnalysis.xlsx", data_format="xlsx")
  db.export_data_frame(output_file="MyAnalysis.csv", data_format="csv")
```

The dataframe can be also accessed using the **db.config.input_df** attribute.

## Usage

### Exact match

The **exact_match** module is used to find exact matches between sequences in the query database and data databases. It allows you to perform exact match searches in single databases or across all configured databases. Users can also take advantage of multiprocessing to speed up the search process.

#### Example of exact match search in single database (first database in the configuration file):

```python
from Database_comparator import db_compare

if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file, log_project="MyProject", log_tag="MyTag")

  db.exact_match.exact_match_search_in_single_database(database_index=0, parallel=True) # Multiprocessing enabled (parallel=True)
```

#### Example of exact match search in all databases:

```python
from Database_comparator import db_compare

if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file, log_project="MyProject", log_tag="MyTag")

  db.exact_match.exact_match_search_in_all_databases(parallel=True) # Multiprocessing enabled (parallel=True)
```

### Aligner

The **aligner** module is based on the Smith-Waterman/Needleman-Wunsch algorithm for sequence alignment. It provides the capability to execute single-core or multiprocessing-based match searches. Algorithm complexity is **O(n*m)**, where n is the length of the first sequence and m is the length of the second sequence. Tolernace parameter is used to determine the minimum score that the alignment must achieve to be considered a **hit**. The gap score, mismatch score, and match score are used to calculate the alignment score. The scoring matrix is used to determine the score for each pair of aligned residues. The alignment mode can be set to either global or local.

#### Example of Smith-Waterman algorithm search in single database (first database in the configuration file):

```python
from Database_comparator import db_compare

if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file, log_project="MyProject", log_tag="MyTag")

  db.aligner.aligner_search_in_single_database(database_index=0, parallel=True) # Multiprocessing enabled (parallel=True)
```

#### Example of Smith-Waterman algorithm search in all databases:

```python
from Database_comparator import db_compare

if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file, log_project="MyProject", log_tag="MyTag")

  db.aligner.aligner_search_in_all_databases(parallel=True) # Multiprocessing enabled (parallel=True)
```

### BLAST

The **blast** module enables users to create BLAST databases, perform BLAST searches for matches. The E-value threshold is used to determine the significance of the match. The database name is used to specify the name of the BLAST database that will be created if needed. The output name is used to specify the name of the output file. 

In future versions, the user will be able to specify if they want to use aligner or hammer distance for analyzing the BLAST results.

#### Example of BLAST search in database:

```python
from Database_comparator import db_compare

if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file, log_project="MyProject", log_tag="MyTag")

  db.blast.blast_database_info() # Provides information about the BLAST database
  
  db.blast.blast_make_database() # Creates BLAST database
  db.blast.blast_search_for_match_in_database() #Query is input database
  db.blast.analyze_matches_in_database() #BLAST output will be analyzed with aligner

  """
  User can also use this function.
  db.blast.blast_search_and_analyze_matches_in_database() - This function will perform both BLAST search and write the results to the ouput dataframe.
  """
```

### Hamming distances

The **hamming_distances** module calculates Hamming distances between sequences. Users can explore Hamming distances in single databases or across all databases. The maximum allowed Hamming distance is used to determine the maximum number of mismatches allowed between two sequences. The user can also analyze the Hamming distance matrices to identify patterns in the data. The Hamming distances can be calculated using standard hamming distance function, that will return matrix with hamming distances between all sequences in the database. This matrix can be analyzed using the **analyze_single_hamming_matrix** function. The user can also calculate Hamming distances for all databases and analyze them using the **analyze_all_hamming_matrices** function.

#### Example of Hamming distance search in single database (first database in the configuration file):

```python
from Database_comparator import db_compare

if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file)
  
  # Hamming distances will be analyzed - The hits under the maximum allowed Hamming distance will be stored in the output file
  db.hamming_distances.find_hamming_distances_for_single_database(database_index=0, analyze=True) 

  # Hamming matrices are stored in >hamming_matrices_for_all_databases<
  db_matrices = db.hamming_distances.hamming_matrices_for_all_databases
```

This aproach is very space consuming, so the user can also calculate Hamming distances **without generating the matrix**. This aproach is much faster and uses less memory. **Sequnces that are under the maximum allowed Hamming distance will be stored in the output file.** No further analysis of the matrix is possible, because was never generated.

#### Example of fast Hamming distance search in single database (first database in the configuration file):

```python
if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file,log_project="MyProject", log_tag="MyTag")

  db.fast_hamming_distances.find_hamming_distances_for_single_database(0, parallel=True) # Multiprocessing enabled (parallel=True)
```

#### Example of fast Hamming distance search in all databases:

```python
if __name__ == "__main__":
  cfg_file = 'path_to_config_file.txt'
  db = db_compare.DB_comparator(cfg_file, log_project="MyProject", log_tag="MyTag")

  db.fast_hamming_distances.find_hamming_distances_for_all_databases(parallel=True) # Multiprocessing enabled (parallel=True)
```

**We highly recommend using the fast Hamming distance search.**
