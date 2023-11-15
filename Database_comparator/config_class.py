import pandas as pd
import pyreadr as pr
import os
import warnings
from pathlib import Path
import multiprocessing as mp
import numpy as np
from Bio import Align

warnings.simplefilter(action='ignore', category=FutureWarning)

# TODO udělat kontrolu formátu pro databaze - zatím je podporován jen .csv formát => rozšířit podporované formáty
# TODO updatovat README o nové parametry v config file (swa_mode, matice)

StopCodon = "#"
class cfg:
    """
    Initialize the configuration class for a bioinformatics sequence analysis program.

    Args:
        config_file (str): Path to the configuration file.

    Note:
        This constructor initializes various parameters and loads settings from the
        specified configuration file to customize the behavior of the program.
    """
    def __init__(self, config_file=None) -> None:
        """
        Initialize the configuration class for a bioinformatics sequence analysis program.

        Args:
            config_file (str): Path to the configuration file.

        Note:
            This constructor initializes various parameters and loads settings from the
            specified configuration file to customize the behavior of the program.
        """
        # Paths
        self.input_file_path = None
        self.blast_database_path = "BLAST_database"
        self.blast_query_files_path = "Query_files"

        # Input data_info 
        self.input_file_info = None
        
        # Databases info list
        self.data_info = []

        # Report switches
        self.show_report_while_inserting_match_to_input_df = False
        self.show_alignments = False

        # Smith–Waterman algorithm
        self.aligner = Align.PairwiseAligner()
        self.tolerance = 0.93

        # Blastp Algorithm
        self.e_value = 0.05
        self.blast_database_name = "clip_seq_db"
        self.blast_database_full_name = self.blast_database_path + "//" + self.blast_database_name
        self.blast_output_name = "blastp_output.txt"
        self.blast_default_input_query = self.blast_query_files_path + "//QUERY.fasta"

        self.in_blast_database = self.data_info
        self.blast_outfmt = "6 qseqid sseqid qseq sseq bitscore score"

        # Hamming distance
        self.max_hamming_distance = 1

        # Multiprocessing
        self.number_of_processors = 1
        self.ns = None
        try:
            self.ns = mp.Manager().Namespace()
        except:
            raise Exception("Multiprocessing is not supported on this system. Exception was raised do to >>>mp.Manager().Namespace()<<<. Comparator must be used under if __name__ == '__main__':")

        self.__load_config(config_file)

        # dataframe of input file
        self.repair_input_df = True
        self.input_df = None
        self.__load_input_df()

    def __str__(self) -> str:
        temp = vars(self)
        for item in temp:
            print(item, ' : ', temp[item])

        return ""
    
    
    def __load_config(self, config_file):
        """
        Load configuration settings from the specified file.

        Args:
            config_file (str): Path to the configuration file.

        Note:
            This method reads and interprets settings from the provided configuration file
            and populates the class properties accordingly.
        """
        if config_file is None:
            print("A configuration file was not provided. Please provide the configuration file and restart the program")
            print("See the documentation for more information: https://pypi.org/project/Database-comparator")
            raise Exception("Configuration file was not provided")
        
        file = open(config_file, 'r')
        for line in file:
            line = line.split()
            if len(line) == 0:
                continue
            

            if line[0].upper() == "DB":
                try:
                    data = {
                        "path": line[1],
                        "sequence_column_name": line[2],
                        "results_column": line[1][:-4],
                        "identifier_of_seq": "".join(line[3:]).strip('][').split(',')
                    }
                    self.data_info.append(data)
                except: raise Exception(f"line: {line}... Database path, sequence column name or results column name is missing. Please check your config file.")

            elif line[0].upper() == "QUERY":
                try:
                    self.input_file_path = line[1]
                    self.input_file_info = {
                        "path": self.input_file_path,
                        "sequence_column_name": line[2],
                        "starting_row": 0
                    }
                except:
                    raise Exception(f"line: {line}... Input file path or sequence column name is missing. Please check your config file.")
                
            elif line[0].upper() == "SWA_tolerance".upper(): 
                try: self.tolerance = float(line[1])
                except: raise Exception(f"line: {line}... Tolerance must be float")

            elif line[0].upper() == "SWA_gap_score".upper():
                try: 
                    self.aligner.open_gap_score = float(line[1])
                    self.aligner.extend_gap_score = float(line[1])

                except: raise Exception(f"line: {line}... Gap score must be float")
            elif line[0].upper() == "SWA_mismatch_score".upper():
                try: self.aligner.mismatch_score = float(line[1])
                except: raise Exception(f"line: {line}... Mismatch score must be float")

            elif line[0].upper() == "SWA_match_score".upper():
                try: self.aligner.match_score = float(line[1])
                except: raise Exception(f"line: {line}... Match score must be float")

            elif line[0].upper() == "BLAST_e_value".upper():
                try: self.e_value = float(line[1])
                except: raise Exception(f"line: {line}... E-value must be float")

            elif line[0].upper() == "BLAST_database_name".upper():
                try: self.blast_database_name = line[1]
                except: raise Exception(f"line: {line}... BLAST database name must be string")

            elif line[0].upper() == "BLAST_output_name".upper():
                try: self.blast_output_name = line[1]
                except: raise Exception(f"line: {line}... BLAST output name must be string")

            elif line[0].upper() == "HD_max_distance".upper():
                try: self.max_hamming_distance = int(line[1])
                except: raise Exception(f"line: {line}... Max hamming distance must be integer")

            elif line[0].upper() == "number_of_processors".upper():
                try: self.number_of_processors = int(line[1])
                except: raise Exception(f"line: {line}... Number of processors must be integer")

            elif line[0].upper() == "SWA_matrix".upper():
                if not line[1] in Align.substitution_matrices.load():
                    err = f"Substitution matrix not found. Substitution matrices: {Align.substitution_matrices.load()}"
                    raise Exception(err) 
                self.aligner.substitution_matrix = Align.substitution_matrices.load(line[1])
            elif line[0].upper == "SWA_mode".upper():
                if not line[1].lower in ["local", "global"]:
                    err = "Mode not found. Please use only global/local"
                    raise Exception(err)
                
                self.aligner.mode = line[1].lower()

            elif line[0].upper() == "#".upper():
                pass
            else:
                print(f"Command not recognized: {line[0]}. Check your config file for possible typos.")
                raise Exception(f"line: {line}... Error in config file.")

        file.close()

    def __load_input_df(self):
        """
        Load input data from a file and prepare it for processing.

        Note:
            This method loads data from the input file, performs data cleaning and
            preprocessing tasks, and stores the resulting DataFrame in the class.
        """
        supported_formats = [".csv", ".tsv" ".xlsx", ".xls", ".RData", ".Rbin", ".RDATA"]
        path = self.input_file_info["path"]
        if Path(path).suffix == ".csv":
            self.input_df = pd.DataFrame(pd.read_csv(self.input_file_info["path"]))
        elif Path(path).suffix in [".xlsx", ".xls"]:
            self.input_df = pd.DataFrame(pd.read_excel(self.input_file_info["path"]))
        elif Path(path).suffix in [".RData", ".Rbin", ".RDATA"]:
            data = pr.read_r(path)
            self.input_df = data[os.path.splitext(path)[0]]
        elif Path(path).suffix == ".tsv":
            self.input_df =  pd.DataFrame(pd.read_csv(self.input_file_info["path"], sep="\t"))
        else:
            print(f"File format is not supported. Supported formats: {supported_formats}")
        
        if self.repair_input_df:
            self.__repair_input_df()

    # Resets the class to default values before running a new analysis
    # TODO Add reset_before_analysis infront of every analysis - Decorator propably
    def reset_before_analysis(self):
        if self.__check_if_input_df_changed():
            print("#"*200)
            print("Resetting the class to default values before running a new analysis")
            self.input_df = None
            self.__load_input_df()
            print("Reset was successfuly done")
            print("Analysing the data...")
            print("#"*200)

    def __check_if_input_df_changed(self) -> bool:
        """
        Check if the input DataFrame has been changed.

        Returns:
            bool: True if the input DataFrame has been changed, False otherwise.

        Note:
            This method is used to check if the input DataFrame has been changed since
            the last analysis.
        """

        # TODO: check if the input_df has been 
        return not all([pd.isnull(self.input_df[self.data_info[i]["results_column"]]).all() for i in range(len(self.data_info))])
        # Not sure if this is the best way to check if the input_df has been changed

    def __repair_input_df(self):
        """
        Perform data cleaning and preprocessing on the input DataFrame.

        Note:
            This method cleans the input data, removing unwanted characters and filtering
            out sequences with stop codons, and ensures that the DataFrame is in a suitable
            format for analysis.
        """
        # deleting string (pre-)filtered from seq column
        self.input_df[self.input_file_info["sequence_column_name"]] = self.input_df[self.input_file_info["sequence_column_name"]].str.replace("(pre-)filtered", "")
        self.input_df[self.input_file_info["sequence_column_name"]] = self.input_df[self.input_file_info["sequence_column_name"]].str.replace(" ", "")
        # deleting sequences with stop codon
        self.input_df = self.input_df[~self.input_df[self.input_file_info["sequence_column_name"]].astype(str).str.contains(StopCodon)]
        self.input_df = self.input_df.reset_index(drop=True)
        self.input_df[self.input_file_info["sequence_column_name"]].fillna(value="********", inplace=True)

        for db in self.data_info:
            self.input_df[db["results_column"]] = np.nan
        
    def fill_Nans(self, database_index: int):
        """
        Fill missing values in the input DataFrame with "False" for a specified database.

        Args:
            database_index (int): Index of the database to fill missing values for.

        Note:
            This method is used to fill missing values in the DataFrame when analyzing
            data from a specific database.
        """
        result_column_name = self.data_info[database_index]["results_column"]
        for i in range(len(self.input_df)):
            if pd.isnull(self.input_df.loc[i, result_column_name]):
                self.input_df.loc[i, result_column_name] = "False"

    @staticmethod
    def merge_all_identifiers(data_df: pd.DataFrame, identifier_column_names: list, output_sequence_index: int) -> str:
        """
        Merge all identifiers for a sequence in the output data.

        Args:
            data_df (pd.DataFrame): The data DataFrame.
            identifier_column_names (list): List of column names containing identifiers.
            output_sequence_index (int): Index of the sequence in the output data.

        Returns:
            str: A merged string of identifiers.

        Note:
            This static method is used to merge multiple identifiers associated with a
            sequence in the output data.
        """
        full_identifier = ""
        for identifier_column_name in identifier_column_names:
            identifier = data_df[identifier_column_name][output_sequence_index]
            full_identifier += f" ({identifier_column_name}: {identifier})"

        return full_identifier
    
    def insert_match_to_input_df(self, data_df: pd.DataFrame, database_index: int, input_sequence_index: int, output_sequence_index: int, mp_input_df=None):
        """
        Insert a match from the output data into the input DataFrame.

        Args:
            data_df (pd.DataFrame): The output data DataFrame.
            database_index (int): Index of the database being analyzed.
            input_sequence_index (int): Index of the sequence in the input DataFrame.
            output_sequence_index (int): Index of the sequence in the output data.
            mp_input_df (Optional[mp_input_df]): A multiprocessing input DataFrame (if applicable).

        Note:
            This method is used to insert a matching result from the output data into the
            input DataFrame for further analysis and reporting.
        """
        if mp_input_df is not None:
            input_df = mp_input_df

        else: input_df = self.input_df

        filename, file_extension = os.path.splitext(self.data_info[database_index]["path"])
        identifier_column_names = self.data_info[database_index]["identifier_of_seq"]
        sseq = data_df[self.data_info[database_index]["sequence_column_name"]][output_sequence_index]

        if identifier_column_names is None:
            identifier = os.path.basename(filename)
            s = "file_name"

        else: identifier = self.merge_all_identifiers(data_df=data_df, identifier_column_names=identifier_column_names,output_sequence_index=output_sequence_index)

        identifier = ";".join(set(identifier.split(sep=";")))
        if pd.isnull(input_df.loc[input_sequence_index, self.data_info[database_index]["results_column"]]):
            input_df.loc[input_sequence_index, self.data_info[database_index]["results_column"]] = f"[seq: {sseq}{identifier}]\n"
        else:
            input_df.loc[input_sequence_index, self.data_info[database_index]["results_column"]] = str(input_df[self.data_info[database_index]["results_column"]][input_sequence_index]) \
                                                                                                        + f"[seq: {sseq}{identifier}]\n"

        if self.show_report_while_inserting_match_to_input_df:
            input_sequence = input_df[self.input_file_info["sequence_column_name"]][input_sequence_index]
            output_sequence = data_df[self.data_info[database_index]["sequence_column_name"]][output_sequence_index]
            print(f"Match found in {os.path.basename(filename)}", flush=True)
            print(f"inp sequence: {input_sequence}", flush=True)
            print(f"out sequence: {output_sequence}", flush=True)
            print(f"Identifier: {identifier}", flush=True)
            print("-"*200)
            
    def load_datafiles_names_from_stored_path(self, database_index: int) -> list:
        """
        Load the names of data files from the specified database path.

        Args:
            database_index (int): Index of the database being analyzed.

        Returns:
            list: A list of file names from the specified database path.

        Note:
            This method is used to retrieve file names from the database path for further
            analysis.
        """
        front = []
        path_to_data = self.data_info[database_index]["path"]
        if os.path.isdir(path_to_data):
            for file in os.listdir(path_to_data):
                front.append(os.path.join(path_to_data, file))
        else:
            front.append(path_to_data)
        return front

    def find_database_index(self, filename: str) -> int:
        """
        Find the index of a database based on a filename.

        Args:
            filename (str): Name of the file to search for.

        Returns:
            int: Index of the database (if found), or -1 if not found.

        Note:
            This method is used to locate a database by searching for a specific filename.
        """
        for i in range(len(self.data_info)):
            if filename in self.data_info[i]["path"]:
                return i

    @staticmethod
    def load_database(path: str, engine=None) -> pd.DataFrame:
        suffix = Path(path).suffix
        supported_formats = [".csv", ".tsv" ".xlsx", ".xls", ".RData", ".Rbin", ".RDATA"]

        if suffix == ".csv":
            return pd.DataFrame(pd.read_csv(path, engine=engine))
        elif suffix in [".xlsx", ".xls"]:
            return pd.DataFrame(pd.read_excel(path, engine=engine))
        elif suffix in [".RData", ".Rbin", ".RDATA"]:
            data = pr.read_r(path)
            return data[os.path.splitext(path)[0]]
        elif suffix == ".tsv":
            return pd.DataFrame(pd.read_csv(path, sep="\t", engine=engine))
        else:
            err = f"File format is not supported for database. Supported formats: {supported_formats}"
            raise Exception(err)
