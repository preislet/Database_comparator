import db_exact_match as db_exact_match
import db_aligner as db_aligner
import db_blast as db_blast
import db_hamming as db_hamming
import config_class as config_class
import db_fast_hamming as db_fast

import warnings
from Bio import BiopythonDeprecationWarning, BiopythonWarning

warnings.simplefilter(action='ignore', category=BiopythonDeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=BiopythonWarning, module="Bio")

from typing import Literal


class DB_comparator:
    """
    The DB_comparator class is responsible for comparing and analyzing databases using various methods.

    It utilizes the provided configuration to perform exact matching, sequence alignment, BLAST searches,
    and calculates Hamming distances between sequences. The class allows for exporting the results to
    different file formats, such as Excel, CSV, and Markdown.
    """
    def __init__(self, config_file) -> None:
        """
        Initialize the DB_comparator class to compare databases based on the provided configuration.

        Args:
            config_file (str): Path to the configuration file.
        Note:
            This constructor initializes various database comparison components based on the
            configuration settings and provides the ability to perform exact matching, alignment,
            and BLAST-based comparisons.
        """
        self.config = config_class.cfg(config_file)
        self.exact_match = db_exact_match.ExactMatch(self.config)  # Done
        self.aligner = db_aligner.Aligner(self.config)   # Done
        self.blast = db_blast.Blast(self.config, self.aligner)  # Done
        self.hamming_distances = db_hamming.hamming_distance(self.config)  # Deprecated - use fast_hamming_distances instead
        self.fast_hamming_distances = db_fast.FastHammingDistance(self.config)  # Hamming distances without distance matrix
        # Place for new modules...
        # self.new_module = new_module.NewModule(self.config)
        # TODO: Add fuzzy matching module (e.g., Levenshtein distance)


    def __str__(self) -> str:
        return str(self.config)

    def export_data_frame(self, output_file: str="Results.xlsx", data_format: Literal["xlsx", "csv", "tsv", "md"] = "xlsx",  control: bool = True):
        """
        Export the data frame to a file in the specified format.

        Args:
            output_file (str): Name of the target file for exporting the data frame.
            data_format (str): The data format to which you want to convert the data frame (e.g., "xlsx", "csv", "tsv", "md").
            control (bool): Flag to control the data format and handle long cells.

        Note:
            This method allows for exporting the data frame to a file in various formats (Excel, CSV, Markdown).
            It can also handle cases where the data frame contains cells with excessive string lengths.
        """
        excel_max_cell_string_len: int = 32767 - 17

        if control:
            is_longer = self.config.input_df.applymap(lambda x: len(str(x)) > excel_max_cell_string_len)
            if is_longer.any().any() and data_format == "xlsx":
                print("The dataframe has a cell that cannot be saved to an .xlsx file. The dataframe will be also exported as backup_save_ExcelCellLengthError.csv")
                print(f"Max len of cell: {excel_max_cell_string_len}")
                self.config.input_df.to_csv("backup_save_ExcelCellLengthError.csv")

        if data_format == "xlsx":
            try:
                self.config.input_df.to_excel(output_file, index=False)
            except Exception as e:
                self.config.input_df.to_csv("Backup_save_EXCEPCTION_WHILE_EXPORTING.csv", index=False)
                raise(f"Exception while exporting to Excel: {e}. Backup file was created.")

        elif data_format == "csv":
            try:
                self.config.input_df.to_csv(output_file, index=False)
            except Exception as e:
                raise(f"Exception while exporting to CSV: {e}. Backup cannot be created.")
        

        elif data_format == "tsv":
            try:
                self.config.input_df.to_csv(output_file, sep="\t", index=False)
            except Exception as e:
                self.config.input_df.to_csv("Backup_save_EXCEPCTION_WHILE_EXPORTING.csv", index=False)
                raise(f"Exception while exporting to TSV: {e}. Backup file was created.")
            
        elif data_format == "md":
            try:
                self.config.input_df.to_markdown(output_file, index=False)
            except Exception as e:
                self.config.input_df.to_csv("Backup_save_EXCEPCTION_WHILE_EXPORTING.csv", index=False)
                raise(f"Exception while exporting to Markdown: {e}. Backup file was created.")
        else:
            print("Unknown error while exporting the data frame. Please check the provided data format. Exporting to Backup_save_EXCEPCTION_WHILE_EXPORTING.csv")
            self.config.input_df.to_csv("Backup_save_EXCEPCTION_WHILE_EXPORTING.csv", index=False)

        

