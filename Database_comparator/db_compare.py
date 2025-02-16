import db_exact_match as db_exact_match
import db_aligner as db_aligner
import db_blast as db_blast
import db_hamming as db_hamming
import config_class as config_class
import db_fast_hamming as db_fast

import warnings
from Bio import BiopythonDeprecationWarning

warnings.simplefilter(action='ignore', category=BiopythonDeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)


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
        self.aligner = db_aligner.aligner(self.config)   # Done
        self.blast = db_blast.blast(self.config, self.aligner)  # Done
        self.hamming_distances = db_hamming.hamming_distance(self.config)  # Done

        # Hamming distances without distance matrix - Experimental
        self.fast_hamming_distances = db_fast.FastHammingDistance(self.config)  # Done
        # Place for new modules...

    def __str__(self) -> str:
        return str(self.config)

    def export_data_frame(self, output_file="Results.xlsx", data_format="xlsx", control = True):
        """
        Export the data frame to a file in the specified format.

        Args:
            output_file (str): Name of the target file for exporting the data frame.
            data_format (str): The data format to which you want to convert the data frame (e.g., "xlsx", "csv", "md").
            control (bool): Flag to control the data format and handle long cells.

        Note:
            This method allows for exporting the data frame to a file in various formats (Excel, CSV, Markdown).
            It can also handle cases where the data frame contains cells with excessive string lengths.
        """
        excel_max_cell_string_len = 32767 - 17

        if control:
            is_longer = self.config.input_df.applymap(lambda x: len(str(x)) > excel_max_cell_string_len)
            if is_longer.any().any() and data_format == "xlsx":
                print("The dataframe has a cell that cannot be saved to an .xlsx file. The dataframe will be also exported as >Backup_save.csv<")
                print(f"Max len of cell: {excel_max_cell_string_len}")
                self.config.input_df.to_csv("Backup_save.csv", index=False)

        if data_format == "xlsx":
            self.config.input_df.to_excel(output_file, index=False)
        elif data_format == "csv":
            self.config.input_df.to_csv(output_file, index=False)
        elif data_format == "md":
            self.config.input_df.to_markdown(output_file, index=False)
        else:
            print("Unsupported data format. Please choose from: xlsx, csv, md")

        

