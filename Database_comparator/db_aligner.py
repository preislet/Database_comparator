from Database_comparator.config_class import cfg
import pandas as pd
import multiprocessing as mp
import numpy as np
from tqdm import tqdm

class aligner:
    """
    The aligner class provides methods for performing sequence alignments using the Smith-Waterman algorithm.
    It allows for single-core and multi-core parallel processing for aligning sequences from the input DataFrame
    with sequences from databases.

    Additionally, this class offers methods for performing sequence alignments, searching for matches, and analyzing
    results within the provided configuration settings.
    """
    def __init__(self, config: cfg) -> None:
        self.config = config

# ------------------------------------Smith–Waterman algorithm------------------------------------
    # PRIVATE Smith–Waterman algorithm

    def __align_all_sequences_in_input_df_with_data_df(self, data_df: pd.DataFrame, database_index: int):
        """
        Align all sequences in the input DataFrame with sequences from the specified data DataFrame for a
        specific database.

        Args:
            data_df (pd.DataFrame): Data frame containing sequences from a specific database.
            database_index (int): Index of the database for alignment.

        Note:
            This method performs sequence alignments using the Smith-Waterman algorithm for all sequences
            in the input DataFrame with sequences from the specified database.
        """
        for i in tqdm(range(self.config.input_file_info["starting_row"], len(self.config.input_df)), desc=f"Making alignments in db: {database_index}", bar_format='{l_bar}{bar:100}{r_bar}'):
            input_seq = self.config.input_df[self.config.input_file_info["sequence_column_name"]][i]
            if "*" in input_seq: continue
            for j in range(len(data_df)):
                output_seq = data_df[self.config.data_info[database_index]["sequence_column_name"]][j]
                if self.align_sequences(input_seq, output_seq):
                    self.config.insert_match_to_input_df(
                        data_df=data_df,
                        database_index=database_index,
                        input_sequence_index=i,
                        output_sequence_index=j)

    def _align_all_sequences_in_input_df_with_data_df_MULTIPROCESSING(self, mp_input_df: pd.DataFrame, data_df: pd.DataFrame, database_index: int) -> pd.DataFrame:
        """
        Align sequences from the input DataFrame with sequences from a data DataFrame for a specific database
        using multiprocessing.

        Args:
            mp_input_df (pd.DataFrame): Data frame with sequences from the input DataFrame.
            data_df (pd.DataFrame): Data frame with sequences from a specific database.
            database_index (int): Index of the database for alignment.

        Returns:
            pd.DataFrame: Data frame with matched sequences.

        Note:
            This private method aligns sequences using the Smith-Waterman algorithm in parallel for sequences
            in the input DataFrame with sequences from a specific database.
        """
        mp_input_df = mp_input_df.reset_index(drop=True)
        for i in range(len(mp_input_df)):
            for j in range(len(data_df)):
                input_seq = mp_input_df[self.config.input_file_info["sequence_column_name"]][i]
                output_seq = data_df[self.config.data_info[database_index]["sequence_column_name"]][j]
                if self.align_sequences(input_seq, output_seq):
                    self.config.insert_match_to_input_df(
                        data_df=data_df,
                        database_index=database_index,
                        input_sequence_index=i,
                        output_sequence_index=j,
                        mp_input_df=mp_input_df)

        return mp_input_df

    # PUBLIC Smith–Waterman algorithm
    def align_sequences(self, seqA: str, seqB: str) -> bool:
        """
        Align two sequences and determine whether they match.

        Args:
            seqA (str): The first sequence.
            seqB (str): The second sequence.

        Returns:
            bool: True if the sequences match, False otherwise.

        Note:
            This method performs sequence alignment using the Smith-Waterman algorithm and returns True if
            the sequences match, based on a specified tolerance.
        """
        score = self.config.aligner.score(seqA, seqB)
        if self.config.aligner.substitution_matrix is None:
            max_score = max(len(seqA), len(seqB))
        else:
            max_score = max(self.config.aligner.score(seqA, seqA), self.config.aligner.score(seqB, seqB))

        if self.config.show_alignments:
            alignments = self.config.aligner.align(seqA, seqB)
            print(f"Score: {score}")
            print(f"Max score: {max_score}")
            print(alignments[0])

        if score/max_score >= self.config.tolerance: return True
        return False

    def smithWatermanAlgorithm_match_search_in_single_database(self, database_index: int, parallel=False) -> pd.DataFrame:
        """
        Perform Smith-Waterman algorithm-based match search in a single database.

        Args:
            database_index (int): Index of the database to perform the search on.
            parallel (bool): Whether to use parallel processing (default is False).

        Note:
            This method performs a Smith-Waterman algorithm-based match search in a single database and updates
            the input DataFrame with match results.
        """
        if parallel:
            self.smithWatermanAlgorithm_match_search_in_single_database_MULTIPROCESSING(database_index=database_index)
            return
        
        data_df = self.config.load_database(database_index=database_index, engine="python")
        self.__align_all_sequences_in_input_df_with_data_df(data_df=data_df, database_index=database_index)

        self.config.fill_Nans(database_index)
        return self.config.input_df.copy(deep=True)

    def smithWatermanAlgorithm_match_search_in_all_databases(self, parallel=False) -> pd.DataFrame:
        """
        Perform Smith-Waterman algorithm-based match search in all databases.

        Args:
            parallel (bool): Whether to use parallel processing (default is False).

        Note:
            This method performs Smith-Waterman algorithm-based match searches in all databases and updates
            the input DataFrame with match results.
        """
        self.config.reset_before_analysis()
        if parallel:
            self.smithWatermanAlgorithm_match_search_in_all_databases_MULTIPROCESSING()
            return

        for i in range(len(self.config.data_info)):
            self.smithWatermanAlgorithm_match_search_in_single_database(database_index=i)

        return self.config.input_df.copy(deep=True)

    def smithWatermanAlgorithm_match_search_in_single_database_MULTIPROCESSING(self, database_index: int) -> pd.DataFrame:
        """
        Perform Smith-Waterman algorithm-based match search in a single database using multiprocessing.

        Args:
            database_index (int): Index of the database to perform the search on.

        Note:
            This method performs Smith-Waterman algorithm-based match searches in a single database using
            multiprocessing and updates the input DataFrame with match results.
        """

        multiprocessing_input_dfs = mp.Manager().list(np.array_split(self.config.input_df, self.config.number_of_processors))

        multiprocessing_result_dfs = None

        self.config.ns.df = self.config.load_database(database_index=database_index, engine="python")
        pool = mp.Pool(processes=self.config.number_of_processors)
        multiprocessing_result_dfs = pool.starmap(self._align_all_sequences_in_input_df_with_data_df_MULTIPROCESSING,
                                                    [(multiprocessing_input_dfs[i], self.config.ns.df, database_index) for i
                                                    in range(self.config.number_of_processors)])
        pool.close()
        pool.join()

        self.config.input_df = pd.concat(multiprocessing_result_dfs, ignore_index=True)
        self.config.input_df[[self.config.data_info[database_index]["results_column"]]] = \
            self.config.input_df[[self.config.data_info[database_index]["results_column"]]].fillna(value="False")

        return self.config.input_df.copy(deep=True)

    def smithWatermanAlgorithm_match_search_in_all_databases_MULTIPROCESSING(self) -> pd.DataFrame:
        """
        Perform Smith-Waterman algorithm-based match search in all databases using multiprocessing.

        Note:
        This method performs Smith-Waterman algorithm-based match searches in all databases using
        multiprocessing and updates the input DataFrame with match results.
        """
        self.config.reset_before_analysis()
        for i in range(len(self.config.data_info)):
            self.smithWatermanAlgorithm_match_search_in_single_database_MULTIPROCESSING(database_index=i)

        return self.config.input_df.copy(deep=True)

# ------------------------------------------------------------------------------------------------
