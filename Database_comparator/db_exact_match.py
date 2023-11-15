from Database_comparator.config_class import cfg
import pandas as pd
import numpy as np
import multiprocessing as mp
from tqdm import tqdm

class exact_match:
    """
    The exact_match class provides methods for performing exact match searches in databases.

    It allows for exact match searches in a single database or across all databases, using both single-core and
    multiprocessing techniques. The class handles the insertion of matching results into the input DataFrame
    based on the provided configuration.
    """
    def __init__(self, config: cfg) -> None:
        """
        Initialize the exact_match class with a configuration.

        Args:
            config (cfg): An instance of the configuration class.

        Note:
            This constructor sets up the exact_match class with the provided configuration.
        """
        self.config = config

    # ------------------------------------Exact search single core------------------------------------------------
    # PRIVATE Exact search

    def __search_for_single_sequence_in_input_df(self, data_df: pd.DataFrame, database_index: int,sequence_index: int):
        """
        Search for a single sequence in the input DataFrame.

        Args:
            data_df (pd.DataFrame): The data DataFrame for a specific database.
            database_index (int): Index of the database being searched.
            sequence_index (int): Index of the sequence in the input DataFrame.

        Note:
            This private method performs an exact match search for a single sequence in the input DataFrame
            and inserts matching results into the input DataFrame using the configuration settings.
        """
        output_indexes = data_df.index[data_df[self.config.data_info[database_index]["sequence_column_name"]] == self.config.input_df[self.config.input_file_info["sequence_column_name"]][sequence_index]]
        for j in output_indexes:
            self.config.insert_match_to_input_df(
                data_df=data_df,
                database_index=database_index,
                input_sequence_index=sequence_index,
                output_sequence_index=j)

    def __search_for_all_sequences_in_input_df(self, data_df: pd.DataFrame, database_index: int):
        """
        Search for all sequences in the input DataFrame for a specific database.

        Args:
            data_df (pd.DataFrame): The data DataFrame for a specific database.
            database_index (int): Index of the database being searched.

        Note:
            This private method iterates over all sequences in the input DataFrame and searches for matches
            in the specified database. It inserts matching results into the input DataFrame using the configuration settings.
        """
        for i in tqdm(range(self.config.input_file_info["starting_row"], len(self.config.input_df)), desc=f"Searching for exact match in database with index: {database_index}"):
            self.__search_for_single_sequence_in_input_df(data_df, database_index, i)

    # PUBLIC Exact search
    def exact_match_search_in_single_database(self, database_index: int, parallel=False) -> pd.DataFrame:
        """
        Perform an exact match search in a single database.

        Args:
            database_index (int): Index of the database to perform the search on.
            parallel (bool): Whether to use parallel processing (default is False).

        Returns:
            pd.DataFrame: The input DataFrame with matching results.

        Note:
            This method performs an exact match search in a single database and can optionally use parallel processing.
            Matching results are inserted into the input DataFrame using the configuration settings.
        """
        if parallel:
            self.exact_match_search_in_single_database_MULTIPROCESSING(database_index=database_index)
            return

        front = self.config.load_datafiles_names_from_stored_path(database_index=database_index)

        for data_file in front:
            data_df = self.config.load_database(path = data_file, engine="python")
            self.__search_for_all_sequences_in_input_df(data_df, database_index)

        self.config.fill_Nans(database_index)
        return self.config.input_df.copy(deep=True)

    def exact_match_search_in_all_databases(self, parallel=False) -> pd.DataFrame:
        """
        Perform an exact match search in all databases.

        Args:
            parallel (bool): Whether to use parallel processing (default is False).

        Note:
            This method performs an exact match search in all databases and can optionally use parallel processing.
        """
        self.config.reset_before_analysis()
        if parallel:
            self.exact_match_search_in_all_databases_MULTIPROCESSING()
            return
        for i in range(len(self.config.data_info)):
            self.exact_match_search_in_single_database(database_index=i)

        return self.config.input_df.copy(deep=True)


    # ------------------------------------Exact search multiprocessing------------------------------------------------

    def _search_for_all_sequences_in_input_df_MULTIPROCESSING(self, mp_input_df: pd.DataFrame,data_df: pd.DataFrame, database_index: int):
        """
        Search for all sequences in the input DataFrame for a specific database using multiprocessing.

        Args:
            mp_input_df (pd.DataFrame): Input DataFrame for a specific processor.
            data_df (pd.DataFrame): The data DataFrame for a specific database.
            database_index (int): Index of the database being searched.

        Returns:
            pd.DataFrame: The input DataFrame with matching results.

        Note:
            This private method searches for all sequences in the input DataFrame using multiprocessing.
            Matching results are inserted into the input DataFrame using the configuration settings.
        """
        mp_input_df.reset_index(drop=True, inplace=True)
        for i in range(len(mp_input_df)):
            output_indexes = data_df.index[data_df[self.config.data_info[database_index]["sequence_column_name"]] == mp_input_df[self.config.input_file_info["sequence_column_name"]][i]]
            for j in output_indexes:
                self.config.insert_match_to_input_df(
                    data_df=data_df,
                    database_index=database_index,
                    input_sequence_index=i,
                    output_sequence_index=j,
                    mp_input_df=mp_input_df)

        return mp_input_df

    def exact_match_search_in_single_database_MULTIPROCESSING(self, database_index: int) -> pd.DataFrame:
        """
        Perform an exact match search in a single database using multiprocessing.

        Args:
            database_index (int): Index of the database to perform the search on.

        Returns:
            pd.DataFrame: The input DataFrame with matching results.

        Note:
            This method performs an exact match search in a single database using multiprocessing.
            Matching results are inserted into the input DataFrame using the configuration settings.
        """
        front = self.config.load_datafiles_names_from_stored_path(database_index=database_index)
        multiprocessing_input_dfs = mp.Manager().list(np.array_split(self.config.input_df, self.config.number_of_processors))

        multiprocessing_result_dfs = None
        for data_file in front:
            self.config.ns.df = self.config.load_database(path=data_file, engine="python")
            pool = mp.Pool(processes=self.config.number_of_processors)
            multiprocessing_result_dfs = pool.starmap(self._search_for_all_sequences_in_input_df_MULTIPROCESSING,
                                                      [(multiprocessing_input_dfs[i], self.config.ns.df, database_index)
                                                       for i in range(self.config.number_of_processors)])
            pool.close()
            pool.join()

        self.config.input_df = pd.concat(multiprocessing_result_dfs, ignore_index=True)
        self.config.input_df[[self.config.data_info[database_index]["results_column"]]] = \
            self.config.input_df[[self.config.data_info[database_index]["results_column"]]].fillna(value="False")

        return self.config.input_df.copy(deep=True)

    def exact_match_search_in_all_databases_MULTIPROCESSING(self) -> pd.DataFrame:
        """
        Perform an exact match search in all databases using multiprocessing.

        Returns:
            pd.DataFrame: The input DataFrame with matching results.

        Note:
            This method performs an exact match search in all databases using multiprocessing.
        """
        self.config.reset_before_analysis()
        for i in range(len(self.config.data_info)):
            self.exact_match_search_in_single_database_MULTIPROCESSING(database_index=i)

        return self.config.copy(deep=True)
