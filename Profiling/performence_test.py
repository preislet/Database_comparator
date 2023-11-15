import sys
import cProfile

sys.path.insert(1,"../Database_comparator")
from Database_comparator import DB_comparator

def reload():
    return DB_comparator(config_file="DEFAULT_config_file.txt")


def profile(function, output_file):
    with open(output_file, 'w') as f:
        sys.stdout = f

        # Use cProfile.run() to profile the function
        cProfile.run(function, sort='cumulative')

if __name__ == "__main__":
    """
    db = reload()
    profile("db.exact_match.exact_match_search_in_all_databases()", "exact_match_search_in_all_databases.txt")
    db = reload()
    profile("db.aligner.smithWatermanAlgorithm_match_search_in_all_databases()", "smithWatermanAlgorithm_match_search_in_all_databases.txt")
    db = reload()
    profile("db.hamming_distances.find_hamming_distances_for_all_databases()", "hamming_distances_for_all_databases.txt")
    print("Done")
    """
    db = reload()
    db.exact_match.exact_match_search_in_all_databases()
    db.exact_match.exact_match_search_in_all_databases()
    db.exact_match.exact_match_search_in_all_databases()
    db.export_data_frame()

    

    