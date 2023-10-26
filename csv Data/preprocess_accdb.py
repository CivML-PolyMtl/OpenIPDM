from accdbtools import accdb2csv
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: preprocess_mtq.py <accdb_filename>")
        sys.exit(1)
    
    accdb_filename = sys.argv[1]
    accdb2csv(accdb_filename)