import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import monica.genomes.fetcher as fetcher
import monica.genomes.tables as tables
import monica.genomes.database as database
import monica.genomes.aligner as aligner