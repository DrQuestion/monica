import os

from context import database

def test_database_builder():
    database.builder()

# def test_database_builder_alt():
#     os.chdir(database.GENOMES_PATH)
#     os.remove(database.DATABASE_NAME)
#     print(os.listdir('.'))
#     os.chdir(database.CWD)
#     database.builder_alt()

if __name__=='__main__':
    test_database_builder()