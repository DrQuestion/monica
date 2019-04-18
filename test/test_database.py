import os

from context import database


def test_database_builder(oldies, database_name=None, keep_genomes=None, oldies_path=None):
    if database_name:
        db = database.builder(oldies, database_name=database_name, keep_genomes=keep_genomes, oldies_path=oldies_path)
    else:
        db = database.builder(oldies, keep_genomes=keep_genomes, oldies_path=oldies_path)
    return db
