import os

from context import database


def test_database_builder(oldies, database_name=None):
    if database_name:
        db = database.builder(oldies, database_name=database_name)
    else:
        db=database.builder(oldies)
    return db


if __name__=='__main__':
    test_database_builder()