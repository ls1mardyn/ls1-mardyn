#!/usr/bin/python3

"""@package PerformanceTesting

This code aims to automate and standardize performance testing for ls1-mardyn.
The implementation should be able to be called by Jenkins whenever new commits are incoming.
To simplify setup and make the system more portable, it relies on tinyDB for storage. This is a filebased JSON database
system which does not require a server to run in the background, like MongoDB etc.

"""

from Commit import Commit
import ujson
from tinydb import TinyDB, Query
import time


""" tinyDB is super lightweight both in installation and during runtime. Not optimal for performance or multi-process
accesses, but these limitations make it very portable and the database file is still human-readable. It doesnt need a
server to run in the background like most other database setups.
"""

""" Available Dimensions """
dMPI = [True, False]
dOpenMP = [True, False]
dVec = ["AVX", "AVX2", "SSE"]
dRMM = [True, False]
dSize = [50]


if __name__ == "__main__":

    start = time.time()

    print("Starting main tests")

    # Commits to test
    # TODO: Get from github
    commits = ["abcdef", "zxywv"]

    # Databse file
    db = TinyDB("results.json")
    # Clean the entire DB
    db.purge_tables()

    # Iterate over given commits
    for commit in commits:
        c = Commit(commit)
        c.fullRun(db)
        c.singleDimension(db, {"commit": commit,
                               "mpi": False,
                               "openMP": False,
                               "vec": "AVX",
                               "RMM": True,
                               "size": 50})

    # Query to test functions
    q = Query()
    print(db.all())
    avx2 = db.search((q.vec == "AVX2") & q.mpi)
    for result in avx2:
        print(result["commit"], result["MMUPS"])

    end = time.time()
    print("DURATION", end - start, "seconds")
