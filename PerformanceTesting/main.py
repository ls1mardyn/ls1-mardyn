#!/usr/bin/python3

"""@package PerformanceTesting

This code aims to automate and standardize performance testing for ls1-mardyn.
The implementation should be able to be called by Jenkins whenever new commits are incoming.
To simplify setup and make the system more portable, it relies on tinyDB for storage. This is a filebased JSON database
system which does not require a server to run in the background, like MongoDB etc.

"""

import ujson
from tinydb import TinyDB, Query
import time

""" tinyDB is super lightweight both in installation and during runtime. Not optimal for performance or multi-process
accesses, but these limitations make it very portable and the database file is still human-readable.
"""

""" Available Dimensions """
dMPI = [True, False]
dOpenMP = [True, False]
dVec = ["AVX", "AVX2", "SSE"]
dRMM = [True, False]
dSize = [50]

class SingleTest:
    """SingleTest class

    Associated with a commit and all important test parameters. Builds and runs the configuration and saves the result.
    Params: - Commit
            - Vectorisation
            - MPI
            - OpenMP
            - RMM
            - Domain size
    """

    def __init__(self, commit="abcdef", mpi=False, openMP=False, vec="AVX2", RMM=True, size=50):
        """Constructor

        Init the Test instance with all necessary params.
        """
        self.commit = commit
        self.mpi = mpi
        self.openMP = openMP
        self.vec = vec
        self.RMM = RMM
        self.size = size

    def test(self):
        """Test

        Build the config. Run on Cluster or local machine. Receive Performance Measure.
        Saving the data is not part of running the test, so running multiple test concurrently somewhere is ok for the
        database file.
        """
        self.MMUPS = 10.50

    def save(self, db):
        """Saving

        Save the result to the database. This MUST only be called sequentially for all tests, as behaviours otherwise is
        not defined. This is due to the file based / non-server appraoch of tinyDB.
        """
        string = {"commit": self.commit,
                  "mpi": self.mpi,
                  "openMP": self.openMP,
                  "vec": self.vec,
                  "RMM": self.RMM,
                  "size": self.size,
                  "MMUPS": self.MMUPS}
        print("Saving: ", string)
        q = Query()

        # Upsert checks if there already exists an entry with the config, if so it gets updated with the MMUPS
        # If no entry exists, it just inserts.
        # TODO: Define behaviour for full run after partial run. Redo every entry or keep already given test results?
        db.upsert(string, (q.commit == self.commit)
                            & (q.mpi == self.mpi)
                            & (q.openMP == self.openMP)
                            & (q.vec == self.vec)
                            & (q.RMM == self.RMM)
                            & (q.size == self.size))

        # TODO: OR just to insert and deal with doubles later. MUCH faster than upsert, as no query needed
        #db.insert(string)


class Commit:
    """Commit

    Running a full or partial test suite for a given commit
    """
    def __init__(self, commit):
        self.commit = commit
        print("New commit:", commit)

    def singleDimension(self, db, config):
        """Run tests for a commit only on a single config"""
        print("Single config:", self.commit, config)
        t = SingleTest(commit=commit, mpi=config["mpi"], openMP=config["openMP"], vec=config["vec"], RMM=config["RMM"], size=config["size"])
        t.test()
        t.save(db)

    def partialRun(self):
        """Partial run for test configurations. Test Dimensions..."""
        print("partial run")

    def fullRun(self, db):
        """Run through the entire configuration space for a given commit."""
        print("Full Run ", self.commit)
        for mpi in dMPI:
            for openMP in dOpenMP:
                for vec in dVec:
                    for RMM in dRMM:
                        for size in dSize:
                            t = SingleTest(commit=commit, mpi=mpi, openMP=openMP, vec=vec, RMM=RMM, size=size)
                            t.test()
                            t.save(db)


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
    avx2 = db.search((q.vec == "AVX2") & (q.mpi == True))
    for result in avx2:
        print(result["commit"], result["MMUPS"])

    end = time.time()
    print("DURATION", end-start)