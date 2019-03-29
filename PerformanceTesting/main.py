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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Modes import FullHistory

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

    gitPath = "/var/lib/jenkins/workspace/Testing/"

    # Databse file
    db = TinyDB(gitPath + "results.json")

    # TODO: this is for debugging only!
    # Clean the entire DB
    db.purge_tables()

    print("Getting repo history")

    # Commits to test
    # TODO: Get from github
    #commits = ["abcdef", "zxywv"]

    

    fullH = FullHistory(gitPath)

    print("Starting main tests")

    fullH.run(db)

    # Iterate over given commits
    '''for commit in commits:
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
    '''

    q = Query()
    res = db.search(q.vec == "AVX")
    Cs = []
    Ms = []
    csv = open(gitPath + "data.csv", "w+")
    csv.write("commit,MMUPS\n")
    for r in res:
    	print(r)
    	csv.write(r["commit"] + "," + str(r["MMUPS"]) + "\n")
    	Cs.append(r["commit"])
    	Ms.append(r["MMUPS"])
    csv.close()

    plt.plot(Ms)
    plt.xticks([i for i in range(0,len(Cs))], Cs, rotation="vertical")
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.savefig(gitPath + "perf.png")

    end = time.time()
    print("DURATION", end - start, "seconds")
