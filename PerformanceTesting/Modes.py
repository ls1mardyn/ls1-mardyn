import os
from git import Repo
from Commit import Commit


# TODO: Make optional if full set of tests or single dimension
class FullHistory:
    """FullHistory

    runs a single set of tests on all commits in the repo
    """
    def __init__(self, gitPath, branch="master"):
        # Repo Object of already cloned Repo, expects a pulled and clean repo
        self.repo = Repo(gitPath)
        # get current head to reset later
        self.initialHead = self.repo.head.commit
        # Check for proper Repo
        if not self.repo.bare:
            self.commits = list(self.repo.iter_commits())
        else:
            print("empty repo")
            exit(-1)

    def run(self, db):
        # TODO: max commits
        MAX_COMMITS = 1000
        INTERVAL = 50
        for ci in self.commits[0:MAX_COMMITS:INTERVAL]:
            c = Commit(self.repo, ci.hexsha)
            #c.fullRun(db)
            c.singleDimension(db, {"commit": ci,
                                   "mpi": False,
                                   "openMP": False,
                                   "vec": "AVX",
                                   "RMM": True,
                                   "size": 50})

        # reset to previous state
        self.repo.head.reset(self.initialHead, index=True, working_tree=True)

# TODO: Implement
class UpdateHistory:
    """UpdateHistory

    runs a full set of tests on only the newest commit and adds it to the database
    """
    def __init__(self):
        pass

    def run(self):
        pass

