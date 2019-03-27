from tinydb import Query
import os
from subprocess import run, PIPE
import re

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

    def __init__(self, path, commit, mpi=False, openMP=False, vec="AVX2", RMM=True, size=50):
        """Constructor

        Init the Test instance with all necessary params.
        """
        self.path = path
        self.commit = commit
        self.mpi = mpi
        self.openMP = openMP
        self.vec = vec
        self.RMM = RMM
        self.size = size
        self.MMUPS = -1

    def test(self):
        """Test

        Build the config. Run on Cluster or local machine. Receive Performance Measure.
        Saving the data is not part of running the test, so running multiple test concurrently somewhere is ok for the
        database file.
        """
        old_dir = os.getcwd()
        source_dir = self.path.strip(".git") + "src"
        # Run make in source directory
        os.chdir(source_dir)
        print(os.getcwd())
        #output = run(["make", "-f","../makefile/Makefile", "-j", "2", "-B"])
        output = run(["make", "-f","../makefile/Makefile", "-j", "2"])
        if output.returncode == 0:
            # Run tests
            print("success")
            # TODO: run real MarDyn with a proper test scenario
            os.chdir("..")
            test = run(["./src/MarDyn", "examples/general-plugins/example-plugin.xml"], stdout=PIPE, stderr=PIPE)
            if test.returncode == 0:
                # Parse results
                # TODO: Scientific notation
                try:
                    performance = re.search("([0-9]*\.[0-9]*e\+[0-9]*) Molecule-updates per second", str(test.stdout)).group(1)
                except:
                    performance = -1
                self.MMUPS = performance
            else:
                print("error", test.stderr)
                # dont break on error, but set MMUPS to negative to indicate failure
                self.MMUPS = -1

        # change back
        os.chdir(old_dir)

        # TODO: remove exit
        #exit(1)

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
                  &  (q.openMP == self.openMP)
                  & (q.vec == self.vec)
                  & (q.RMM == self.RMM)
                  & (q.size == self.size))

        # TODO: OR just to insert and deal with doubles later. MUCH faster than upsert, as no query needed
        # db.insert(string)
