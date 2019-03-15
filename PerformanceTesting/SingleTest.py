from tinydb import Query


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
        self.MMUPS = -1

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
                  &  (q.openMP == self.openMP)
                  & (q.vec == self.vec)
                  & (q.RMM == self.RMM)
                  & (q.size == self.size))

        # TODO: OR just to insert and deal with doubles later. MUCH faster than upsert, as no query needed
        # db.insert(string)
