#!/usr/bin/groovy

def printVariation(def it) {
  println ("Variation: " +
    "VECTORIZE_CODE is " + it[0] + " \n" +
    "TARGET is " + it[1] + " \n" +
    "OPENMP is " + it[2] + " \n" +
    "PARTYPE is " + it[3] + " \n" +
    "PRECISION is " + it[4] + " \n" +
    "REDUCED_MEMORY_MODE is " + it[5] + " \n")
}

def ciMatrix = [
  ["SSE","AOS","AVX","AVX2","SOA","KNL_MASK","KNL_G_S"], // VECTORIZE_CODE
  ["DEBUG","RELEASE"],                                   // TARGET
  ["0","1"],                                             // OPENMP
  ["PAR","SEQ"],                                         // PARTYPE
  ["DOUBLE","SINGLE","MIXED"],                           // PRECISION
  ["0","1"]                                              // REDUCED_MEMORY_MODE
]

def combinationFilter(def it) {
  def VECTORIZE_CODE      = it[0]
  def TARGET              = it[1]
  def OPENMP              = it[2]
  def PARTYPE             = it[3]
  def PRECISION           = it[4]
  def REDUCED_MEMORY_MODE = it[5]

  return (
    ( // Mandatory restrictions
      ( // REDUCED_MEMORY_MODE can't be w/ AOS
        (VECTORIZE_CODE=="AOS").implies(REDUCED_MEMORY_MODE=="0")
      )
    ) &&
    (
      ( // Combination filters
        ( // 42
          TARGET=="RELEASE" && OPENMP=="1" && PARTYPE=="PAR"
        ) ||
        ( // 16
          VECTORIZE_CODE=="AVX2" && PRECISION=="DOUBLE"
        )
      )
    )
  )
}

// Holds the build results
def results = [:]
// Holds the id of the allocated slurm job and state
def knl_jobid
def knl_jobstate

pipeline {
  agent none
  options {
    buildDiscarder logRotator(artifactDaysToKeepStr: '7', artifactNumToKeepStr: '1', daysToKeepStr: '7', numToKeepStr: '2')
    skipDefaultCheckout true
    skipStagesAfterUnstable()
    timeout(activity: true, time: 3, unit: 'HOURS')
  }
  stages{
    stage('fetch git repo') {
      agent { label 'atsccs11_prio' }
      steps {
        checkout scm
        stash 'repo'
      }
    }
    stage('check AutoPas integration') {
      agent { label 'atsccs11' }
      stages {
        stage('build with autopas') {
          steps {
            unstash 'repo'
            dir ("build"){
              sh """
                cmake -DENABLE_AUTOPAS=ON -DOPENMP=ON -DENABLE_UNIT_TESTS=1 ..
                make -j8
              """
            }
            stash includes: "build/src/MarDyn", name: "autopas_exec"
          }
        }
        stage('unit test with autopas') {
          steps {
            unstash 'repo'
            unstash 'autopas_exec'
            dir ("build"){
              sh """
                ./src/MarDyn -t -d ../test_input/
              """
            }
          }
        }
        stage('validation tests with autopas') {
          parallel {
            stage('run with autopas aos') {
              steps {
                dir('aostest'){
                  unstash 'repo'
                  unstash 'autopas_exec'
                  dir ("build"){
                    sh """
                      ./src/MarDyn ../examples/Argon/200K_18mol_l/config_autopas_aos.xml --steps=20 | tee autopas_run_log.txt
                      grep "Simstep = 20" autopas_run_log.txt > simstep20.txt
                      grep "T = 0.000633975" simstep20.txt
                      grep "U_pot = -2.14161" simstep20.txt
                      grep "p = 5.34057e-07" simstep20.txt
                    """
                  }
                }
              }
            }
            stage('run with autopas soa') {
              steps {
                dir('soatest'){
                  unstash 'repo'
                  unstash 'autopas_exec'
                  dir ("build"){
                    sh """
                      ./src/MarDyn ../examples/Argon/200K_18mol_l/config_autopas_soa.xml --steps=20 | tee autopas_run_log.txt
                      grep "Simstep = 20" autopas_run_log.txt > simstep20.txt
                      grep "T = 0.000633975" simstep20.txt
                      grep "U_pot = -2.14161" simstep20.txt
                      grep "p = 5.34057e-07" simstep20.txt
                    """
                  }
                }
              }
            }
          }
        }
      }
    }
    stage('run ci-matrix') {
      steps {
        println "setting up and running ci-matrix..."
        script {
          def variations = [:]
          def matrixEntry = []
          def matrixBuilder // needs to be defined beforehand to overload

          // construct and return jobs
          def variation = {
            def VECTORIZE_CODE      = it[0]
            def TARGET              = it[1]
            def OPENMP              = it[2]
            def PARTYPE             = it[3]
            def PRECISION           = it[4]
            def REDUCED_MEMORY_MODE = it[5]

            // return the job
            if (combinationFilter(it)) {
              return {
                node(VECTORIZE_CODE) {
                  def ARCH = (NODE_NAME == "KNL-Cluster-login") ? "KNL" : "HSW"
                  def build_result = "not run"
                  def unit_test_result = "not run"
                  def validation_test_result = "not run"
                  try {
                    stage("${it.join('-')}") {
                      sh "rm -rf ${it.join('-')} || echo ''"
                      ws("${WORKSPACE}/${it.join('-')}") {
                        unstash 'repo'
                        stage("build/${it.join('-')}") {
                          try {
                            printVariation(it)
                            dir("src") {
                              if(ARCH == "HSW") {
                                sh "export VTK_INCDIR=/usr/include/vtk-6.3"
                                sh "export VTK_LIBDIR=/usr/lib"
                              } else {
                                sh "source /etc/profile.d/modules.sh"
                              }
                              sh """
                                make TARGET=$TARGET \
                                  PARTYPE=$PARTYPE \
                                  PRECISION=$PRECISION \
                                  OPENMP=$OPENMP \
                                  REDUCED_MEMORY_MODE=$REDUCED_MEMORY_MODE \
                                  UNIT_TESTS=1 \
                                  VECTORIZE_CODE=$VECTORIZE_CODE \
                                  """ +
                                  ((ARCH=="KNL") ?
                                    "VTK=0 CFG=" + (PARTYPE == "PAR" ? "icc-impi" : "icc") :
                                    "VTK=1"
                                  ) + " -j4"
                              sh "mv `readlink MarDyn` ${it.join('-')}"
                              if (it.join('-') == "AVX2-DEBUG-0-PAR-DOUBLE-0") {
                                stash includes: "AVX2-DEBUG-0-PAR-DOUBLE-0", name: "build"
                              }
                              build_result = "success"
                            }
                          } catch (err) {
                            build_result = "failure"
                            error err
                          }
                        }
                        // FIXME: Mixed precision unit-tests fail with rmm
                        if ((PRECISION=="MIXED").implies(REDUCED_MEMORY_MODE=="0")) {
                          stage("unit-test/${it.join('-')}") {
                            try {
                              printVariation(it)
                              // Wait for allocation if necessary
                              while (ARCH=="KNL" && knl_jobstate=="PENDING") {
                                sleep 150
                              }
                              if (ARCH=="HSW" && PARTYPE=="PAR") {
                                sh "mpirun -n 4 ./src/${it.join('-')} -t -d ./test_input/"
                              } else if (ARCH=="HSW" && PARTYPE=="SEQ") {
                                sh "./src/${it.join('-')} -t -d ./test_input/"
                              } else if (ARCH=="KNL" && PARTYPE=="PAR") {
                                sh """
                                  source /etc/profile.d/modules.sh
                                  export OMP_NUM_THREADS=64
                                  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
                                  export SLURM_CONF=$HOME/slurm.conf
                                  while : ; do
                                    set +e
                                    output=`srun --jobid=$knl_jobid -n 4 ./src/${it.join('-')} -t -d ./test_input/ 2>&1`
                                    rc=\$?
                                    if [[ \$rc == 1 && (\$output == *"Job violates accounting/QOS policy"* || \$output == *"Socket timed out on send/recv"* || \$output == *"Communication connection failure"*)]] ; then
                                      echo "srun submit limit reached or socket timed out error, trying again in 60s"
                                      sleep 60
                                      continue
                                    fi
                                    set -e
                                    echo \$output
                                    exit \$rc
                                  done
                                """
                              } else if (ARCH=="KNL" && PARTYPE=="SEQ") {
                                sh """
                                  source /etc/profile.d/modules.sh
                                  export OMP_NUM_THREADS=64
                                  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
                                  export SLURM_CONF=$HOME/slurm.conf
                                  while : ; do
                                    set +e
                                    output=`srun --jobid=$knl_jobid -n 1 ./src/${it.join('-')} -t -d ./test_input/ 2>&1`
                                    rc=\$?
                                    if [[ \$rc == 1 && (\$output == *"Job violates accounting/QOS policy"* || \$output == *"Socket timed out on send/recv"*)]] ; then
                                      echo "srun submit limit reached or socket timed out error, trying again in 60s"
                                      sleep 60
                                      continue
                                    fi
                                    set -e
                                    echo \$output
                                    exit \$rc
                                  done
                                """
                              }
                              unit_test_result = "success"
                            } catch (err) {
                              unit_test_result = "failure"
                              error err
                            }
                            xunit([CppUnit(deleteOutputFiles: true, failIfNotNew: false, pattern: 'results.xml', skipNoTestFiles: false, stopProcessingIfError: true)])
                          }
                        }
                        if (PRECISION=="DOUBLE" && REDUCED_MEMORY_MODE=="0") {
                          stage("validation-test/${it.join('-')}") {
                            try {
                              copyArtifacts filter: '**/*', fingerprintArtifacts: true, projectName: 'MardynUpdateValidationBase', selector: lastSuccessful()
                              sh "pwd && ls"
                              def configDirs = [
                                "simple-lj",
                                "simple-lj-mp",
                                "surface-tension_LRC_CO2_Merker_280",
                                "simple-lj_kdd",
                                "simple-lj-direct-mp",
                                "simple-lj-direct"
                              ]
                              for (configDirVar in configDirs) {
                                println configDirVar
                                dir("validation/validationBase") {
                                  def sameParTypeOption = (fileExists ('../validationInput/' + configDirVar + '/compareSameParType')) ? "" : "-b"
                                  def mpicmd = (PARTYPE=="PAR" ? "-m 4" : "-m 1")
                                  def mpiextra = (ARCH=="KNL" ? "-M \"srun --jobid=" + knl_jobid + "\"" : "")
                                  def allmpi = (ARCH=="KNL" ? "--allMPI" : "")
                                  def srunfix = (ARCH=="KNL" ? "--srunFix" : "")
                                  def icount = (ARCH=="KNL" ? "-I 5" : "")
                                  def oldexec = "MarDyn_ValidationBase" + (ARCH=="KNL" ? "_KNL" : "")
                                  def plugins = (VECTORIZE_CODE=="AOS" ? "" : "--disablePlugin RDF") +
                                                ((fileExists ('../validationInput/' + configDirVar + '/noGamma') ) ?
                                                " --disablePlugin GammaWriter" : "")
                                  for (partypeToCompare in ["PAR", "SEQ"]) {
                                    if ((sameParTypeOption=="").implies(partypeToCompare == PARTYPE)) {
                                      dir(partypeToCompare) {
                                        printVariation(it)
                                        if (!fileExists ("./${oldexec}")) {
                                          sh "ls -Ra"
                                          error 'Validation base not found!'
                                        } else if (!fileExists ("${WORKSPACE}/src/${it.join('-')}")) {
                                          sh "ls -Ra ../../.."
                                          error 'New build not found!'
                                        }
                                        sh """
                                          export OMP_NUM_THREADS=64
                                          export SLURM_CONF=$HOME/slurm.conf
                                          echo "Running validation tests"
                                          rm ${it.join('-')} || echo ""
                                          cp ${WORKSPACE}/src/${it.join('-')} .
                                          pwd && ls
                                          ../../validationRun/validationRun.py \
                                            $srunfix \
                                            -o ./$oldexec \
                                            -n ./${it.join('-')} \
                                            $allmpi \
                                            -c \"\$(realpath \$(find ../../validationInput/$configDirVar/ -type f -name *.xml) )\" \
                                            -i \"\$(realpath \$(find ../../validationInput/$configDirVar/ -type f \\( -iname \\*.dat -o -iname \\*.inp -o -iname \\*.xdr \\)) )\" \
                                            $plugins $icount $sameParTypeOption $mpicmd $mpiextra
                                        """
                                      }
                                    }
                                  }
                                }
                              }
                              validation_test_result = "success"
                            } catch (err) {
                              validation_test_result = "failure"
                              error err
                            }
                          }
                        }
                      }
                    }
                    results.put(it.join('-'), [:])
                    results[it.join('-')].put("runOn", ARCH)
                    results[it.join('-')].put("build", build_result)
                    results[it.join('-')].put("unit-test", unit_test_result)
                    results[it.join('-')].put("validation-test", validation_test_result)
                  }
                  catch (err) {
                    results.put(it.join('-'), [:])
                    results[it.join('-')].put("runOn", ARCH)
                    results[it.join('-')].put("build", build_result)
                    results[it.join('-')].put("unit-test", unit_test_result)
                    results[it.join('-')].put("validation-test", validation_test_result)
                    error err
                 }
                }
              }
            }
          }

          // return all combinations
          // FIXME can't be defined globally due to a bug in Jenkins:
          // https://issues.jenkins-ci.org/browse/JENKINS-49826
          matrixBuilder = { def matrix, int level ->
            variations.failFast = true
            // HACK Jobs to manage resource allocation on the knl cluster
            variations["slurm"] = {
              try {
                node("KNL_PRIO") { // Executor on the CoolMUC3 login node reserved for slurm allocation and management
                  parallel "allocation": {
                    try {
                      timeout(time: 6, unit: 'HOURS') {
                        // Allocate a new interactive job with up to three nodes
                        // and two hours maximum run time. The ci-matrix will
                        // attach subjobs to this via srun and the slurm.slurmcontrol
                        // stage will revoke the allocation as soon as we are done.
                        // Until that, the shell needs to stay open for the job
                        // to survive. Hence the second "sleep 7200".
                        sh """
                          export SLURM_CONF=$HOME/slurm.conf
                          salloc --job-name=mardyn-test --nodes=1-4 --partition=mpp3_batch\
                            --tasks-per-node=3 --time=00:45:00 --begin=now+150\
                            sleep 45m || echo 0
                          sleep 6h
                        """
                      }
                    }
                    // Stop the slurm.allocation on any error, but do not
                    // propatate. Needed to allow slurm.slurmcontrol (below)
                    // to kill slurm.* by throwing an error.
                    catch (err) {
                      println err
                    }
                  }, "slurmcontrol": {
                    // Wait for slurm.allocation to work its magic
                    sleep 10
                    // Store jobid
                    knl_jobid = sh(
                      returnStdout: true,
                      script: 'export SLURM_CONF=$HOME/slurm.conf && squeue -O jobid | sed -n 2p'
                    ).replace("\n", "")
                    println "Scheduled job " + knl_jobid
                    // Wait for all KNL jobs to finish by comparing the list
                    // of scheduled jobs with the list of results
                    while (results.count { key, value -> key.contains("KNL") } < variations.count { key, value -> key.contains("KNL") }) {
                      sleep 150
                      knl_jobstate = sh(
                        returnStdout: true,
                        script: 'export SLURM_CONF=$HOME/slurm.conf && squeue -j $knl_jobid -O state | sed -n 2p'
                      ).replace("\n", "")
                      println knl_jobstate
                    }
                    // Revoke slurm job allocation
                    sh "scancel $knl_jobid -f --user=ga38cor3"
                    // Throw an error (that will not be propagated) to kill the
                    // slurm.* stages. Without this, slurm.allocation would keep
                    // running even after the allocation is revoked.
                    error "all finished"
                  },
                  failFast: true // If one of the slurm.* stages throws an error, slurm.* is cancelled
                }
              }
              // Only print all error signals. See above.
              catch (err) {
                println err
              }
            }
            // Assemble the ci-matrix in a map
            for ( entry in matrix[0] ) {
              matrixEntry[level] = entry
              if (matrix.size() > 1 ) {
                matrixBuilder(matrix.drop(1), level + 1)
              }
              else {
                def _variation = variation(matrixEntry.collect())
                if (_variation != null)
                  variations[matrixEntry.join("-")] = _variation
              }
            }
          }

          node { matrixBuilder(ciMatrix, 0) }
          // Lock the KNL cluster while running the ci-matrix. This is needed
          // because we are not allowed to submit more than one interactive job
          // at a time. As soon as we're in post-build, the lock is released
          // and another running job can access KNL.
          // https://doku.lrz.de/display/PUBLIC/Resource+limits+for+parallel+jobs+on+Linux+Cluster
          lock('KNL_slurm_allocation') {
            parallel variations
          }
        }
      }
    }
    stage('post-build'){
      parallel {
        stage('documentation') {
          agent { label 'atsccs11' }
          stages {
            stage('build documentation') {
              steps {
                unstash 'repo'
                sh "mkdir doxygen_doc || echo 'doxygen_doc Folder exists already'"
                sh "doxygen"
              }
            }
            stage('release documentation') {
              when { branch 'master' }
              steps {
                sh "rm -rf /home/wwwsccs/html/mardyn/doc /home/wwwsccs/html/mardyn/doxys_docs"
                sh "cp -r doc /home/wwwsccs/html/mardyn"
                sh "cp -r doxygen_doc /home/wwwsccs/html/mardyn"
                sh "chmod -R 775 /home/wwwsccs/html/mardyn/doc"
                sh "chmod -R 775 /home/wwwsccs/html/mardyn/doxygen_doc"
              }
            }
          }
        }
        stage('libs and generators') {
          agent { label 'atsccs11' }
          steps {
            // Build lib
            sh "ls"
            unstash 'repo'
            dir('src') {
              sh """
                pwd
                make -f ../makefile/Makefile.lib UNIT_TESTS=0 VTK=0 PARTYPE=SEQ TARGET=RELEASE -j3 clean
                make -f ../makefile/Makefile.lib UNIT_TESTS=0 VTK=0 PARTYPE=SEQ TARGET=RELEASE -j3 lib
              """
            }

            // Build generators
            dir ('tools/gui') {
              sh """
                export VTKINCLUDEPATH=/usr/include/vtk-6.3
                tar xfz ScenarioGenerator.tar.gz

                qmake DEFINES+="MARDYN_DPDP" DropletGenerator.pro -o Makefile.droplet
                qmake DEFINES+="MARDYN_DPDP" CubicGridGenerator.pro -o Makefile.cubic
                qmake DEFINES+="MARDYN_DPDP" AqueousNaClGenerator.pro -o Makefile.aqueous
                qmake DEFINES+="MARDYN_DPDP" CrystalLatticeGenerator.pro -o Makefile.crystal
                qmake DEFINES+="MARDYN_DPDP" MS2RSTGenerator.pro -o Makefile.ms2
                qmake DEFINES+="MARDYN_DPDP" RayleighTaylorGenerator.pro -o Makefile.rayleigh

                rm -r libMardyn* || :

                make -f Makefile.droplet -j2
                make -f Makefile.cubic -j2
                make -f Makefile.aqueous -j2
                make -f Makefile.crystal -j2
                make -f Makefile.ms2 -j2
                make -f Makefile.rayleigh -j2
              """
            }
            dir ('executablerun') {
                // Mktcts generator
                unstash "build"
                sh """
                  chmod 770 AVX2-DEBUG-0-PAR-DOUBLE-0
                  mpirun -n 2 ./AVX2-DEBUG-0-PAR-DOUBLE-0 ../examples/Generators/mkTcTS/config.xml --steps 100 --final-checkpoint=0
                """
            }
          }
        }
        stage('export-src') {
          agent { label 'atsccs11' }
          steps {
            unstash 'repo'
            sh './export-src.sh'
            dir('unarchive') {
              sh 'tar -xf ../Mardyn-src.tar.gz'
              dir('src') {
                sh 'make'
              }
            }
          }
        }
      }
    }
  }
}
