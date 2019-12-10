#!/usr/bin/groovy
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
      when { expression { return null } }
      agent { label 'atsccs11' }
      stages {
        stage('build with autopas') {
          parallel {
            stage('build sequential') {
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
            stage('build MPI') {
              steps {
                unstash 'repo'
                sh "rm -rf libs/ALL/ALL"
                sh "cp -r /work/jenkins/ALL libs/ALL/ALL"
                dir ("build-mpi"){
                  sh """
                    CC=mpicc CXX=mpicxx cmake -DENABLE_ALLLBL=ON -DENABLE_MPI=ON -DENABLE_AUTOPAS=ON -DOPENMP=ON -DENABLE_UNIT_TESTS=1 ..
                    make -j8
                  """
                }
                stash includes: "build-mpi/src/MarDyn", name: "autopas_mpi_exec"
              }
            }
          }
        }
        stage('unit test with autopas') {
          parallel {
            stage('test sequential') {
              steps {
                dir("seq"){
                  unstash 'repo'
                  unstash 'autopas_exec'
                  dir ("build"){
                    sh """
                      ./src/MarDyn -t -d ../test_input/
                    """
                  }
                }
              }
            }
            stage('test mpi') {
              steps {
                dir("mpi"){
                  unstash 'repo'
                  unstash 'autopas_mpi_exec'
                  dir ("build-mpi"){
                    sh """
                      mpirun -n 1 ./src/MarDyn -t -d ../test_input/
                      mpirun -n 4 ./src/MarDyn -t -d ../test_input/
                    """
                  }
                }
              }
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
            stage('run with autopas DD') {
              steps {
                dir('ddtest'){
                  unstash 'repo'
                  unstash 'autopas_mpi_exec'
                  dir ("build-mpi"){
                    sh """
                      mpirun -n 4 ./src/MarDyn ../examples/Argon/200K_18mol_l/config_autopas_aos.xml --steps=20 | tee autopas_run_log.txt
                      grep "Simstep = 20" autopas_run_log.txt > simstep20.txt
                      grep "T = 0.000633975" simstep20.txt
                      grep "U_pot = -2.14161" simstep20.txt
                      grep "p = 5.34057e-07" simstep20.txt
                    """
                  }
                }
              }
            }
            stage('run with autopas ALLLB') {
              steps {
                dir('alltest'){
                  unstash 'repo'
                  unstash 'autopas_mpi_exec'
                  dir ("build-mpi"){
                    sh """
                      mpirun -n 2 ./src/MarDyn ../examples/Argon/200K_18mol_l/config_autopas_lc_ALL.xml --steps=20 | tee autopas_run_log.txt
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
      matrix {
        agent none
        axes {
          axis {
            name 'VECTORIZE_CODE'
            values "SSE", "NOVEC", "AVX", "AVX2", "KNL_MASK", "KNL_G_S"
          }
          axis {
            name 'TARGET'
            values 'DEBUG', 'RELEASE'
          }
          axis {
            name 'OPENMP'
            values '0', '1'
          }
          axis {
            name 'PARTYPE'
            values 'PAR', 'SEQ'
          }
          axis {
            name 'PRECISION'
            values 'SINGLE', 'DOUBLE', 'MIXED'
          }
        }
        excludes {
          exclude {
            axis {
              name 'VECTORIZE_CODE'
              notValues 'AVX2'
            }
            axis {
              name 'PRECISION'
              notValues 'DOUBLE'
            }
          }
        }
        stages {
          stage('Build') {
            agent { label VECTORIZE_CODE }
            steps {
              echo "Seems to work"
            }
          }
          stage('Test') {
            agent { label VECTORIZE_CODE }
            steps {
              echo "Yep"
            }
          }
        }
      }
    }
    stage('post-build'){
      when { expression { return null } }
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
            unstash 'repo'
            dir('src') {
              sh label: "Build libs", script: """
                pwd
                make -f ../makefile/Makefile.lib UNIT_TESTS=0 VTK=0 PARTYPE=SEQ TARGET=RELEASE -j3 clean
                make -f ../makefile/Makefile.lib UNIT_TESTS=0 VTK=0 PARTYPE=SEQ TARGET=RELEASE -j3 lib
              """
            }

            dir ('tools/gui') {
              sh label: "Build generators", script: """
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

            dir ('src') {
              sh label: "Build again normally", script: """
                export VTK_INCDIR=/usr/include/vtk-6.3
                export VTK_LIBDIR=/usr/lib7
                make cleanall
                make TARGET=DEBUG PARTYPE=PAR PRECISION=DOUBLE OPENMP=0 \
                  REDUCED_MEMORY_MODE=0 UNIT_TESTS=1 VECTORIZE_CODE=AVX2 VTK=1 -j4
              """
            }
            dir ('executablerun') {
                sh label: "Mktcts generator", script: """
                  mv ../src/MarDyn*.PAR_DEBUG_AVX2 ./MarDyn
                  chmod 770 ./MarDyn
                  mpirun -n 2 ./MarDyn ../examples/Generators/mkTcTS/config.xml \
                    --steps 100 --final-checkpoint=0
                """
            }

            dir ('tools/standalone-generators/build') {
                sh label: "Build package", script: """
                    cmake ..
                    make package
                    make
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
                sh 'make -j4'
              }
              dir('build') {
                sh 'cmake .. && make -j4'
              }
            }
          }
        }
      }
    }
  }
}
