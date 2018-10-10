#!groovy
pipeline {
    agent any
    options {
        gitLabConnection('ls1-mardyn')
        disableConcurrentBuilds()
        buildDiscarder logRotator(artifactDaysToKeepStr: '', artifactNumToKeepStr: '', daysToKeepStr: '60', numToKeepStr: '5')
        skipStagesAfterUnstable()
        timeout(time: 12, unit: 'HOURS')
    }
    stages{
        stage('build and test') {
            parallel {
                stage('x86'){
                    agent { label 'atsccs11' }
                    stages {
                        stage('build (x86)') {
                            steps {
                                updateGitlabCommitStatus name: 'build (x86)', state: 'pending'
                                build job: 'mardyn-build-git-x86', parameters: [string(name: 'GIT_BRANCH', value: BRANCH_NAME)]
                            }
                            post {
                                success { updateGitlabCommitStatus name: 'build (x86)', state: 'success' }
                                failure { updateGitlabCommitStatus name: 'build (x86)', state: 'failed' }
                            }
                        }
                        stage('unit-test (x86)') {
                            steps {
                                updateGitlabCommitStatus name: 'unit-test (x86)', state: 'pending'
                                build job: 'mardyn-unit-test-git-x86', parameters: [string(name: 'GIT_BRANCH', value: BRANCH_NAME)]
                            }
                            post {
                                success { updateGitlabCommitStatus name: 'unit-test (x86)', state: 'success' }
                                failure { updateGitlabCommitStatus name: 'unit-test (x86)', state: 'failed' }
                            }
                        }
                        stage('validation-test (x86)') {
                            steps {
                                updateGitlabCommitStatus name: 'validation-test (x86)', state: 'pending'
                                build job: 'mardyn-validation-test-git-x86', parameters: [string(name: 'GIT_BRANCH', value: BRANCH_NAME)]
                            }
                            post {
                                success { updateGitlabCommitStatus name: 'validation-test (x86)', state: 'success' }
                                failure { updateGitlabCommitStatus name: 'validation-test (x86)', state: 'failed' }
                            }
                        }
                    }
                }
                stage('knl'){
                    agent { label 'KNL' }
                    stages {
                        stage('build (knl)') {
                            steps {
                                updateGitlabCommitStatus name: 'build (knl)', state: 'pending'
                                build job: 'mardyn-build-git-knl', parameters: [string(name: 'GIT_BRANCH', value: BRANCH_NAME)]
                            }
                            post {
                                success { updateGitlabCommitStatus name: 'build (knl)', state: 'success' }
                                failure { updateGitlabCommitStatus name: 'build (knl)', state: 'failed' }
                            }
                        }
                        stage('unit-test (knl)') {
                            steps {
                                updateGitlabCommitStatus name: 'unit-test (knl)', state: 'pending'
                                build job: 'mardyn-unit-test-git-knl', parameters: [string(name: 'GIT_BRANCH', value: BRANCH_NAME)]
                            }
                            post {
                                success { updateGitlabCommitStatus name: 'unit-test (knl)', state: 'success' }
                                failure { updateGitlabCommitStatus name: 'unit-test (knl)', state: 'failed' }
                            }
                        }
                        stage('validation-test (knl)') {
                            steps {
                                updateGitlabCommitStatus name: 'validation-test (knl)', state: 'pending'
                                build job: 'mardyn-validation-test-git-knl', parameters: [string(name: 'GIT_BRANCH', value: BRANCH_NAME)]
                            }
                            post {
                                success { updateGitlabCommitStatus name: 'validation-test (knl)', state: 'success' }
                                failure { updateGitlabCommitStatus name: 'validation-test (knl)', state: 'failed' }
                            }
                        }
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
                                sh "mkdir doxygen_doc"
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
                        sh """cd src
                            pwd
                            make -f ../makefile/Makefile.lib UNIT_TESTS=0 VTK=0 PARTYPE=SEQ TARGET=RELEASE -j3 clean
                            make -f ../makefile/Makefile.lib UNIT_TESTS=0 VTK=0 PARTYPE=SEQ TARGET=RELEASE -j3 lib

                            mkdir -p /home/wwwsccs/html/mardyn/lastSuccessfulBuild/Git/$BRANCH_NAME/lib/
                            cp libMardyn.so.1.0 /home/wwwsccs/html/mardyn/lastSuccessfulBuild/Git/$BRANCH_NAME/lib/ || :"""
                        // Build generators
                        sh """export VTKINCLUDEPATH=/usr/include/vtk-6.3
                            cd tools/gui
                            tar xfz ScenarioGenerator.tar.gz

                            qmake DEFINES+="MARDYN_DPDP" DropletGenerator.pro -o Makefile.droplet
                            qmake DEFINES+="MARDYN_DPDP" CubicGridGenerator.pro -o Makefile.cubic
                            qmake DEFINES+="MARDYN_DPDP" AqueousNaClGenerator.pro -o Makefile.aqueous
                            qmake DEFINES+="MARDYN_DPDP" CrystalLatticeGenerator.pro -o Makefile.crystal
                            qmake DEFINES+="MARDYN_DPDP" MS2RSTGenerator.pro -o Makefile.ms2
                            qmake DEFINES+="MARDYN_DPDP" RayleighTaylorGenerator.pro -o Makefile.rayleigh


                            if [ -e libMardyn* ]; then
                            rm -r libMardyn*
                            fi

                            #cp /home/wwwsccs/html/mardyn/lastSuccessfulBuild/lib/libMardyn.so.1.0 .
                            #ln -s libMardyn.so.1.0 libMardyn.so

                            make -f Makefile.droplet -j2
                            make -f Makefile.cubic -j2
                            make -f Makefile.aqueous -j2
                            make -f Makefile.crystal -j2
                            make -f Makefile.ms2 -j2
                            make -f Makefile.rayleigh -j2"""
                        // Build mktcts generator
                        sh """if [ -e MarDyn ]; then
                            rm MarDyn*
                            fi
                            wget https://www5.in.tum.de/mardyn/lastSuccessfulBuild/Git/$BRANCH_NAME/PAR/DEBUG/AVX2/WR_0/DOUBLE/OPENMP_0/MarDyn
                            chmod 770 MarDyn
                            pwd
                            ls

                            #mpdboot
                            #mpirun -n 2 ./MarDyn mkTcTS -c 0.06482 -C 0.6223 -N 10240 -R 2.5 -T 0.95 --steps 500 --final-checkpoint=0
                            mpirun -n 2 ./MarDyn examples/Generators/mkTcTS/config.xml --steps 100 --final-checkpoint=0 """
                    }
                }
            }
            post {
                success { updateGitlabCommitStatus name: 'post-build', state: 'success' }
                failure { updateGitlabCommitStatus name: 'post-build', state: 'failed' }
            }
        }
    }
}
