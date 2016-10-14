#!groovy
echo 'Starting Mardyn Pipeline'
def SVN_REVISION_GROOVY="HEAD"
node('atsccs11') {
    checkout([$class: 'SubversionSCM', additionalCredentials: [], excludedCommitMessages: '', excludedRegions: '', excludedRevprop: '', excludedUsers: '', filterChangelog: false, ignoreDirPropChanges: false, includedRegions: '', locations: [[credentialsId: '59cd9845-6050-481b-b42f-17142ea72f2f', depthOption: 'infinity', ignoreExternalsOption: true, local: '.', remote: 'svn+ssh://seckler@atsccs11.informatik.tu-muenchen.de:/home_local/repositories/svn/mardyn-mirror/MarDyn/trunk']], workspaceUpdater: [$class: 'UpdateUpdater']])
    SVN_REVISION_GROOVY = sh([script: 'svnversion .', returnStdout: true]).trim()
    stage('Build normally'){
        print SVN_REVISION_GROOVY
        build job: 'MardynBuild', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
    }
    stage('Unittests (normally)'){
        build job: 'MardynTest', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
    }
    stage('Validation Tests (normally)'){
        build job: 'MardynValidationTest', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
    }
    build job: 'mardyn-mirror-trunk-supermic'
    stage('Build SuperMIC'){
        print SVN_REVISION_GROOVY
        build job: 'MardynBuild-SuperMIC', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
    }
    stage('Unittests SuperMIC'){
        build job: 'MardynTest-supermic', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
    }
    stage('Validation Tests SuperMIC'){
        build job: 'MardynValidationTest-SuperMIC', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
    }
    stage('Docu+Lib'){
        build job: 'Mardyn Documentation', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
        build job: 'Mardyn-BuildLib', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
    }
    stage('generators'){//BuildMardynGenerators requires Lib to be build!
        build job: 'BuildMardynGenerators', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
        build job: 'mkTcTS generator (internal)', parameters: [string(name: 'SVN_REVISION_INPUT', value: "$SVN_REVISION_GROOVY")]
    }
}

