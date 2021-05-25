pipeline {
  agent any
  stages {
    stage('prep') {
        steps{
            sh 'printenv'
            sh 'ls -a'
            sh "git branch -a"
            sh '''
                for versionFilePath in $(git --no-pager diff --name-only origin/${CHANGE_TARGET} | grep VERSION);
                do
                    echo ${versionFilePath};
                done;
               '''
        }
    }
    stage('build docker container') {
        when {
            anyOf {
                branch 'master'
                branch 'develop'
            }
        }
        steps{
            sh 'printenv'
        }
       }
    }
  }
