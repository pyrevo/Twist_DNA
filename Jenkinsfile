pipeline {
  agent any

  stages {
    stage('build') {
      steps {
        sh 'pip3.8 install -r requirements.txt'
      }
    }
    stage('test') {
      steps {
        sh 'python test.py'
      }
    }
  }
}
