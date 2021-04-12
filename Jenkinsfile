pipeline {
  agent any

  stages {
    stage('build') {
      steps {
        sh '''virtualenv -p python3.8 venv
              source venv/bin/activate
              pip3.8 install -r requirements.txt'''
      }
    }
    stage('test') {
      steps {
        sh '''
              source venv/bin/activate
              pytest src/lib/python/test_utils.py
              echo "Test"
           '''
         sh '''
              source venv/bin/activate
              snakemake -n -s demultiplex.smk --directory tests/workflow_dry_run/demultiplex/              
            '''

      }
    }
  }
}
