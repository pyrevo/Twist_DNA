pipeline {
  agent { dockerfile {
      filename 'tests/dockerfiles/twist_dna.dockerfile'
      dir './'
    }
  }
  stages {
    stage('Dry run tests') {
      steps {
        sh 'snakemake --version && echo "TEST 3"'
      }
    }
  }
}
