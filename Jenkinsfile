pipeline {
  agent { dockerfile {
      filename 'tests/dockerfiles/twist_dna.dockerfile'
      dir .
      label "test-env"
    }
  }

  stages {
    stage('test') {
      steps {
        snakemake --version
      }
    }
  }
}
