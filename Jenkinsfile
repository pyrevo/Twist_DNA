pipeline {
  agent { dockerfile {
      filename 'tests/dockerfiles/twist_dna.dockerfile'
      dir './'
    }
  }
  stages {
    stage('Dry run tests') {
      steps {
        sh 'run: snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_fgbio '
        sh 'run: snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_markdup'
        sh 'run: snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_gpu'
        sh 'run: snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_cutadapt'
        sh 'run: snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_fastp'
        sh 'run: snakemake -n -s gms_somatic.smk --directory /data_gms_somatic/'
        sh 'run: snakemake -n -s demultiplex.smk --directory /data_demultiplex/'
      }
    }
  }
}
