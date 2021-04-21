pipeline {
  agent { dockerfile {
      filename 'tests/dockerfiles/twist_dna.dockerfile'
      dir './'
    }
  }
  stages {
    stage('Dry run tests') {
      steps {
        sh 'snakemake -n -s /Twist_DNA/Twist_DNA.smk --directory /data_twist_dna_fgbio '
        sh 'snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_markdup'
        sh 'snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_gpu'
        sh 'snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_cutadapt'
        sh 'snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_fastp'
        sh 'snakemake -n -s gms_somatic.smk --directory /data_gms_somatic/'
        sh 'snakemake -n -s demultiplex.smk --directory /data_demultiplex/'
      }
    }
  }
}
