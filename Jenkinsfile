def isPullRequest = false

pipeline {
  agent any
  stages {
    stage('First stage') {
      steps {
        script {
          isPullRequest = env.CHANGE_ID != null
        }
        sh 'printenv'
      }
    }
    stage('Build - dependent software container') {
        when {
            anyOf {
                   expression { isPullRequest == false && env.BRANCH_NAME == 'master' }
                   expression { isPullRequest == false && env.BRANCH_NAME == 'develop' }
            }
        }
        environment {
            DOCKERHUB_CREDS = credentials('dokcerhhub')
            DOCKER_REGISTRY_URL = credentials('dockerhub_registry_url')
            CGU_CREDS = credentials('cgu-registry')
            CGU_REGISTRY_URL = credentials('cgu_registry_url')
            REPOSITORY = credentials('REPOSITORY')
        }
        steps{
            sh '''
                for versionFilePath in $( git --no-pager diff --name-only @~..@ | grep VERSION);
                do
                    folder=${versionFilePath%"/VERSION"};
                    IMAGE_NAME=${folder##*/};

                    tmpName="image-$RANDOM";
                    docker build $folder --file $folder/Dockerfile --tag $tmpName --no-cache;
                    VERSION=${BRANCH_NAME};


                    echo "${DOCKERHUB_CREDS_PSW}" | docker login ${DOCKER_REGISTRY_URL} -u ${DOCKERHUB_CREDS_USR} --password-stdin;
                    IMAGE_ID=${DOCKER_REGISTRY_URL}/${REPOSITORY}/$IMAGE_NAME;
                    docker tag $tmpName $IMAGE_ID:$VERSION;
                    docker push $IMAGE_ID:$VERSION;
                    docker logout;

                    echo "${CGU_CREDS_PSW}" | docker login ${CGU_REGISTRY_URL} -u ${CGU_CREDS_USR} --password-stdin
                    IMAGE_ID="docker-registry.cgu10.igp.uu.se/gmsuppsala/$IMAGE_NAME";
                    docker tag $tmpName $IMAGE_ID:$VERSION;
                    docker push $IMAGE_ID:$VERSION;
                    docker logout;
                done;
               '''
        }
    }
    stage('Build - test container') {
        when {
            anyOf {
                   expression { isPullRequest == false && env.BRANCH_NAME == 'master' }
                   expression { isPullRequest == false &&  env.BRANCH_NAME == 'develop' }
            }
        }
        environment {
            CGU_CREDS = credentials('cgu-registry')
            CGU_REGISTRY_URL = credentials('cgu_registry_url')
            REPOSITORY = credentials('REPOSITORY')
        }
        steps{
            sh '''
                tmpName1="image1-$RANDOM";
                docker build . --file tests/dockerfiles/test_env.dockerfile --tag $tmpName1 --no-cache;
                VERSION=${BRANCH_NAME};
                IMAGE_NAME1="test_env_somatic"
                echo "${CGU_CREDS_PSW}" | docker login ${CGU_REGISTRY_URL} -u ${CGU_CREDS_USR} --password-stdin;
                IMAGE_ID1="docker-registry.cgu10.igp.uu.se/gmsuppsala/$IMAGE_NAME1";
                docker tag $tmpName1 $IMAGE_ID1:$VERSION;
                docker push $IMAGE_ID1:$VERSION;
                
                tmpName2="image2-$RANDOM";
                docker build . --file tests/dockerfiles/test_env_full.dockerfile --tag $tmpName2  --no-cache;
                VERSION=${BRANCH_NAME};
                IMAGE_NAME2="test_env_somatic_full"
                echo "${CGU_CREDS_PSW}" | docker login ${CGU_REGISTRY_URL} -u ${CGU_CREDS_USR} --password-stdin;
                IMAGE_ID2="docker-registry.cgu10.igp.uu.se/gmsuppsala/$IMAGE_NAME2";
                docker tag $tmpName2 $IMAGE_ID2:$VERSION;
                docker push $IMAGE_ID2:$VERSION;
                docker logout;
               '''
        }
    }
    stage('Dry run tests') {
        when {
            expression { isPullRequest == true }
        }
        agent {
            dockerfile {
                    filename 'tests/dockerfiles/twist_dna_working.dockerfile'
                    dir './'
                    args '-u 0 --privileged -v $HOME/.m2:/home/jenkins/.m2:ro -v /beegfs-storage:/beegfs-storage:ro'
                    registryUrl 'https://docker-registry.cgu10.igp.uu.se'
                    registryCredentialsId 'cgu-registry'
           }
        }
        steps {
            sh 'cp -r /data_twist_dna_fgbio . && snakemake -n -s /Twist_DNA/Twist_DNA.smk --directory /data_twist_dna_fgbio'
            sh 'cp -r /data_twist_dna_markdup . && snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_markdup'
            sh 'cp -r /data_twist_dna_gpu . && snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_gpu'
            sh 'cp -r /data_twist_dna_cutadapt . && snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_cutadapt'
            sh 'cp -r /data_twist_dna_fastp . && snakemake -n -s Twist_DNA.smk --directory /data_twist_dna_fastp'
            sh 'cp -r /data_gms_somatic . && snakemake -n -s gms_somatic.smk --directory /data_gms_somatic/'
            sh 'cp -r /data_demultiplex . && snakemake -n -s demultiplex.smk --directory /data_demultiplex/'
        }
     }
     stage('Small dataset 1') {
     when {
         anyOf {
                expression { env.BRANCH_NAME == 'master' || isPullRequest == true }
                expression { env.BRANCH_NAME == 'develop' || isPullRequest == true }
            }
        }
        agent {
            dockerfile {
                 filename 'tests/dockerfiles/twist_dna_working_full.dockerfile'
                 dir './'
                 args '-u 0 --privileged -v $HOME/.m2:/home/jenkins/.m2:ro -v /beegfs-storage:/beegfs-storage:ro'
                 registryUrl 'https://docker-registry.cgu10.igp.uu.se'
                 registryCredentialsId 'cgu-registry'
            }
        }
        steps {
            sh 'snakemake -j 4 -s /Twist_DNA/gms_somatic.smk --directory /data_gms_somatic --use-singularity --singularity-prefix /Twist_DNA --singularity-args  " --cleanenv --bind /beegfs-storage  --bind /projects --bind /data --bind /Twist_DNA "'
            sh 'snakemake -j 4 -s /Twist_DNA/Twist_DNA.smk --directory /data_twist_dna_markdup --use-singularity --singularity-prefix /Twist_DNA --singularity-args  " --cleanenv --bind /beegfs-storage  --bind /projects --bind /data --bind /Twist_DNA "'
        }
    }
    stage('Small dataset 2') {
       when {
           anyOf {
                   expression { env.BRANCH_NAME == 'master' && isPullRequest == false }
                   expression { env.BRANCH_NAME == 'develop' && isPullRequest == false }
           }
       }
       agent {
           dockerfile {
                filename 'tests/dockerfiles/twist_dna_working_full.dockerfile'
                dir './'
                args '-u 0 --privileged -v $HOME/.m2:/home/jenkins/.m2:ro -v /beegfs-storage:/beegfs-storage:ro'
                registryUrl 'https://docker-registry.cgu10.igp.uu.se'
                registryCredentialsId 'cgu-registry'
           }
       }
       steps {
           sh 'snakemake -j 4 -s /Twist_DNA/Twist_DNA.smk --directory /data_twist_dna_fgbio --use-singularity --singularity-prefix /Twist_DNA --singularity-args  " --cleanenv --bind /beegfs-storage  --bind /projects --bind /data --bind /Twist_DNA "'
           sh 'snakemake -j 4 -s /Twist_DNA/Twist_DNA.smk --directory /data_twist_dna_cutadapt --use-singularity --singularity-prefix /Twist_DNA --singularity-args  " --cleanenv --bind /beegfs-storage  --bind /projects --bind /data --bind /Twist_DNA "'
           sh 'snakemake -j 4 -s /Twist_DNA/Twist_DNA.smk --directory /data_twist_dna_fastp --use-singularity --singularity-prefix /Twist_DNA --singularity-args  " --cleanenv  --bind /beegfs-storage  --bind /projects --bind /data --bind /Twist_DNA "'
       }
     }
  }
  post {
     always {
         sh 'docker image prune -fa'
     }
  }
}
