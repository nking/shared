// if jenkins SCM polling is set to trigger on post-commit,
// can use this for a trigger:
// curl http://jenkins:8080/git/notifyCommit?url=git@gitserver:tools/common.git
// or
// curl -G --data-urlencode url=local_git_project_path --data-urlencode token=token http://localhost:8080/git/notifyCommit
//
pipeline {
    // agent can be any, a label, or none
    agent any
    options {
        skipStagesAfterUnstable()
    }
    stages {
        stage('Build') {
            steps {
               echo 'Building..'
               sh 'ant compileTests'
            }
        }
        stage('Test') {
            steps {
               echo 'Testing..'
               sh 'ant runTests'
            }
        }
        stage('Deploy') {
            steps {
               // increment version in shared.mf too
               echo 'Deploying..'
               //sh 'ant runCoverage'
               //sh 'ant javadoc'
               //sh 'ant package'
            }
        }
    }
}
