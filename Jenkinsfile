pipeline {
  agent { docker { image 'wolfram-docker-11.3.0' } }

  options {
    timeout(time: 10, unit: 'MINUTES')
    skipDefaultCheckout(true)
  }

  stages {
    stage('Run tests') {
      steps {
        dir('KerrGeodesics') {
          checkout scm
          sh 'Tests/AllTests.wls'
          junit 'TestReport.xml'
        }
      }
    }
  }
}