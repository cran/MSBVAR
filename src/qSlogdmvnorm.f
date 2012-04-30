      SUBROUTINE qSlogdmvnorm(nbeta,nvar,bigK,S,
     &     betadraw,cpm,cpinv,logdet,qSmvn)

c     External function calls
      EXTERNAL detabs, inv
      EXTERNAL MatrixMultiply, MatrixTranspose

cccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
c     INPUTS

      INTEGER S, nbeta, nvar, bigK
      DOUBLE PRECISION  betadraw(nbeta*nvar,bigK)
      DOUBLE PRECISION  cpm(S,nbeta*nvar,bigK)
      DOUBLE PRECISION  cpinv(S,nbeta*nvar,nbeta*nvar,bigK)
      DOUBLE PRECISION  logdet(S,1,bigK)

cccccc
c     OUTPUTS

c     Log-likelihood value
      DOUBLE PRECISION  qSmvn(S,1)


cccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER its, ntotbeta
      DOUBLE PRECISION  ymu(nbeta*nvar,1), ymut(1,nbeta*nvar)
      DOUBLE PRECISION  matone(1,nbeta*nvar), mattwo(1,1)

c     pi, obtained using vpa(pi,40) in Matlab
      DOUBLE PRECISION pi
      PARAMETER(pi = 3.141592653589793238462643383279502884197)


cccccccccccccccccccccccccccccccccccccccccccccccc
c

      ntotbeta = nbeta*nvar

      qSmvn = 0.0

      its=1
      itk=1

      DO itk = 1,bigK
      DO its = 1,S
         ymu(:,1) = betadraw(:,itk) - cpm(its,:,itk)
         CALL MatrixTranspose(ymu,ntotbeta,1,ymut)
         CALL MatrixMultiply(ymut, 1, ntotbeta,
     &        cpinv(its,:,:,itk), ntotbeta, ntotbeta, matone) 
         CALL MatrixMultiply(matone, 1, ntotbeta,
     &        ymu, ntotbeta, 1, mattwo) 
         qSmvn(its,1) = qSmvn(its,1) - (REAL(ntotbeta)/2.0)*LOG(2.0*pi) 
     &        - (1.0/2.0)*logdet(its,1,itk) - (1.0/2.0)*mattwo(1,1)
      END DO
      END DO
      
      RETURN


      END SUBROUTINE qSlogdmvnorm

