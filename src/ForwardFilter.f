      SUBROUTINE ForwardFilter(nvar, bigRk, bigK, bigT, 
     &     nbeta, sigdraw, xidraw, llh, pfilt)

      Implicit None

c     External function calls
      EXTERNAL detabs, inv
      EXTERNAL MatrixMultiply, MatrixTranspose

cccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
c     INPUTS

      INTEGER nvar, bigK, bigT, nbeta

c     sigdraw = residual VCVM
      DOUBLE PRECISION  sigdraw(nvar,nvar,bigK)
      DOUBLE PRECISION  xidraw(bigK,bigK)

cccccc
c     OUTPUTS

c     Log-likelihood value
      DOUBLE PRECISION  llh

c     Filtered probabilities (including initial), (bigT+1) by bigK
      DOUBLE PRECISION  pfilt(bigT+1,bigK)
   

cccccccccccccccccccccccccccccccccccccccccccccccc

c     Residuals
      DOUBLE PRECISION  bigRk(bigT,nvar,bigK)

c     Reshaped beta
      DOUBLE PRECISION  betak(nbeta, nvar)

      DOUBLE PRECISION  bigXbetak(bigT, nvar)

c     matrix of log likelihood by regime at each observation
      DOUBLE PRECISION  llht(bigT,bigK)

c     maximum of likelihoods, ma = max-adjusted
      DOUBLE PRECISION  lhma(bigT,bigK)

c     filtering
      DOUBLE PRECISION  sp, st, pinit(1,bigK), ptmp(bigK,1)
      DOUBLE PRECISION  transpfilt(bigK,1), transptmp(1,bigK)
      DOUBLE PRECISION  transxi(bigK,bigK), translhma(bigK,1)
      DOUBLE PRECISION  matmulfilt(bigK,1)
      DOUBLE PRECISION  summaxllh

c     loop iteration counters
      INTEGER itk,itt

c     temporary storage
      DOUBLE PRECISION  tmpe(1,nvar), transtmpe(nvar,1)
      DOUBLE PRECISION  tmpsig(nvar,nvar), tmpsigld(1,1)
      DOUBLE PRECISION  tmpsiginv(nvar,nvar)
      DOUBLE PRECISION  matmulone(1,nvar), matmultwo(1,1)
      DOUBLE PRECISION  tmpmax

c     pi, obtained using vpa(pi,40) in Matlab
      DOUBLE PRECISION pi
      PARAMETER(pi = 3.141592653589793238462643383279502884197)


cccccccccccccccccccccccccccccccccccccccccccccccc
c

      llh = 0.0

      itk=1

      DO itk = 1,bigK

cccccc
c     Setup sigma matrices
         tmpsig = sigdraw(:,:,itk)
         tmpsigld = tmpsig(1,1)
         tmpsiginv(1,1) = 1/tmpsig(1,1)

c     calculate determinant of Sigma if m > 1
         IF (nvar > 1) THEN 
            CALL detabs(nvar,tmpsig(:,:),tmpsigld)
            CALL inv(nvar,tmpsig(:,:),tmpsiginv)
         END IF

c     Log the determinant
         tmpsigld = LOG(tmpsigld)

cccccc

         itt = 1

         DO itt = 1,bigT

c           tmpe = 1 by nvar, with residuals at time itt
            tmpe = bigRk(itt:itt,1:nvar,itk)

c           transtmpe = nvar by 1, transpose of tmpe
            CALL MatrixTranspose(tmpe, 1, nvar, transtmpe)

            CALL MatrixMultiply(tmpe, 1, nvar, tmpsiginv, nvar, nvar, 
     &           matmulone, 1)

            CALL MatrixMultiply(matmulone, 1, nvar, transtmpe, nvar, 1, 
     &                          matmultwo, 1)

            llht(itt:itt,itk:itk) = -(REAL(nvar)/2.0)*LOG(2.0*pi) 
     &           - (1.0/2.0)*tmpsigld - (1.0/2.0)*matmultwo
         END DO

      END DO


cccccccccccccccccccccccccccccccccccccccccccccccc
c     max-adjusted llht for numerical stabilization

      summaxllh = 0.0
      DO itt = 1,bigT
c        set temporary max to first regime
         tmpmax = llht(itt,1)
         IF (bigK > 1) THEN 
            DO itk = 2,bigK
               IF (llht(itt,itk) > tmpmax) THEN
                  tmpmax = llht(itt,itk)
               END IF
            END DO
         END IF
         summaxllh = summaxllh + tmpmax
         DO itk = 1,bigK
            lhma(itt,itk) = EXP(REAL(llht(itt,itk)) 
     &                              - REAL(tmpmax))
         END DO
      END DO

cccccccccccccccccccccccccccccccccccccccccccccccc
c     Forward filtering

      sp = 0.0
      pfilt = 0.0

      IF (bigK > 1) THEN 
c        initialize
         pinit(1,:) = 1.0/bigK
         pfilt(1,:) = pinit(1,:)

c        get transpose
         CALL MatrixTranspose(xidraw, bigK, bigK, transxi)

         DO itt = 1,bigT

            CALL MatrixTranspose(pfilt(itt,:), 1, bigK, transpfilt)

            CALL MatrixMultiply(transxi, bigK, bigK, 
     &           transpfilt, bigK, 1, 
     &           matmulfilt, 1)

            CALL MatrixTranspose(lhma(itt,:), 1, bigK, translhma)
            ptmp = matmulfilt*translhma

c           now to normalize, sum up and divide.  these probabilities become
c           pfilt values at next iteration, t+1
            st = 0.0
            DO itk = 1,bigK
               st = st + REAL(ptmp(itk,1))
            END DO
            CALL MatrixTranspose(ptmp, bigK, 1, transptmp)
            pfilt(itt+1,:) = transptmp(1,:) / st

            sp = sp + LOG(st)
            
         END DO
      END IF

c     Add back max adjustment
      llh = sp + summaxllh
      
      RETURN


      END SUBROUTINE ForwardFilter
