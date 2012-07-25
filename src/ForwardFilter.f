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

c$$$      summaxllh = 0.0
c$$$      DO itt = 1,bigT
c$$$c        set temporary max to first regime
c$$$         tmpmax = llht(itt,1)
c$$$         IF (bigK > 1) THEN 
c$$$            DO itk = 2,bigK
c$$$               IF (llht(itt,itk) > tmpmax) THEN
c$$$                  tmpmax = llht(itt,itk)
c$$$               END IF
c$$$            END DO
c$$$         END IF
c$$$         summaxllh = summaxllh + tmpmax
c$$$         DO itk = 1,bigK
c$$$            lhma(itt,itk) = EXP(REAL(llht(itt,itk)) 
c$$$     &                              - REAL(tmpmax))
c$$$         END DO
c$$$      END DO

cccccccccccccccccccccccccccccccccccccccccccccccc
c     Forward filtering

      sp = 0.0
      pfilt = 0.0

c$$$      IF (bigK > 1) THEN 
c$$$c        initialize
c$$$         pinit(1,:) = 1.0/bigK
c$$$         pfilt(1,:) = pinit(1,:)
c$$$
c$$$c        get transpose
c$$$         CALL MatrixTranspose(xidraw, bigK, bigK, transxi)
c$$$
c$$$         DO itt = 1,bigT
c$$$
c$$$            CALL MatrixTranspose(pfilt(itt,:), 1, bigK, transpfilt)
c$$$
c$$$            CALL MatrixMultiply(transxi, bigK, bigK, 
c$$$     &           transpfilt, bigK, 1, 
c$$$     &           matmulfilt, 1)
c$$$
c$$$            CALL MatrixTranspose(lhma(itt,:), 1, bigK, translhma)
c$$$            ptmp = matmulfilt*translhma
c$$$
c$$$c           now to normalize, sum up and divide.  these probabilities become
c$$$c           pfilt values at next iteration, t+1
c$$$            st = 0.0
c$$$            DO itk = 1,bigK
c$$$               st = st + REAL(ptmp(itk,1))
c$$$            END DO
c$$$            CALL MatrixTranspose(ptmp, bigK, 1, transptmp)
c$$$            pfilt(itt+1,:) = transptmp(1,:) / st
c$$$
c$$$            sp = sp + LOG(st)
c$$$            
c$$$         END DO
c$$$      END IF
c$$$
c$$$c     Add back max adjustment
c$$$      llh = sp + summaxllh
      
      RETURN


      END SUBROUTINE ForwardFilter
