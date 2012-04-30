      SUBROUTINE BackwardSampler(bigK,bigT,pfilt,xidraw,rvu,
     &     backsampbigS, transmat)

      EXTERNAL  genSt

cccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
c     INPUTS
      INTEGER bigK, bigT
c     Filtered probabilities (including initial), (bigT+1) by bigK
      DOUBLE PRECISION  pfilt(bigT+1,bigK)
c     Xi draw
      DOUBLE PRECISION  xidraw(bigK,bigK)
c     Random uniforms
      DOUBLE PRECISION  rvu(bigT+1)

cccccc
c     OUTPUTS
c     Sampled states (including initial)
      INTEGER  backsampbigS(bigT+1,1)
      INTEGER  transmat(bigK,bigK)

cccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER  Imc(bigT+1,1)
      INTEGER  STT(bigT,1)
      INTEGER  genSt
      INTEGER  itt, itk
      DOUBLE PRECISION  p(bigK,1), pfiltrans(bigK,1), st


cccccccccccccccccccccccccccccccccccccccccccccccc
c

c     Initialize Imc
      Imc(:,1) = 1

c     Set last t probability
      IF (bigK > 2) THEN
         Imc(bigT+1,1) = genSt(bigK,pfilt(bigT+1,1:(bigK-1)),
     &        rvu(bigT+1))
      ELSE 
         Imc(bigT+1,1) = genSt(bigK,pfilt(bigT+1,1),rvu(bigT+1))
      END IF

      DO itt=bigT,1,-1
         CALL MatrixTranspose(pfilt(itt,:),1,bigK,pfiltrans)
         p(:,1) = pfiltrans(:,1) * xidraw(:,Imc(itt+1,1))
         st = 0.0
         DO itk = 1,bigK
            st = st + REAL(p(itk,1))
         END DO
         p(:,1) = p(:,1) / st
         IF (bigK > 2) THEN
            Imc(itt,1) = genSt(bigK,p(1:(bigK-1),1),rvu(itt))
         ELSE
            Imc(itt,1) = genSt(bigK,p(1,1),rvu(itt))
         END IF
      END DO

      backsampbigS = Imc

c     Construct transition matrix
      STT(:,1) = backsampbigS(2:(bigT+1),1)
      transmat(:,:) = 0
      DO itt = 1,(bigT-1)
         transmat(STT(itt+1,1), STT(itt,1)) = 
     &        transmat(STT(itt+1,1), STT(itt,1)) + 1
      END DO

      RETURN
      
      END SUBROUTINE BackwardSampler


      INTEGER FUNCTION genSt(bigK,prob,rvu)

      INTEGER  bigK,itk
      DOUBLE PRECISION  prob(bigK-1,1),tmpprob,tmprv,rvu

      IF (bigK > 2) THEN
         tmpprob = 0.0
         tmprv = rvu
         genSt = 1
         DO itk = 1,(bigK-1)
            tmpprob = tmpprob + prob(itk,1)
            IF ((tmpprob < tmprv) .EQV. .TRUE.) THEN
               genSt = genSt + 1
            END IF
         END DO
      ELSE
         IF ((prob(1,1) < rvu) .EQV. .TRUE.) THEN
            genSt = 1 + 1
         ELSE
            genSt = 1
         END IF
      END IF

      RETURN
      END FUNCTION genSt
