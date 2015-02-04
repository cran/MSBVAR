      SUBROUTINE FFBS(bigt,m,p,h,e,sig2,Q,f,filtprSt,SS,transmat)

      EXTERNAL HamiltonFilter

c     bigt = total number of observations in time series
c     m    = number of dependent variables.  i.e., VAR if m > 1
c     p    = number of lags: AR(p) / VAR(p)
c     h    = number of regimes

      INTEGER bigt,m,p,h,st1
      DOUBLE PRECISION e(bigt-p,m,h)
      DOUBLE PRECISION sig2(m,m,h), Q(h,h)
      DOUBLE PRECISION filtprSt(bigt-p,h)
      DOUBLE PRECISION f
      INTEGER SS(bigt-p, h), STT(bigt-p, 1)
      INTEGER transmat(h,h)
      INTEGER SScurr(h,1)
      INTEGER t


c     Forward-Filtering
      CALL HamiltonFilter(bigt,m,p,h,e,sig2,Q,f,filtprSt)

c     SS is (bigt-p)x(h) matrix of zeros currently
c     (SS was initialized in R code).
c     Using last observation from the filter output,
c     generate state
c     SScurr = current SS based on t
      st1 = 1
      CALL bingen(filtprSt(bigt-p,:), Q, st1, h, SScurr)
      SS(bigt-p,1:h) = SScurr(1:h,1)
      STT(bigt-p,1) = st1

      DO t = (bigt-p-1),1,-1
         CALL bingen(filtprSt(t,:), Q, st1, h, SScurr)
         SS(t,1:h) = SScurr(1:h,1)
         STT(t,1) = st1
      END DO

c     Construct transition matrix
      transmat = 0.0
      DO t = 1,(bigt-p-1)
         transmat(STT(t+1,1), STT(t,1)) =
     &        transmat(STT(t+1,1), STT(t,1)) + 1
      END DO

      RETURN

      END SUBROUTINE FFBS


c     see Kim & Nelson page 212
c     st1 is changed to for the next call
      SUBROUTINE bingen(prob,Q,st1,h,SScurr)

      EXTERNAL diag

      INTEGER st1, h
      DOUBLE PRECISION prob(h,1), Q(h,h)
      INTEGER SScurr(h,1)
      DOUBLE PRECISION pr0, tmpsum, MVUNI

      INTEGER i

c     identity matrix used below
      DOUBLE PRECISION Imath(h,h)

c     store identity matrix in Imath
      CALL diag(h,Imath)
c
      i = 1
      DO WHILE (i < h)
         tmpsum = 0.0
         DO j = i,h
            tmpsum = tmpsum + prob(j,1)*Q(st1,j)
         END DO
         pr0 = (prob(i,1)*Q(st1,i)) / tmpsum
         IF (MVUNI(1.0) .LE. pr0) THEN
            SScurr(1:h,1) = Imath(i,:)
            st1 = i
            GO TO 3
         ELSE
            i = i + 1
            st1 = i
         END IF
      END DO

c     return h'th row if have not already returned
      SScurr(1:h,1) = Imath(h,:)

 3    RETURN

      END SUBROUTINE bingen
