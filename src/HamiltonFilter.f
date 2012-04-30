      SUBROUTINE HamiltonFilter(bigt,m,p,h,e,sig2,Q,f,filtprSt)

c     bigt = total number of observations in time series
c     m    = number of dependent variables.  i.e., VAR if m > 1
c     p    = number of lags: AR(p) / VAR(p)
c     h    = number of regimes
      INTEGER bigt,m,p,h
      INTEGER its

c     regime is in 3rd dimension
c     sig2 is scalar if m=1, or var-cov matrix if m>1
      DOUBLE PRECISION  sig2(m,m,h),Q(h,h), transQ(h,h)

c     array to store residuals
      DOUBLE PRECISION  e(bigt-p,m,h), tmpfit(1,m), transtmpfit(m,1)

c     likelihoods
      DOUBLE PRECISION  ylik(bigt-p,1,h), ypwlik(bigt-p,h,h)
      DOUBLE PRECISION  detsig2, invsig2(m,m), tmpsig2(m,m)
      DOUBLE PRECISION  matmulone(1,m), matmultwo(1,1)
      DOUBLE PRECISION  matmulthree(h,h)

      DOUBLE PRECISION  ssAmat(h+1,h), transssAmat(h,h+1)

c     identity matrices used below
      DOUBLE PRECISION  Imath(h,h), Invmatss(h,h), ssEvec(h+1,1)
      DOUBLE PRECISION  rho(h,1), rhomatmulone(h,h+1)

      DOUBLE PRECISION  pSt1_t1(h,1), pytSt1St_t1_itert(h,h)
      DOUBLE PRECISION  filtprSt1St(bigt-p,h,h), filt_llfval(bigt-p)
      DOUBLE PRECISION  filtprSt(bigt-p,h)
      DOUBLE PRECISION  f, tmpdist(1,1), tmpsum

      DOUBLE PRECISION pi
      PARAMETER(pi = 3.1415926535897932384626433)

c     Loop over regimes (1 to h) to calculate univariate Normal
c     or multivariate Normal density given parameter values
      DO iterh = 1,h
         tmpsig2 = sig2(:,:,iterh)
         detsig2 = tmpsig2(1,1)
         invsig2(1,1) = 1/tmpsig2(1,1)

c     calculate determinant of Sigma if m > 1
         IF (m > 1) THEN 
            CALL detabs(m,tmpsig2(:,:),detsig2)
            CALL inv(m,tmpsig2(:,:),invsig2)
         END IF
         
         DO itert = 1,(bigt-p)
c           Recall dimensions of tmpfit(1,m)
            tmpfit(1,:) = e(itert,:,iterh)
            CALL MatrixTranspose(tmpfit, 1, m, transtmpfit)
            CALL MatrixMultiply(tmpfit, 1, m, invsig2, m, m, 
     &           matmulone, 1)
            CALL MatrixMultiply(matmulone, 1, m, transtmpfit, m, 1, 
     &                          matmultwo, 1)
            ylik(itert:itert,:,iterh) = EXP(-(REAL(m)/2.0)*LOG(2.0*pi) 
     &           - (1.0/2.0)*LOG(detsig2) - (1.0/2.0)*matmultwo)
         END DO
      END DO

      DO iterh1 = 1,h
         DO iterh2 = 1,h
            ypwlik(:,iterh1,iterh2) = Q(iterh1,iterh2)
     &           * ylik(:,1,iterh2)
         END DO
      END DO

      ssAmat = 1.0
      ssEvec = 0.0; ssEvec(h+1,1) = 1.0
      CALL diag(h,Imath)

      CALL MatrixTranspose(Q, h, h, transQ)
      ssAmat(1:h,:) = Imath-transQ

c     Recall ssAmat(h+1,h)
      CALL MatrixTranspose(ssAmat, h+1, h, transssAmat)

      CALL MatrixMultiply(transssAmat, h, h+1, ssAmat, h+1, h, 
     &                    matmulthree, 1)

      CALL inv(h, matmulthree, Invmatss)

      CALL MatrixMultiply(Invmatss, h, h, transssAmat, h, h+1, 
     &                    rhomatmulone, 1)

      CALL MatrixMultiply(rhomatmulone, h, h+1, ssEvec, h+1, 1, 
     &                    rho, 1)

      pSt1_t1 = rho

      its = 1
      f = 0.0

      DO WHILE (its .LE. (bigt-p))

         tmpsum = 0.0
         DO i = 1,h
            DO j = 1,h
               pytSt1St_t1_itert(i,j) = pSt1_t1(i,1)*ypwlik(its,i,j)
               tmpsum = tmpsum + pytSt1St_t1_itert(i,j)
            END DO
         END DO

         filt_llfval(its) = tmpsum
         f = f + LOG(filt_llfval(its))

         pytSt1St_t1_itert = pytSt1St_t1_itert/filt_llfval(its)

         filtprSt1St(its,:,:) = pytSt1St_t1_itert
         DO i = 1,h
            tmpsum = 0.0
            DO j = 1,h
               tmpsum = tmpsum + pytSt1St_t1_itert(j,i)
            END DO
            pSt1_t1(i,1) = tmpsum
         END DO

         its = its + 1

      END DO

c     integrate over St-1
      filtprSt = 0.0
      DO i = 1,h
         filtprSt = filtprSt + filtprSt1St(:,i,:)
      END DO
      
      RETURN

      END SUBROUTINE HamiltonFilter

