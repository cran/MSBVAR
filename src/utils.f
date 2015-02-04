      SUBROUTINE chol(n, matIn, matOut)

      INTEGER n, INFO
      DOUBLE PRECISION matIn(n,n), matOut(n,n)
      INTEGER i,j

c      EXTERNAL DPOTRF

      matOut = matIn
      CALL DPOTRF('U', n, matOut, n, INFO)

c     Now set lower triangular off-diagonal elements
c     to zero.
      IF (n > 1) THEN
         DO i = 2,n
            DO j = 1,(i-1)
               matOut(i,j) = 0.0
            END DO
         END DO
      END IF

      RETURN

      END SUBROUTINE chol

ccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine to create identity matrix
c     eyemat is an n x n identity matrix
      SUBROUTINE diag(n, eyemat)

      INTEGER j, n
      DOUBLE PRECISION eyemat (n,n)

c     Set each element to zero
      eyemat = 0.0

c     Change value of each element along diagonal to zero
      DO j = 1,n
         eyemat(j,j) = 1.0
      END DO

      END SUBROUTINE diag


ccccccccccccccccccccccccccccccccccccccccccccc
c     Return inverse of square matrices
      SUBROUTINE inv(n,matIn,matRet)

c      EXTERNAL DGETRF, DGETRI

      INTEGER n, LWORK, INFO
      DOUBLE PRECISION matRet(n,n), matIn(n,n)
      DOUBLE PRECISION WORK(n,n) ! not really sure why this is, but...
      INTEGER IPIV(n)

      matRet = matIn
      CALL DGETRF(n,n,matRet,n,IPIV,INFO)

      LWORK=-1
      CALL DGETRI(n,matRet,n,IPIV,WORK,LWORK,INFO)
      LWORK=WORK(1,1)
      CALL DGETRI(n,matRet,n,IPIV,WORK,LWORK,INFO)

      END SUBROUTINE inv


ccccccccccccccccccccccccccccccccccccccccccccc
c     Return absolute value of determinant
c     for square matrices
      SUBROUTINE detabs(n,matIn,det)

c      EXTERNAL DGETRF

      INTEGER n, INFO
      DOUBLE PRECISION matIn(n,n), matOut(n,n), det
      INTEGER IPIV(n)

      matOut=matIn
      CALL DGETRF(n,n,matOut,n,IPIV,INFO)

c     the unit diagonal elements of L are not stored, so
c     lower elements are L, while diagonal and upper are U
c     so, just multiply together diagonal elements of mat
c     to obtain determinant

      det = 1.0
      DO i = 1,n
         DO j = 1,n
            IF (i==j) det = det*matOut(i,j)
         END DO
      END DO

      det = ABS(det)

      RETURN

      END SUBROUTINE detabs


ccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine for matrix matrix multiplication
      SUBROUTINE MatrixMultiply(matA, nrowA, ncolA,
     &                          matB, nrowB, ncolB,
     &                          matC, tmpint)

c      EXTERNAL DGEMM

      INTEGER nrowA, ncolA, nrowB, ncolB, tmpint
      DOUBLE PRECISION matA(nrowA,ncolA), matB(nrowB,ncolB)
      DOUBLE PRECISION matC(nrowA,ncolB)

c     Note the d0's in the call.
      CALL DGEMM('N','N',nrowA,ncolB,ncolA,1.0d0,
     &           matA,nrowA,matB,nrowB,0.0d0,matC,nrowA)

      RETURN

      END SUBROUTINE MatrixMultiply



ccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine to transpose matrix
      SUBROUTINE MatrixTranspose(matIn, nrow, ncol, matOut)

      INTEGER nrow, ncol, irow, icol
      DOUBLE PRECISION matIn(nrow,ncol), matOut(ncol,nrow)

c     Change value of each element along diagonal to zero
      DO irow = 1,nrow
         DO icol = 1,ncol
            matOut(icol,irow) = matIn(irow,icol)
         END DO
      END DO

      RETURN

      END SUBROUTINE MatrixTranspose

ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Uniform (0,1) random number generator
c
c     use R's random number generator directly
c     the way `Writing R extentions' advertises.
c

      SUBROUTINE MVUNI(y)
      DOUBLE PRECISION unifrnd, y

      CALL rndstart()
      y = unifrnd()
      CALL rndend()
      RETURN

      END SUBROUTINE MVUNI

