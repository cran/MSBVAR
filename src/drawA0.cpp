#include "MSBVARcpp.h"
#include "A0_W.h"

//extern "C"
SEXP drawA0(SEXP A0gbs, SEXP UT, SEXP dfR, SEXP n0, SEXP Wout)
{
  GetRNGstate(); //Ensures seed consistency

  int *dA0=INTEGER(getAttrib(A0gbs,R_DimSymbol));	//Pointer for A0gbs
  int i, j;						//Index ints
  int m=dA0[0], n=dA0[1];				//Values from A0
  int df=INTEGER(dfR)[0];     				//Extract df int from df R expression
  double tol=1E-12;					//Set tolerance level

  // Initialize input objects
  PROTECT(Wout = allocVector(VECSXP, m+8));		//Protect the return expression, Wout
  Matrix A0=R2Cmat(A0gbs,m,n);				//A0 is a matrix composed from A0gbs

  // Debug
  SEXP UTs;
  PROTECT(UTs=allocVector(VECSXP,m));

  // Main loop
  for(i=1;i<=m;i++)
  {
    //Variable assignment etc
    Matrix X=A0, Xt, R, dR, w(m,1);

    w=0.0;
    X.Column(i)=0.0;
    Xt=X.t();

    QRD qr(Xt, tol);

    R=qr.R();

    dR=diag(R);
    dR=absmat(dR);
    dR/=dR.Maximum();

    double *sing=dR.Store();
    double stol=pow(2.0,-16);
    int jidx=0, sidx=0;

    for(j=0;j<m;j++)
	if(sing[j]<=stol)
	{
		jidx=qr.pivot(j);
		sidx=j+1;
		break;
	}

    w.Row(jidx)=1.0;

    if(jidx>1 && jidx<=m)
    {
      Matrix jA = R.SubMatrix(1,jidx-1,1,jidx-1);
      Matrix jb = R.SubMatrix(1,jidx-1,sidx,sidx);
      QRD jy((-jA), tol);
      w.Rows(1,jidx-1) = jy.Solve(jb);
    }

    SEXP UTR;
    PROTECT(UTR = VECTOR_ELT(UT,i-1));

    int *UTdim = INTEGER(getAttrib(UTR,R_DimSymbol));
    Matrix UTtmp = R2Cmat(UTR,UTdim[0],UTdim[1]);

    UNPROTECT(1);

    Matrix W, W0=UTtmp.t()*w, W1=W0/W0.NormFrobenius();
    W=0.0;
    QRD tmp(W1,tol);
    W=tmp.Q();

    //Handle bad Q
    if(ISNAN(W.Sum()))
    {
      Rprintf("\nInvalid QR.Q\n");
      SEXP WoutR;
      PROTECT(WoutR=VECTOR_ELT(Wout,i-1));
      W = R2Cmat(WoutR,W.Nrows(),W.Ncols()).AsColumn();
      UNPROTECT(1);
      Rprintf("W:");
      printMatrix(W);
    }

    int *n0vec=INTEGER(coerceVector(n0,INTSXP));

    Real jstd=sqrt(1.0/df);
    Matrix gkb(n0vec[i-1],1), jr = jstd*rnorms(df+1);
    gkb=0.0;

    if(n0vec[i-1]>1)
      gkb.SubMatrix(2,n0vec[i-1],1,1) = jstd*rnorms(n0vec[i-1]-1);

    Real jrx=sqrt((jr.t()*jr).AsScalar()), rnd=unif_rand();

    if(rnd<0.5)
	gkb(1,1)=jrx;
    else
	gkb(1,1)=-jrx;

    A0.Column(i) = UTtmp*(W*gkb);

    SET_VECTOR_ELT(Wout, i-1, C2Rmat(W));
    SET_VECTOR_ELT(UTs, i-1, C2Rmat(UTtmp));
  }
  PutRNGstate();

  // R output stuff
  SEXP output, names, A0out;

  PROTECT(A0out = C2Rmat(A0));
  PROTECT(output = allocVector(VECSXP, 3));
  PROTECT(names = allocVector(STRSXP, 3));

  SET_VECTOR_ELT(output,0,A0out);
  SET_STRING_ELT(names,0,mkChar("A0gbs"));
  SET_VECTOR_ELT(output,1,Wout);
  SET_STRING_ELT(names,1,mkChar("W"));
  SET_VECTOR_ELT(output,2,UTs);
  SET_STRING_ELT(names,2,mkChar("UT"));

  setAttrib(output, R_NamesSymbol, names);

  UNPROTECT(5);

  return output;
}

// Function to draw posterior samples (A0/W) on the C side
//extern "C"
ReturnMatrix drawA0cpp(const Matrix &A0gbs, const class UTobj &UT,
				  const int df, const int *n0, class Wobj &W)
{
  Matrix A0, X, R, Wmat, W0, W1, UTi;
  int i, j, Wct, m, *Widx, jidx, sidx;
  double  tmp, maxR, jstd, tol, stol, *sing, *pR, *pW, *Wval;

  A0=A0gbs;
  m=A0.Nrows();
  Wct=0;
  tol=1E-12;
  stol=pow(2.0,-16);

  // Allocate space for to store indices/values of W's
  Widx=(int*)R_alloc(m,sizeof(int));
  Wval=(double*)R_alloc(m*m*m,sizeof(double));

  // Main loop
  GetRNGstate();//Seed consistency
  for(i=1;i<=m;i++)
  {
	X=A0;
	X.Column(i)=0.0;
	QRD qr(X.t(),tol);
	R=qr.R(); pR=R.Store();
	maxR=pR[0];

    	sing=(double*)R_alloc(m,sizeof(double));
	jidx=0;
	sidx=0;
	maxR=0;

    	for(j=0;j<m;j++)
	{
		sing[j]=fabs(pR[j*m+j]);
		if(sing[j]>maxR)
			maxR=sing[j];
	}
    	for(j=0;j<m;j++)
		if(sing[j]/maxR<=stol)
		{
			jidx=qr.pivot(j);
			sidx=j+1;
			break;
		}

    	ColumnVector w(m);
	w=0.0;
	w(jidx)=1.0;

   	if(jidx>1 && jidx<=m)
	{
      		Matrix jA = R.SubMatrix(1,jidx-1,1,jidx-1);
      		Matrix jb = R.SubMatrix(1,jidx-1,sidx,sidx);
      		QRD jy(-jA, tol);
		w.Rows(1,jidx-1) = jy.Solve(jb);
    	}

	UTi=UT.getUTelt(i);
	W0=UTi.t()*w;
	W1=W0/W0.NormFrobenius();

    	// Store W object values/indices
    	QRD qrtmp(W1,tol);
	Wmat=qrtmp.Q();
	Widx[i-1]=Wmat.Storage();
    	pW=C2F(Wmat);

	for(j=0;j<Widx[i-1];j++)
		Wval[Wct++]=pW[j];

    	ColumnVector gkb(n0[i-1]), jr;
	gkb=0.0;
	jstd=sqrt(1.0/df);
    	jr=jstd*rnorms(df+1);
	tmp=sqrt((jr.t()*jr).AsScalar());

    	if(unif_rand()<0.5)
		gkb(1)=tmp;
	else
		gkb(1)=-tmp;

    	if(n0[i-1]>1)
		gkb.Rows(2,n0[i-1])=jstd*rnorms(n0[i-1]-1);

    	A0.Column(i) = UTi*(Wmat*gkb);
  }
  PutRNGstate();

  // Assign W and return.
  W.setWflat(m, Wval, Widx);

  return A0.ForReturn();
}
