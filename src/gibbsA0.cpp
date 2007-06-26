#include "MSBVARcpp.h"
#include "A0_W.h" 

extern "C" SEXP gibbsA0(SEXP varobj, SEXP N1R, SEXP N2R, SEXP thinR,
			SEXP method, SEXP UTR) 
{
  int i, j, k, m, df, norm_method=INTEGER(method)[0], pctct=0;
  int N1=INTEGER(N1R)[0], N2=INTEGER(N2R)[0], thin=INTEGER(thinR)[0], switchct=0;   
  SEXP n0R, dfR, mode, identR; double *pA0;
  Matrix A0ml, A0gbs, ident;

//   Rprintf("N1:\t%d\nN2:\t%d\nthin:\t%d\nmethod:\t%d\n", N1, N2, thin, norm_method);

  // Initialize varobj$A0.mode/m/df/n0/ident
  PROTECT(mode=listElt(varobj,"A0.mode")); 
  m=getdims(mode)[0]; A0ml=R2Cmat(mode,m,m); 
  UNPROTECT(1);  
  
  PROTECT(dfR=coerceVector(listElt(varobj,"df"),INTSXP));
  df=INTEGER(dfR)[0]; 
  UNPROTECT(1);  
  
  PROTECT(identR=listElt(varobj,"ident")); 
  ident=R2Cmat(identR,m,m); double *pid=ident.Store(); 
  UNPROTECT(1); 

  PROTECT(n0R=listElt(varobj,"n0")); 
  int len=length(n0R), *pn0=INTEGER(n0R), *n0=(int*)Calloc(len,int); 
  for(i=0;i<len;i++) n0[i]=pn0[i]; 
  UNPROTECT(1); 

  // Normalize A0's mode/Set A0gbs
  A0ml=norm_svar(A0ml,A0ml,norm_method,&switchct); A0gbs=A0ml; 

  // Initialize storage/loop variables
  A0obj A0all(ident, N2); Wlist Wall(N2,m); Wobj W; UTobj UT(UTR); 

  // Find Matrix pointer locations of restrictions 
  int numx=m*m-A0all.numfp(), *restx=(int*)Calloc(numx,int), xct=0; 
  for(i=0;i<m*m;i++) if(!pid[i]) restx[xct++]=i;

  //////////////////
  // Burn-in Loop //
  //////////////////
    
  for(i=1;i<=N1;i++){
    // Draw posterior sample and immediately free unused allocated space 
    A0gbs=drawA0cpp(A0gbs, UT, df, n0, W); pA0=A0gbs.Store(); W.clear(); 

    // Clean up A0 to minimize roundoff error and normalize 
    for(j=0;j<numx;j++) pA0[restx[j]]=0; 
    A0gbs=norm_svar(A0gbs, A0ml, norm_method, &switchct); 

    if(!(i%(N1/10))) Rprintf("Gibbs Burn-in %d %% Complete\n",++pctct*10); 
  } 

  // Reset counters 
  switchct=0; pctct=0;

  ///////////////////
  // Sampling Loop //
  ///////////////////

  for(i=1;i<=N2;i++){
    // Take 'thin' number of draws 
    for(j=1; j<=thin; j++){
      // Draw posterior sample: Assign A0gbs explicitly, Wobj by reference 
      A0gbs=drawA0cpp(A0gbs, UT, df, n0, W); pA0=A0gbs.Store();

      // Clean up A0 to minimize roundoff error and normalize 
      for(k=0;k<numx;k++) pA0[restx[k]]=0; 
      A0gbs=norm_svar(A0gbs, A0ml, norm_method, &switchct);  
    }

    // Store A0/W in compact form, free Wobj memory for reuse  
    A0all.setA0(A0gbs,i); Wall.setWobj(W,i); W.clear(); 

    // Report iteration status and current A0's log determinant 
    if(!(i%(N2/10))){
      Rprintf("Gibbs Sampling %d %% Complete (%d draws)\n", ++pctct*10, i);
      LogAndSign ld = A0gbs.LogDeterminant(); 
      Rprintf("A0 log-det \t = %f \n", ld.Sign()*ld.LogValue()); 
    }
  }

  // Free allocated memory 
  Free(restx); Free(n0); 

  // Create output list w/ A0.posterior, W.posterior, ident, thin, N2 
  SEXP out, names;   
  PROTECT(out=allocVector(VECSXP,5)); 
  PROTECT(names=allocVector(STRSXP,5));
  SET_VECTOR_ELT(out, 0, A0all.toR()); 
  SET_STRING_ELT(names, 0, mkChar("A0.posterior")); 
  SET_VECTOR_ELT(out, 1, Wall.toR()); 
  SET_STRING_ELT(names, 1, mkChar("W.posterior")); 
  SET_VECTOR_ELT(out, 2, C2Rmat(ident)); 
  SET_STRING_ELT(names, 2, mkChar("ident")); 
  SET_VECTOR_ELT(out, 3, thinR); 
  SET_STRING_ELT(names, 3, mkChar("thin")); 
  SET_VECTOR_ELT(out, 4, N2R); 
  SET_STRING_ELT(names, 4, mkChar("N2")); 
  setAttrib(out, R_NamesSymbol, names); 
  UNPROTECT(2); 

  // Return output list
  return out;
}

// extern "C" SEXP gibbsA0(SEXP varobj, SEXP N1R, SEXP N2R, SEXP thinR,
// 			SEXP method, SEXP UTR) 
// {
//   SEXP modeR, dfR, n0R, draw;
//   int i, j, k, m, *dA0, df, *n0, *pn0R, norm_method=INTEGER(method)[0], pctct=0;
//   int N1=INTEGER(N1R)[0], N2=INTEGER(N2R)[0], thin=INTEGER(thinR)[0], switchct=0;   
//   double *pA0gbs, *pident; 

// //   Rprintf("N1:\t%d\nN2:\t%d\nthin:\t%d\nmethod:\t%d\n", N1, N2, thin, norm_method);

//   // Get A0.mode, ident, n0, df, UT, m from varobj/gibbs.setup

//   // varobj$A0.mode and varobj$m
//   PROTECT(modeR=listElt(varobj, "A0.mode")); 
//   dA0=getdims(modeR); m=dA0[0];      // varobj$m
//   Matrix A0mode=R2Cmat(modeR,m,m);   // varobj$A0.mode
//   UNPROTECT(1); 

//   // varobj$df
//   PROTECT(dfR=coerceVector(listElt(varobj,"df"),INTSXP));
//   df=INTEGER(dfR)[0]; 
//   UNPROTECT(1); 

//   // varobj$n0 (note: 'Free(n0);' must be called before return  
//   PROTECT(n0R=listElt(varobj,"n0")); 
//   n0=Calloc(length(n0R),int); pn0R=INTEGER(n0R); 
//   for(i=0;i<length(n0R);i++) n0[i]=pn0R[i]; 
//   UNPROTECT(1); 

//   // varobj$ident
//   Matrix ident=R2Cmat(listElt(varobj,"ident"),m,m); pident=ident.Store(); 
//   for(i=0;i<ident.Storage();i++) pident[i] = (!pident[i]) ? 0 : 1; 

//   // gibbs.setup$UT
//   UTobj UT(UTR); 

//   // Normalize A0
//   Matrix A0ml=norm_svar(A0mode,A0mode,norm_method,&switchct), A0gbs=A0ml;

//   //////////////////
//   // Burn-in Loop //
//   //////////////////

// //   SEXP Wtmp; int mtmp=m; PROTECT(Wtmp=allocVector(VECSXP,m)); 
// //   while(mtmp--) SET_VECTOR_ELT(Wtmp,m-mtmp,R_NilValue);
    
//   Rprintf("\n\nStarting Burn-in...\n");
//   Wobj W;
//   for(i=1;i<=N1;i++){
//     PROTECT(draw=drawA0cpp(A0gbs, UT, df, n0, W));  
//     A0gbs=R2Cmat(listElt(draw, "A0gbs"),m,m); pA0gbs=A0gbs.Store();
//     PROTECT(Wtmp=listElt(draw,"W")); W.setW(Wtmp); UNPROTECT(2);
    
//     // Clean up A0 to minimize roundoff error and normalize 
//     for(j=0;j<m*m;j++) pA0gbs[j]*=pident[j]; 
//     A0gbs=norm_svar(A0gbs, A0ml, norm_method, &switchct); 
//     if(!(i%(N1/10))) Rprintf("Gibbs Burn-in %d %% Complete\n",++pctct*10); 
//   }

//   ///////////////////
//   // Sampling Loop //
//   ///////////////////

//   // Initialize input and storage variables 
//   A0obj A0list(A0gbs, N2); switchct=0; pctct=0;
//   int *pWidx, *tmpidx, idxct=0, nidx=m*N2, lidx;
//   SEXP Widx; PROTECT(Widx=allocVector(INTSXP, nidx)); pWidx=INTEGER(Widx);
//   //  for(i=0;i<nidx;i++) pWidx[i]=0;
//   RowVector Wvals; Wvals=0.0; 

//   Rprintf("\n\n=============\nSampling...\n=============\n\n"); 
//   for(i=1;i<=N2;i++){
//     // Take thin number of draws 
//     for(j=1; j<=thin; j++){

//       // Draw posterior sample 
//       PROTECT(draw=drawA0cpp(A0gbs,UT,df,n0,&W));
//       A0gbs=R2Cmat(listElt(draw,"A0gbs"),m,m); pA0gbs=A0gbs.Store();
//       PROTECT(Wtmp=listElt(draw,"W")); W.setW(Wtmp); UNPROTECT(2); 

//       for(k=0;k<m*m;k++) pA0gbs[k]*=pident[k]; 
//       A0gbs=norm_svar(A0gbs, A0ml, norm_method, &switchct);  
//     }

//     // Store A0/W of thin'th draw in compact form 
//     A0list.setA0(A0gbs,i); Wvals|=W.getWvals(); 
//     tmpidx=W.getWindex(); lidx=0;
//     for(j=0;j<m;j++)
//       {idxct+=tmpidx[j]-lidx; lidx=tmpidx[j]; pWidx[(i-1)*m+j]=idxct;}

//     // Report iteration status and current A0's log determinant 
//     if(!(i%(N2/10))){
//       Rprintf("Gibbs Sampling %d %% Complete (%d draws)\n", ++pctct*10, i);
//       LogAndSign ld = A0gbs.LogDeterminant(); 
//       Rprintf("A0 log-det \t = %f \n", ld.Sign()*ld.LogValue()); 
// //       Matrix A0tmp=A0list.getA0(i); Rprintf("\nA0:"); printMatrix(A0tmp); 
// //       int tmpct=0; 
// //       PROTECT(Wtmp=W.getW());
// //       Matrix Welt; SEXP WeltR;
// //       int *dWelt; double *pWeltR, *pWelt=Welt.Store(); 
// //       for(j=0;j<m;j++){
// // 	PROTECT(WeltR=VECTOR_ELT(Wtmp,j)); pWeltR=REAL(WeltR); 
// // 	dWelt=getdims(WeltR); Welt=R2Cmat(WeltR,dWelt[0],dWelt[1]);
// // 	//Welt.ReSize(dWelt[0],dWelt[1]); Welt<<pWeltR; 
// // 	Rprintf("W[[%d]] (%dx%d)\n",j,dWelt[0],dWelt[1]); printMatrix(Welt); 
// // 	UNPROTECT(1); 
// //       }
// //       UNPROTECT(1);  
//     }
//   }
  
//   // Construct return objects/assign to R list for return
  
//   // Convert C types to SEXP 
//   SEXP Wpost, WR, mR, names, out; 
//   int Wlen=Wvals.Storage(); double *pWval=Wvals.Store(), *pWR; 
//   PROTECT(mR=allocVector(INTSXP,1)); INTEGER(mR)[0]=m; 
//   PROTECT(WR = allocVector(REALSXP, Wlen)); pWR=REAL(WR);  
//   for(i=0;i<Wlen;i++) *pWR++ = pWval[i]; 
  
//   // Create W.posterior object 
//   PROTECT(Wpost = allocVector(VECSXP, 3)); 
//   PROTECT(names = allocVector(STRSXP, 3));
//   SET_VECTOR_ELT(Wpost, 0, WR); 
//   SET_STRING_ELT(names, 0, mkChar("W")); 
//   SET_VECTOR_ELT(Wpost, 1, Widx); 
//   SET_STRING_ELT(names, 1, mkChar("W.index")); 
//   SET_VECTOR_ELT(Wpost, 2, mR); 
//   SET_STRING_ELT(names, 2, mkChar("m")); 
//   setAttrib(Wpost, R_NamesSymbol, names); 
//   UNPROTECT(1); 

//   // Create output list w/ A0.posterior, W.posterior, ident, thin, N2 
//   PROTECT(out=allocVector(VECSXP,5)); 
//   PROTECT(names=allocVector(STRSXP,5));
//   SET_VECTOR_ELT(out, 0, A0list.toR()); 
//   SET_STRING_ELT(names, 0, mkChar("A0.posterior")); 
//   SET_VECTOR_ELT(out, 1, Wpost); 
//   SET_STRING_ELT(names, 1, mkChar("W.posterior")); 
//   SET_VECTOR_ELT(out, 2, C2Rmat(ident)); 
//   SET_STRING_ELT(names, 2, mkChar("ident")); 
//   SET_VECTOR_ELT(out, 3, thinR); 
//   SET_STRING_ELT(names, 3, mkChar("thin")); 
//   SET_VECTOR_ELT(out, 4, N2R); 
//   SET_STRING_ELT(names, 4, mkChar("N2")); 
//   setAttrib(out, R_NamesSymbol, names); 
//   UNPROTECT(6); 

//   return out;
// }


// // void burnin(const int n, const Matrix& A0, const Matrix& UT, 
// // 	    const int df, const int *n0, const int A0idx, Matrix& W)
// //

// extern "C" SEXP burnin(SEXP varobj, SEXP A0gbs, SEXP A0idxR, SEXP UTR,
// 		       SEXP method, SEXP N1R)
// {
//   SEXP dfR, n0R, A0mlR, Wtmp, draw, A0R;
//   int i, j, N1=INTEGER(N1)[0], pctct=0, *n0, numrestrict=0, switchct=0;
  
//   // Get A0gbs(m,n) and assign A0ml
//   int *dA0=getdims(A0gbs), m=dA0[0], norm_method=INTEGER(method)[0]; 
//   Matrix A0=R2Cmat(A0gbs,m,m), A0ml=A0, A0str=A0; 
//   double *pA0=A0.Store(), *pA0str=A0str.Store(); 

//   for(i=0;i<A0.Storage();i++)
//     if(!pA0[i]){ numrestrict++; A0str(i+1)=0; } else { A0str(i+1)=1; }  

//   //  int *A0gbsd=INTEGER(getAttrib(A0gbs,R_DimSymbols)), m=A0gbsd[0], n=A0gbsd[1];

//   // Get int varobj$df 
//   int df=INTEGER(listElt(varobj,"df"))[0]; 
//   //  PROTECT(dfR=coerceVector(listElt(varobj,"df"), INTSXP)); 
//   //  int df=INTEGER(getAttrib(dfR,R_DimSymbol))[0]; UNPROTECT(1); 
  
//   // Get int array n0 
//   PROTECT(n0R=coerceVector(listElt(varobj,"n0"),INTSXP)); n0=INTEGER(n0R); 
  
//   // Create R list obj WR
//   //  PROTECT(WR=allocVector(VECSXP,m)); 

//   for(i=1;i<=N1;i++){
//     // Make R copy of A0 every iteration
//     A0gbs=C2Rmat(A0); 

//     // Call drawA0 fn and assign results of A0/W
//     PROTECT(draw=drawA0(A0gbs, UTR, n0R, WR)); 
//     A0=R2Cmat(listElt(draw, "A0gbs"),m,m); pA0=A0.Store();
//     PROTECT(Wtmp=listElt(draw,"W")); 
    
//     // Clean up A0 to minimize roundoff error and normalize 
//     for(j=0;j<m*m;j++) pA0[j]*=pA0str[j]; 
//     A0=svar_norm(A0, A0ml, norm_method, switchct); 

//     if(!(i%(N1/10))) Rprintf("Gibbs Burn-in %d Percent Complete\n",++pctct*10); 
//   }
  
//   // Make output object 
//   SEXP output, names; 
//   PROTECT(output=allocVector(VECSXP, 2)); PROTECT(names=allocVector(STRSXP, 2)); 
//   SET_VECTOR_ELT(output, 0, C2Rmat(A0)); SET_STRING_ELT(names, 0, mkChar("A0gbs")); 
//   SET_VECTOR_ELT(output, 1, Wtmp); SET_STRING_ELT(names, 1, mkChar("W")); 
//   setAttrib(output, R_NamesSymbol, names);

//   // Free memory and return 
//   UNPROTECT(6); return output; 
// }

// extern "C" SEXP gibbsdraw(SEXP varobj, SEXP N2R, SEXP thinR, 
// 			  SEXP A0gbs, SEXP A0ml, SEXP W, SEXP method,
// 			  SEXP UT)  
// {
//   // Initialize loop variables
//   int *dA0=getdims(A0gbs), m=dA0[0], switchct=0, pctct=0, i;
//   int N2=INTEGER(N2R)[0], thin=INTEGER(thinR)[0], norm_method=INTEGER(method)[0]; 

//   SEXP draw, n0, Wtmp, Widx; Wobj W(m); 
//   PROTECT(n0=coerceVector(listElt(varobj,"n0"),INTSXP)); 
//   PROTECT(Widx=allocVector(INTSXP, m*N2/thin)); 
//   int *pWidx=INTEGER(Widx), *tmpidx, idxct=0;
//   RowVector Wvals(N2*m*(m-1)*(m-1)); Wvals=0.0;

//   A0obj A0list(N2, R2Cmat(A0gbs,m,m));
//   Matrix A0(m,m), A0norm(m,m); A0=0.0; A0norm=0.0; 
//   Matrix A0str=A0list.structure(), A0mode=R2Cmat(A0ml,m,m); 

//   for(i=1;i<=N2*thin;i++){
//     PROTECT(draw=drawA0(A0gbs,UT,n0,W));
//     A0=R2Cmat(listElt(draw,"A0gbs"),m,m); 
//     double *pA0=A0.Store(), *pA0str=A0str.Store();
//     for(j=1;j<=m*m;j++) if(!pA0str[j]) pA0[j]=0;
//     A0=norm_svar(A0, A0mode, norm_method, &switchct);  
    
//     if(!(i%thin)){
//       PROTECT(Wtmp=listElt(draw,"W")); 
//       W.setW(Wtmp); tmpidx = W.getWindex(); 
//       for(j=0;j<m;j++)
// 	{ pWidx[(i/thin)*m+j]=idxct+tmpidx[j]; idxct+=tmpidx[j]; }
//       Wvals(pWidx[i/thin*m],pWidx[(i/thin)*m+j]) << W.getWvals(); 
//       A0list.setA0(i/thin,A0);
//     }
//     A0gbs=C2Rmat(A0); 
//     UNPROTECT(2); 

//     if(!(i%((N2*thin)/10))){
//       Rprintf("Gibbs Sampling %d%s Complete (%d draws)\n", ++pctct*10, '%', i/thin);
//       LogAndSign ld = A0gbs.LogDeterminant(); 
//       Rprintf("A0 log-det \t = %f \n", ld.Sign()*ld.LogValue()); 
//     }
//   }

//   // Construct return objects/assign to R list for return

//   SEXP Wpost; PROTECT(Wpost = allocVector(VECSXP, 3)); 
//   SEXP names; PROTECT(names = allocVector(STRSXP, 3));
//   SET_VECTOR_ELT(Wpost, 0, C2Rdouble(Wvals)); 
//   SET_STRING_ELT(names, 0, mkChar("W")); 
//   SET_VECTOR_ELT(Wpost, 1, Widx); 
//   SET_STRING_ELT(names, 1, mkChar("W.index")); 
//   SET_VECTOR_ELT(Wpost, 2, C2Rint(&m)); 
//   SET_STRING_ELT(names, 2, mkChar("m")); 
//   setAttrib(Wpost, R_NamesSymbol, names); 
//   UNPROTECT(1); 

//   SEXP A0ident; int dA0ident[]={m,m}; 
//   PROTECT(A0ident=coerceVector(listElt(varobj,"ident"),INTSXP));
//   setdims(A0ident, 2, dA0ident); 

//   SEXP out; PROTECT(out=allocVector(VECSXP,5)); 
//   PROTECT(names=allocVector(STRSXP,5));
//   SET_VECTOR_ELT(out, 0, A0list.toR());  
//   SET_STRING_ELT(names, 0, mkChar("A0.posterior"));  
//   SET_VECTOR_ELT(out, 1, Wpost); 
//   SET_STRING_ELT(names, 1, mkChar("W.posterior")); 
//   SET_VECTOR_ELT(out, 2, A0ident); 
//   SET_STRING_ELT(names, 2, mkChar("ident")); 
//   SET_VECTOR_ELT(out, 3, thinR); 
//   SET_STRING_ELT(names, 3, mkChar("thin")); 
//   SET_VECTOR_ELT(out, 4, N2R); 
//   SET_STRING_ELT(names, 4, mkChar("N2")); 
//   setAttrib(out, R_NamesSymbol, names); 

//   UNPROTECT(6); 
//   return out;
// }
