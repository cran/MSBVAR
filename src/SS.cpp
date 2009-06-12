#include "SS.h"
#include "MSBVARcpp.h"
#include "Rdefines.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
////                                                  ////
////  Functions for class SSobj: an ANSI C compatible ////
////  implementation of state-spaces object storage   ////
////  and processing.                                 ////
////                                                  ////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// Constructors that populate integer representation array _int
SSobj::SSobj(UINT **ss, int **tmap, int h, int T)
{ SSbuild(ss,tmap,h,T); }

// Constructor to initialize SSobj and allocate memory
SSobj::SSobj(int h, int T){ SSinit(h,T); }

// Initialize and allocate SSobj
void SSobj::SSinit(int h, int T)
{
  // Initialize number of regimes, number of time periods, number of
  // ints needed, object size, integer rep of SS, and the trans map

  _h=h; _T=T; _size=_h*_T; _st0=0; _nint=ceil(_T/NBITSF); 
  _int = (UINT **)Calloc(_h, UINT*);
  _tmap = (int **)Calloc(_h, int*); 

  // Initialize arrays to zero
  int i, j;
  for(i=0;i<_h;i++){ 
    _int[i] = (UINT *)Calloc(_nint+1, UINT); 
    for(j=0;j<_nint;j++) _int[i][j] = 0;
    _tmap[i] = (int *)Calloc(_h, int);  
    for(j=0;j<_h;j++) _tmap[i][j] = 0;
  } 

}

void SSobj::SSbuild(UINT **ss, int **tmap, int h, int T)
{
  int i, j;  _h=h; _T=T; _size=_h*_T; _st0=0; _nint=ceil(_T/NBITSF); 
  _int = (UINT**)Calloc(_h, UINT*); _tmap = (int**)Calloc(_h, int*); 

  for(i=0;i<_h;i++){
    _int[i]=(UINT*)Calloc(_T, UINT);  _tmap[i]=(int*)Calloc(_h, int); 
    for(j=0;j<_T;j++) _int[i][j]=ss[i][j];
    for(j=0;j<_h;j++) _tmap[i][j]=tmap[j][i];
  }
}

// Function to generate draws of the SS

void SSobj::SSgenerate(Matrix &fp, Matrix &Q)
{
  GetRNGstate();
  // Draw regime for _T -- this one has to be done separately, since
  // it does not depend on Q.
  Matrix tmp=fp.Row(_T); 
  int i=1, h=Q.Nrows(), state=-1;
  for(i=1; i<=h-1; i++) 
    {
      Real pp = tmp.Column(i).AsScalar()/tmp.Columns(i,h).Sum();
      double u = runif(0,1);
      if(pp<=u) state=i-1; break;    // sampled value
    }
  if(state>-1) {setregime(_T, state); }else{ setregime(_T, h-1); }

  //  setregime(_T, state);

  // Now back-sample from T-1:1
  for(int t=_T; t>0; t--){ 
    // Get the filter probability and smooth it and sample it
    tmp=fp.Row(t); 
    setregime(t-1, bingen(tmp, Q, getregime(t)-1)); 
  }
  PutRNGstate();
}

// Set the bit for the regime
void SSobj::setregime(int t, int st)
{
  _int[st][t/NBITS] |= (1<<(NBITS-t%NBITS));   // Set regime at time t 
  _tmap[_st0][st]++;                  // Update transition map 
  _st0=st;                            // Update previous regime 
}

// Return the regime at time T as an integer
int SSobj::getregime(int t)
{
  int i, intidx=t/NBITS; UINT mask=(1<<(NBITS-t%NBITS));
  for(i=0;i<_h-1;i++) 
    if(_int[i][intidx] & mask){ return i+1; }
  return _h;
}


////////////////////////////////////////////////////////////////////////
// SSobj::normalizeSS() 
// Function to normalize objects generated in the SS class.
// Applied to the objects for the newly generated state-space, this
// function ranks the states by prevalence.  So the state that most of
// the observations are in is ranked first, the second most common
// state is second, etc.  The function re-ranks the states and re-maps
// the transition matrix in the SSobj.
//
////////////////////////////////////////////////////////////////////////
void SSobj::normalizeSS()
{

  // alloc a variable to hold the counts of the observations in each
  // state
  int i, j, *ct; ct=(int*)Calloc(_h, int);   

  // Count the time in the states from the _tmap
  for(i=0;i<_h;i++) for(j=0;j<_h;j++) ct[i]+=_tmap[j][i]; 

  // allocate / declare variable to hold the temporary objects.
  bool exchanges; int h=_h;
  UINT *SStmp; SStmp=(UINT*)Calloc(_nint, UINT);
  int **maptmp = (int**)Calloc(_h,int*);
  int *one2h; one2h=(int*)Calloc(_h,int);

  // Create 0->h-1 index and allocate memory for a copy of _tmap
  for(i=0;i<_h;i++){ 
    one2h[i]=i;  
    maptmp[i]=(int*)Calloc(_h,int);
  }

  //  below is a bubble sort that is then fed into a swap of the
  //  elements in the SS matrix.
  
  do {
    h--;               // make loop smaller each time.
    exchanges = false; // assume this is last pass over array
    for (int i=0; i<h; i++) {
      if (ct[i] < ct[i+1]) {

	// Swap the values in count array
        int cttmp = ct[i];  ct[i] = ct[i+1];  ct[i+1] = cttmp; 
	int tmp = one2h[i]; one2h[i] = one2h[i+1]; one2h[i+1] = tmp;

	// Swap the pointers to the columns of the state space 
	memcpy(SStmp, _int[i], _nint); 
	memcpy(_int[i], _int[i+1], _nint);
	memcpy(_int[i+1], SStmp, _nint);

        exchanges = true;  // after an exchange, must look again 
      }
    }
  } while (exchanges);
  
  // No re-map the transition matrix in _tmap into a temporary
  // matrix.  This is so we do not clobber the original in memory.
  for(i=0;i<_h;i++){
    for(j=0;j<_h;j++){
      maptmp[i][j] = _tmap[one2h[i]][one2h[j]];
    }
  }

  // Point _tmap to the memory address of maptmp -- this lets us
  // control the memory used in _tmap.

  for(i=0;i<_h;i++) for(j=0;j<_h;j++) _tmap[i][j]=maptmp[i][j];  

  // Clean up after yourself
  Free(maptmp);  Free(one2h);  Free(SStmp);  Free(ct);
}



// // Basic population of one SSobj from one slimmed state space.
// void SSobj::SSpopulate(UINT **ss, int **tmap, int h, int T)
// {
//   int i, j; _h=h; _T=T; _size=_h*_T; _st0=0; _nint=ceil(_T/NBITSF); 
//   _int = (UINT**)Calloc(_h, UINT*); _tmap = (int**)Calloc(_h, int*); 

//   for(i=0;i<_h;i++){
//     _int[i]=(UINT*)Calloc(_T, UINT);  _tmap[i]=(int*)Calloc(_h, int); 
//     for(j=0;j<_T;j++) _int[i][j]=ss[i][j];
//     for(j=0;j<_h;j++) _tmap[i][j]=tmap[j][i];
//   }
// }

//void SSobj::SSclean(){ Free(_int); Free(_tmap); }

// int* SSobj::SSinfo()
// { 
//   int *info; info=(int*)Calloc(4,int); 
//   info[0]=_T; info[1]=_h; info[2]=_size; info[3]=_nint;
//   return info;
// }

// int SSobj::SS_h(){return _h;}
// int SSobj::SS_T(){return _T;}
// int SSobj::SS_size(){return _size;}
// int SSobj::SS_nint(){return _nint;}

// // Utility functions that encode a 0/1 sequence SS[_size] of into a
// // vector of integer representations _int[_size/NBITS]
// void SSobj::encodeSS(SEXP SS)
// {
//   // Determine number of states to encode and number of integers
//   // necessary to store the states. Data type 'UINT' from stdint.h
//   // is used to ensure a NBITS-bit allocation on all platforms, without
//   // losing bit twiddling capabilities. Populate private data members
//   // and allocate/populate 2D array _int[h][ceil(T/NBITSF)].

//   int i, *dSS = getdims(SS); _T=dSS[0]; _h=dSS[1]; _size = _T*_h; 
//   _int = (UINT**)Calloc(_h, UINT*); _nint=ceil(_T/NBITSF); 
//   for(i=0;i<_h;i++){ _int[i] = (UINT*)Calloc(_nint, UINT); } 
  
//   // Loop over every regime in the state-space in blocks of NBITS. For
//   // each block, tabulate the integer representation of the sequence
//   // of 0,1 as if it were in base 2 (value+=2^index(sequence==1)).

//   UINT *ss=(UINT*)REAL(SS), j, k=0; 
//   for(i=0;i<_h;i++)
//     for(j=0;j<_T;j++){
//       if(j%NBITS==0) k++;
//       if(ss[i*_T+j]){_int[i][k]|=(1<<(j-NBITS*k));}
//     }
// }

// SEXP SSobj::getregimeR(int T)
// {
//   SEXP out; 
//   PROTECT(out=allocVector(INTSXP,1));
//   INTEGER(out)[0]=getregime(T); 
//   UNPROTECT(1); 
//   return out;
// }

SEXP SSobj::getSSmapR()
{
  SEXP out; PROTECT(out=allocMatrix(INTSXP,_h,_h));
  int i, j, *pout=INTEGER(out);
  for(i=0;i<_h;i++) for(j=0;j<_h;j++) pout[j+i*_h]=_tmap[j][i]; 
  UNPROTECT(1); return out;
}

// int** SSobj::getSSmap()
// {
//   int i, **tmap; tmap=(int**)Calloc(_h, int*);
//   for(i=0;i<_h;i++) memcpy(tmap[i],_tmap[i],_h);
//   return tmap;
// //   int i, j; tmap=(int**)Calloc(_h,int*);
// //   for(i=0;i<_h;i++){ 
// //     tmap[i]=(int*)Calloc(_h,int);
// //     for(j=0;j<_h;j++) tmap[i][j]=_tmap[i][j];
// //   }
// }

// Return integer representation array _int[1:_size/NBITS]  This gives us
// back the most compact representation of the state-space.
// SEXP SSobj::SSslim()
// { 
//   SEXP SS, tmap, out, names; 
//   PROTECT(SS=allocMatrix(REALSXP, _nint, _h)); 
//   PROTECT(tmap=getSSmapR()); 

//   int i, j; double *ss=REAL(SS);
//   for(i=0;i<_h;i++) for(j=0;j<_nint;j++) ss[i*_nint+j] = _int[i][j]; 
//   //  Rprintf("In SSobj::SSslim()...  ss(%dx%d) assigned\n", _nint, _h);

//   PROTECT(out=allocVector(VECSXP, 2)); 
//   PROTECT(names=allocVector(STRSXP, 2)); 

//   SET_VECTOR_ELT(out, 0, SS); 
//   SET_STRING_ELT(names, 0, mkChar("intSS"));

//   SET_VECTOR_ELT(out, 1, tmap); 
//   SET_STRING_ELT(names, 1, mkChar("transitions"));

//   setAttrib(out, R_NamesSymbol, names); 
//   UNPROTECT(4); 
//   return out; 
// }

// SEXP SSobj::SSfull()
// { 
//   int i; Matrix ss(_h,_T); ss=0; 
//   for(i=0;i<_T;i++) ss(i+1,getregime(i))=1; 

//   Rprintf("Populated 0-1!!\n");

//   SEXP out; PROTECT(out=allocVector(LISTSXP,2)); 
//   SEXP names; PROTECT(names=allocVector(STRSXP,2));

//   SET_VECTOR_ELT(out, 0, C2Rmat(ss)); 
//   SET_STRING_ELT(names, 0, mkChar("SSfull"));
  
//   SET_VECTOR_ELT(out, 1, getSSmapR()); 
//   SET_STRING_ELT(names, 1, mkChar("transitions"));
  
//   setAttrib(out, R_NamesSymbol, names); 
//   UNPROTECT(2); return out;
// }


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
////                                                  ////
////  Functions for class SSlist: a class-based list  ////
////  of SSobj (state-space objects)                  //// 
////                                                  ////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// SSlist::SSlist(int n){ _ndraws=n; SSlistinit(_ndraws); }

// SSlist::SSlist(SEXP SSl)
// { SEXP SSlR; PROTECT(SSlR=SSl); SSlistR(SSlR); UNPROTECT(1); }

// void SSlist::SSlistinit(int n){ 
//   //  int i, h, T;  
//   _ndraws=n; 
//   _SS=(SSobj*)Calloc(_ndraws,SSobj); 
// }

// void SSlist::SSlistR(SEXP SSsample)
// {
//   int i, j, k, nint, **tmap, *dtmap, *ptmap; 
//   UINT **ssint, *pSS; SEXP SSiR, SSintR, tmapR;  

//   // Get number of draws and initialize SSlist obj
//   _ndraws = length(SSsample); SSlistinit(_ndraws); 

//   // Get dimension information (h, nint, T)
//   PROTECT(SSiR = VECTOR_ELT(SSsample,0));
//   PROTECT(tmapR = VECTOR_ELT(SSiR,1)); 
//   dtmap = getdims(tmapR);  
//   _h = dtmap[0];  
//   _T = INTEGER(VECTOR_ELT(SSiR,2))[0];  
//   nint = ceil(_T/NBITSF);
//   UNPROTECT(2);
  
//   // Pre-allocate memory for temporary objects
//   SSobj SStmp(_h,_T);  
//   ssint=(UINT**)Calloc(_h, UINT*);  
//   tmap=(int**)Calloc(_h, int*);

//   for(i=0;i<_h;i++){ 
//     ssint[i]=(UINT*)Calloc(nint, UINT); 
//     tmap[i]=(int*)Calloc(_h,int); 
//   }

//   for(i=0; i<_ndraws; i++){
//     // Get i'th SSobj from SSsample list
//     PROTECT(SSiR = VECTOR_ELT(SSsample,i));

//     // Get SSint/tmap from SSiR (SS[[i]])
//     PROTECT(SSintR = coerceVector(VECTOR_ELT(SSiR,0), INTSXP));  
//     PROTECT(tmapR = coerceVector(VECTOR_ELT(SSiR,1), INTSXP));  
//     pSS = (UINT*)INTEGER(SSintR);  ptmap = INTEGER(tmapR);
// //     Rprintf("Length of SSintR = %d\n", length(SSintR));
//     for(j=0; j<_h; j++){
//       for(k=0; k<nint; k++) { ssint[j][k]=pSS[j*nint+k]; 
// //       Rprintf("ssint[%d][%d] = %d\n", i, j, ssint[j][k]);
// //       Rprintf("pSS[%d] = %d\n", j*nint+k, pSS[j*nint+k]);
//       }
//       for(k=0; k<_h; k++) tmap[j][k]=ptmap[j*_h+k];
//     }

// //     SStmp.SSbuild(ssint, tmap, _h, _T);  
//     _SS[i].SSpopulate(ssint, tmap, _h, _T); 

//     UNPROTECT(3); 
//   }

//   // Pack out what we created
//   Free(ssint);  Free(tmap);  
//   Free(dtmap);
//   Free(pSS);
// }

// void SSlist::setSSobj(const class SSobj &SS, int i){ _SS[i]=SS; }


// int** SSlist::SSsum()
// {
//   int i, j, k, **sums;  sums=(int**)Calloc(_h, int*);
//   for(i=0;i<_h*_T;i++){ 
//     sums[i]=(int*)Calloc(_T,int);
//     for(j=0;j<_T;j++) sums[i][j]=0;
//   }

//   //  Rprintf("sums[%d][%d] allocated\n",_h,_T);
//   for(i=0;i<_ndraws;i++) {
//       for(k=0;k<_T;k++) { 
// // 	tmp=_SS[i].getregime(k)-1;
// // 	Rprintf("_SS[%d] index = %d\n", i, tmp);
//  	sums[_SS[i].getregime(k)-1][k]+=1.0;
// // 	Rprintf("[i=%d][k=%d]\n", i, k);
//       }
//   }

// //   Rprintf("sums[%d][%d] populated\n",_h,_T);

//   return sums;
// }

// double** SSlist::SSmean(const int method)
// {
//   // Declare variables, get num regimes (h), allocate memory for
//   // cumulative sums and mean scores and initialize to zero
  
//   int i, j, **sums;  //double nh, *out; nh=(double)_ndraws*(double)_h;
//   double **out; out=(double**)Calloc(_h,double*);
//   for(i=0;i<_h;i++) {
//     out[i]=(double*)Calloc(_T,double); 
//     for(j=0;j<_T;j++) out[i][j]=0.0;
//   }

//   // Choose method based on input
//   switch(method){

//   default:
//   case 1: {
//     // Compute the mean probabilities from the sums of the states over
//     // the state space draws  
//     Rprintf("Default case selected\n");
//     sums = SSsum(); 
//     Rprintf("sums = SSsum() returned\n");
//     for(i=0;i<_h;i++) for(j=0;j<_T;j++) out[i][j]=sums[i][j]/(double)_ndraws;
//     Rprintf("out[%d][%d]/=%d\n",i,j,_ndraws);
//     break;
//   }

//   case 2: {
	
//     // Allocate out[h], the long-run mean probs of being in each
//     // state, and csums[h], the cumulative sum of counts within each
//     // state, and initialize them both to zero

// //     out=(double*)Calloc(_h, double);  csums=(int*)Calloc(_h, int);  
// //     for(i=0;i<_h;i++){ csums[i]=0; out[i]=0.0; }
    
//     // Compute cumulative sums for each state over all _ndraws draws by
//     // finding the column sums of the state transition map
    
// //     for(i=0;i<_ndraws-1;i++)
// //       for(j=0;j<_h;j++) 
// // 	for(k=0;k<_h;k++) 
// // 	  csums[j]+=_SS[i]._tmap[j][k]; 
	  
//     // Compute mean scores by dividing the cumulative sums for each
//     // state by the number of draws (csums[i]/_ndraws) and free memory 

// //     for(i=0;i<_h;i++) out[i]=(double)csums[i]/_ndraws;  
// //     Free(csums); 

//     break;
//   }}

//   return out;
// }

// double** SSlist::SSvar(const int method)
// {
//   // Declare variables, assign num regimes (_h), num periods (_T), and
//   // expected values (SSmean())

//   //  int i, T=_SS[0].SS_T(), h=_SS[0].SS_h(); 
//   double **var=SSmean(1);  

//   // Find the probability of being in a state (p=mean[i]/T), compute
//   // the variance of each state (var[i]=p*(1-p)), and return the
//   // computed variances (var)

// //   for(i=0;i<h;i++){ 
// //     var[i]/=T;  
// //     var[i]*=(1-var[i]); 
// //   } 

//   return var;
// }
