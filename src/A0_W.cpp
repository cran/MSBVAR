#include <R.h>
#include "A0_W.h"

///////////////////////////////////////////////////////////////////////
// A0obj constructor takes any A0 representation (as a 0/1           //
// ident matrix or an example filled with 0/[double value]) and a    //
// scalar indicating the number of gibbs sampler draws to be         //
// saved. Constructor allocates all necessary memory and initializes //
// private data members.                                             //
///////////////////////////////////////////////////////////////////////

// Empty Constructor initialization
A0obj::A0obj(){ _m=0; _nident=0; _ndraw=0; _cdraw=0; }

// Free allocated memory on destruction of object
A0obj::~A0obj(){Free(_str); Free(_val); Free(_fpidx); Free(_xidx);}

A0obj::A0obj(SEXP A0)
{
  int i, ct=0, *str; 
  SEXP mR, strR, valsR; 
  double *vals; 

//   Rprintf("--- Setting _m...\n"); 

  PROTECT(mR=coerceVector(listElt(A0,"m"),INTSXP));
  _m=INTEGER(mR)[0]; 
  UNPROTECT(1); 
  
//   Rprintf("--- Setting _fpidx, strR...\n"); 

  PROTECT(strR=coerceVector(listElt(A0,"struct"),INTSXP)); 
  _nident=length(strR); str=INTEGER(strR);
  _fpidx=(int*)Calloc(_nident,int); 
  for(i=0;i<_nident;i++) _fpidx[i]=str[i]-1;
  UNPROTECT(1); 

//   Rprintf("--- Setting _str...\n"); 

  _str=(int*)Calloc(_m*_m,int); 
  _xidx=(int*)Calloc(_m*_m-_nident,int);
  for(i=0;i<_m*_m;i++) 
    if(_fpidx[ct]==i){_str[i]=1; ct++;}else{_str[i]=0; _xidx[i-ct]=i;}

//   Rprintf("--- Setting valsR...\n"); 

  PROTECT(valsR=coerceVector(listElt(A0,"A0"),REALSXP)); 
  _ndraw=length(valsR)/_nident; _cdraw=0; vals=REAL(valsR);

//   Rprintf("--- Setting _val...\n"); 

  _val=(double*)Calloc(_ndraw*_nident,double); 
  for(i=0;i<_ndraw*_nident;i++) _val[i]=vals[i]; 
  UNPROTECT(1); 
}

A0obj::A0obj(const Matrix &ident, const int ndraw)
{
  _m=ident.Nrows(); _ndraw=ndraw; _cdraw=0; _str=(int*)Calloc(_m*_m,int); 

  int i, ct=0; double *str=C2F(ident); _nident=0; 
  for(i=0;i<_m*_m;i++) if(str[i]){_str[i]=1; _nident++;}else{_str[i]=0;} 
  Free(str); 

  // Allocate memory to store values/indices
  _val=(double*)Calloc(_nident*_ndraw,double); 
  _fpidx=(int*)Calloc(_nident,int); 
  _xidx=(int*)Calloc(_m*_m-_nident,int); 
  for(i=0;i<_m*_m;i++) if(_str[i]){_fpidx[ct++]=i;}else{_xidx[i-ct]=i;}
}

///////////////////////////////////////
// Store A0 matrix at supplied index //
///////////////////////////////////////

void A0obj::setA0(const Matrix &mat, const int index)
{
  int i, st=(index-1)*_nident; double *A0=C2F(mat); 
  for(i=0;i<_nident;i++) _val[st++]=A0[_fpidx[i]]; Free(A0);
}

int A0obj::fpidx(const int i) const { return _fpidx[i]; }
int A0obj::xidx(const int i) const { return _xidx[i]; }
int A0obj::numfp() const { return _nident; }

//////////////////////////////////////
// Store A0 matrix at current index //
//////////////////////////////////////

void A0obj::setA0c(const Matrix &mat) 
{
  double *A0=mat.Store(); int i, st=_cdraw*_nident;
  for(i=0;i<_m*_m;i++) if(_str[i]!=0) _val[st++]=A0[i]; 
}

/////////////////////////////////////////////////////////////////////
// Returns A0's structure matrix as a 0/1 identification matrix 
/////////////////////////////////////////////////////////////////////

ReturnMatrix A0obj::structure() const
{
  Matrix out(_m,_m); out=0; double *pout=out.Store();
  for(int i=0; i<_m*_m; i++) pout[i]=_str[i]; 
  out.Release(); return out.ForReturn(); 
}

/////////////////////////////////////////
// Returns a COPY of the values stored //
/////////////////////////////////////////

double* A0obj::getA0vals() const
{
  double* out; out=(double*)Calloc(_nident*_ndraw,double);  
  for(int i=0;i<_nident*_ndraw;i++) out[i]=_val[i]; return out; 
}

////////////////////////////////////////////////////////////
// Returns the A0 matrix at the supplied index [1:_ndraw] //
////////////////////////////////////////////////////////////

ReturnMatrix A0obj::getA0(const int idx) const
{
  int i, st=idx*_nident; double *pmat=(double*)Calloc(_m*_m,double); 
  for(i=0;i<_m*_m;i++) pmat[i]=(_str[i])?_val[st++]:0; 
  Matrix mat=F2C(pmat,_m,_m); Free(pmat); 
  mat.Release(); return mat.ForReturn(); 
}

///////////////////////////////////////////////
// Returns the A0 matrix at the current draw //
///////////////////////////////////////////////

ReturnMatrix A0obj::getA0c() const
{
  int i, st=_cdraw*_nident; 
  Matrix mat(_m,_m); double *pmat=mat.Store(); 
  for(i=0;i<_m*_m;i++) pmat[i]=(_str[i]!=0)?_val[st++]:0; 
  mat.Release(); return mat.ForReturn(); 
}

//////////////////////////////////////////////////////
// Returns SEXP list containing A0flat, A0struct, m //
//////////////////////////////////////////////////////

SEXP A0obj::toR() const 
{
  int i, A0len = _nident*_ndraw; 
  SEXP out; PROTECT(out=allocVector(VECSXP,3)); 
  SEXP names; PROTECT(names=allocVector(STRSXP,3));
  SEXP A0; PROTECT(A0=allocVector(REALSXP, A0len));
  SEXP m; PROTECT(m=allocVector(INTSXP,1)); INTEGER(m)[0]=_m; 
  SEXP str; PROTECT(str=allocVector(INTSXP,_nident)); 

  // Copy member data to SEXP objects
  double *pA0=REAL(A0); for(i=0;i<A0len;i++) pA0[i]=_val[i]; 
  int *pstr=INTEGER(str); for(i=0;i<_nident;i++) pstr[i]=_fpidx[i]+1; 

  SET_VECTOR_ELT(out, 0, A0); SET_STRING_ELT(names, 0, mkChar("A0"));
  SET_VECTOR_ELT(out, 1, str); SET_STRING_ELT(names, 1, mkChar("struct"));
  SET_VECTOR_ELT(out, 2, m); SET_STRING_ELT(names, 2, mkChar("m"));
  setAttrib(out, R_NamesSymbol, names); UNPROTECT(5); 
  return out;
}

////////////////////
// Wobj functions //
////////////////////

Wobj::Wobj(){ _m=0; _nval=0; }

// Allocate memory and assign values from SEXP W
Wobj::Wobj(SEXP W)
{
  int i, j, ct=0, len, *dW; double *pW; SEXP Welt; 

  // Allocate and populate _idx
  _m=length(W); _nval=0; _idx=(int*)Calloc(_m,int);
  for(i=0;i<_m;i++){
    PROTECT(Welt=VECTOR_ELT(W,i)); 
    dW=getdims(Welt); _idx[i]=_nval+dW[0]*dW[1]; _nval=_idx[i];
    UNPROTECT(1); 
  }

  // Allocate and populate _val
  _val=(double*)Calloc(_nval,double); ct=0;
  for(i=0;i<_m;i++){
    PROTECT(Welt=VECTOR_ELT(W,i)); pW=REAL(Welt); len=_idx[i]-ct; 
    for(j=ct;j<len+ct;j++) _val[j]=pW[j-ct]; ct=_idx[i];  
    UNPROTECT(1);
  }
  UNPROTECT(1);
}

// Populate Wobj from flat vector, indicies, and number of equations. 
Wobj::Wobj(int m, double *vals, int *idx)
{ 
  // Set _idx (W's _m indices) and _val (W's vals) 
  _m=m; _nval=0; _idx=(int*)Calloc(_m,int); 
  for(int i=0;i<_m;i++){_idx[i]=idx[i]; _nval+=idx[i];}
  _val=(double*)Calloc(_nval,double); 
  for(int i=0;i<_nval;i++) _val[i]=vals[i]; 
}

// Free allocated memory on destruction of object
Wobj::~Wobj(){Free(_idx); Free(_val);}

// Fill Wobj with values from SEXP W
void Wobj::setW(SEXP W)
{
  int i, j, len, ct=0, *dW; double *pW; SEXP Welt; 
  _m=length(W); _nval=0;

  // Allocate and populate _idx
  _idx=(int*)Calloc(_m,int);
  for(i=0;i<_m;i++){
    PROTECT(Welt=VECTOR_ELT(W,i)); dW=getdims(Welt); 
    _idx[i]=_nval+dW[0]*dW[1]; _nval=_idx[i];
    UNPROTECT(1); 
  }
  
  // Allocate and populate _val
  _val=(double*)Calloc(_nval,double); ct=0;
  for(i=0;i<_m;i++){
    PROTECT(Welt=VECTOR_ELT(W,i)); pW=REAL(Welt); len=_idx[i]-ct; 
    for(j=ct;j<len+ct;j++) _val[j]=pW[j-ct]; ct=_idx[i];  
    UNPROTECT(1);
  }
}

// Populate Wobj from flat vector, indicies, and number of equations. 
void Wobj::setWflat(int m, double *vals, int *idx)
{
  // Set _idx (W's _m indices) and _val (W's vals) 
  _m=m; _nval=0; _idx=(int*)Calloc(_m,int); 
  for(int i=0;i<_m;i++){_idx[i]=idx[i]; _nval+=idx[i];}
  _val=(double*)Calloc(_nval,double); 
  for(int i=0;i<_nval;i++) _val[i]=vals[i]; 
}

// Fill W[[idx]] with input Matrix
void Wobj::setWelt(const Matrix &Wmat, const int idx)
{
  SEXP tmp; PROTECT(tmp=C2Rmat(Wmat)); double *ptmp=REAL(tmp);

  // Get/Update size/index information 
  int i, len=Wmat.Storage(); _idx[idx-1]=_nval+len; 

  // alloc valtmp, _val==>valtmp, realloc _val, valtmp==>_val, free valtmp
  if(_nval){
    double *valtmp=(double*)Calloc(_nval,double); for(i=0;i<_nval;i++) valtmp[i]=_val[i]; 
    _val=(double*)Realloc(_val,_nval+len,double); for(i=0;i<_nval;i++) _val[i]=valtmp[i]; 
    Free(valtmp); 
  } else { _val=(double*)Calloc(len,double); }

  // Copy new element values to _val[_nval:_nval+len]
  for(i=_nval;i<_nval+len;i++) _val[i]=ptmp[i-_nval]; _nval+=len; 
  UNPROTECT(1); 
}

// Return W as a fully typed, compatible R object 
SEXP Wobj::getW() const
{
  Rprintf("Wobj::getW() called...\n");
  int i, j, ct=0, len;  double *pWelt;
  SEXP W, Welt; PROTECT(W=allocVector(VECSXP,_m));
  
  Rprintf("Looping over %d W elements\n",_m); 
  for(i=0;i<_m;i++){
    len=_idx[i]-ct; 

    PROTECT(Welt=allocVector(REALSXP,len)); pWelt=REAL(Welt); 
    for(j=ct;j<ct+len;j++) *pWelt++=_val[j]; ct=_idx[i]; 
    int n=(int)sqrt(len), dims[]={n,n}; setdims(Welt,2,dims); 
    SET_VECTOR_ELT(W,i,Welt); UNPROTECT(1);
  }
  UNPROTECT(1); 
  return W;
}

// Return W as flattened vector 
double* Wobj::getWvals() const
{
  int i; double *out=(double*)Calloc(_nval,double); 
  for(i=0;i<_nval;i++) out[i]=_val[i]; return out; 
}

// Return W indicies for reconstruction from the flattened vector
int* Wobj::getWindex() const
{ int i, *idx=(int*)Calloc(_m,int); for(i=0;i<_m;i++) idx[i]=_idx[i]; return idx; }

// Return W[[idx]] as a Matrix
ReturnMatrix Wobj::getWelt(int idx) const
{
  // Get start/end of index range
  idx-=1; int i, st=0, end, n; 
  if(idx) for(i=0;i<idx;i++) st+=_idx[i]; end=st+_idx[idx]; n=(int)sqrt(end-st); 
//   Rprintf("Wobj::getWelt ---- W[%d](%dx%d)==>_val[%d:%d] \n",idx,n,n,st,end);

  double *tmp=(double*)Calloc(end-st,double); for(i=st;i<end;i++) tmp[i-st]=_val[i]; 
  Matrix W=F2C(tmp,n,n); Free(tmp); W.Release(); return W.ForReturn(); 
}

void Wobj::clear(){_m=0; _nval=0; Free(_val); Free(_idx);}

/////////////////////
// Wlist functions //
/////////////////////

Wlist::Wlist(){_m=0; _n=0; _cdraw=0; _nval=0;}

Wlist::Wlist(SEXP Wpost, int num_draws)
{
//  Rprintf("Creating Wlist object with %d draws\n",num_draws); 
  // Initialize draw count variables
  _n=num_draws; _cdraw=0;

  // Allocate memory and populate index and value arrays 
  SEXP W, Widx, m; int i; 

  // Put W$W.posterior values in _val
//  Rprintf("W$W.posterior...\n"); 
  PROTECT(W=listElt(Wpost,"W")); double *pW=REAL(W); 
  _nval=length(W); _val=(double*)Calloc(_nval,double); 
  for(i=0;i<_nval;i++) _val[i]=pW[i];
  UNPROTECT(1); 

  // Put W$W.index values in _idx
//  Rprintf("W$W.index...\n"); 
  PROTECT(Widx=listElt(Wpost,"W.index")); int *pWidx=INTEGER(Widx); 
  _nidx=length(Widx);  _idx=(int*)Calloc(_nidx,int);
  for(i=0;i<_nidx;i++) _idx[i]=pWidx[i]; 
  UNPROTECT(1); 

  // Put W$m (number of equations) in _m
  PROTECT(m=listElt(Wpost,"m")); _m=INTEGER(m)[0]; UNPROTECT(1); 
}

Wlist::Wlist(const int num_draws, const int num_eq)
{
  _n=num_draws; _m=num_eq; _cdraw=0; _nidx=_n*_m; _nval=0;
  _idx=(int*)Calloc(_nidx,int); _val=(double*)Calloc(_n*_m*_m*(_m-1),double); 
}

Wlist::~Wlist(){Free(_idx); Free(_val);}

void Wlist::setWobj(Wobj &W, const int idx)
{
  int i, st=_nval, *Widx=W.getWindex(); double *pW=W.getWvals(); 
  for(i=0;i<_m;i++){_nval+=Widx[i]; _idx[(idx-1)*_m+i]=_nval;} 
  for(i=st;i<_nval;i++) _val[i]=pW[i-st]; Free(pW); Free(Widx); 
}

void Wlist::getWobj(Wobj& W, int idx) const 
{
  idx-=1; int i;

  // Find start/end index, # vals 
  int stidx=idx*_m, endidx=idx*_m+_m-1;
  int stval=(idx>0)?_idx[stidx-1]:0, endval=_idx[endidx];

  // Get index array
  int *idxs=(int*)Calloc(_m,int), ctidx=stval; 
  for(i=0;i<_m;i++){idxs[i]=_idx[i+stidx]-ctidx; ctidx+=idxs[i];}

  // Get values array
  double *vals=(double*)Calloc(endval-stval,double); 
  for(i=stval;i<endval;i++) vals[i-stval]=_val[i];

//   Rprintf("Wlist::getWobj(%d)--vals[%d:%d]\n",idx,stval,endval); 
//   Rprintf("--idxs[1:%d]:    ",_m);  for(i=0;i<_m;i++) Rprintf("  %d  ",idxs[i]); Rprintf("\n");

  // Set Wobj W, free allocated memory 
  W.setWflat(_m,vals,idxs); Free(vals); Free(idxs); 
}

int Wlist::num_eq() const {return _m;} 
int Wlist::num_draws() const {return _n;}

SEXP Wlist::toR() const 
{
  SEXP out, names, W, Widx, m; int i, *pidx; double *pW; 
  PROTECT(out=allocVector(VECSXP,3)); 
  PROTECT(names=allocVector(STRSXP,3)); 

  PROTECT(W=allocVector(REALSXP,_nval)); 
  pW=REAL(W); for(i=0;i<_nval;i++) pW[i]=_val[i]; 
  SET_VECTOR_ELT(out,0,W); SET_STRING_ELT(names,0,mkChar("W")); 
  UNPROTECT(1); 

  PROTECT(Widx=allocVector(INTSXP,_nidx)); 
  pidx=INTEGER(Widx); for(i=0;i<_nidx;i++) pidx[i]=_idx[i]; 
  SET_VECTOR_ELT(out,1,Widx); SET_STRING_ELT(names,1,mkChar("W.index")); 
  UNPROTECT(1); 
  
  PROTECT(m=allocVector(INTSXP,1)); INTEGER(m)[0]=_m; 
  SET_VECTOR_ELT(out,2,m); SET_STRING_ELT(names,2,mkChar("m")); 
  UNPROTECT(1); 

  setAttrib(out, R_NamesSymbol, names); 
  UNPROTECT(2); return out; 
}


/////////////////////
// UTobj functions //
/////////////////////

UTobj::UTobj(){_m=0;}

UTobj::UTobj(SEXP UT)
{
  int i, j, ct; double *pUTi; SEXP UTi;

  _m=length(UT); _dim=(int*)Calloc(_m,int); _nval=0; 
  for(i=0;i<_m;i++) {
    PROTECT(UTi=VECTOR_ELT(UT,i)); 
    _dim[i]=INTEGER(getAttrib(UTi,R_DimSymbol))[1]; _nval+=_m*_dim[i]; 
    UNPROTECT(1); 
  }

  _val=(double*)Calloc(_nval,double); ct=0;
  for(i=0;i<_m;i++){
    PROTECT(UTi=VECTOR_ELT(UT,i)); pUTi=REAL(UTi); 
    for(j=0;j<_m*_dim[i];j++) _val[ct++]=pUTi[j]; 
    UNPROTECT(1); 
  }
}

UTobj::~UTobj(){Free(_dim); Free(_val);}

void UTobj::setUTelt(Matrix &UT, const int idx)
{
  int i;
  if(!_m && idx-1){ 
    Rprintf("Error in UTobj::setUTelt - index must be zero for uninitialized objects");
    return;
  } 
  else if(!_m) {
    _m=UT.Nrows();  _dim=(int*)Calloc(_m,int); 
    _nval=0; _val=(double*)Calloc(_m*_m*_m,double); 
  }

  // Set dimension/value information 
  _dim[idx-1]=UT.Ncols(); double *pUT=UT.Store(); 
  for(i=_nval;i<_nval+UT.Storage();i++) _val[i]=pUT[i-_nval]; 
  _nval+=UT.Storage(); 
}

ReturnMatrix UTobj::getUTelt(int idx) const 
{
  idx-=1; int i, st=0; for(i=0;i<idx;i++) st+=_m*_dim[i]; 
  double *UT=(double*)Calloc(_m*_dim[idx],double);
  for(i=st;i<st+_m*_dim[idx];i++) UT[i-st]=_val[i]; 
  Matrix UTi=F2C(UT,_m,_dim[idx]); Free(UT); 
//   Matrix UTi(_m,_dim[idx]); UTi<<UT; 
  UTi.Release(); return UTi.ForReturn(); 
}

SEXP UTobj::toR() const
{
  int i, j, len, ct=0; SEXP out, tmp; double *ptmp; 
  PROTECT(out=allocVector(VECSXP,_m)); 
  for(i=0;i<_m;i++){
    len=_m*_dim[i]; 
    PROTECT(tmp=allocVector(REALSXP,len)); ptmp=REAL(tmp); 
    for(j=ct;j<ct+len;j++){ptmp[j-ct]=_val[j];} ct+=len; 
    int dUT[]={_m,_dim[i]}; setdims(tmp,2,dUT); SET_VECTOR_ELT(out,i,tmp); 
    UNPROTECT(1); 
  }
  UNPROTECT(1); return out;
}
