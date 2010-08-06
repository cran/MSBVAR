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
A0obj::A0obj(){ eqnCount=0; idRltnshpCount=0; numDraws=0; currentDraw=0; }

//Deconstruction handled by R

A0obj::A0obj(SEXP A0)
{
  int i, ct=0, *structPtr;//*str; 
  SEXP eqnListR, structR, valsR;//mR 
  double *vals; 

//   Rprintf("--- Setting eqnCount...\n"); 

  PROTECT(eqnListR=coerceVector(listElt(A0,"m"),INTSXP));
  eqnCount=INTEGER(eqnListR)[0]; 
  UNPROTECT(1); 
  
//   Rprintf("--- Setting freeParamIDX, strR...\n"); 

  PROTECT(structR=coerceVector(listElt(A0,"struct"),INTSXP)); 

  idRltnshpCount=length(structR); 
  structPtr=INTEGER(structR);
  freeParamIDX=(int*)R_alloc(idRltnshpCount,sizeof(int)); 

  for(i=0;i<idRltnshpCount;i++) 
	freeParamIDX[i]=structPtr[i]-1;

  UNPROTECT(1); 

//   Rprintf("--- Setting A0Struct...\n"); 

  A0Struct=(int*)R_alloc(eqnCount*eqnCount,sizeof(int)); 
  restrictionIDX=(int*)R_alloc(eqnCount*eqnCount-idRltnshpCount,sizeof(int));

  for(i=0;i<eqnCount*eqnCount;i++) 
  {  
    if(freeParamIDX[ct]==i)
    {
	A0Struct[i]=1; 
	ct++;
    }
    else
    {
	A0Struct[i]=0; 
	restrictionIDX[i-ct]=i;
    }
  }

//   Rprintf("--- Setting valsR...\n"); 

  PROTECT(valsR=coerceVector(listElt(A0,"A0"),REALSXP)); 

  numDraws=length(valsR)/idRltnshpCount; 
  currentDraw=0; 
  vals=REAL(valsR);

//   Rprintf("--- Setting A0Values...\n"); 

  A0Values=(double*)R_alloc(numDraws*idRltnshpCount,sizeof(double)); 

  for(i=0;i<numDraws*idRltnshpCount;i++) 
	A0Values[i]=vals[i]; 

  UNPROTECT(1); 
}

A0obj::A0obj(const Matrix &ident, const int ndraw)
{
  eqnCount=ident.Nrows(); 
  numDraws=ndraw; 
  currentDraw=0; 
  A0Struct=(int*)R_alloc(eqnCount*eqnCount,sizeof(int)); 

  int i, ct=0; 
  double *structPtr=C2F(ident); 
  idRltnshpCount=0; 

  for(i=0;i<eqnCount*eqnCount;i++) 
  {
    if(structPtr[i])
    {
	A0Struct[i]=1; 
	idRltnshpCount++;
    }
    else
    	A0Struct[i]=0; 
  }

  // Allocate memory to store values/indices
  A0Values=(double*)R_alloc(idRltnshpCount*numDraws,sizeof(double)); 
  freeParamIDX=(int*)R_alloc(idRltnshpCount,sizeof(int)); 
  restrictionIDX=(int*)R_alloc(eqnCount*eqnCount-idRltnshpCount,sizeof(int)); 
  
  for(i=0;i<eqnCount*eqnCount;i++)
  { 
    if(A0Struct[i])
	freeParamIDX[ct++]=i;
    else
	restrictionIDX[i-ct]=i;
  }
}

///////////////////////////////////////
// Store A0 matrix at supplied index //
///////////////////////////////////////

void A0obj::setA0(const Matrix &mat, const int index)
{
  int i, st=(index-1)*idRltnshpCount; 
  double *A0=C2F(mat); 

  for(i=0;i<idRltnshpCount;i++) 
	A0Values[st++]=A0[freeParamIDX[i]]; 
  //Free(A0);
}

int A0obj::fpidx(const int i) const { return freeParamIDX[i]; }
int A0obj::xidx(const int i) const { return restrictionIDX[i]; }
int A0obj::numfp() const { return idRltnshpCount; }

//////////////////////////////////////
// Store A0 matrix at current index //
//////////////////////////////////////

void A0obj::setA0c(const Matrix &mat) 
{
  double *A0=mat.Store(); 
  int i, st=currentDraw*idRltnshpCount;

  for(i=0;i<eqnCount*eqnCount;i++) 
    if(A0Struct[i]!=0) 
	A0Values[st++]=A0[i]; 
}

/////////////////////////////////////////////////////////////////////
// Returns A0's structure matrix as a 0/1 identification matrix 
/////////////////////////////////////////////////////////////////////

ReturnMatrix A0obj::structure() const
{
  Matrix out(eqnCount,eqnCount); 
  out=0; 
  double *pout=out.Store();

  for(int i=0; i<eqnCount*eqnCount; i++) 
	pout[i]=A0Struct[i]; 

  return out.ForReturn(); 
}

/////////////////////////////////////////
// Returns a COPY of the values stored //
/////////////////////////////////////////

double* A0obj::getA0vals() const
{
  double* out; 
  out=(double*)R_alloc(idRltnshpCount*numDraws,sizeof(double));  

  for(int i=0;i<idRltnshpCount*numDraws;i++) 
	out[i]=A0Values[i]; 

  return out; 
}

////////////////////////////////////////////////////////////
// Returns the A0 matrix at the supplied index [1:numDraws] //
////////////////////////////////////////////////////////////

ReturnMatrix A0obj::getA0(const int idx) const
{
  int i, st=idx*idRltnshpCount; 
  double *pmat=(double*)R_alloc(eqnCount*eqnCount,sizeof(double)); 

  for(i=0;i<eqnCount*eqnCount;i++) 
	pmat[i]=(A0Struct[i])?A0Values[st++]:0; 
  Matrix mat=F2C(pmat,eqnCount,eqnCount); 

  return mat.ForReturn(); 
}

///////////////////////////////////////////////
// Returns the A0 matrix at the current draw //
///////////////////////////////////////////////

ReturnMatrix A0obj::getA0c() const
{
  int i, st=currentDraw*idRltnshpCount; 
  Matrix mat(eqnCount,eqnCount); 
  double *pmat=mat.Store(); 

  for(i=0;i<eqnCount*eqnCount;i++) 
	pmat[i]=(A0Struct[i]!=0)?A0Values[st++]:0;

  return mat.ForReturn(); 
}

//////////////////////////////////////////////////////
// Returns SEXP list containing A0flat, A0struct, m //
//////////////////////////////////////////////////////

SEXP A0obj::toR() const 
{
  int i, A0len = idRltnshpCount*numDraws; 

  SEXP out; 
  PROTECT(out=allocVector(VECSXP,3)); 
  SEXP names; 
  PROTECT(names=allocVector(STRSXP,3));
  SEXP A0; 
  PROTECT(A0=allocVector(REALSXP, A0len));
  SEXP m; 
  PROTECT(m=allocVector(INTSXP,1)); 
  INTEGER(m)[0]=eqnCount; 
  SEXP str; 
  PROTECT(str=allocVector(INTSXP,idRltnshpCount)); 

  // Copy member data to SEXP objects
  double *pA0=REAL(A0); 
  for(i=0;i<A0len;i++) 
	pA0[i]=A0Values[i]; 
  
  int *pstr=INTEGER(str); 
  for(i=0;i<idRltnshpCount;i++) 
	pstr[i]=freeParamIDX[i]+1; 

  SET_VECTOR_ELT(out, 0, A0); SET_STRING_ELT(names, 0, mkChar("A0"));
  SET_VECTOR_ELT(out, 1, str); SET_STRING_ELT(names, 1, mkChar("struct"));
  SET_VECTOR_ELT(out, 2, m); SET_STRING_ELT(names, 2, mkChar("m"));
  setAttrib(out, R_NamesSymbol, names); UNPROTECT(5); 

  return out;
}

////////////////////
// Wobj functions //
////////////////////

Wobj::Wobj(){ eqnCount=0; numVals=0; }

// Allocate memory and assign values from SEXP W
Wobj::Wobj(SEXP W)
{
  int i, j, ct=0, len, *dimsW; 
  double *pW; 
  SEXP Welt; 

  // Allocate and populate WMatIDX
  eqnCount=length(W); 
  numVals=0; 
  WMatIDX=(int*)R_alloc(eqnCount,sizeof(int));

  for(i=0;i<eqnCount;i++)
  {
    PROTECT(Welt=VECTOR_ELT(W,i)); 

    dimsW=getdims(Welt); 
    WMatIDX[i]=numVals+dimsW[0]*dimsW[1]; 
    numVals=WMatIDX[i];

    UNPROTECT(1); 
  }

  // Allocate and populate flatWVals
  flatWVals=(double*)R_alloc(numVals,sizeof(double)); 
  ct=0;

  for(i=0;i<eqnCount;i++)
  {
    PROTECT(Welt=VECTOR_ELT(W,i)); 

    pW=REAL(Welt); 
    len=WMatIDX[i]-ct; 

    for(j=ct;j<len+ct;j++) 
	flatWVals[j]=pW[j-ct]; 

    ct=WMatIDX[i];  

    UNPROTECT(1);
  }
  //I'm pretty sure this is an extra unprotect.
  //As in, it's unprotecting something it should not.
  //If something breaks around A0 this is probably the culprit.
  //UNPROTECT(1);
}

// Populate Wobj from flat vector, indicies, and number of equations. 
Wobj::Wobj(int m, double *vals, int *idx)
{ 
  // Set WMatIDX (W's eqnCount indices) and flatWVals (W's vals) 
  eqnCount=m; 
  numVals=0; 
  WMatIDX=(int*)R_alloc(eqnCount,sizeof(int)); 

  for(int i=0;i<eqnCount;i++)
  {
	WMatIDX[i]=idx[i]; 
	numVals+=idx[i];
  }

  flatWVals=(double*)R_alloc(numVals,sizeof(double)); 

  for(int i=0;i<numVals;i++) 
	flatWVals[i]=vals[i]; 
}

// Fill Wobj with values from SEXP W
void Wobj::setW(SEXP W)
{
  int i, j, len, ct=0, *dimsW; 
  double *pW; 
  SEXP Welt; 
  eqnCount=length(W); 
  numVals=0;

  // Allocate and populate WMatIDX
  WMatIDX=(int*)R_alloc(eqnCount,sizeof(int));

  for(i=0;i<eqnCount;i++)
  {
    PROTECT(Welt=VECTOR_ELT(W,i)); 

    dimsW=getdims(Welt); 
    WMatIDX[i]=numVals+dimsW[0]*dimsW[1]; 
    numVals=WMatIDX[i];

    UNPROTECT(1); 
  }
  
  // Allocate and populate flatWVals
  flatWVals=(double*)R_alloc(numVals,sizeof(double)); ct=0;

  for(i=0;i<eqnCount;i++)
  {
    PROTECT(Welt=VECTOR_ELT(W,i)); 

    pW=REAL(Welt); 
    len=WMatIDX[i]-ct; 

    for(j=ct;j<len+ct;j++) 
	flatWVals[j]=pW[j-ct]; 

    ct=WMatIDX[i];  

    UNPROTECT(1);
  }
}

// Populate Wobj from flat vector, indicies, and number of equations. 
void Wobj::setWflat(int m, double *vals, int *idx)
{
  // Set WMatIDX (W's eqnCount indices) and flatWVals (W's vals) 
  eqnCount=m; 
  numVals=0; 
  WMatIDX=(int*)R_alloc(eqnCount,sizeof(int)); 

  for(int i=0;i<eqnCount;i++)
  {
	WMatIDX[i]=idx[i]; 
	numVals+=idx[i];
  }

  flatWVals=(double*)R_alloc(numVals,sizeof(double)); 

  for(int i=0;i<numVals;i++) 
	flatWVals[i]=vals[i]; 
}

// Fill W[[idx]] with input Matrix
void Wobj::setWelt(const Matrix &Wmat, const int idx)
{
  SEXP tmp; 
  PROTECT(tmp=C2Rmat(Wmat)); 

  double *ptmp=REAL(tmp);

  // Get/Update size/index information 
  int i, len=Wmat.Storage(); 
  WMatIDX[idx-1]=numVals+len; 

  // alloc valtmp, flatWVals==>valtmp, realloc flatWVals, valtmp==>flatWVals, free valtmp
  if(numVals)
  {
    double *valtmp=(double*)R_alloc(numVals,sizeof(double)); 

    for(i=0;i<numVals;i++) 
	valtmp[i]=flatWVals[i]; 

    flatWVals=(double*)R_alloc(numVals+len,sizeof(double)); 
    for(i=0;i<numVals;i++) 
	flatWVals[i]=valtmp[i];
  } 
  else 
  	flatWVals=(double*)R_alloc(len,sizeof(double));

  // Copy new element values to flatWVals[numVals:numVals+len]
  for(i=numVals;i<numVals+len;i++) 
	flatWVals[i]=ptmp[i-numVals]; 
  numVals+=len; 

  UNPROTECT(1); 
}

// Return W as a fully typed, compatible R object 
SEXP Wobj::getW() const
{
  Rprintf("Wobj::getW() called...\n");

  int i, j, ct=0, len;  
  double *pWelt;
  SEXP W, Welt; 
  PROTECT(W=allocVector(VECSXP,eqnCount));

  Rprintf("Looping over %d W elements\n",eqnCount); 

  for(i=0;i<eqnCount;i++)
  {
    len=WMatIDX[i]-ct; 

    PROTECT(Welt=allocVector(REALSXP,len)); 
 
    pWelt=REAL(Welt); 

    for(j=ct;j<ct+len;j++) 
	*pWelt++=flatWVals[j]; 

    ct=WMatIDX[i]; 
    int n=(int)sqrt((double)len);
    int dims[]={n,n}; 
    setdims(Welt,2,dims); 

    SET_VECTOR_ELT(W,i,Welt); 
    UNPROTECT(1);
  }
  UNPROTECT(1); 

  return W;
}

// Return W as flattened vector 
double* Wobj::getWvals() const
{
  int i; 
  double *out=(double*)R_alloc(numVals,sizeof(double)); 

  for(i=0;i<numVals;i++) 
	out[i]=flatWVals[i]; 

  return out; 
}

// Return W indicies for reconstruction from the flattened vector
int* Wobj::getWindex() const
{ 
  int i, *idx=(int*)R_alloc(eqnCount,sizeof(int)); 

  for(i=0;i<eqnCount;i++) 
	idx[i]=WMatIDX[i]; 

  return idx; 
}

// Return W[[idx]] as a Matrix
ReturnMatrix Wobj::getWelt(int idx) const
{
  // Get start/end of index range
  idx-=1; 
  int i, st=0, end, n; 

  if(idx) 
	for(i=0;i<idx;i++) 
		st+=WMatIDX[i]; 

  end=st+WMatIDX[idx]; 
  n=(int)sqrt((double)end-st); 
//   Rprintf("Wobj::getWelt ---- W[%d](%dx%d)==>flatWVals[%d:%d] \n",idx,n,n,st,end);

  double *tmp=(double*)R_alloc(end-st,sizeof(double)); 

  for(i=st;i<end;i++) 
	tmp[i-st]=flatWVals[i]; 

  Matrix W=F2C(tmp,n,n); 

  return W.ForReturn(); 
}

void Wobj::clear()
{
	eqnCount=0; 
	numVals=0; 
}

/////////////////////
// Wlist functions //
/////////////////////

Wlist::Wlist(){eqnCount=0; numDraws=0; currentDraw=0; numVals=0;}

Wlist::Wlist(SEXP Wpost, int num_draws)
{
  // Initialize draw count variables
  numDraws=num_draws; 
  currentDraw=0;

  // Allocate memory and populate index and value arrays 
  SEXP W, Widx, m; 
  int i; 

  // Put W$W.posterior values in flatWValList
  PROTECT(W=listElt(Wpost,"W")); 

  double *pW=REAL(W); 
  numVals=length(W); 
  flatWValList=(double*)R_alloc(numVals,sizeof(double)); 

  for(i=0;i<numVals;i++) 
	flatWValList[i]=pW[i];

  UNPROTECT(1); 

  // Put W$W.index values in arrayIDX
  PROTECT(Widx=listElt(Wpost,"W.index")); 

  int *pWidx=INTEGER(Widx); 
  IDXcount=length(Widx);  
  arrayIDX=(int*)R_alloc(IDXcount,sizeof(int));

  for(i=0;i<IDXcount;i++) 
	arrayIDX[i]=pWidx[i]; 

  UNPROTECT(1); 

  // Put W$m (number of equations) in eqnCount
  PROTECT(m=listElt(Wpost,"m")); 
  eqnCount=INTEGER(m)[0]; 
  UNPROTECT(1); 
}

Wlist::Wlist(const int num_draws, const int num_eq)
{
  numDraws=num_draws; 
  eqnCount=num_eq; 
  currentDraw=0; 
  IDXcount=numDraws*eqnCount; 
  numVals=0;

  arrayIDX=(int*)R_alloc(IDXcount+8,sizeof(int)); 
  flatWValList=(double*)R_alloc(numDraws*eqnCount*eqnCount*eqnCount,sizeof(double)); 
}

//R handles deconstruction

void Wlist::setWobj(Wobj &W, const int idx)
{
  int i, st=numVals, *Widx=W.getWindex(); 
  double *pW=W.getWvals(); 

  for(i=0;i<eqnCount;i++)
  {
	numVals+=Widx[i]; 
	arrayIDX[(idx-1)*eqnCount+i]=numVals;
  } 

  for(i=st;i<numVals;i++) 
	flatWValList[i]=pW[i-st]; 
}

void Wlist::getWobj(Wobj& W, int idx) const 
{
  idx-=1; 
  int i;

  // Find start/end index, # vals 
  int stidx=idx*eqnCount;
  int endidx=idx*eqnCount+eqnCount-1;
  int stval; 
  int endval=arrayIDX[endidx];
  if(idx>0)
	stval=arrayIDX[stidx-1];
  else
	stval=0;

  // Get index array
  int *idxs=(int*)R_alloc(eqnCount,sizeof(int));
  int ctidx=stval; 

  for(i=0;i<eqnCount;i++)
  {
	idxs[i]=arrayIDX[i+stidx]-ctidx; 
	ctidx+=idxs[i];
  }

  // Get values array
  double *vals=(double*)R_alloc(endval-stval,sizeof(double)); 

  for(i=stval;i<endval;i++) 
	vals[i-stval]=flatWValList[i];

  // Set Wobj W
  W.setWflat(eqnCount,vals,idxs); 
}

int Wlist::num_eq() const {return eqnCount;} 
int Wlist::num_draws() const {return numDraws;}

SEXP Wlist::toR() const 
{
  SEXP out, names, W, Widx, m; int i, *pidx; 
  double *pW; 
  
  PROTECT(out=allocVector(VECSXP,3)); 
  PROTECT(names=allocVector(STRSXP,3)); 

  PROTECT(W=allocVector(REALSXP,numVals)); 
  pW=REAL(W); 
  for(i=0;i<numVals;i++) 
	pW[i]=flatWValList[i]; 
  SET_VECTOR_ELT(out,0,W); 
  SET_STRING_ELT(names,0,mkChar("W")); 
  UNPROTECT(1); 

  PROTECT(Widx=allocVector(INTSXP,IDXcount)); 
  pidx=INTEGER(Widx); 
  for(i=0;i<IDXcount;i++) 
	pidx[i]=arrayIDX[i]; 
  SET_VECTOR_ELT(out,1,Widx); 
  SET_STRING_ELT(names,1,mkChar("W.index")); 
  UNPROTECT(1); 
  
  PROTECT(m=allocVector(INTSXP,1)); INTEGER(m)[0]=eqnCount; 
  SET_VECTOR_ELT(out,2,m); 
  SET_STRING_ELT(names,2,mkChar("m")); 
  UNPROTECT(1); 

  setAttrib(out, R_NamesSymbol, names); 
  UNPROTECT(2); 

  return out; 
}


/////////////////////
// UTobj functions //
/////////////////////

UTobj::UTobj(){eleCount=0;}

UTobj::UTobj(SEXP UT)
{
  int i, j, ct; 
  double *pUTi; SEXP UTi;

  eleCount=length(UT); 
  dimensionList=(int*)R_alloc(eleCount,sizeof(int)); numVals=0; 
  
  for(i=0;i<eleCount;i++) 
  {
    PROTECT(UTi=VECTOR_ELT(UT,i)); 
    dimensionList[i]=INTEGER(getAttrib(UTi,R_DimSymbol))[1]; numVals+=eleCount*dimensionList[i]; 
    UNPROTECT(1); 
  }

  objVals=(double*)R_alloc(numVals,sizeof(double)); ct=0;

  for(i=0;i<eleCount;i++)
  {
    PROTECT(UTi=VECTOR_ELT(UT,i)); 
    pUTi=REAL(UTi); 

    for(j=0;j<eleCount*dimensionList[i];j++) 
	objVals[ct++]=pUTi[j]; 
    UNPROTECT(1); 
  }
}

void UTobj::setUTelt(Matrix &UT, const int idx)
{
  int i;

  if(!eleCount && idx-1)
  { 
    Rprintf("Error in UTobj::setUTelt - index must be zero for uninitialized objects");
    return;
  } 
  else if(!eleCount) 
  {
    eleCount=UT.Nrows();  
    dimensionList=(int*)R_alloc(eleCount,sizeof(int)); 
    numVals=0;
    objVals=(double*)R_alloc(eleCount*eleCount*eleCount,sizeof(double)); 
  }

  // Set dimension/value information 
  dimensionList[idx-1]=UT.Ncols(); 
  double *pUT=UT.Store(); 

  for(i=numVals;i<numVals+UT.Storage();i++) 
	objVals[i]=pUT[i-numVals]; 
  numVals+=UT.Storage(); 
}

ReturnMatrix UTobj::getUTelt(int idx) const 
{
  idx-=1; 
  int i, st=0; 

  for(i=0;i<idx;i++) 
	st+=eleCount*dimensionList[i]; 

  double *UT=(double*)R_alloc(eleCount*dimensionList[idx],sizeof(double));

  for(i=st;i<st+eleCount*dimensionList[idx];i++) 
	UT[i-st]=objVals[i]; 
  Matrix UTi=F2C(UT,eleCount,dimensionList[idx]); 

  return UTi.ForReturn(); 
}

SEXP UTobj::toR() const
{
  int i, j, len, ct=0; 
  SEXP out, tmp; 
  double *ptmp; 

  PROTECT(out=allocVector(VECSXP,eleCount)); 

  for(i=0;i<eleCount;i++)
  {
    len=eleCount*dimensionList[i]; 

    PROTECT(tmp=allocVector(REALSXP,len)); 

    ptmp=REAL(tmp); 

    for(j=ct;j<ct+len;j++)
	ptmp[j-ct]=objVals[j]; 
    ct+=len; 
    int dUT[]={eleCount,dimensionList[i]}; 
    setdims(tmp,2,dUT); 
    SET_VECTOR_ELT(out,i,tmp); 

    UNPROTECT(1); 
  }

  UNPROTECT(1); 

  return out;
}
