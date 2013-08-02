#include <iostream>
#include <math.h>
#include <float.h>
#include "active.h"

factor::factor() {

}

factor::factor(int nVar, int c[], int ns, int s[]) {
  NumVar = nVar;
  Card = new int[NumVar];
  for (int i=0; i<NumVar; i++)
    Card[i] = c[i];

  NumScope = ns;
  std::sort(s,s+ns);
  Scope = new int[NumScope];
  int nVal = 1;
  for (int i=0; i<NumScope; i++) {
    Scope[i] = s[i];
    nVal = nVal * Card[s[i]];
  }
  
  
  NumValues = nVal;

  makeStrides();

  Values = new double[NumValues];
}

factor::factor (int nVar, int c[], int ns, int s[], double val[]) {
  
  factor(nVar,c,ns,s);
  setValues(val);

}

void factor::makeStrides() {
  
  Stride = new int[NumVar];
  for (int i=0; i< NumVar; i++)
    Stride[i] = 0;

  int j;
  int sLast = 1;
  for (int i=0; i< NumScope; i++){
    j = Scope[i];
    Stride[j] = sLast;
    sLast = sLast * Card[j];
  }
}

double factor::entropy() {
  double ent = 0;
  for (int i=0; i<NumValues; i++) {
    if (Values[i] > DBL_EPSILON)
	ent = ent - Values[i] * log2(Values[i]);
  }
  return ent;
}

// Replaces P.Values with values of PJ with variable y marginalized. 
// Assumes that P.Scope = (PJ.Scope - {y}). 
void factor::marginalize(factor PJ, int y) {

  int assignment [NumScope];
  
  for (int i=0; i<NumScope; i++) {
    assignment[i] = 0;
  }
  
  // Main summation step
  int StrideY = PJ.Stride[y];
  int CardY = Card[y];

  int k = 0;
  double curVal;

  for (int i=0; i < NumValues; i++) {
    curVal = 0;
    for (int j=0; j < CardY; j++) 
      curVal += PJ.Values[k+j*StrideY];
    
    Values[i] = curVal;
    
    for (int l=0; l < NumScope; l++) {
      assignment[l] = assignment[l]+1;
      if (assignment[l] == Card[Scope[l]]) {
	assignment[l] = 0;
	k = k - (Card[Scope[l]]-1)*PJ.Stride[Scope[l]];
      }
      else {
	k = k + PJ.Stride[Scope[l]];
	break;
      }
    }
  }
  
}





// Assumes that P.Scope = (PS.Scope - {y});
void factor::substitute(factor PS, int y, int yInd) {
  int assignment [NumScope];
  
  for (int i=0; i<NumScope; i++) {
    assignment[i] = 0;
  }

  int StrideY = PS.Stride[y];
  

  int k = 0;

  for (int i=0; i<NumValues; i++) {
    Values[i] = PS.Values[k+yInd*StrideY];

    for (int l=0; l<NumScope; l++) {
      assignment[l] = assignment[l]+1;
      if (assignment[l] == Card[Scope[l]]) {
	assignment[l] = 0;
	k -= (Card[Scope[l]] - 1)*PS.Stride[Scope[l]];
      }
      else {
	k += PS.Stride[Scope[l]];
	break;
      }
    }
  }
}

// Replaces the P.product(P1,P2) replaces the values of P with
// the values of P1 * P2; Assumes that P.Scope is the union of the 
// of P1.Scope and P2.Scope.
void factor::product(factor P1, factor P2) {
  
  int assignment[NumScope];

  for (int i=0; i<NumScope; i++) {
    assignment[i] = 0;
  }

  int j = 0;
  int k = 0;
  for (int i=0; i<NumValues; i++) {
    Values[i] = P1.Values[j] * P2.Values[k];
    for (int l = 0; l < NumScope; l++) {
      assignment[l] = assignment[l]+1;
      if (assignment[l] == Card[Scope[l]]) {
	assignment[l] = 0;
	j -= (Card[Scope[l]]-1)*P1.Stride[Scope[l]];
	k -= (Card[Scope[l]]-1)*P2.Stride[Scope[l]];
      }
      else {
	j += P1.Stride[Scope[l]];
	k += P2.Stride[Scope[l]];
	break;
      }
    } 
  }
}

void factor::divide(factor P1, factor P2) {
  int assignment[NumScope];

  for (int i=0; i<NumScope; i++) {
    assignment[i] = 0;
  }

  int j = 0;
  int k = 0;
  for (int i=0; i<NumValues; i++) {
    if (P2.Values[k] > DBL_EPSILON) 
      Values[i] = P1.Values[j] / P2.Values[k];
    else
      Values[i] = 1.0;

    for (int l = 0; l < NumScope; l++) {
      assignment[l] = assignment[l]+1;
      if (assignment[l] == Card[Scope[l]]) {
	assignment[l] = 0;
	j -= (Card[Scope[l]]-1)*P1.Stride[Scope[l]];
	k -= (Card[Scope[l]]-1)*P2.Stride[Scope[l]];
      }
      else {
	j += P1.Stride[Scope[l]];
	k += P2.Stride[Scope[l]];
	break;
      }
    } 
  }
}



void factor::setValues(double* val) {
  for (int i=0; i < NumValues; i++) {
    Values[i] = val[i];
  }
}



void factor::normalize() {
  double tot = 0.0;
  for (int i=0; i< NumValues; i++) 
    tot += Values[i];
  for (int i=0; i< NumValues; i++)
    Values[i] = Values[i] / tot;
}

// Assumes that P and PS have the same scope
void factor::surprise( factor PS ) {
  for (int i=0; i < NumValues; i++) {
    if (PS.Values[i] > DBL_EPSILON) 
      Values[i] = -PS.Values[i] * log2(PS.Values[i]);
    else
      Values[i] = 0;
  }
}

void factor::print() {
  for (int i=0; i<NumValues; i++) 
    std::cout << Values[i] << " ";

  std::cout << std::endl;
}

// Active Learning 
int maxIndex(int NumX, double* X) {
  double maxValue = DBL_MIN;
  int maxInd = 0;
  for (int i=0; i<NumX; i++) {
    if (X[i] > maxValue) { 
      maxValue = X[i];
      maxInd = i;
    }
  }
  return maxInd;
} 

double* linspace(int nX,double minX,double maxX) {
  double dX = (maxX - minX) / double(nX-1); 

  double x = minX;

  double *Xarr = new double [nX];

  for (int i=0; i<nX; i++) {
    Xarr[i] = x;
    x = x+dX;
  }
  return Xarr;
}

activeLearningModel::activeLearningModel(int nX, double minX, double maxX, int nT, double minT, double maxT, int nS, double minS, double maxS) {

  numX = nX;
  X = new double[nX];
  X = linspace(nX,minX,maxX);

  numT = nT;
  T = new double[nT];
  T = linspace(nT,minT,maxT);

  numS = nS;
  S = new double[nS];
  S = linspace(nS,minS,maxS);

  int card [] = {2,numT,numS,numX};
  int scopeTot [] = {0,1,2,3};
  int nV = 4;

  Pc_tsx = factor(nV,card,4,scopeTot);
  
  int scopeTS [] = {1,2};
  Pts = factor(nV,card,2,scopeTS);

  Pcts_x = factor(nV,card,4,scopeTot);

  int scopeCX [] = {0,3};
  Pc_x = factor(nV,card,2,scopeCX);
  Sc_x = factor(nV,card,2,scopeCX);
  int scopeX [] = {3};
  HC_x = factor(nV,card,1,scopeX);

  int scopeT [] = {1};
  Pt = factor(nV,card,1,scopeT);
  int scopeCTX [] = {0,1,3};
  Pct_x = factor(nV,card,3,scopeCTX);
  Pc_tx = factor(nV,card,3,scopeCTX);
  Sc_tx = factor(nV,card,3,scopeCTX);
  int scopeTX [] = {1,3};
  HC_tx = factor(nV,card,2,scopeTX);
  SC_tx = factor(nV,card,2,scopeTX);
  HC_Tx = factor(nV,card,1,scopeX);
  IC_T_x = new double [numX];

  int scopeS [] = {2};
  Ps = factor(nV,card,1,scopeS);
  int scopeCSX [] = {0,2,3};
  Pcs_x = factor(nV,card,3,scopeCSX);
  Pc_sx = factor(nV,card,3,scopeCSX);
  Sc_sx = factor(nV,card,3,scopeCSX);
  int scopeSX [] = {2,3};
  HC_sx = factor(nV,card,2,scopeSX);
  SC_sx = factor(nV,card,2,scopeSX);
  HC_Sx = factor(nV,card,1,scopeX);
  IC_S_x = new double [numX];

  Sc_tsx = factor(nV,card,4,scopeTot);
  int scopeTSX [] = {1,2,3};
  HC_tsx = factor(nV,card,3,scopeTSX);
  SC_tsx = factor(nV,card,3,scopeTSX);
  SC_Tsx = factor(nV,card,2,scopeSX);
  HC_TSx = factor(nV,card,1,scopeX);	       
  IC_TS_x = new double [numX];

  int scopeCTS [] = {0,1,2};
  likelihood_X = factor(nV,card,3,scopeCTS);
  likelihood = factor(nV,card,2,scopeTS);
  posterior = factor(nV,card,2,scopeTS);

  defaultPrior();
  modelDistribution();
}

void activeLearningModel::defaultPrior() {
  double Tmean = (T[numT-1] + T[0]) / 2.0;
  double Tstd = (T[numT-1] - T[0]) / 2.0;
  double Smean = (S[numS-1] + S[0]) / 2.0;
  double Sstd = (S[numS-1] - S[0]) / 2.0;


  double *Tvalues = new double [numT];

  int card [] = {2,numT,numS,numX};

  for (int i=0; i<numT; i++) {
    Tvalues[i] = exp(-.5*pow(T[i] - Tmean,2) / pow(Tstd,2));
  }

  Pt.setValues(Tvalues);
  Pt.normalize();

  delete [] Tvalues;

  double *Svalues = new double [numS]; 

  for (int i=0; i<numS; i++) {
    Svalues[i] = exp(-.5*pow(S[i]-Smean,2) / pow(Sstd,2));
  }

  Ps.setValues(Svalues);
  Ps.normalize();

  delete[] Svalues;

  Pts.product(Pt,Ps);
}

void activeLearningModel::modelDistribution() {
  int nVal = numX * numT * numS;

  int strideT = 1;
  int strideS = numT;
  int strideX = numT * numS;
  
  int iT = 0;
  int iS = 0;
  int iX = 0;

  double t;
  double s;
  double x;

  double* modValues = new double [2*nVal];
  double y;

  for (int i=0; i<nVal; i++) {
   
    iT = i % numT;
    iS = (i / strideS) % numS;
    iX = (i / strideX) % numX;
   
    t = T[iT];
    s = S[iS];
    x = X[iX];
    
    
    y = 1./(1.+exp(-s*(x-t)));

    modValues[2*i] = 1.0 - y;
    modValues[2*i + 1] = y;
    
  }

  Pc_tsx.setValues(modValues);
  
  delete[] modValues;
}

int activeLearningModel::optimalInput(const int objectiveType) {
  
  int xInd = 0;
  
  Pcts_x.product(Pc_tsx, Pts);
  
  Pct_x.marginalize(Pcts_x,2);
  Pc_x.marginalize(Pct_x,1);

  Sc_x.surprise(Pc_x);
  HC_x.marginalize(Sc_x,0);

  int maxH = maxIndex(numX,HC_x.Values);

  switch (objectiveType) {
  case THRESHOLD: {
    
    Pt.marginalize(Pts,2);
    Pct_x.marginalize(Pcts_x,2);
    Pc_tx.divide(Pct_x, Pt);
    Sc_tx.surprise(Pc_tx);
    HC_tx.marginalize(Sc_tx,0);
    SC_tx.product(Pt,HC_tx);
    HC_Tx.marginalize(SC_tx,1);
    for (int i=0; i< numX; i++) {
      IC_T_x[i] = HC_x.Values[i] - HC_Tx.Values[i];
    }
    xInd = maxIndex(numX,IC_T_x);
    
    return xInd;
  }
  case SLOPE: {
    Ps.marginalize(Pts,1);
    Pcs_x.marginalize(Pcts_x,1);
    Pc_sx.divide(Pcs_x,Ps);
    Sc_sx.surprise(Pc_sx);
    HC_sx.marginalize(Sc_sx,0);
    SC_sx.product(Ps,HC_sx);
    HC_Sx.marginalize(SC_sx,2);
    for (int i=0; i<numX; i++) 
      IC_S_x[i] = HC_x.Values[i] - HC_Sx.Values[i];

    xInd = maxIndex(numX,IC_S_x);
   
    return xInd;
  }
  default: {
    Sc_tsx.surprise(Pc_tsx);
    HC_tsx.marginalize(Sc_tsx,0);
    SC_tsx.product(Pts,HC_tsx);
    SC_Tsx.marginalize(SC_tsx,1);
    HC_TSx.marginalize(SC_Tsx,2);
 
    for (int i=0; i < numX; i++) {
      IC_TS_x[i] = HC_x.Values[i] - HC_TSx.Values[i];
    }   

    xInd = maxIndex(numX,IC_TS_x);
    return xInd;
  }
  }
}

void activeLearningModel::updatePosterior(int xInd, int choice) {
  likelihood_X.substitute(Pc_tsx,3,xInd);
  likelihood.substitute(likelihood_X,0,choice);
  posterior.product( Pts , likelihood);
  posterior.normalize();
  Pts.setValues(posterior.Values);

  updateMap();

}

void activeLearningModel::updateMap() {
  int maxInd = maxIndex(Pts.NumValues,Pts.Values);
  int iT = maxInd % numT;
  int iS = (maxInd / numT) % numS;
  mapT = T[iT];
  mapS = S[iS];
}
