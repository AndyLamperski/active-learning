#ifndef ACTIVE_H
#define ACTIVE_H

class factor {
public:
  int NumVar;
  int* Card;
  int NumScope;
  int* Scope;
  int NumValues;
  int* Stride;
  double* Values;
  
  factor();
  factor(int, int [], int, int []);
  factor(int, int [], int, int [], double []);
  void makeStrides();
  double entropy();
  void marginalize(factor, int);
  void product(factor,factor);
  void divide(factor, factor);
  void substitute(factor, int , int);
  void setValues(double*);
  void normalize();
  void surprise( factor );
  void print();


};

int maxIndex(int, double*);
double* linspace(int, double, double);

#define BOTH 0
#define THRESHOLD 1
#define SLOPE 2

class activeLearningModel {
public:

  int numX;
  int numT;
  int numS;

  double* X;
  double* S;
  double* T;

  // Main factors
  factor::factor Pc_tsx;
  factor::factor Pts;

  // Auxiallary Factors for learning step
  factor::factor Pcts_x;
  
  factor::factor Pc_x;
  factor::factor Sc_x;
  factor::factor HC_x;

  factor::factor Pt;
  factor::factor Pct_x;
  factor::factor Pc_tx;
  factor::factor Sc_tx;
  factor::factor HC_tx;
  factor::factor SC_tx;
  factor::factor HC_Tx;
  double *IC_T_x;

  factor::factor Ps;
  factor::factor Pcs_x;
  factor::factor Pc_sx;
  factor::factor Sc_sx;
  factor::factor HC_sx;
  factor::factor SC_sx;
  factor::factor HC_Sx;
  double *IC_S_x;

  factor::factor Sc_tsx;
  factor::factor HC_tsx;
  factor::factor SC_tsx;
  factor::factor SC_Tsx;
  factor::factor HC_TSx;
  double *IC_TS_x;

  // Auxilliary update factors
  factor::factor likelihood_X;
  factor::factor likelihood;
  factor::factor posterior;

  // MAP Values
  double mapT;
  double mapS;

  activeLearningModel(int, double, double, int, double, double, int, double, double);
  void defaultPrior(); 
  void modelDistribution();
  int optimalInput(int objectiveType);
  void updatePosterior(int, int);
  void updateMap();
};


#endif
