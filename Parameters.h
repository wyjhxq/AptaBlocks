#ifndef Parameters_h
#define Parameters_h

#include  "fold_vars.h"
#include  "fold.h"
#include  "part_func.h"
#include  "inverse.h"
#include  "RNAstruct.h"
#include  "treedist.h"
#include  "stringdist.h"
#include  "profiledist.h"
#include  "part_func.h"
#include  "params.h"



struct Parameters
{
    pf_paramT *pf_parameters;
    double KT;
    double Energy_Threshold;
    double lambda1;
    double lambda2;
    char NcBasePair[4][2];
    char NuclBase[4];
    double T_start;
    double T_end;
    double Cooling;
};

//golbal variable
extern double StackingEnerogyDNARNA[4][4];
extern double StackingEnerogyRNARNA[6][6];
extern double StackingEnerogyRNARNA1999[6][6];

void Init_Param(struct Parameters *Param, double betaScale, double temperature, double EnergyThreshold, double lambda1, double lambda2, double T_start, double T_end, double Cooling);


#endif


