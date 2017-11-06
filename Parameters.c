#include  <math.h>
#include "Parameters.h"

// DNA <-> RNA:
//>Stacking DNA 5' X1 Y1 3'
//>         RNA 3' X2 Y2 5'
//>X1X2 = TA GC CG AU (row headings)
//>Y1Y2 = TA GC CG AU (column headings)
double StackingEnerogyDNARNA[4][4]={
    // TA    GC    CG     AU
    /*TA*/ { -1.0, -0.9, -1.3, -0.6},
    /*GC*/ { -2.1, -2.1, -2.7, -1.5},
    /*CG*/ { -1.8, -1.7, -2.9, -1.6},
    /*AU*/ { -0.9, -0.9, -1.1, -0.2}
};

//Turner_2004
//>Stacking 5' X1 Y1 3'
//>         3' X2 Y2 5'
//>X1X2 = AU CG GC UA GU UG (row headings)
//>Y1Y2 = AU CG GC UA GU UG (column headings)
double StackingEnerogyRNARNA[6][6]={
    //         AU    CG    GC     UA     GU    UG
    /*AU*/ { -0.9, -2.2, -2.1,  -1.1,  -0.6, -1.4},
    /*CG*/ {-2.1, -3.3, -2.4, -2.1, -1.4, -2.1},
    /*GC*/ {-2.4, -3.4, -3.3, -2.2, -1.5, -2.5},
    /*UA*/ {-1.3, -2.4, -2.1,  -0.9,  -1.0, -1.3},
    /*GU*/ {-1.3, -2.5, -2.1, -1.4,  -0.5,  1.3},
    /*UG*/ { -1.0, -1.5, -1.4,  -0.6,  0.3,  -0.5}
};




//nupack rna_1999
//>Stacking 5' X1 Y1 3'
//>         3' X2 Y2 5'
//>X1X2 = AU CG GC UA GU UG (row headings)
//>Y1Y2 = AU CG GC UA GU UG (column headings)
double StackingEnerogyRNARNA1999[6][6]={
    //         AU    CG    GC     UA     GU    UG
    /*AU*/ { -0.9, -2.1, -1.7,  -0.9,  -0.5, -1.0},
    /*CG*/ {-1.8, -2.9, -2.0, -1.7, -1.2, -1.9},
    /*GC*/ {-2.3, -3.4, -2.9, -2.1, -1.4, -2.1},
    /*UA*/ {-1.1, -2.3, -1.8,  -0.9,  -0.8, -1.1},
    /*GU*/ {-1.1, -2.1, -1.9, -1.0,  -0.4,  1.5},
    /*UG*/ { -0.8, -1.4, -1.2,  -0.5,  -0.2,  -0.4}
};


void Init_Param(struct Parameters *Param, double betaScale, double temperature, double EnergyThreshold, double lambda1, double lambda2, double T_start, double T_end, double Cooling){
    
    Param->Energy_Threshold = EnergyThreshold;
    read_parameter_file("rna_turner2004.par");// rna_turner1999.par
    Param->pf_parameters = get_scaled_pf_parameters();
    Param->pf_parameters->temperature = temperature;
    Param->KT = (betaScale*((temperature+K0)*GASCONST))/1000.; // in Kcal
    Param->pf_parameters->pf_scale = -1;//exp(-(-185+(Param->pf_parameters->temperature-37.)*7.27)/Param->pf_parameters->kT);
    Param->lambda1 = lambda1;
    Param->lambda2 = lambda2;

    Param->NcBasePair[0][0] = 'A';Param->NcBasePair[0][1] = 'U';
    Param->NcBasePair[1][0] = 'U';Param->NcBasePair[1][1] = 'A';
    Param->NcBasePair[2][0] = 'G';Param->NcBasePair[2][1] = 'C';
    Param->NcBasePair[3][0] = 'C';Param->NcBasePair[3][1] = 'G';
    //Param->NcBasePair[4][0] = 'U';Param->NcBasePair[4][1] = 'G';
    //Param->NcBasePair[5][0] = 'G';Param->NcBasePair[5][1] = 'U';
    
    Param->NuclBase[0] = 'A';
    Param->NuclBase[1] = 'U';
    Param->NuclBase[2] = 'G';
    Param->NuclBase[3] = 'C';
    
    
    Param->T_start = T_start;
    Param->T_end = T_end;
    Param->Cooling = Cooling;

}









