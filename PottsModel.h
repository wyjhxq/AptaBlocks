#ifndef PottsModel_h
#define PottsModel_h

double ProbWithStruct(struct SE_Design SED, struct Parameters Param);
double GetStickyEnerogyRNARNA(struct SE_Design SED);
double Prob4Interact(struct SE_Design SED, struct Parameters Param);
double ObjctiveFunction(struct SE_Design SED, struct Parameters Param);
void PrintSingleSEDesign(struct SE_Design SED, struct Parameters Param);
void Potts_ModelM(struct SE_Design SED[], struct Parameters Param, int NumPairs, FILE *fh);
void Get_Optimal_Spacer_M(struct SE_Design SED[], struct Parameters Param, int NumPairs);
void PrintSingleSEDesignFile(struct SE_Design SED, struct Parameters Param, FILE *fh);

#endif
