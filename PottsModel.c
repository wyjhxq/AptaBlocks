#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
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
#include "Initialization.h"
#include "Parameters.h"


int NT = 4;


double ProbWithStruct(struct SE_Design SED, struct Parameters Param){
    double Obj, enc, ec, pws_rna, pws_apta;
    //char *structtmp;
    int is_circular, is_constrained, calculate_bppm, length;
    
    // for aptamer_SE
    /* compute partition function with constrained structure */
    calculate_bppm = 0;
    is_constrained = 1;
    is_circular    = 0;
    ec = pf_fold_par(SED.AptamerInfo.SeqAll, SED.AptamerInfo.SpacerSEConst, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    /* compute partition function without constrained structure */
    is_constrained = 0;
    enc = pf_fold_par(SED.AptamerInfo.SeqAll, NULL, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    pws_apta = exp(-ec/Param.KT)/exp(-enc/Param.KT);
    
    // for RNA_SE
    calculate_bppm = 0;
    is_constrained = 1;
    is_circular    = 0;
    ec = pf_fold_par(SED.RNAInfo.SeqAll, SED.RNAInfo.SpacerSEConst, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    is_constrained = 0;
    enc = pf_fold_par(SED.RNAInfo.SeqAll, NULL, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    pws_rna = exp(-ec/Param.KT)/exp(-enc/Param.KT);
    
    
    Obj = pws_apta*pws_rna;// + lambda*pcg;
    
    //printf("Inside: %f %f %le\n", pws_apta, pws_rna, Obj);
    
    return Obj;
}

double ProbWithStruct_Detail(struct SE_Design SED, struct Parameters Param, double prob[2]){
    double Obj, enc, ec, pws_rna, pws_apta;
    //char *structtmp;
    int is_circular, is_constrained, calculate_bppm, length;
    
    // for aptamer_SE
    /* compute partition function with constrained structure */
    calculate_bppm = 0;
    is_constrained = 1;
    is_circular    = 0;
    ec = pf_fold_par(SED.AptamerInfo.SeqAll, SED.AptamerInfo.SpacerSEConst, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    /* compute partition function without constrained structure */
    is_constrained = 0;
    enc = pf_fold_par(SED.AptamerInfo.SeqAll, NULL, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    pws_apta = exp(-ec/Param.KT)/exp(-enc/Param.KT);
    
    // for RNA_SE
    calculate_bppm = 0;
    is_constrained = 1;
    is_circular    = 0;
    ec = pf_fold_par(SED.RNAInfo.SeqAll, SED.RNAInfo.SpacerSEConst, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    is_constrained = 0;
    enc = pf_fold_par(SED.RNAInfo.SeqAll, NULL, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    pws_rna = exp(-ec/Param.KT)/exp(-enc/Param.KT);
    
    
    Obj = pws_apta*pws_rna;// + lambda*pcg;
    
    prob[0] = pws_apta;
    prob[1] = pws_rna;
    
    //printf("Inside: %f %f %le\n", pws_apta, pws_rna, Obj);
    
    return Obj;
}

double ProbWithStruct_Print(struct SE_Design SED, struct Parameters Param){
    double Obj, enc, ec, pws_rna, pws_apta;
    //char *structtmp;
    int is_circular, is_constrained, calculate_bppm, length;
    
    // for aptamer_SE
    /* compute partition function with constrained structure */
    //printf("Aptamer: %s\n", SED.AptamerInfo.SeqAll);
    //printf("Aptamer: %s\n", SED.AptamerInfo.SpacerSEConst);
    calculate_bppm = 0;
    is_constrained = 1;
    is_circular    = 0;
    ec = pf_fold_par(SED.AptamerInfo.SeqAll, SED.AptamerInfo.SpacerSEConst, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    /* compute partition function without constrained structure */
    is_constrained = 0;
    enc = pf_fold_par(SED.AptamerInfo.SeqAll, NULL, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    pws_apta = exp(-ec/Param.KT)/exp(-enc/Param.KT);
    
    //printf("pws_apta %f\n", pws_apta);
    
    // for RNA_SE
    //printf("RNA: %s\n", SED.RNAInfo.SeqAll);
    //printf("RNA: %s\n", SED.RNAInfo.SpacerSEConst);
    calculate_bppm = 0;
    is_constrained = 1;
    is_circular    = 0;
    ec = pf_fold_par(SED.RNAInfo.SeqAll, SED.RNAInfo.SpacerSEConst, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    //printf("RNA: %s\n", SED.RNAInfo.SeqAll);
    //printf("RNA: %s\n", SED.RNAInfo.SpacerSEConst);
    
    
    is_constrained = 0;
    enc = pf_fold_par(SED.RNAInfo.SeqAll, NULL, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    pws_rna = exp(-ec/Param.KT)/exp(-enc/Param.KT);
    
    //printf("pws_rna %f\n", pws_rna);
    
    
    Obj = pws_apta*pws_rna;// + lambda*pcg;
    
    printf("Probability with structure constraints: %0.5f * %0.5f = %0.5f\n", pws_apta, pws_rna, Obj);
    
    return Obj;
}

// mapping base pair to enerogy index
// Ncp1 in seq 5'--3'
// Ncp2 in seq 3'--5'
int EnerogyIndexMappingRNARNA(char Ncp1, char Ncp2){
    int index;
    
    if ((Ncp1 == 'A' && Ncp2 == 'U') || Ncp1 == 'A' && Ncp2 == 'T') {
        index = 0;
        return index;
    }
    
    if (Ncp1 == 'C' && Ncp2 == 'G') {
        index = 1;
        return index;
    }
    
    if (Ncp1 == 'G' && Ncp2 == 'C') {
        index = 2;
        return index;
    }
    
    if ((Ncp1 == 'U' && Ncp2 == 'A') || (Ncp1 == 'T' && Ncp2 == 'A')) {
        index = 3;
        return index;
    }
    
    if ((Ncp1 == 'G' && Ncp2 == 'U') || (Ncp1 == 'G' && Ncp2 == 'T')) {
        index = 4;
        return index;
    }
    
    if ((Ncp1 == 'U' && Ncp2 == 'G') || (Ncp1 == 'T' && Ncp2 == 'G')) {
        index = 5;
        return index;
    }
    
    return -1; // no such binding
}


// compute the enerogy of the sticky end
double GetStickyEnerogyRNARNA(struct SE_Design SED){
    int length_sticky = SED.SEInfo.SELen;
    int i, basepairIndex0, basepairIndex1;
    char ncp1, ncp2;
    double Ebp, Sbp;
    
    Ebp = 0;
    Sbp = 4.09;//initfactor;
    //printf("len(SE): %d\n", length_sticky);
    //printf("sbp: %f\n", Sbp);
    for (i = 0; i < length_sticky; i++) {
        if (i > 0) {
            ncp1 = SED.AptamerInfo.SeqAll[SED.SEInfo.SEIndex[i].InApta];// AptaSeq->seq[GIndex[i].ind_apta];
            ncp2 = SED.RNAInfo.SeqAll[SED.SEInfo.SEIndex[i].InRNA]; //RNASeqs[0].seq[GIndex[i].ind_rna[0]];
            basepairIndex1 = EnerogyIndexMappingRNARNA(ncp1, ncp2);
            Sbp += StackingEnerogyRNARNA[basepairIndex0][basepairIndex1];
            //printf("2:%c|%c -- %d\n", ncp1, ncp2, basepairIndex1);
            //printf("e[%d][%d] %f, sbp: %f\n", basepairIndex0, basepairIndex1, StackingEnerogyRNARNA[basepairIndex0][basepairIndex1], Sbp);
            basepairIndex0 = basepairIndex1;
            
        }else{
            ncp1 = SED.AptamerInfo.SeqAll[SED.SEInfo.SEIndex[i].InApta]; //AptaSeq->seq[GIndex[i].ind_apta];
            ncp2 = SED.RNAInfo.SeqAll[SED.SEInfo.SEIndex[i].InRNA]; // RNASeqs[0].seq[GIndex[i].ind_rna[0]];
            basepairIndex0 = EnerogyIndexMappingRNARNA(ncp1, ncp2);
            //printf("1:%c|%c -- %d\n", ncp1, ncp2, basepairIndex0);
            //Ebp += BasepairEnerogy[basepairIndex0]; // base pair enerogy
        }
    }
    
    // compute the best possible energy
    // use turner 2004 model
    double LargestEnergy;
    LargestEnergy = 4.09 - 3.4 - 3.3*(length_sticky-2);
    //printf("largest energy is %f\n", LargestEnergy);
    
    if(Sbp > 0){
        return 0.000000001;
    }else{
        return Sbp/LargestEnergy;
    }

    
    
}


// seqd53 is a sequence in direction 5'--3' and seqd35 is a sequence in opposite direction 3'--5'
// the probability of not inter-interaction of the following case:
// 5'-seqd53-3' this is the aptamer or siRNA
//     \\
// 3'-seqd35-5' this is the sticky end
// only considering stacking binding at one position
double Prob_InterInteraction(char seqd53[], char seqd35[], double **PUseq53, int window_size, double kT, double energy_threshold){
    int len1, len2;
    int i, j, t;
    int pstacklen;
    int index0, index1;
    double *Zs; // partition function for stacking in [k, k+t, k*, k*+t]
    double Zi, Z, Zmax, *E;
    double interInitiG = 4.09;
    double energyt = 0;
    
    len1 = strlen(seqd53);
    len2 = strlen(seqd35);
    
    Zs = malloc(sizeof(double)*window_size);
    E = malloc(sizeof(double)*window_size);
    
    //printf("5'-%s-3'\n", seqd53);
    //printf("3'-%s-5'\n", seqd35);
    Zi = 0;
    Z  = 0;
    Zmax = 0;
    for (i = 0; i < len1; i++) {
        for (j = 0; j < len2; j++) {
            index0 = EnerogyIndexMappingRNARNA(seqd53[i], seqd35[j]);
            if (index0 != -10000) {
                Zs[0] = exp(-interInitiG / kT);
                E[0] = interInitiG;
            }else{
                Zs[0] = 1;
                E[0] = 0;
            }
            //printf("begin %d(%c) %d(%c): Zs[0]=%f E[0] = %f\n", i, seqd53[i], j, seqd35[j], Zs[0], E[0]);
            
            pstacklen = (len1-1 - i < window_size) ? len1-1 - i : window_size;
            pstacklen = (pstacklen < len2 -1 - j) ? pstacklen : len2-1 - j;
            //printf("pstacklen is %d\n", pstacklen);
            for (t = 1; t < pstacklen; t++) {
                index1 = EnerogyIndexMappingRNARNA(seqd53[i+t], seqd35[j+t]);
                //printf("index0 %d, index1 %d\n", index0, index1);
                if (index1 != -1 && index0 != -1) {
                    Zs[t] = Zs[t-1]*exp(-StackingEnerogyRNARNA[index0][index1] / kT)*PUseq53[i][t];
                    E[t] = E[t-1] + StackingEnerogyRNARNA[index0][index1];
                    energyt = -StackingEnerogyRNARNA[index0][index1];
                    index0 = index1;
                    if (E[t] > energy_threshold) {
                        Zi += Zs[t];
                    }
                    if (E[t] < Zmax) {
                        Zmax = E[t];
                    }
                    Z  += Zs[t];
                    //printf("stacking %d(%c) %d(%c): Zs[%d]=%f = %f * %f * %f  E[%d]=%f\n", i+t, seqd53[i+t], j+t, seqd35[j+t], t, Zs[t], Zs[t-1], exp(-StackingEnerogyRNARNA[index0][index1] / kT), PUseq53[i][t], t, E[t]);
                }else{
                    //printf("t is %d\n", t);
                    //Zs[t] = 1;
                    //E[t] = 0;
                    //Z += Zs[t];
                    //index0 = -1;
                    //printf("stacking %d(%c) %d(%c): Zs[%d]=%f E[%d]=%f\n", i+t, seqd53[i+t], j+t, seqd35[j+t], t, Zs[t], t, E[t]);
                    break;
                }
            }
            //printf("end!!!\n");
        }
    }
    Z = Z; //11/22 add
    
    free(Zs);
    free(E);
    //printf("Probability of not inter-interaction is %f=  1/%f\n", 1/Z, Z);
    
    //return 1/Z;
    
    if(Zi == Z ){
        return 1.0;
    }else{
        return Zi/Z;
    }
}



// seqd53 is a sequence in direction 5'--3' and seqd35 is a sequence in opposite direction 3'--5'
// the probability of not inter-interaction of the following case:
// 5'-seqd53-3' this is the aptamer or siRNA
//     \\
// 3'-seqd35-5' this is the sticky end
// only considering stacking binding at one position
double ProbHomoDimerSE(char seqd53[], char seqd35[],  double kT){
    int len1, len2;
    int i, j, t;
    int pstacklen;
    int index0, index1;
    double *Zs; // partition function for stacking in [k, k+t, k*, k*+t]
    double Zi, Z, Zmax, *E;
    double interInitiG = 4.09;
    double energyt = 0;
    int window_size = 0;
    
    len1 = strlen(seqd53);
    len2 = strlen(seqd35);
    window_size = len1;
    
    
    Zs = malloc(sizeof(double)*window_size);
    E = malloc(sizeof(double)*window_size);
    
    //printf("5'-%s-3'\n", seqd53);
    //printf("3'-%s-5'\n", seqd35);
    Zi = 0;
    Z  = 0;
    Zmax = 0;
    for (i = 0; i < len1; i++) {
        for (j = 0; j < len2; j++) {
            index0 = EnerogyIndexMappingRNARNA(seqd53[i], seqd35[j]);
            if (index0 != -10000) {
                Zs[0] = exp(-interInitiG / kT);
                E[0] = interInitiG;
            }else{
                Zs[0] = 1;
                E[0] = 0;
            }
            //printf("begin %d(%c) %d(%c): Zs[0]=%f E[0] = %f\n", i, seqd53[i], j, seqd35[j], Zs[0], E[0]);
            
            pstacklen = (len1-1 - i < window_size) ? len1-1 - i : window_size;
            pstacklen = (pstacklen < len2 -1 - j) ? pstacklen : len2-1 - j;
            //printf("pstacklen is %d\n", pstacklen);
            for (t = 1; t < pstacklen; t++) {
                index1 = EnerogyIndexMappingRNARNA(seqd53[i+t], seqd35[j+t]);
                //printf("index0 %d, index1 %d\n", index0, index1);
                if (index1 != -1 && index0 != -1) {
                    Zs[t] = Zs[t-1]*exp(-StackingEnerogyRNARNA[index0][index1] / kT);//*PUseq53[i][t];
                    E[t] = E[t-1] + StackingEnerogyRNARNA[index0][index1];
                    energyt = -StackingEnerogyRNARNA[index0][index1];
                    index0 = index1;
                    if (Zs[t] > Zmax) {
                        Zmax = Zs[t];
                    }
                    Z  += Zs[t];
                    //printf("stacking %d(%c) %d(%c): Zs[%d]=%f = %f * %f * %f  E[%d]=%f\n", i+t, seqd53[i+t], j+t, seqd35[j+t], t, Zs[t], Zs[t-1], exp(-StackingEnerogyRNARNA[index0][index1] / kT), PUseq53[i][t], t, E[t]);
                }else{
                    //printf("t is %d\n", t);
                    //Zs[t] = 1;
                    //E[t] = 0;
                    //Z += Zs[t];
                    //index0 = -1;
                    //printf("stacking %d(%c) %d(%c): Zs[%d]=%f E[%d]=%f\n", i+t, seqd53[i+t], j+t, seqd35[j+t], t, Zs[t], t, E[t]);
                    break;
                }
            }
            //printf("end!!!\n");
        }
    }
    Z = Z; //11/22 add
    
    free(Zs);
    free(E);
    //printf("Probability of not inter-interaction is %f=  1/%f\n", 1/Z, Z);
    
    return Zmax/Z;
}





// reverse the given null-terminated string
void strrev(char *p)
{
    char *q = p;
    while(q && *q) ++q;
    for(--q; p < q; ++p, --q)
    *p = *p ^ *q,
    *q = *p ^ *q,
    *p = *p ^ *q;
}

double Prob4Interact(struct SE_Design SED, struct Parameters Param){
    
    double prob = 1;
    char *SEseq;
    
    if (SED.AttachType == 1) {
        SEseq = malloc(sizeof(char)*(SED.SEInfo.SELen+1));
        strncpy(SEseq, SED.AptamerInfo.SeqAll+SED.SEInfo.SEIndex[0].InApta, SED.SEInfo.SELen);
        SEseq[SED.SEInfo.SELen] = '\0';
        strrev(SEseq);
        prob = Prob_InterInteraction(SED.RNAInfo.Seq, SEseq, SED.RNAInfo.ProbContextSeq, SED.kmer, Param.KT, Param.Energy_Threshold);
        prob = prob*Prob_InterInteraction(SED.AptamerInfo.AptaSeq, SEseq, SED.AptamerInfo.ProbContext, SED.kmer, Param.KT, Param.Energy_Threshold);
        //printf("prob in1: %f\n", prob);
        
        strncpy(SEseq, SED.RNAInfo.SeqAll+SED.SEInfo.SEIndex[SED.SEInfo.SELen-1].InRNA, SED.SEInfo.SELen);
        SEseq[SED.SEInfo.SELen] = '\0';
        strrev(SEseq);
        prob = prob*Prob_InterInteraction(SED.AptamerInfo.AptaSeq, SEseq, SED.AptamerInfo.ProbContext, SED.kmer, Param.KT, Param.Energy_Threshold);
        prob = prob*Prob_InterInteraction(SED.RNAInfo.Seq, SEseq, SED.RNAInfo.ProbContextSeq, SED.kmer, Param.KT, Param.Energy_Threshold);
        //printf("prob in2: %f\n", Prob_InterInteraction(SED.AptamerInfo.AptaSeq, SEseq, SED.AptamerInfo.ProbContext, SED.kmer, Param.KT, Param.Energy_Threshold));
        free(SEseq);
        
        
        return prob;
    }
    
    if (SED.AttachType == 2) {
        SEseq = malloc(sizeof(char)*(SED.SEInfo.SELen+1));
        strncpy(SEseq, SED.RNAInfo.SeqAll+SED.SEInfo.SEIndex[SED.SEInfo.SELen-1].InRNA, SED.SEInfo.SELen);
        SEseq[SED.SEInfo.SELen] = '\0';
        strrev(SEseq);
        
        prob = Prob_InterInteraction(SED.AptamerInfo.AptaSeq, SEseq, SED.AptamerInfo.ProbContext, SED.kmer, Param.KT, Param.Energy_Threshold);
        //printf("prob in1: %f\n", prob);
        prob = prob*Prob_InterInteraction(SED.RNAInfo.AntiSeq, SEseq, SED.RNAInfo.ProbContextAntiSeq, SED.kmer, Param.KT, Param.Energy_Threshold);
        //printf("prob in2: %f\n", Prob_InterInteraction(SED.RNAInfo.AntiSeq, SEseq, SED.RNAInfo.ProbContextAntiSeq, SED.kmer, Param.KT, Param.Energy_Threshold));
        free(SEseq);
        
        return prob;
    }
    
    prob  = -100;
    return prob;
}

double ProbHomoDimerBetweenSE(struct SE_Design SED, struct Parameters Param){
    double prob = 1;
    char *SEseq1, *SEseq2;
    
    SEseq1 = malloc(sizeof(char)*(SED.SEInfo.SELen+1));
    strncpy(SEseq1, SED.AptamerInfo.SeqAll+SED.SEInfo.SEIndex[0].InApta, SED.SEInfo.SELen);
    SEseq1[SED.SEInfo.SELen] = '\0';
    
    SEseq2 = malloc(sizeof(char)*(SED.SEInfo.SELen+1));
    strncpy(SEseq2, SED.RNAInfo.SeqAll+SED.SEInfo.SEIndex[SED.SEInfo.SELen-1].InRNA, SED.SEInfo.SELen);
    SEseq2[SED.SEInfo.SELen] = '\0';
    strrev(SEseq2);
    
    prob = ProbHomoDimerSE(SEseq1, SEseq2, Param.KT);

    return prob;
}


double Num_Consecutive(struct SE_Design SED, int ConsecutiveLen){
    int i, j;
    int count;
    double NumC=0.0;
    
    //aptamter
    for(i = 0; i < SED.AptamerInfo.SpacerSELen - ConsecutiveLen + 1; i++){
        count = 1;
        for(j = 1; j < ConsecutiveLen; j++){
            if(SED.AptamerInfo.SeqAll[SED.AptamerInfo.SpacerSEIndex[i+j]] == SED.AptamerInfo.SeqAll[SED.AptamerInfo.SpacerSEIndex[i]] ){
                count++;
            }else{
                break;
            }
        }
        if(count == ConsecutiveLen){
            NumC = NumC + 1.0;
        }
    }
    
    //rna
    for(i = 0; i < SED.RNAInfo.SpacerSELen - ConsecutiveLen + 1; i++){
        count = 1;
        for(j = 1; j < ConsecutiveLen; j++){
            if(SED.RNAInfo.SeqAll[SED.RNAInfo.SpacerSEIndex[i+j]] == SED.RNAInfo.SeqAll[SED.RNAInfo.SpacerSEIndex[i]] ){
                count++;
            }else{
                break;
            }
        }
        if(count == ConsecutiveLen){
            NumC = NumC + 1.0;
        }
    }
    
    if(NumC > 0){
        return -10000000;
    }else{
        return 0;
    }
    //return NumC;
}






double Num_Consecutive_ALL(struct SE_Design SED, int ConsecutiveLen){
    int i, j;
    int count;
    double NumC=0.0;
    
    //aptamter
    for(i = 0; i < SED.AptamerInfo.TotalLen - ConsecutiveLen + 1; i++){
        if(SED.AptamerInfo.SeqAll[i]!='A' && SED.AptamerInfo.SeqAll[i]!='U' && SED.AptamerInfo.SeqAll[i]!='G' && SED.AptamerInfo.SeqAll[i]!='C'){
            continue;
        }
        count = 1;
        for(j = 1; j < ConsecutiveLen; j++){
            if(SED.AptamerInfo.SeqAll[i+j] == SED.AptamerInfo.SeqAll[i] ){
                count++;
            }else{
                break;
            }
        }
        if(count == ConsecutiveLen){
            NumC = NumC + 1.0;
        }
    }
    
    //rna
    for(i = 0; i < SED.RNAInfo.TotalLen - ConsecutiveLen + 1; i++){
        count = 1;
        for(j = 1; j < ConsecutiveLen; j++){
            if(SED.RNAInfo.SeqAll[i+j] == SED.RNAInfo.SeqAll[i] ){
                count++;
            }else{
                break;
            }
        }
        if(count == ConsecutiveLen){
            NumC = NumC + 1.0;
        }
    }
    
    if(NumC > 0){
        return -NumC*10000000;
    }else{
        return 0;
    }
    //return NumC;
}


double GC_Content(struct SE_Design SED, double probgc[2]){
    int i;
    double pa, pr;
    double count;
    double gccv;
    
    //aptamer
    count = 0.0;
    for(i = 0; i < SED.AptamerInfo.SpacerSELen; i++){
        if(SED.AptamerInfo.SeqAll[SED.AptamerInfo.SpacerSEIndex[i]] == 'G' || SED.AptamerInfo.SeqAll[SED.AptamerInfo.SpacerSEIndex[i]] == 'C'){
            count  = count + 1.0;
        }
    }
    pa = count / (double)SED.AptamerInfo.SpacerSELen;
    
    //rna
    count = 0.0;
    for(i = 0; i < SED.RNAInfo.SpacerSELen; i++){
        if(SED.RNAInfo.SeqAll[SED.RNAInfo.SpacerSEIndex[i]] == 'G' || SED.RNAInfo.SeqAll[SED.RNAInfo.SpacerSEIndex[i]] == 'C'){
            count  = count + 1.0;
        }
    }
    pr = count / (double)SED.RNAInfo.SpacerSELen;
    
    probgc[0] = pa;
    probgc[1] = pr;
    
    if(probgc[0] == 0 || probgc[1] == 0){
        gccv = -100.0;
    }else{
        gccv = log(probgc[0])/log(2) + log(probgc[1])/log(2);
    }
    
    return gccv;
}

double AUGC_Content(struct SE_Design SED){
    int i;
    double pa, pr;
    double count;
    double augc;
    
    //aptamer
    double countA = 0.0;
    for(i = 0; i < SED.AptamerInfo.SpacerSELen; i++){
        if(SED.AptamerInfo.SeqAll[SED.AptamerInfo.SpacerSEIndex[i]] == 'A'){
            countA  = 1.0;
            break;
        }
    }
    
    double countU = 0.0;
    for(i = 0; i < SED.AptamerInfo.SpacerSELen; i++){
        if(SED.AptamerInfo.SeqAll[SED.AptamerInfo.SpacerSEIndex[i]] == 'U'){
            countU  = 1.0;
            break;
        }
    }
    
    double countG = 0.0;
    for(i = 0; i < SED.AptamerInfo.SpacerSELen; i++){
        if(SED.AptamerInfo.SeqAll[SED.AptamerInfo.SpacerSEIndex[i]] == 'G'){
            countG  = 1.0;
            break;
        }
    }
    
    double countC = 0.0;
    for(i = 0; i < SED.AptamerInfo.SpacerSELen; i++){
        if(SED.AptamerInfo.SeqAll[SED.AptamerInfo.SpacerSEIndex[i]] == 'C'){
            countC  = 1.0;
            break;
        }
    }
    
    if(countA + countU + countG + countC >=3){
        return 0;
    }else{
        return -100000000;
    }
    
    
}

double ObjctiveFunction(struct SE_Design SED, struct Parameters Param){
    double obj,objstruct, objenergy, objinter, numc, gcc[2], gccv, augc, objhomo;
    
    
    // probality with structure constraints
    objstruct = log(ProbWithStruct(SED, Param)) / log(2);
    
    // energy of the SE
    objenergy = log(GetStickyEnerogyRNARNA(SED))/log(2);
                            
    // porbaility of inter-interactions
    objinter = log(Prob4Interact(SED, Param))/log(2);
    
    // probability of homo dimer between SE
    //objhomo = log(ProbHomoDimerBetweenSE(SED, Param))/log(2);
    
    // number of consecutive identical
    numc  = Num_Consecutive_ALL(SED, 4);//Num_Consecutive(SED, 4);
    //printf("number of consective identical: %f\n", numc);
    
    // GC content
    //gccv = GC_Content(SED, gcc);
    //printf("gc content: %f\n", gccv);
    
    // AGUC content
    //augc = AUGC_Content(SED);
    
    // overall objective
    obj = objstruct + objinter + Param.lambda2*objenergy + numc; //+ numc + augc;// + 10*objhomo;// - numc + Param.lambda1*gccv;
    
    //printf("Inside: %f + %f + %f = %f\n", objstruct, objenergy, objinter, obj);
    
    return obj;
}


double ObjctiveFunctionM(struct SE_Design SED[], struct Parameters Param, int NumPairs){
    
    int i;
    double objm;
    double weight[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    
    objm = 0.0;
    for(i = 0; i < NumPairs; i++){
        objm += weight[i]*ObjctiveFunction(SED[i], Param);
    }
    
    return objm;
}


void PrintSingleSEDesign(struct SE_Design SED, struct Parameters Param){
    double obj,objstruct, objenergy, objinter, tmp, numc, gcc[2], gccv, augc;
    
    // print sequence
    printf("Aptamer: %s\n", SED.AptamerInfo.SeqAll);
    printf("RNA: %s\n", SED.RNAInfo.SeqAll);
    if(SED.AttachType == 2){
        printf("AntiRNA: %s\n", SED.RNAInfo.AntiSeq);
    }
    
    // probality with structure constraints
    tmp = ProbWithStruct_Print(SED, Param);
    objstruct = log(tmp) / log(2);
    printf("Probability with structure constraints: %0.3f log_2: %f\n", tmp, objstruct);
    
    // energy of the SE
    tmp = GetStickyEnerogyRNARNA(SED);
    double LargestEnergy = 4.09 - 3.4 - 3.3*(SED.SEInfo.SELen-2);
    objenergy = log(tmp)/log(2);
    printf("Energy of SE (converted in [0,1]): %f %f log_2: %f\n", tmp, tmp*LargestEnergy, objenergy);
    
    // probility of homo dimer between SE
    //tmp = ProbHomoDimerBetweenSE(SED, Param);
    //printf("Probability of homo dimer between SE is %f\n", tmp);
    
    // porbaility of inter-interactions
    tmp = Prob4Interact(SED, Param);
    objinter = log(tmp)/log(2);
    printf("Probability of inter-interactions: %f log_2: %f\n", tmp, objinter);
    
    // number consecutive
    numc  = Num_Consecutive_ALL(SED,4);//Num_Consecutive(SED, 4);
    printf("number of consective identical: %f\n", numc);
    
    // GC content
    gccv = GC_Content(SED, gcc);
    printf("GC content: gcc_apta=%f, gcc_rna=%f, log_2(gcc_apta*gcc_rna)=%f\n", gcc[0], gcc[1], gccv);
    
    // AUGC content
    //augc  = AUGC_Content(SED);
    //printf("number of AUGC appeared: %f\n", augc);
    
    // overall objective
    obj = objstruct + objinter + Param.lambda2*objenergy;// - numc + Param.lambda1*gccv;
    printf("Overall objectives: %f\n", obj);
    
}

void PrintSEDeisgnM(struct SE_Design SED[], struct Parameters Param, int NumPairs){
    int i;
    
    for(i = 0; i< NumPairs; i++){
        PrintSingleSEDesign(SED[i], Param);
    }
}


void PrintSingleSEDesignFile(struct SE_Design SED, struct Parameters Param, FILE *fh){
    double obj,objstruct, objenergy, objinter, tmp;
    double probdetail[2];
    double numc, gccv, gcc[2];
    
    // print sequence
    fprintf(fh, "Aptamer: %s\n", SED.AptamerInfo.SeqAll);
    fprintf(fh, "RNA: %s\n", SED.RNAInfo.SeqAll);
    if(SED.AttachType == 2){
        fprintf(fh, "AntiRNA: %s\n", SED.RNAInfo.AntiSeq);
    }
    
    // probality with structure constraints
    tmp = ProbWithStruct_Detail(SED, Param, probdetail);
    objstruct = log(tmp) / log(2);
    fprintf(fh, "Probability with structure constraints: %0.5f * %0.5f = %0.5f log_2: %f\n", probdetail[0], probdetail[1], tmp, objstruct);
    
    // energy of the SE
    tmp = GetStickyEnerogyRNARNA(SED);
    double LargestEnergy = 4.09 - 3.4 - 3.3*(SED.SEInfo.SELen-2);
    objenergy = log(tmp)/log(2);
    fprintf(fh, "Energy of SE (converted in [0,1]): %f %f log_2: %f\n", tmp, tmp*LargestEnergy, objenergy);
    
    // porbaility of inter-interactions
    tmp = Prob4Interact(SED, Param);
    objinter = log(tmp)/log(2);
    fprintf(fh, "Probability of inter-interactions: %f log_2: %f\n", tmp, objinter);
    
    // number consecutive
    numc  = Num_Consecutive_ALL(SED,4);//Num_Consecutive(SED, 4);
    fprintf(fh, "number of consective identical: %f\n", numc);
    
    // GC content
    gccv = GC_Content(SED, gcc);
    fprintf(fh, "GC content: gcc_apta=%f, gcc_rna=%f, log_2(gcc_apta*gcc_rna)=%f\n", gcc[0], gcc[1], gccv);
    
    // overall objective
    obj = objstruct + Param.lambda1*objinter + Param.lambda2*objenergy;// - numc + Param.lambda1*gccv;
    fprintf(fh, "Overall objectives: %f\n", obj);
    
}


void PrintSEDeisgnMFile(struct SE_Design SED[], struct Parameters Param, int NumPairs, FILE *fh){
    int i;
    
    for(i = 0; i< NumPairs; i++){
        PrintSingleSEDesignFile(SED[i], Param, fh);
    }
}


int GoodPair(char s, char t){
    if(s == 'A' && t == 'U') return 1;
    if(s == 'U' && t == 'A') return 1;
    
    if(s == 'G' && t == 'U') return 1;
    if(s == 'U' && t == 'G') return 1;
    
    if(s == 'C' && t == 'G') return 1;
    if(s == 'G' && t == 'C') return 1;
    
    return 0;
}

int validSEComplement(struct SE_Design SED[], int NumPairs){
    int i, j;
    
    for(j = 0; j < NumPairs; j++){
        for (i = 0; i < SED[i].SEInfo.SELen; i++) {
            if (GoodPair(SED[j].AptamerInfo.SeqAll[SED[j].SEInfo.SEIndex[i].InApta], SED[j].RNAInfo.SeqAll[SED[j].SEInfo.SEIndex[i].InRNA]) == 0) {
                return 0;
            }
        }
    }
    
    
    return 1;
}

void Reshuffle(int *index, int len_index, int *seed){
    int *array;
    int i;
    int temp;
    int randomIndex;
    array = malloc(sizeof(int)*len_index);
    
    
    for (i = 0; i < len_index; i++) {     // fill array
        array[i] = index[i];
    }
    
    for (i = 0; i < len_index; i++) {    // shuffle array
        temp = array[i];
        randomIndex = rand_r(seed) % len_index;
        
        array[i]           = array[randomIndex];
        array[randomIndex] = temp;
    }
    
    
    for (i = 0; i < len_index; i++) {    // print array
        index[i] = array[i];
    }
    
    free(array);
}


void SpinSeqM2(struct SE_Design SED[], int NumPairs, int loc, char aptaNc, char rnaNc){
    
    int i;
    
    if (aptaNc == 'A' && rnaNc == 'U') {
        for(i = 0; i < NumPairs; i++){
            SED[i].AptamerInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InApta] = 'A';
            SED[i].RNAInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InRNA] = 'U';
        }
    }
    else if (aptaNc == 'U' && rnaNc == 'A') {
        for(i = 0; i < NumPairs; i++){
            SED[i].AptamerInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InApta] = 'U';
            SED[i].RNAInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InRNA] = 'A';
        }
    }
    else if (aptaNc == 'U' && rnaNc == 'G') {
        for(i = 0; i < NumPairs; i++){
            SED[i].AptamerInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InApta] = 'U';
            SED[i].RNAInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InRNA] = 'G';
        }
    }
    else if (aptaNc == 'C' && rnaNc == 'G') {
        for(i = 0; i < NumPairs; i++){
            SED[i].AptamerInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InApta] = 'C';
            SED[i].RNAInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InRNA] = 'G';
        }
    }
    else if (aptaNc == 'G' && rnaNc == 'C') {
        for(i = 0; i < NumPairs; i++){
            SED[i].AptamerInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InApta] = 'G';
            SED[i].RNAInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InRNA] = 'C';
        }
    }
    else if (aptaNc == 'G' && rnaNc == 'U') {
        for(i = 0; i < NumPairs; i++){
            SED[i].AptamerInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InApta] = 'G';
            SED[i].RNAInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InRNA] = 'U';
        }
    }
    else if (aptaNc == 'T' && rnaNc == 'A'){
        for(i = 0; i < NumPairs; i++){
            SED[i].AptamerInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InApta] = 'T';
            SED[i].RNAInfo.SeqAll[SED[i].SEInfo.SEIndex[loc].InRNA] = 'A';
        }
    }
    else {
        printf("Problem!!! in SpinSeqM2!!!\n");
        exit(1);
    }
    
}



// index to NcBasePair: 0: A-U 1: U-A 2: U-G 3: C-G 4: G-C 5:G-U
int Roulette_wheel(double *Prob, int Num){
    double rndNumber = rand() /( (double) RAND_MAX + 1.0);
    double offset = 0.0;
    int pick = 0;
    int i;
    
    for (i = 0; i < Num; i++) {
        offset += Prob[i];
        if (rndNumber < offset) {
            pick = i;
            break;
        }
    }
    
    return pick;
}


void Get_wheel_probability_SE_M(struct SE_Design SED[], struct Parameters Param, int NumPairs, int loc, double Temperature, double Prob[NT]) {
    double OBJ, OBJ0;
    int i;
    
    OBJ0 = ObjctiveFunctionM(SED, Param, NumPairs);
    //printf("%f\n", OBJ0);
    
    char NctC, NctC1;
    double DeltaO[NT];
    double SumD=0;
    NctC = SED[0].AptamerInfo.SeqAll[SED[0].SEInfo.SEIndex[loc].InApta];//SeqApta->seq[GroupIndex[loc].ind_apta];// sequence1[PairIndex[loc].ind1];
    NctC1 = SED[0].RNAInfo.SeqAll[SED[0].SEInfo.SEIndex[loc].InRNA];//SeqRNAs[0].seq[GroupIndex[loc].ind_rna[0]];
    for (i = 0; i < NT; i++) {
        //printf("*************LOCATION:%d (%d %d)**********\n", loc, SED[0].SEInfo.SEIndex[loc].InApta, SED[0].SEInfo.SEIndex[loc].InRNA);
        //printf("%s\n%s\n", SED[0].AptamerInfo.SeqAll, SED[0].RNAInfo.SeqAll);
        if (NctC == Param.NcBasePair[i][0] && NctC1 == Param.NcBasePair[i][1]) {
            DeltaO[i] = exp(0);
            //printf("%c -- %c\n", NctC, NctC1);
            //printf("Delta0[%d]=%f\n", i, DeltaO[i]);
            //SpinSeqM(SeqApta, SeqRNAs, GroupIndex, NumComb, loc, NcBasePair[i][0], NcBasePair[i][1]);
        } else{
            SpinSeqM2(SED, NumPairs, loc, Param.NcBasePair[i][0], Param.NcBasePair[i][1]);
            OBJ = ObjctiveFunctionM(SED, Param, NumPairs);
            DeltaO[i] = exp((OBJ - OBJ0)/Temperature);
            //printf("Delta0[%d]=%f, obj = %f\n", i, DeltaO[i], OBJ);
            //printf("NcPair: %c %c\n", NcBasePair[i][0], NcBasePair[i][1]);
            //PrintSeqs(SeqApta, SeqRNAs, NumComb);
        }
        SumD += DeltaO[i];
    }
    
    for (i = 0; i < NT; i++) {
        Prob[i] = DeltaO[i]/SumD;
        //printf("%f ", Prob[i]);
    }
    //printf("\n");
    
    
}




double ProbWithStruct_Single(char *seq, char *structure, struct Parameters Param){
    double Obj, enc, ec, pws_rna, pws_apta;
    int is_circular, is_constrained, calculate_bppm, length;
    double prob;
    
    /* compute partition function with constrained structure */
    calculate_bppm = 0;
    is_constrained = 1;
    is_circular    = 0;
    ec = pf_fold_par(seq, structure, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    /* compute partition function without constrained structure */
    is_constrained = 0;
    enc = pf_fold_par(seq, NULL, Param.pf_parameters, calculate_bppm, is_constrained, is_circular);
    free_pf_arrays();
    
    prob = exp(-ec/Param.KT)/exp(-enc/Param.KT);
    
    return prob;
}


double Num_Consecutive_Seq(char *Seq, int *Index, int Len, int ConsecutiveLen){
    int i, j;
    int count;
    double NumC=0.0;
    
    for(i = 0; i < Len - ConsecutiveLen + 1; i++){
        count = 1;
        for(j = 1; j < ConsecutiveLen; j++){
            if(Seq[Index[i+j]] == Seq[Index[i]] ){
                count++;
            }else{
                break;
            }
        }
        if(count == ConsecutiveLen){
            NumC = NumC + 1.0;
        }
    }
    
    return NumC;
}


double GC_Content_Seq(char *Seq, int *Index, int Len){
    int i;
    double pa, pr;
    double count;
    double gccv;
    
    count = 0.0;
    for(i = 0; i < Len; i++){
        if(Seq[Index[i]] == 'G' || Seq[Index[i]] == 'C'){
            count  = count + 1.0;
        }
    }
    pa = count / (double)Len;
    
    if(pa == 0){
        gccv = -100.0;
    }else{
        gccv = log(pa)/log(2);
    }
    
    return gccv;
}


void Get_wheel_probability_Spacer_Aptamer(struct Aptamer_Info *AptamerInfo, struct Parameters Param, int loc, double Temperature, double *Prob) {
    double OBJ, OBJ0;
    int i;
    double numc, gccv, gcc[2];
    
    // number of consecutive identical
    numc  = Num_Consecutive_Seq(AptamerInfo->SeqAll, AptamerInfo->SpacerIndex, AptamerInfo->SpacerLen, 4);
    //printf("number of consective identical: %f\n", numc);
    
    // GC content
    gccv = GC_Content_Seq(AptamerInfo->SeqAll, AptamerInfo->SpacerIndex, AptamerInfo->SpacerLen);
    
    OBJ0 = ProbWithStruct_Single(AptamerInfo->SeqAll, AptamerInfo->SpacerSEConst, Param) - numc + Param.lambda1*gccv;
    //printf("%f\n", OBJ0);
    
    char NctC, NctC1;
    double DeltaO[NT];
    double SumD=0;
    NctC = AptamerInfo->SeqAll[AptamerInfo->SpacerIndex[loc]];
    //printf("%s\n", AptamerInfo->SeqAll);
    for(i = 0; i < 4; i++){
        if(NctC == Param.NuclBase[i]){
            DeltaO[i] = exp(0);
        }else{
            AptamerInfo->SeqAll[AptamerInfo->SpacerIndex[loc]] = Param.NuclBase[i];
            // number of consecutive identical
            numc  = Num_Consecutive_Seq(AptamerInfo->SeqAll, AptamerInfo->SpacerIndex, AptamerInfo->SpacerLen, 4);
            // GC content
            gccv = GC_Content_Seq(AptamerInfo->SeqAll, AptamerInfo->SpacerIndex, AptamerInfo->SpacerLen);
            //printf("%s\n", AptamerInfo->SeqAll);
            OBJ = ProbWithStruct_Single(AptamerInfo->SeqAll, AptamerInfo->SpacerSEConst, Param)- numc + Param.lambda1*gccv;
            DeltaO[i] = exp((OBJ - OBJ0)/Temperature);
            //printf("%c %f\n", Param.NuclBase[i], DeltaO[i]);
        }
        SumD += DeltaO[i];
    }
    
    for (i = 0; i < 4; i++) {
        Prob[i] = DeltaO[i]/SumD;
        //printf("%f ", Prob[i]);
    }
    //printf("\n");
    
}

void Get_wheel_probability_Spacer_RNA(struct RNA_Info *RNAInfo, struct Parameters Param, int loc, double Temperature, double *Prob) {
    double OBJ, OBJ0;
    int i;
    double numc, gccv, gcc[2];
    
    // number of consecutive identical
    numc  = Num_Consecutive_Seq(RNAInfo->SeqAll, RNAInfo->SpacerIndex, RNAInfo->SpacerLen, 4);
    // GC content
    gccv = GC_Content_Seq(RNAInfo->SeqAll, RNAInfo->SpacerIndex, RNAInfo->SpacerLen);
    
    OBJ0 = ProbWithStruct_Single(RNAInfo->SeqAll, RNAInfo->SpacerSEConst, Param) - numc + Param.lambda1*gccv;
    //printf("%f\n", OBJ0);
    
    char NctC, NctC1;
    double DeltaO[NT];
    double SumD=0;
    NctC = RNAInfo->SeqAll[RNAInfo->SpacerIndex[loc]];
    //printf("%s\n", RNAInfo->SeqAll);
    for(i = 0; i < 4; i++){
        if(NctC == Param.NuclBase[i]){
            DeltaO[i] = exp(0);
        }else{
            RNAInfo->SeqAll[RNAInfo->SpacerIndex[loc]] = Param.NuclBase[i];
            // number of consecutive identical
            numc  = Num_Consecutive_Seq(RNAInfo->SeqAll, RNAInfo->SpacerIndex, RNAInfo->SpacerLen, 4);
            // GC content
            gccv = GC_Content_Seq(RNAInfo->SeqAll, RNAInfo->SpacerIndex, RNAInfo->SpacerLen);
            OBJ = ProbWithStruct_Single(RNAInfo->SeqAll, RNAInfo->SpacerSEConst, Param)- numc + Param.lambda1*gccv;
            DeltaO[i] = exp((OBJ - OBJ0)/Temperature);
            //printf("%c %f\n", Param.NuclBase[i], DeltaO[i]);
        }
        SumD += DeltaO[i];
    }
    
    for (i = 0; i < 4; i++) {
        Prob[i] = DeltaO[i]/SumD;
        //printf("%f ", Prob[i]);
    }
    //printf("\n");
    
}


// need to test!!!!!
void Get_Optimal_Apta_Spacer(struct Aptamer_Info *AptamerInfo, struct Parameters Param){
    int i, k;
    char t[3];
    double Best, ProbTemp;
    unsigned long long x, y, bestx;
    
    k = AptamerInfo->SpacerLen;
    Best = 0;
    for (x = 0; x < 1ULL<<(2*k); ++x) {
        //printf("%s\n", AptamerInfo->SeqAll);
        //printf("%s\n", AptamerInfo->SpacerConst);
        for (i = 0, y = x; i < k; ++i, y >>= 2){
            AptamerInfo->SeqAll[AptamerInfo->SpacerIndex[i]]=(char)"ACGT"[y&3];
        }
        if(ProbTemp > Best){
            Best = ProbTemp;
            bestx = x;
            //printf("BestX: %d, BestScore: %f\n", bestx, Best);
            //printf("%s\n", AptamerInfo->SeqAll);
            //printf("%s\n", AptamerInfo->SpacerConst);
            ProbTemp = ProbWithStruct_Single(AptamerInfo->SeqAll, AptamerInfo->SpacerSEConst, Param);
        }
    }
    
    for (i = 0, y = bestx; i < k; ++i, y >>= 2){
        AptamerInfo->SeqAll[AptamerInfo->SpacerIndex[i]]=(char)"ACGT"[y&3];
    }
    
}

// need to test!!!!!
void Get_Optimal_RNA_Spacer(struct RNA_Info *RNAInfo, struct Parameters Param){
    int i, k;
    char t[3];
    double Best, ProbTemp;
    unsigned long long x, y, bestx;
    
    k = RNAInfo->SpacerLen;
    Best = 0;
    for (x = 0; x < 1ULL<<(2*k); ++x) {
        for (i = 0, y = x; i < k; ++i, y >>= 2){
            RNAInfo->SeqAll[RNAInfo->SpacerIndex[i]]=(char)"ACGT"[y&3];
        }
        ProbTemp = ProbWithStruct_Single(RNAInfo->SeqAll, RNAInfo->SpacerSEConst, Param);
        if(ProbTemp > Best){
            Best = ProbTemp;
            bestx = x;
            //printf("BestX: %d, BestScore: %f\n", bestx, Best);
            //printf("%s\n", RNAInfo->SeqAll);
            //printf("%s\n", RNAInfo->SpacerConst);
        }
    }
    
    for (i = 0, y = bestx; i < k; ++i, y >>= 2){
        RNAInfo->SeqAll[RNAInfo->SpacerIndex[i]]=(char)"ACGT"[y&3];
    }
    
}

// *** need to be finish
void Get_Optimal_Spacer_M(struct SE_Design SED[], struct Parameters Param, int NumPairs){
    int i;
    
    // spacer for aptamer for different SEDesign
    for(i = 0; i < NumPairs; i++){
        if(SED[i].AptamerInfo.SpacerLen != 0){
            Get_Optimal_Apta_Spacer(&SED[i].AptamerInfo, Param);
        }
    }
    
    // spacer for RNA for different SEDesign
    for(i = 0; i < NumPairs; i++){
        if(SED[i].RNAInfo.SpacerLen != 0){
            Get_Optimal_RNA_Spacer(&SED[i].RNAInfo, Param);
        }
    }
}


void Sweep_SE_M(struct SE_Design SED[], struct Parameters Param, int NumPairs, double Tc){
    int i, j, pick, loc;
    double Prob[NT];
    int selen = SED[0].SEInfo.SELen;
    int *IINDEX;
    int s, pid, seed;
    
    s           = time(NULL);
    seed        = abs(((s*181)*((pid-83)*359))%104729);
    
    IINDEX = (int *)malloc(sizeof(int)*selen);
    for (i = 0; i < selen; i++) {
        IINDEX[i] = i;
    }
    
    Reshuffle(IINDEX, selen, &seed);
    for (j = 0; j < selen; j++) {
        loc=IINDEX[j];
        //printf("loc:%d\n", loc);
        Get_wheel_probability_SE_M(SED, Param, NumPairs, loc, Tc, Prob);
        pick = Roulette_wheel(Prob, NT);
        //printf("pick:%d\n",pick);
        SpinSeqM2(SED, NumPairs, loc, Param.NcBasePair[pick][0], Param.NcBasePair[pick][1]);
        //printf("%s|%s\n", seqt1, seqt2);
    }
    
    free(IINDEX);
}


void Sweep_Spacer_Aptamer(struct Aptamer_Info *AptamerInfo, struct Parameters Param, double Tc){
    int i, j, pick, loc;
    double Prob[4];
    int selen = AptamerInfo->SpacerLen;
    int *IINDEX;
    int s, pid, seed;
    
    s           = time(NULL);
    seed        = abs(((s*181)*((pid-83)*359))%104729);
    
    IINDEX = (int *)malloc(sizeof(int)*selen);
    for (i = 0; i < selen; i++) {
        IINDEX[i] = i;
    }
    
    Reshuffle(IINDEX, selen, &seed);
    for (j = 0; j < selen; j++) {
        loc=IINDEX[j];
        //printf("loc:%d\n", loc);
        Get_wheel_probability_Spacer_Aptamer(AptamerInfo, Param, loc, Tc, Prob);
        pick = Roulette_wheel(Prob, 4);
        AptamerInfo->SeqAll[AptamerInfo->SpacerIndex[loc]] = Param.NuclBase[pick];
    }
    
    free(IINDEX);
}

void Sweep_Spacer_RNA(struct RNA_Info *RNAInfo, struct Parameters Param, double Tc){
    int i, j, pick, loc;
    double Prob[4];
    int selen = RNAInfo->SpacerLen;
    int *IINDEX;
    int s, pid, seed;
    
    s           = time(NULL);
    seed        = abs(((s*181)*((pid-83)*359))%104729);
    
    IINDEX = (int *)malloc(sizeof(int)*selen);
    for (i = 0; i < selen; i++) {
        IINDEX[i] = i;
    }
    
    Reshuffle(IINDEX, selen, &seed);
    for (j = 0; j < selen; j++) {
        loc=IINDEX[j];
        //printf("loc:%d\n", loc);
        Get_wheel_probability_Spacer_RNA(RNAInfo, Param, loc, Tc, Prob);
        pick = Roulette_wheel(Prob, 4);
        RNAInfo->SeqAll[RNAInfo->SpacerIndex[loc]] = Param.NuclBase[pick];
    }
    
    free(IINDEX);
}



void Sweep_Spacer_M(struct SE_Design SED[], struct Parameters Param, int NumPairs, double Tc){
    int i;
    
    // spacer for aptamer for different SEDesign
    for(i = 0; i < NumPairs; i++){
        if(SED[i].AptamerInfo.SpacerLen != 0){
            Sweep_Spacer_Aptamer(&SED[i].AptamerInfo, Param, Tc);
        }
    }
    
    // spacer for RNA for different SEDesign
    for(i = 0; i < NumPairs; i++){
        if(SED[i].RNAInfo.SpacerLen != 0){
            Sweep_Spacer_RNA(&SED[i].RNAInfo, Param, Tc);
        }
    }
}


void Potts_ModelM(struct SE_Design SED[], struct Parameters Param, int NumPairs, FILE *fh){
    double Tc, ObjTc, Prob[NT];
    int Time_sweep, loc, length1, length2;
    int i, j, pick;
    int selen;
    int count;
    int seed, pid, s;
    
    srand(time(NULL));
    int rint = rand();
    
    
    
    Time_sweep  =    5; //
    loc         =    0; // current location in the sequence
    pick        =    0;
    selen      = SED[0].SEInfo.SELen;
    s           = time(NULL);
    seed        = abs(((s*181)*((pid-83)*359))%104729);
    //printf("1:%s|%s\n", seqt1, seqt2);
    
    
    // find the length and poitions of the uncertain part and find their correspondences
    if (validSEComplement(SED, NumPairs)) {
        
    }else{
        printf("lengths of the sticky ends are different or SE are not complement!!! (check the input sequences!!!)\n");
        exit(0);
    }
    
    
    double obj00;
    obj00=ObjctiveFunctionM(SED, Param, NumPairs);
    printf("!!!obj00:%f\n", obj00);
    fflush(stdout);
    
    // simulated annealing
    count = 0;
    Tc = Param.T_start;
    while (Tc > Param.T_end) {
        //printf("%s|%s\n", seqt1, seqt2);
        // For Tc, sweep sequence
        for (i = 0; i < Time_sweep; i++) {
            
            rint = rand()%2;
            //printf("rint: %d\n", rint);
            switch (rint) {
                case 0:
                    Sweep_SE_M(SED, Param, NumPairs, Tc);
                    break;
                case 1:
                    Sweep_Spacer_M(SED, Param, NumPairs, Tc);
                    break;
                    
                default:
                    printf("rint error!!!\n");
                    break;
            }
            
        }
        //PrintSeqs(SeqApta, SeqRNAs, NumComb);
        //PrintIndividualOBJ(SeqApta, SeqRNAs, GroupIndex, NumComb, parameters, kT);
        ObjTc = ObjctiveFunctionM(SED, Param, NumPairs);
        printf("T %f Obj %f \n", Tc, ObjTc);
        if(isnan(ObjTc)){
            PrintSEDeisgnM(SED, Param, NumPairs);
        }
        fflush(stdout);
        
        if (count++ == 50) {
            count = 0;
            PrintSEDeisgnM(SED, Param, NumPairs);
        }
        
        Tc = Tc*Param.Cooling;
    }
    
    
    // final results
    PrintSEDeisgnM(SED, Param, NumPairs);
    PrintSEDeisgnMFile(SED, Param, NumPairs, fh);
    
    // print seq result:
    //PrintSeqs(SeqApta, SeqRNAs, NumComb);
    //PrintSeqsFile(SeqApta, SeqRNAs, NumComb, fh);
    
    
}











































