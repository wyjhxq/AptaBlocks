#include "Initialization.h"
#include  "part_func_up.h"

int Length_Seq(char *structure){
    int len;
    int i;
    int count = 0;
    
    len = strlen(structure);
    for(i = 0; i < len; i++){
        if(structure[i] == '*'){
            count++;
        }
    }
    
    return count;
}

int Length_Spacer_SE(char *structure){
    int len;
    int i;
    int count = 0;
    
    len = strlen(structure);
    for(i = 0; i < len; i++){
        if(structure[i] != '*'){
            count++;
        }
    }
    
    return count;
}

int Length_Sapcer_Seq(char *structure){
    int len;
    int i;
    int count = 0;
    
    len = strlen(structure);
    for(i = 0; i < len; i++){
        if(structure[i] == '.'){
            count++;
        }
    }
    
    return count;
}

int GetOSeq(char *seq, char *seqr, char *oseq1){
    int i;
    char seqt[1000];
    int count, ind;
    
    if (strlen(seq) != strlen(seqr)) {
        printf("#Error: seq and seqr do not have the same length!!! %s(%d) vs %s(%d)\n", seq, strlen(seq), seqr, strlen(seqr));
        exit(1);
    }
    
    count = 0;
    for (i = 0; seqr[i] != '\0'; i++) {
        if (seqr[i] == '*' && seqr[i] == seqr[i+1]) {
            seqt[count] = seq[i];
            count++;
        }else{
            if (seqr[i] == '*') {
                seqt[count] = seq[i];
                count++;
                seqt[count] = '\0';
            }
        }
    }
    
    strcpy(oseq1, seqt);
    
    return count;
}

void UnpairProb4InterInteractionComputation(char seq[], double **P, int kmer){
    int i,j;
    pu_contrib *PU;
    PU = get_pu_contrib_struct((unsigned)strlen(seq), (unsigned)kmer);
    pf_fold(seq, NULL);
    PU = pf_unstru(seq, kmer);
    for (i = 0; i < strlen(seq); i++) {
        for (j = 0; j < kmer; j++) {
            P[i][j] = PU->H[i+1][j]+PU->I[i+1][j]+PU->M[i+1][j]+PU->E[i+1][j];
            //printf("PU[%d][%d]=%f ", i, j, P[i][j]);
        }
        //printf("\n");
    }
    
}

void GetConstraint(char *seq, char *structure){
    int len;
    int i;
    
    strcpy(seq, structure);
    len = strlen(structure);
    
    for(i = 0; i < len; i++){
        if(seq[i] != 'x'){
            seq[i] = '.';
        }
    }
}

void GetSapcerConstraint(char *seq, char *structure){
    int len;
    int i;
    
    strcpy(seq, structure);
    len = strlen(structure);
    
    for(i = 0; i < len; i++){
        if(seq[i] != '.'){
            seq[i] = '.';
        }else{
            seq[i] = 'x';
        }
    }
}

void CombineSpacerSEConstraint(char *seq, char *structure1, char *structure2){
    int len;
    int i;
    
    len = strlen(structure1);
    for(i = 0; i < len; i++){
        if(structure1[i] == 'x' || structure2[i] == 'x'){
            seq[i] = 'x';
        }else{
            seq[i] = '.';
        }
    }
}

void GetSpacerIndex(int *Index, char *structure){
    int i;
    int len;
    
    len = strlen(structure);
    
    int count = 0;
    for(i = 0; i < len; i++){
        if(structure[i] == '.'){
            Index[count] = i;
            count++;
        }
    }
}

void GetSpacerSEIndex(int *Index, char *structure){
    int i;
    int len;
    
    len = strlen(structure);
    
    int count = 0;
    for(i = 0; i < len; i++){
        if(structure[i] == 'x'){
            Index[count] = i;
            count++;
        }
    }
}

void Init_Aptamer(struct Aptamer_Info *AptamerInfo, char *AptamerSeq, char *AptamerStruct, int kmer){
    int i;
    
    AptamerInfo->TotalLen = strlen(AptamerSeq);
    AptamerInfo->SeqLen = Length_Seq(AptamerStruct);
    AptamerInfo->SpacerLen = Length_Sapcer_Seq(AptamerStruct);
    AptamerInfo->SpacerSELen = Length_Spacer_SE(AptamerStruct);
    
    printf("AllLEN %d, AptaLEN %d\n", AptamerInfo->TotalLen, AptamerInfo->SeqLen);
    
    // make space for struct Aptamer_Info
    AptamerInfo->SeqAll = (char *)malloc(sizeof(char)*(AptamerInfo->TotalLen+1));
    AptamerInfo->SeqStruct = (char *)malloc(sizeof(char)*(AptamerInfo->TotalLen+1));
    AptamerInfo->SeqConst = (char *)malloc(sizeof(char)*(AptamerInfo->TotalLen+1));
    AptamerInfo->SpacerConst = (char *)malloc(sizeof(char)*(AptamerInfo->TotalLen+1));
    AptamerInfo->SpacerSEConst = (char *)malloc(sizeof(char)*(AptamerInfo->TotalLen+1));
    AptamerInfo->AptaSeq = (char *)malloc(sizeof(char)*(AptamerInfo->SeqLen+1));
    AptamerInfo->SpacerIndex = (int *)malloc(sizeof(int)*AptamerInfo->SpacerLen);
    AptamerInfo->SpacerSEIndex = (int *)malloc(sizeof(int)*AptamerInfo->SpacerSELen);
    
    // close the string
    AptamerInfo->SeqAll[AptamerInfo->TotalLen] = '\0';
    AptamerInfo->SeqStruct[AptamerInfo->TotalLen] = '\0';
    AptamerInfo->SeqConst[AptamerInfo->TotalLen] = '\0';
    AptamerInfo->SpacerConst[AptamerInfo->TotalLen] = '\0';
    AptamerInfo->SpacerSEConst[AptamerInfo->TotalLen] = '\0';
    
    // copy sequences
    strcpy(AptamerInfo->SeqAll, AptamerSeq);
    strcpy(AptamerInfo->SeqStruct, AptamerStruct);
    GetConstraint(AptamerInfo->SeqConst, AptamerStruct);
    GetSapcerConstraint(AptamerInfo->SpacerConst, AptamerStruct);
    GetOSeq(AptamerSeq, AptamerStruct, AptamerInfo->AptaSeq);
    CombineSpacerSEConstraint(AptamerInfo->SpacerSEConst, AptamerInfo->SeqConst, AptamerInfo->SpacerConst);
    printf("%s\n%s\n%s\n%s\n%s\n", AptamerInfo->SeqAll, AptamerInfo->SeqStruct, AptamerInfo->SeqConst, AptamerInfo->SpacerConst, AptamerInfo->SpacerSEConst);
    
    //spacer index
    GetSpacerIndex(AptamerInfo->SpacerIndex, AptamerInfo->SeqStruct);
    for(i = 0; i < AptamerInfo->SpacerLen; i++){
        printf("%d ", AptamerInfo->SpacerIndex[i]);
    }
    printf("\n");
    GetSpacerSEIndex(AptamerInfo->SpacerSEIndex, AptamerInfo->SpacerSEConst);
    for(i = 0; i < AptamerInfo->SpacerSELen; i++){
        printf("%d ", AptamerInfo->SpacerSEIndex[i]);
    }
    printf("\n");
    
    
    // probability context
    AptamerInfo->ProbContext = malloc(AptamerInfo->SeqLen*sizeof(double*));
    for (i = 0; i < AptamerInfo->SeqLen; i++) {
        AptamerInfo->ProbContext[i] = malloc(kmer*sizeof(double));
    }
    UnpairProb4InterInteractionComputation(AptamerInfo->AptaSeq, AptamerInfo->ProbContext, kmer);
}

void FreeAptamerInfo(struct Aptamer_Info *AptamerInfo){
    free(AptamerInfo->SeqAll);
    free(AptamerInfo->SeqStruct);
    free(AptamerInfo->SeqConst);
    free(AptamerInfo->SpacerConst);
    free(AptamerInfo->SpacerSEConst);
    free(AptamerInfo->AptaSeq);
    free(AptamerInfo->SpacerIndex);
    free(AptamerInfo->SpacerSEIndex);
    
    int i;
    for(i = 0; i < AptamerInfo->SeqLen; i++){
        free(AptamerInfo->ProbContext[i]);
    }
    free(AptamerInfo->ProbContext);
}


void Init_RNASingle(struct RNA_Info *RNAInfo, char *RNASeq, char *RNAStruct, int kmer){
    int i;
    
    RNAInfo->TotalLen = strlen(RNASeq);
    RNAInfo->SeqLen = Length_Seq(RNAStruct);
    RNAInfo->SpacerLen = Length_Sapcer_Seq(RNAStruct);
    RNAInfo->SpacerSELen = Length_Spacer_SE(RNAStruct);
    
    printf("AllLEN %d, RNALEN %d\n", RNAInfo->TotalLen, RNAInfo->SeqLen);
    
    // make space
    RNAInfo->SeqAll = (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->SeqStruct = (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->SeqConst= (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->SpacerConst= (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->SpacerSEConst= (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->Seq = (char *)malloc(sizeof(char)*(RNAInfo->SeqLen+1));
    RNAInfo->SpacerIndex = (int *)malloc(sizeof(int)*RNAInfo->SpacerLen);
    RNAInfo->SpacerSEIndex = (int *)malloc(sizeof(int)*RNAInfo->SpacerSELen);
    
    // close the string
    RNAInfo->SeqAll[RNAInfo->TotalLen] = '\0';
    RNAInfo->SeqStruct[RNAInfo->TotalLen] = '\0';
    RNAInfo->SeqConst[RNAInfo->TotalLen] = '\0';
    RNAInfo->SpacerConst[RNAInfo->TotalLen] = '\0';
    RNAInfo->SpacerSEConst[RNAInfo->TotalLen] = '\0';
    
    // copy sequence
    strcpy(RNAInfo->SeqAll, RNASeq);
    strcpy(RNAInfo->SeqStruct, RNAStruct);
    GetOSeq(RNASeq, RNAStruct, RNAInfo->Seq);
    GetConstraint(RNAInfo->SeqConst, RNAStruct);
    GetSapcerConstraint(RNAInfo->SpacerConst, RNAStruct);
    CombineSpacerSEConstraint(RNAInfo->SpacerSEConst, RNAInfo->SeqConst, RNAInfo->SpacerConst);
    printf("%s\n%s\n%s\n", RNAInfo->SeqConst, RNAInfo->SpacerConst, RNAInfo->SpacerSEConst);
    
    //spacer index
    GetSpacerIndex(RNAInfo->SpacerIndex, RNAInfo->SeqStruct);
    GetSpacerSEIndex(RNAInfo->SpacerSEIndex, RNAInfo->SpacerSEConst);
    
    // probability context
    RNAInfo->ProbContextSeq = malloc(RNAInfo->SeqLen*sizeof(double*));
    for (i = 0; i < RNAInfo->SeqLen; i++) {
        RNAInfo->ProbContextSeq[i] = malloc(kmer*sizeof(double));
    }
    UnpairProb4InterInteractionComputation(RNAInfo->Seq, RNAInfo->ProbContextSeq, kmer);
    RNAInfo->ProbContextAntiSeq = NULL;
    
}

void Init_RNADouble(struct RNA_Info *RNAInfo, char *RNASeq, char *RNAAntiSeq, char *RNAStruct, int kmer){
    int i;
    
    RNAInfo->TotalLen = strlen(RNASeq);
    RNAInfo->SeqLen = Length_Seq(RNAStruct);
    RNAInfo->SpacerLen = Length_Sapcer_Seq(RNAStruct);
    RNAInfo->AntiSeqLen = strlen(RNAAntiSeq);
    RNAInfo->SpacerSELen = Length_Spacer_SE(RNAStruct);
    
    printf("AllLEN %d, RNALEN %d, AntiRNALen %d\n", RNAInfo->TotalLen, RNAInfo->SeqLen, RNAInfo->AntiSeqLen);
    
    // make space
    RNAInfo->SeqAll = (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->SeqStruct = (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->SeqConst= (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->SpacerConst= (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->SpacerSEConst= (char *)malloc(sizeof(char)*(RNAInfo->TotalLen+1));
    RNAInfo->Seq = (char *)malloc(sizeof(char)*(RNAInfo->SeqLen+1));
    RNAInfo->AntiSeq = (char *)malloc(sizeof(char)*(RNAInfo->AntiSeqLen+1));
    RNAInfo->SpacerIndex = (int *)malloc(sizeof(int)*RNAInfo->SpacerLen);
    RNAInfo->SpacerSEIndex = (int *)malloc(sizeof(int)*RNAInfo->SpacerSELen);
    
    // close the string
    RNAInfo->SeqAll[RNAInfo->TotalLen] = '\0';
    RNAInfo->SeqStruct[RNAInfo->TotalLen] = '\0';
    RNAInfo->SeqConst[RNAInfo->TotalLen] = '\0';
    RNAInfo->SpacerConst[RNAInfo->TotalLen] = '\0';
    RNAInfo->SpacerSEConst[RNAInfo->TotalLen] = '\0';
    RNAInfo->Seq[RNAInfo->SeqLen] = '\0';
    RNAInfo->AntiSeq[RNAInfo->AntiSeqLen] = '\0';
    
    
    // copy sequence
    strcpy(RNAInfo->SeqAll, RNASeq);
    strcpy(RNAInfo->SeqStruct, RNAStruct);
    GetOSeq(RNASeq, RNAStruct, RNAInfo->Seq);
    GetConstraint(RNAInfo->SeqConst, RNAStruct);
    GetSapcerConstraint(RNAInfo->SpacerConst, RNAStruct);
    strcpy(RNAInfo->AntiSeq, RNAAntiSeq);
    CombineSpacerSEConstraint(RNAInfo->SpacerSEConst, RNAInfo->SeqConst, RNAInfo->SpacerConst);
    
    
    //spacer index
    GetSpacerIndex(RNAInfo->SpacerIndex, RNAInfo->SeqStruct);
    GetSpacerSEIndex(RNAInfo->SpacerSEIndex, RNAInfo->SpacerSEConst);
    
    // probability context
    RNAInfo->ProbContextAntiSeq = malloc(RNAInfo->AntiSeqLen*sizeof(double*));
    for (i = 0; i < RNAInfo->AntiSeqLen; i++) {
        RNAInfo->ProbContextAntiSeq[i] = malloc(kmer*sizeof(double));
    }
    UnpairProb4InterInteractionComputation(RNAInfo->AntiSeq, RNAInfo->ProbContextAntiSeq, kmer);
    RNAInfo->ProbContextSeq = NULL;
    
}

void FreeRNAInfo(struct RNA_Info *RNAInfo, int type){
    int i;
    
    if(type == 1){
        free(RNAInfo->SeqAll);
        free(RNAInfo->SeqStruct);
        free(RNAInfo->SeqConst);
        free(RNAInfo->SpacerConst);
        free(RNAInfo->SpacerSEConst);
        free(RNAInfo->Seq);
        free(RNAInfo->SpacerIndex);
        free(RNAInfo->SpacerSEIndex);
        
        for(i = 0; i < RNAInfo->SeqLen; i++){
            free(RNAInfo->ProbContextSeq[i]);
        }
        free(RNAInfo->ProbContextSeq);
    }
    
    if(type == 2){
        free(RNAInfo->SeqAll);
        free(RNAInfo->SeqStruct);
        free(RNAInfo->SeqConst);
        free(RNAInfo->SpacerConst);
        free(RNAInfo->SpacerSEConst);
        free(RNAInfo->Seq);
        free(RNAInfo->AntiSeq);
        free(RNAInfo->SpacerIndex);
        free(RNAInfo->SpacerSEIndex);
        
        for(i = 0; i < RNAInfo->AntiSeqLen; i++){
            free(RNAInfo->ProbContextAntiSeq[i]);
        }
        free(RNAInfo->ProbContextAntiSeq);
    }
}

int Length_SE(char *structure){
    int len;
    int i;
    int count = 0;
    
    len = strlen(structure);
    for(i = 0; i < len; i++){
        if(structure[i] == 'x'){
            count++;
        }
    }
    
    return count;
}

void Index_SE(char *str, int *index){
    int len, i, count;
    
    count = 0;
    len = (int)strlen(str);
    for (i = 0; i < len; i++) {
        if (str[i] == 'x') {
            index[count++] = i;
        }
    }
    
}

void Init_SE(struct SE_Info *SEInfo, struct Aptamer_Info *AptamerInfo, struct RNA_Info *RNAInfo){
    SEInfo->SELen = Length_SE(AptamerInfo->SeqStruct);
    
    int *index1, *index2;
    index1 = malloc(sizeof(int)*SEInfo->SELen);
    index2 = malloc(sizeof(int)*SEInfo->SELen);
    
    //make space
    SEInfo->SEIndex = malloc(sizeof(struct SE_Index)*SEInfo->SELen);
    
    int i, j;
    Index_SE(AptamerInfo->SeqStruct, index1);
    Index_SE(RNAInfo->SeqStruct, index2);
    for(i = 0; i < SEInfo->SELen; i++){
        //printf("%d %d\n", index1[i], index2[SEInfo->SELen-1-i]);
        SEInfo->SEIndex[i].InApta = index1[i];
        SEInfo->SEIndex[i].InRNA  = index2[SEInfo->SELen-1-i];
    }
    
    free(index1);
    free(index2);
}


void FreeSEInfo(struct SE_Info *SEInfo){
    free(SEInfo->SEIndex);
}

void Init_SEDesign(struct SE_Design *SED, char * AptamerSeq, char *AptamerStruct, char *RNASeq, char *RNAAntiSeq, char *RNAStruct, int type, int kmer){
    SED->AttachType = type;
    SED->kmer = kmer;
    
    // initi Aptamer
    Init_Aptamer(&SED->AptamerInfo, AptamerSeq, AptamerStruct, kmer);
    
    // initi RNA
    if(type == 1){
        Init_RNASingle(&SED->RNAInfo, RNASeq, RNAStruct, kmer);
    }
    else if(type == 2){
        Init_RNADouble(&SED->RNAInfo, RNASeq, RNAAntiSeq, RNAStruct, kmer);
    }
     
    
    // initi SE
    Init_SE(&SED->SEInfo, &SED->AptamerInfo, &SED->RNAInfo);
}


void Free_SEDesign(struct SE_Design *SED){

    FreeAptamerInfo(&SED->AptamerInfo);
    FreeRNAInfo(&SED->RNAInfo, SED->AttachType);
    FreeSEInfo(&SED->SEInfo);
}

















