#ifndef Initialization_h
#define Initialization_h

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <string.h>



struct Aptamer_Info
{
    int TotalLen;           // length of the overall "GGGAGGACGAUGCGGGCCUUCGUUUGUUUCGUCCACAGACGACUCGCCCGAoooooUUUUUUUUUUUUUUUUUUUU";
    char *SeqAll;           // "GGGAGGACGAUGCGGGCCUUCGUUUGUUUCGUCCACAGACGACUCGCCCGAoooooUUUUUUUUUUUUUUUUUUUU" 'o' indicates we cannot change the spacer, 'A', 'U', 'G' or 'C' means we also need to design the spacer
    char *SeqStruct;        // "***************************************************.....xxxxxxxxxxxxxxxxxxxx": * -> aptamer; x -> SE; . -> spacer
    char *SeqConst;         // |.x(){}><
    int SeqLen;
    char *AptaSeq;          // "GGGAGGACGAUGCGGGCCUUCGUUUGUUUCGUCCACAGACGACUCGCCCGA"
    int SpacerLen;
    int *SpacerIndex;
    char *SpacerConst;      // spacer constraint : ".....................................xxxxx...................."
    char *SpacerSEConst;    // spacer constraint + SE constraint : ".....................................xxxxxxxxxxxxxxxxxxxxxxxxx"
    int SpacerSELen;
    int *SpacerSEIndex;
    double **ProbContext;   // probablity context of each kmer
};


struct RNA_Info
{
    int TotalLen;           // length of the overall "GGGAGGACGAUGCGGGCCUUCGUUUGUUUCGUCCACAGACGACUCGCCCGAoooooUUUUUUUUUUUUUUUUUUUU"
    char *SeqAll;           // 5'-SE-ooo-seq-3' 'o' indicates we cannot change the spacer, 'A', 'U', 'G' or 'C' means we also need to design the spacer
    char *SeqStruct;        // "***********.....xxxxxxxxxxxxxxxxxxxx": * -> aptamer; x -> SE; . -> spacer
    char *SeqConst;         //
    
    int SeqLen;
    char *Seq;              // sequence without SE and ooooo
    int SpacerLen;
    int *SpacerIndex;
    char *SpacerConst;      // spacer constraint : "..........xxxxx.........."
    char *SpacerSEConst;    // spacer constraint + SE constraint : "..........xxxxxxxxxxxxxxxxx"
    int SpacerSELen;
    int *SpacerSEIndex;
    int AntiSeqLen;
    char *AntiSeq;          // type 1: null
    double **ProbContextSeq;  // type 1: probability context for Seq; type 2: probability context for AntiSeq
    double **ProbContextAntiSeq;  // type 1: probability context for Seq; type 2: probability context for AntiSeq
};

struct SE_Index
{
    int InApta;
    int InRNA;
};

struct SE_Info
{
    int SELen;
    struct SE_Index *SEIndex;
};

struct Spacer_Info
{
    int type; //
};

struct SE_Design
{
    // note: only these two types (one aptamer vs. one sgRNA(1) or siRNA(2)).
    int AttachType;  // 1: 5'-seq-oooo-SE-3'             |  2: 5'-seq-oooo-SE-3'
                     //             3'-SE-ooo-seq-5'     |              3'-SE-ooo-seq-5'
                     //                                  |                   5'-anti_seq-3'
    int kmer;
    struct Aptamer_Info AptamerInfo;
    struct RNA_Info RNAInfo;
    struct SE_Info SEInfo;
};

void Init_Aptamer(struct Aptamer_Info *AptamerInfo, char *AptamerSeq, char *AptamerStruct, int kmer);
void Init_SEDesign(struct SE_Design *SED, char * AptamerSeq, char *AptamerStruct, char *RNASeq, char *RNAAntiSeq, char *RNAStruct, int type, int kmer);
void Free_SEDesign(struct SE_Design *SED);



#endif


























































