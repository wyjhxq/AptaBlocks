#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <string.h>
#include <unistd.h>
#include  "utils.h"
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
#include "omp.h"
#include "Initialization.h"
#include "Parameters.h"
#include "PottsModel.h"


void print_AptaBlocks_help() {
    printf("Usage: AptaBlocks [OPTIONS] PREFIX\n");
    printf("Design the sequences of sticky ends (SEs) for gluing aptamer and 'cargo' RNA together. \n");
    printf("Example: AptaBlocks -N 1 -l 1 -g 0.5 -i input.txt -o output.txt \n");
    printf("\n");
    printf("\n");
    printf("Initialization:\n");
    printf(" -N [integer]  the number of aptamer-'cargo' paris (needs to correapond to number in input.txt. )\n");
    printf(" -l [positve number]  lambda that balances between structure require and stickiness of the sticky end. The larger, the stickier. \n");
    printf(" -g [postive number]  encourage larger GC content in the designed sticky ends. The larger, the larger the GC content. \n");
    printf(" -i the name of the input file\n");
    printf(" -o the name of the output file\n");
    printf("\n");
    printf("\n");
    printf("Input file example: # all sequences are from 5'- -3'\n");
    printf(">Aptamer #'Aptamer' is a keyword to indicate the following is the sequence of an aptamer+SE");
    printf("AUGCCUUGGCCCAAUACCGGGGGGGGG # the sequence of an aptamer+SE\n");
    printf("***************...xxxxxxxxx # '*': fixed nuclitude, '.': a spacer between aptamer and SE, could be any nuclitude, 'x': a nuclitude in a sticky end\n");
    printf(">RNA # 'RNA' is a keyword to indicate the following is only one strand\n");
    printf("ACCCCGGGGUUAACCGUUCCCCCCCCC\n");
    printf("******************xxxxxxxxx\n");
    printf(">siRNA # 'siRNA' is a keyword to indicate the following contain the sense and anti-sense part of the siRNA\n");
    printf("CUCUGCUUCGGUGUCGAAAUGAGAACCCCCCCCCCC\n");
    printf("***************************xxxxxxxxx\n");
    printf("UUCUCAUUUCGACACCGAAGCAGAG\n");
    printf("\n");
    printf("\n");
    printf("Output file example: # all sequences are from 5'- -3'\n");
    printf("Aptamer: GGGAGGACGAUGCGGGCCUUCGUUUGUUUCGUCCACAGACGACUCGCCCGAAGAAAUAAAUGUAAUGUGUAGG # sequence of aptamer and a designed SE.\n");
    printf("RNA: GAGCUGGACGGCGACGUAAAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCAUAAACCUACACAUUACAUUUA # desinged SE for RNA\n");
    printf(".... # some output that charactrized the desinged SE.\n");
    exit(0);
}


void ReadInput(struct SE_Design *SEDesign, char FileName[], int kmer, int NumPairs){
    char Apta_Seq[5000], Apta_Structure[5000];
    char Temp1_Seq[5000], Temp2_Seq[5000], Temp1_Structure[5000];
    char temp[100];
    int read;
    char *line;
    size_t len=0;
    int i, count;
    
    FILE *fh = fopen(FileName, "r");
    if(fh == NULL){
        printf("Can not open %s\n", FileName);
        exit(0);
    }
    
    count = 0;
    while ((read = getline(&line, &len, fh)) != -1) {
        if(line[0] == '>'){
            memcpy(temp, &line[1], strlen(line)-2);
            temp[strlen(line)-2] = '\0';
            //printf("temp: %s %d\n", temp, strcmp(temp, "Aptamer"));
            if(strcmp(temp, "Aptamer") == 0){
                read = getline(&line, &len, fh);
                if (read != -1) {
                    strncpy(Apta_Seq, line, strlen(line)-1);
                    Apta_Seq[strlen(line)-1] = '\0';
                    printf("Aptamer: %s\n", Apta_Seq);
                }
                read = getline(&line, &len, fh);
                if (read != -1) {
                    strncpy(Apta_Structure, line, strlen(line)-1);
                    Apta_Structure[strlen(line)-1] = '\0';
                    printf("Aptamer: %s\n", Apta_Structure);
                }
            }else if(strcmp(temp, "RNA") == 0){
                read = getline(&line, &len, fh);
                if (read != -1) {
                    strncpy(Temp1_Seq, line, strlen(line)-1);
                    Temp1_Seq[strlen(line)-1] = '\0';
                    printf("RNA: %s\n", Temp1_Seq);
                }
                read = getline(&line, &len, fh);
                if (read != -1 && line[strlen(line)-1] == '\n') {
                    strncpy(Temp1_Structure, line, strlen(line)-1);
                    Temp1_Structure[strlen(line)-1] = '\0';
                    printf("RNA: %s\n", Temp1_Structure);
                }else{
                    strcpy(Temp1_Structure, line);
                }
                Init_SEDesign(&SEDesign[count], Apta_Seq, Apta_Structure, Temp1_Seq, NULL, Temp1_Structure, 1, kmer);
                count++;
            }else if(strcmp(temp, "siRNA") == 0){
                read = getline(&line, &len, fh);
                if (read != -1) {
                    strncpy(Temp1_Seq, line, strlen(line)-1);
                    Temp1_Seq[strlen(line)-1] = '\0';
                }
                read = getline(&line, &len, fh);
                if (read != -1) {
                    strncpy(Temp1_Structure, line, strlen(line)-1);
                    Temp1_Structure[strlen(line)-1] = '\0';
                }
                read = getline(&line, &len, fh);
                if (read != -1 && line[strlen(line)-1] == '\n') {
                    strncpy(Temp2_Seq, line, strlen(line)-1);
                    Temp2_Seq[strlen(line)-1] = '\0';
                }else{
                    strcpy(Temp2_Seq, line);
                }
                Init_SEDesign(&SEDesign[count], Apta_Seq, Apta_Structure, Temp1_Seq, Temp2_Seq, Temp1_Structure, 2, kmer);
                count++;
            }else{
                printf("Wrong keyword in input file %s\n", FileName);
                exit(0);
            }
        }
    }

    if(count != NumPairs){
        printf("Number of pairs is not consistent with the -N input.\n");
        exit(0);
    }
}


int main(int argc, char *argv[]) {
    
    // every input sequences should be ordered as 5'- xxxxx.... -3'
    char
    *seq1        ="GGGAGGACGAUGCGGGCCUUCGUUUGUUUCGUCCACAGACGACUCGCCCGAGGGGGGGGGGGGGGGGGGGGGG",//UUUUUUUUUUUUUUUUUUUU
    *struct1_con ="***************************************************.....xxxxxxxxxxxxxxxxx",
    *seq2        ="AUGAACACCAGUGAGUAGAGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUUUUUUUUUUUU",//CAAAAUGUAAAAAACU
    *struct2_con ="************************************************************************************************xxxxxxxxxxxxxxxxx",
    *seq3        ="GCUGCCGCCCAGUGGGACUUGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCCCCCCCCCCCCCCCCCCCCCC",
    *struct3_con ="************************************************************************************************xxxxxxxxxxxxxxxxxxxxx",
    *seq4        ="CUCUGCUUCGGUGUCGAAAUGAGAACCCCUAAAUAUAUAUUAAAUCU",
    *struct4_con ="***************************xxxxxxxxxxxxxxxxxxxx",
    *seq4_anti   ="UUCUCAUUUCGACACCGAAGCAGAG",
    *seq5        ="UACACGAGCUGCAAUGUCGGCUUUGCUCCUAAAUAUAUAUUAAAUCU",
    *struct5_con ="***************************xxxxxxxxxxxxxxxxxxxx",
    *seq5_anti   ="CAAAGCCGACAUUGCAGCUCGUGUAUU",
    *seq6        ="UCAAGAGACUCCUCAGUGAGAAGAACCUAAAUAUAUAUUAAAUCU",
    *struct6_con ="*************************xxxxxxxxxxxxxxxxxxxx",
    *seq6_anti   ="UUCUUCUCACUGAGGAGUCUCUUGAUU",
    *seq8        ="GAGCUGGACGGCGACGUAAAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUUUUUUUUUUUUUUUUU",
    *struct8_con ="************************************************************************************************.....xxxxxxxxxxxxxxxxx",
    *test1       ="AUGCCUUGGCCCAAUGGGGGGGGG",
    *test1_con   ="***************xxxxxxxxx",
    *test2       ="ACCCCGGGGUUAACCGUUCCCCCCCCC",
    *test2_con   ="******************xxxxxxxxx";
    
    
    int NP=1;
    double lambda=1;
    double eta=1;
    char Input[500], Output[500];
    int c;

    while ((c = getopt (argc, argv, "i:o:")) != -1)
        switch (c)
    {
        case 'i':
            strcpy(Input, optarg);
            break;
        case 'o':
            strcpy(Output, optarg);
            break;
        case '?':
            print_AptaBlocks_help();
            return 1;
        default:
            abort ();
    }

    //printf("Number of pairs: %d, lambda: %f, eta: %f, input: %s, output: %s\n", NP, lambda, eta, Input, Output);
    
    // command line agurments
    double lambda1, lambda2;
    double kmer;
    int NumPairs = NP;
    lambda1 = eta;
    lambda2 = lambda;
    int i;
    kmer = 5;
    
    // initi sequences
    struct SE_Design *SEDesign;
    SEDesign = (struct SE_Design *)malloc(NumPairs*sizeof(struct SE_Design));
    ReadInput(SEDesign, Input, kmer, NumPairs);
    //Init_SEDesign(&SEDesign[0], test1, test1_con, test2, NULL, test2_con, 1, kmer);
    //Init_SEDesign(&SEDesign[0], seq1, struct1_con, seq2, NULL, struct2_con, 1, kmer);
    //Init_SEDesign(&SEDesign[0], seq1, struct1_con, seq8, NULL, struct8_con, 1, kmer);
    //Init_SEDesign(&SEDesign[1], seq1, struct1_con, seq3, NULL, struct3_con, 1, kmer);
    //Init_SEDesign(&SEDesign[2], seq1, struct1_con, seq4, seq4_anti, struct4_con, 2, kmer);
    //Init_SEDesign(&SEDesign[3], seq1, struct1_con, seq5, seq5_anti, struct5_con, 2, kmer);
    //Init_SEDesign(&SEDesign[4], seq1, struct1_con, seq6, seq6_anti, struct6_con, 2, kmer);
    
    // initi parameters
    struct Parameters Param;
    Init_Param(&Param, 1.0, 37.0, -1.0, lambda1, lambda2, 100, 0.001, 0.99);
    printf("pf_parasms->kT is %f %f and pf_params->pf_scale is %f Tempature %f\n", (Param.pf_parameters)->kT, Param.KT, Param.pf_parameters->pf_scale, Param.pf_parameters->temperature);
    
    printf("energy RNA-RNA %f\n", StackingEnerogyRNARNA[0][3]);
    
    
    // show scores for initi SE
    FILE *fp = fopen(Output, "w");
    for(i = 0; i < NumPairs; i++){
        PrintSingleSEDesignFile(SEDesign[i], Param, fp);
    }
    fclose(fp);

    
    // main computations
    //FILE *fh = fopen(Output, "w");
    //Potts_ModelM(SEDesign, Param, NumPairs, fh);
    
    // free spaces
    for(i = 0; i < NumPairs; i++){
        Free_SEDesign(&SEDesign[i]);
    }
    free(SEDesign);
    
    return 1;
}


