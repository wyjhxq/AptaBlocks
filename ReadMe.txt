# AptaBlocks

#Installation
1. Install ViennaRNA package and find its path
2. Change "viennarna_path" in Makefile 
3. Run make  

#Input parameters
./AptaBlocks -h # show the explanations of all input parameters.

-N [integer(mandatory)]  the number of aptamer-'cargo' paris (For one aptamer and one cargo, N = 2)
-l [positve number]  \lambda in the oriningal paper. The larger the lower the binding energy of stikcy bridges is. default l = 1
-g [postive number]  The larger the larger the GC content. default g = 1
-i [mandatory] the name of the input file
-o [mandatory] the name of the output file

#Input file format
1. All sequences are from 5'- -3'
2. 'Aptamer' is a keyword to indicate the following is the sequence of an aptamer+stikcy bridge
3. 'ssRNA' is a keyword to indicate the following is a ssRNA + stikcy bridge
4. 'dsRNA' is a keyword to indicate the following is a dsRNA + stikcy bridge
5. '*' indicates known bases and 'x' denotes unknown sticky bridge sequences
6. sticky bridges need to be complementary in the input file

Designing 5nt stikcy bridge for one aptamer (UACCUGGUAC) and one ssRNA cargo (CUCGGCCUUA). (./Example/Input_ssRNA_cargo.txt)
>Aptamer
UACCUGGUACGGGGG
**********xxxxx
>ssRNA
CUCGGCCUUACCCCC
**********xxxxx

Designing 5nt stikcy bridge for one aptamer (UACCUGGUAC) and one dsRNA cargo (CUCUGCUUCGGUGUCGAAAUGAGAACC || UUCUCAUUUCGACACCGAAGCAGAG). (./Example/Input_dsRNA_cargo.txt)
>Aptamer
UACCUGGUACGGGGG
**********xxxxx
>dsRNA
CUCUGCUUCGGUGUCGAAAUGAGAACCCCCCC
***************************xxxxx
UUCUCAUUUCGACACCGAAGCAGAG


Designing a universal sticky bridge for the above ssRNA and dsRNA. (./Example/Input_universal_ssRNA_dsRNA.txt)
>Aptamer
UACCUGGUACGGGGG
**********xxxxx
>ssRNA
CUCGGCCUUACCCCC
**********xxxxx
>dsRNA
CUCUGCUUCGGUGUCGAAAUGAGAACCCCCCC
***************************xxxxx
UUCUCAUUUCGACACCGAAGCAGAG


# Example of implement ApatBlocks
1. For one aptamer and one ssRNA cargo:
./AptaBlocks -N 1 -l 2.5 -i ./Example/Input_ssRNA_cargo.txt -o ./Example/Output_ssRNA_cargo.txt 

2. For one aptamer and one dsRNA cargo:
./AptaBlocks -N 1 -l 1.0 -i ./Example/Input_dsRNA_cargo.txt -o ./Example/Output_dsRNA_cargo.txt 

3. For one aptamer and one ssRNA and one dsRNA (universal sticky bridge)
./AptaBlocks -N 2 -l 10 -i ./Example/Input_universal_ssRNA_dsRNA.txt.txt -o ./Example/Output_universal_ssRNA_dsRNA.txt.txt 

# Output file
Take one aptamer and one ssRNA cargo for exampel, the output looks like:
Aptamer: UACCUGGUACGCUGU # aptamer sequence + sticky brdige sequence B1
RNA: CUCGGCCUUAACAGC # ssRNA sequence + sticky brdige sequence B2
Probability with structure constraints: 0.90521 (<-probability of B1 unpaired) * 0.84481 (<-probability of B2 unpaired) = 0.76473 
Probability with structure constraints: 0.765  log_2: -0.386975
Energy of SE (converted in [0,1]): 0.619978 -5.710000 (<-hybridzation energy between B1 and B2) log_2: -0.689710
Probability of inter-interactions: 1.000000 (<-probability of not forming into undesired dimers) log_2: 0.000000
number of consective identical: 0.000000
GC content: gcc_apta=0.600000, gcc_rna=0.600000, log_2(gcc_apta*gcc_rna)=-1.473931
Overall objectives: -1.076686
