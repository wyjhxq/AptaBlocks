## Installation
1. Install ViennaRNA package and find its path
2. Change `viennarna_path` in Makefile to reflect the location of ViennaRNA on your system
3. Run `make`  

## Input parameters
```bash
$ ./AptaBlocks -h 

-N [integer(mandatory)]  the number of aptamer-'cargo' paris (For one aptamer and one cargo, N = 1)
-l [positve number]  \lambda in the oriningal paper. The larger the lower the binding energy of stikcy bridges is. default l = 1
-g [postive number]  The larger the larger the GC content. default g = 1
-i [mandatory] the name of the input file
-o [mandatory] the name of the output file
```

## Input file format
1. All sequences are from 5'- -3'
2. 'Aptamer' is a keyword to indicate the following is the sequence of an aptamer+stikcy bridge
3. 'ssRNA' is a keyword to indicate the following is a ssRNA + stikcy bridge
4. 'dsRNA' is a keyword to indicate the following is a dsRNA + stikcy bridge
5. '*' indicates known bases and 'x' denotes unknown sticky bridge sequences
6. sticky bridges need to be complementary in the input file

### Input file examples
Designing 5nt sticky bridge for one aptamer (UACCUGGUAC) and one ssRNA cargo (CUCGGCCUUA). (./Example/Input_ssRNA_cargo.txt)
```
>Aptamer
UACCUGGUACGGGGG
**********xxxxx
>ssRNA
CUCGGCCUUACCCCC
**********xxxxx
```

Designing 5nt sticky bridge for one aptamer (UACCUGGUAC) and one dsRNA cargo (CUCUGCUUCGGUGUCGAAAUGAGAACC || UUCUCAUUUCGACACCGAAGCAGAG). (./Example/Input_dsRNA_cargo.txt)
```
>Aptamer
UACCUGGUACGGGGG
**********xxxxx
>dsRNA
CUCUGCUUCGGUGUCGAAAUGAGAACCCCCCC
***************************xxxxx
UUCUCAUUUCGACACCGAAGCAGAG
```

Designing a universal sticky bridge for the above ssRNA and dsRNA. (./Example/Input_universal_ssRNA_dsRNA.txt)
```
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
```

## Example of executing AptaBlocks
1. For one aptamer and one ssRNA cargo:
```bash
$ ./AptaBlocks -N 1 -l 2.5 -i ./Example/Input_ssRNA_cargo.txt -o ./Example/Output_ssRNA_cargo.txt 
```

2. For one aptamer and one dsRNA cargo:
```bash
$ ./AptaBlocks -N 1 -l 1.0 -i ./Example/Input_dsRNA_cargo.txt -o ./Example/Output_dsRNA_cargo.txt 
```

3. For one aptamer and one ssRNA and one dsRNA (universal sticky bridge)
```bash
$ ./AptaBlocks -N 2 -l 10 -i ./Example/Input_universal_ssRNA_dsRNA.txt.txt -o ./Example/Output_universal_ssRNA_dsRNA.txt.txt 
```
## Output file
Take one aptamer and one ssRNA cargo for example, the output looks like:
```bash
1   Aptamer: UACCUGGUACGCUGU 
2   RNA: CUCGGCCUUAACAGC 
3   Probability with structure constraints: 0.90521 * 0.84481 = 0.76473 
4   Probability with structure constraints: 0.765  log_2: -0.386975
5   Energy of SE (converted in [0,1]): 0.619978 -5.710000 log_2: -0.689710
6   Probability of inter-interactions: 1.000000 log_2: 0.000000
7   number of consective identical: 0.000000
8   GC content: gcc_apta=0.600000, gcc_rna=0.600000, log_2(gcc_apta*gcc_rna)=-1.473931
9   Overall objectives: -1.076686
```
### Explianation of the output file
1. Line 1 gives the aptamer sequence + the first part of the sticky brdige sequences (we use B1 to present it)
2. Line 2 gives ssRNA sequence + the secand part sticky brdige sequence (we use B2 to present it). B1 and B2 are complementary.
3. Line 3 gives the probability that B1 is unparied in aptamer-stick (the first number) and the probability of B2 is unpaired in cargo-stick (the second number)
4. Line 5 gives the hybridization energy of the sticky bridge (the second number (kcal/mol))
5. Line 6 gives the proability of not forming into undesired dimers (the first number)
