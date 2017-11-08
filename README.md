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

For one aptamer and one ssRNA cargo: ./Example/Input_ssRNA_cargo.txt
>Aptamer
UACCUGGUACGGGGG
**********xxxxx
>ssRNA
CUCGGCCUUACCCCC
**********xxxxx

For one aptamer and one dsRNA cargo: ./Example/Input_dsRNA_cargo.txt

For one aptamer and one ssRNA cargo and one dsRNA cargo: ./Example/Input_universal_ssRNA_dsRNA.txt

# Example 
./AptaBlocks -N 1 -i ./Example/Input_example.txt -o ./Example/Output.txt 
