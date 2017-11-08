# AptaBlocks

#Installation
1. Install ViennaRNA package and find its path
2. Change "viennarna_path" in Makefile 
3. Run make  

#Input parameters
./AptaBlocks -h # show the explanations of all input parameters.

-N [integer]  the number of aptamer-'cargo' paris (For one aptamer and one cargo, N = 2)

-l [positve number]  \lambda in the oriningal paper. The larger the lower the binding energy of stikcy bridges is. 

-g [postive number]  The larger the larger the GC content. 

-i the name of the input file

-o the name of the output file

# Example 
./AptaBlocks -N 1 -i ./Example/Input_example.txt -o ./Example/Output.txt 
