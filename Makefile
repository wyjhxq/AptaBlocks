#adopt this path to point to your local ViennaRNA installation
viennarna_path = /usr/local/ViennaRNA/2.1.9/


AptaBlocks: AptaBlocks.o Initialization.o Parameters.o PottsModel.o
	gcc -o AptaBlocks AptaBlocks.o Initialization.o Parameters.o PottsModel.o -fopenmp -L ${viennarna_path}lib -lRNA -lm
AptaBlocks.o: AptaBlocks.c
	gcc -c AptaBlocks.c -fopenmp -I ${viennarna_path}include/ViennaRNA
Initialization.o: Initialization.c
	gcc -c Initialization.c -fopenmp -I ${viennarna_path}include/ViennaRNA
Parameters.o: Parameters.c
	gcc -c Parameters.c -I ${viennarna_path}include/ViennaRNA
PottsModel.o: PottsModel.c
	gcc -c PottsModel.c -I ${viennarna_path}include/ViennaRNA