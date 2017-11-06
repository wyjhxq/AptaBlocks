AptaBlocks: AptaBlocks.o Initialization.o Parameters.o PottsModel.o
	gcc -o AptaBlocks AptaBlocks.o Initialization.o Parameters.o PottsModel.o -fopenmp -L /usr/local/ViennaRNA/2.1.9/lib -lRNA -lm
AptaBlocks.o: AptaBlocks.c
	gcc -c AptaBlocks.c -fopenmp -I /usr/local/ViennaRNA/2.1.9/include/ViennaRNA
Initialization.o: Initialization.c
	gcc -c Initialization.c -fopenmp -I /usr/local/ViennaRNA/2.1.9/include/ViennaRNA
Parameters.o: Parameters.c
	gcc -c Parameters.c -I /usr/local/ViennaRNA/2.1.9/include/ViennaRNA
PottsModel.o: PottsModel.c
	gcc -c PottsModel.c -I /usr/local/ViennaRNA/2.1.9/include/ViennaRNA