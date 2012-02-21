# Makefile
#
CC= g++
#
# header files
#
#INCLUDE= -I/usr/include
#
# build options
#
CFLAGS= -O3 -g #$(INCLUDE)
#
# standard libraries
#
LIB= -L/usr/lib -lm
#
treeCS: treeCS.o tree_analyzer.o genebank.o
	$(CC) -o treeCS treeCS.o tree_analyzer.o genebank.o $(LIB)
#
treeCS.o: treeCS.cpp tree_analyzer.h
	$(CC) $(CFLAGS) -c treeCS.cpp
#
tree_analyzer.o: tree_analyzer.cpp genebank.h genotype.h
	$(CC) $(CFLAGS) -c tree_analyzer.cpp
#
genebank.o: genebank.cpp genotype.h
	$(CC) $(CFLAGS) -c genebank.cpp
#
clean:
	rm -f *.o treeCS












