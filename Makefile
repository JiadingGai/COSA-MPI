CC=gcc
#pgcc
EnableMacroInGDB=-gdwarf-2 -g3
all:cosa2 parallel_cosa2 compare
parallel_cosa2:parallel_cosa2.c data_io.o matrix_util.o
	mpicc -Wall -g -o $@ $^ -lm $(PROFILING) 
	#gcc $(EnableMacroInGDB) -g -o $@ $^ -lm -lmpi
cosa2:cosa2.o data_io.o matrix_util.o
	$(CC) $(EnableMacroInGDB) -Wall -g -o $@ $^ -lm
cosa2.o:cosa2.c
	$(CC) $(EnableMacroInGDB) -Wall -g -c $^ -lm
matrix_util.o:matrix_util.c
	$(CC) $(EnableMacroInGDB) -Wall -g -c matrix_util.c -lm
data_io.o:data_io.c
	$(CC) $(EnableMacroInGDB) -Wall -g -c data_io.c
compare:compare.c
	$(CC) $(EnableMacroInGDB) -Wall -g -o $@ $^ -lm
clean:
	rm -f cosa2 parallel_cosa2 compare *.o 
