To use MpiP profiler, set LD_LIBRARY_PATH properly as follows:

$ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/afs/crc.nd.edu/user/j/jgai/Private/project/mpiP-3.2.1/libunwind-installed/lib
$ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/afs/crc.nd.edu/user/j/jgai/Private/project/mpiP-3.2.1/binutil-2.20.1-installed/lib

Also, don't forget to set MPIP environment variable:

$ export MPIP="-t 10.0 -k 2"

USAGE:

cosa2:

cosa2 [data file in binary format] [# of objects] [# of attributes] [# of outer iterations] [# of inner iterations]

parallel_cosa2:

mpirun -np [# of processors] parallel_cosa2 [data file in binary format] [# of objects] [# of attributes] [px] [py] [pz] [# of outer iterations] [# of inner iterations]

c22_ceu.txt
     an N-by-n data Matrix saved in a row-wise fashion.
     Each entry of the matrix is delimited by the tab character.

     n=13000: # of attributes
     N=1000 : # of objects

e0.dat: 
     an N-by-n data Matrix saved in a row-wise fashion.
     Each entry of the matrix is delimited by the tab character.

     n=13000: # of attributes
     N=1000 : # of objects
     
e4.dat: 
     an N-by-n data Matrix saved in a row-wise fashion.
     Each entry of the matrix is delimited by the tab character.

     n=10: # of attributes
     N=6: # of objects

e6.dat: 
     an N-by-n data Matrix saved in a row-wise fashion.
     Each entry of the matrix is delimited by the tab character.

     n=10: # of attributes
     N=7: # of objects

e7.dat: 
     an N-by-n data Matrix saved in a row-wise fashion.
     Each entry of the matrix is delimited by the tab character.

     n=10: # of attributes
     N=8: # of objects
