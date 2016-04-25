/*****************************************************************************************
Parallel Implementation of Gaussian Elimination algorithm Using P-threads
******************************************************************************************/





/*******COMPILING THE GAUSSIAN ELIMINATION ALGORITHM****************************/

gcc -O2 -pthread -o gauss gaussian_parallel.c

/**********************************************************************************/





/***************** EXECUTING THE GAUSSIAN ELIMINATION  ************************************/

/usr/bin/time ./gauss -n [size] -I [type] -m [maximum number] -p [0/1]

   I - rand/fast
   m - 23/12/12..............
   p - 0/1

Example : /usr/bin/time ./gauss -n 4096 -I rand -m 30 -p 0

/***************************************************************************************************/
