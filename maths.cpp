#include "maths.h"
#include <sstream>
#ifdef _OPENMP
 #include <omp.h>
#endif

Normal NormalGen; 
Uniform UniformGen;
double nRandSeed;



double NormalExt(double nMean, double nStdDev, double nLowBound, double nHighBound) {
	
	double nRet = nStdDev * NormalGen.Next() + nMean;
	return  ( nRet  > nHighBound ) ? nHighBound : ( nRet < nLowBound? nLowBound: nRet );
}


double fnSum( double * pNums, int nSize ) {

	double nRet = 0;
	for (int i=0; i< nSize; i++) {
		nRet += pNums[i];
	 } 
	 
	return nRet;
}

int fnCompare (const void * a, const void * b)
{

  double nComp = *(double*)a - *(double*)b;
  return ( nComp < 0 ) ? -1: ( (nComp == 0)? 0:1 );
}

int fnGetRandIndex(int nSizeOfArray) {
	int nRet;
	
	do {
		
			nRet = (int)floor(UniformGen.Next() * (double)nSizeOfArray - 0.01);
		
		nRet = nRet>=0? nRet:0;
	}
	while(nRet >= nSizeOfArray);
	
	/*
	const unsigned long n = nSizeOfArray;
    const unsigned long divisor = (RAND_MAX + 1) / n;

    unsigned long nRet;
    do { nRet = rand() / divisor; } while (nRet >= n);
	*/
	return nRet;
};

double round(double nNum) {
	double nIntegral = (double)((int)nNum);
	if (nNum - nIntegral < 0.5) {
		return nIntegral; 
	}
	else {
		return nIntegral + 1;
	}
};

string fnIntToString(int nNum ) {
	std::stringstream sstm;
	sstm << nNum;
    return (sstm.str());

};
