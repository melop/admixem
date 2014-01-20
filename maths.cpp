#include "maths.h"
#include <sstream>
#ifdef _OPENMP
 #include <omp.h>
#endif

#ifdef _RANGEN_ARRAYS_
#else

	#ifdef _OPENMP
		int nTotalCPUCore =  omp_get_num_threads();
	#else
		int nTotalCPUCore = 1;
	#endif

Normal * arrNormalGen = new Normal[nTotalCPUCore]; 
Uniform * arrUniformGen = new Uniform[nTotalCPUCore];
#define _RANGEN_ARRAYS_
#endif

double nRandSeed;

double NormalExt(double nMean, double nStdDev, double nLowBound, double nHighBound) {
	#ifdef _OPENMP
		int nCurrProccess = omp_get_thread_num();
	#else
		int nCurrProccess = 0;
	#endif

	double nRet = nStdDev * arrNormalGen[nCurrProccess].Next() + nMean;
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
	#ifdef _OPENMP
		int nCurrProccess =  omp_get_thread_num();
	#else
		int nCurrProccess = 0;
	#endif

	do {
		
			nRet = (int)floor(arrUniformGen[nCurrProccess].Next() * (double)nSizeOfArray - 0.01);
		
		nRet = nRet>=0? nRet:0;
	}
	while(nRet >= nSizeOfArray);

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
