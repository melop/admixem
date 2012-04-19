#include "maths.h"

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
	int nRet = (int)floor(UniformGen.Next() * (double)nSizeOfArray - 0.01);
	return nRet>=0? nRet:0;
}

double round(double nNum) {
	double nIntegral = (double)((int)nNum);
	if (nNum - nIntegral < 0.5) {
		return nIntegral; 
	}
	else {
		return nIntegral + 1;
	}
}
