/*
 *  makemarkerfile.cpp
 *  admixsimul
 *
 *  Created by Rosenthal Lab on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <map>
#include "def.h"
#include "makemarkerfile.h"
#include "newran02/newran.h"
#include "maths.h"


using namespace std;

extern Normal NormalGen; 
extern Uniform UniformGen;
extern double nRandSeed;


void UIMakeNewMarkerFile() {

int nCurrChromosome = 1;
int nHaploidChromosomeNum, nTotalMarkers;
string sOutFile = "";
double nChromosomeSize, nCentromerePos, nPop1AvgAlleleFreq, nPop1AvgAlleleFreqStdDev,  nPop2AvgAlleleFreq, nPop2AvgAlleleFreqStdDev;

	//Random::Set(0.42);

		while(true) { //Get random seed from user
		
			std::cout << "\nEnter a random seed (0.0~1.0). Write down this number to repeat simulation: ";
			std::cin >> nRandSeed;
			
			if (nRandSeed > 0 && nRandSeed < 1.0 ) {
				Random::Set(nRandSeed);
				break; // 
			}
			
			
			
		}

			
			/*
		  ofstream myfile;
			myfile.open("testnormal.txt");
			
			for (int i=0;i<=999;i++)
			{
				myfile << NormalExt(0.8, 0.1 , 0.5, 1.0) << "\n"; // draw from normal distribution, mean=0.8, stdev=0.1, >=0.5, <= 1.0
			
			}
			myfile.close();
	*/

	//std::cout << NormalExt(0.8, 0.1) << "\n";

		while(true) { //Get haploid chromosome number from user
		
			std::cout << "\nEnter the number of chromosomes in the haploid genome (n): ";
			std::cin >> nHaploidChromosomeNum;
			
			if (nHaploidChromosomeNum > 0 ) {
				break; // 
			}
			
			
			
		}
		
		pChromosomeSizes = new double[nHaploidChromosomeNum];
		double * pCentromerePos = new double[nHaploidChromosomeNum];
		
		while(true) { // Get sizes of each chromosome from user
		
			std::cout << "\nEnter the size of chromosome " << nCurrChromosome << "(enter a negative number to finish, enter 0 set it to the same number as the previous one): ";
			std::cin >> nChromosomeSize;

			if (nChromosomeSize == 0 && nCurrChromosome == 1) {
				std::cout << "Can't use this option\n";
				continue;
			}

			if (nChromosomeSize == 0) { //set it to the same as the previous one
				pChromosomeSizes[nCurrChromosome-1]  = pChromosomeSizes[nCurrChromosome-2];
				nChromosomeSize = pChromosomeSizes[nCurrChromosome-1];
			}
			
			if (nChromosomeSize < 0 ) {
				break; // 
			}
			
			if (nChromosomeSize > 0) {
				pChromosomeSizes[nCurrChromosome-1] = nChromosomeSize;
			}
			
			/*
			Now do not allow users to move centromere any more since centromere can be simply a marker, and it's forced 
			to the end of the chromosome.

			std::cout << "\nEnter the centromere position of chromosome " << nCurrChromosome << " (in same unit as chromosome size. -1 default: at 50% of chromosome, enter 0 set it to the same as the previous): ";
			std::cin >> nCentromerePos;
			 
			if  (nCentromerePos == -1) {
				pCentromerePos[nCurrChromosome-1] = nChromosomeSize / 2 ;
				
			}
			
			else if (nCentromerePos < 0 || nCentromerePos > nChromosomeSize ) {
				cout << "Centromere position out of bound, enter chromosome info again" << endl ;
				continue;
			}
			
			else if (nCentromerePos == 0 && nCurrChromosome > 1) {
				pCentromerePos[nCurrChromosome-1] =  pCentromerePos[nCurrChromosome-2];
			}
			
			else {
				pCentromerePos[nCurrChromosome-1] =  nCentromerePos;
			}

			*/

			pCentromerePos[nCurrChromosome-1] = nChromosomeSize;
			
			nCurrChromosome++;
			if (nCurrChromosome > nHaploidChromosomeNum ) {
				break;
			}
		}
		
		while(true) { // Get total number of markers
		
			std::cout << "\nHow many markers in total? (this will be distributed to all chromosomes weighted by their lengths): ";
			std::cin >> nTotalMarkers;
			
			if (nTotalMarkers <= 0 ) {
				continue; // 
			}
			
			break;
			
		}
		
		while(true) { // Get average freq of the more abundant allele in pop 1
		
			std::cout << "\nFor parental population 1, what is the average frequency of the more abundant allele [ 0.5~1.0 ] : ";
			std::cin >> nPop1AvgAlleleFreq;
			
			if (nPop1AvgAlleleFreq < 0.5 || nPop1AvgAlleleFreq > 1.0) {
				continue; // 
			}
			
			break;
		}
		
		while(true) { //  Get std. dev of the more abundant allele in pop 1
		
			std::cout << "\nFor parental population 1, what is the standard deviation of the more abundant allele [ 0~0.5 ] : ";
			std::cin >> nPop1AvgAlleleFreqStdDev;
			
			if (nPop1AvgAlleleFreqStdDev < 0 || nPop1AvgAlleleFreqStdDev > 0.5) {
				continue; // 
			}
			
			break;
		}
		
		while(true) { // Get average freq of the more abundant allele in pop 2
		
			std::cout << "\nFor parental population 2, what is the average frequency of the more abundant allele [ 0.5~1.0 ] : ";
			std::cin >> nPop2AvgAlleleFreq;
			
			if (nPop2AvgAlleleFreq < 0.5 || nPop2AvgAlleleFreq > 1.0) {
				continue; // 
			}
			
			break;
		}
		
		while(true) { //  Get std. dev of the more abundant allele in pop 2
		
			std::cout << "\nFor parental population 2, what is the standard deviation of the frequency of the more abundant allele [ 0~0.5 ] : ";
			std::cin >> nPop2AvgAlleleFreqStdDev;
			
			if (nPop2AvgAlleleFreqStdDev < 0 || nPop2AvgAlleleFreqStdDev > 0.5) {
				continue; // 
			}
			
			break;
		}
		
		cout << "Output File name : ";
		
		while(true) { //  Get std. dev of the more abundant allele in pop 2
		
			
		
			
			getline(cin, sOutFile) ;
			
			
			if (sOutFile.size() > 3 ) {
				break;
			}
			
			
			continue;
		}
		
		
		// now prepare file:
		//Figure out how to distribute the markers
		
		ofstream fOutFile;
		fOutFile.open(sOutFile.c_str());
		fOutFile << "HaploidChromosomeNum = " << nHaploidChromosomeNum << endl; //Write haploid n.
		fOutFile << "RandomSeed = " << nRandSeed << endl; //Write random seed n.
		fOutFile << "Pop1AvgAlleleFreq = " << nPop1AvgAlleleFreq << endl; //
		fOutFile << "Pop1AvgAlleleFreqStdev = " << nPop1AvgAlleleFreqStdDev << endl; //
		fOutFile << "Pop2AvgAlleleFreq = " << nPop2AvgAlleleFreq << endl; //
		fOutFile << "Pop2AvgAlleleFreqStdev = " << nPop2AvgAlleleFreqStdDev << endl; //
			/*	
			for (int i=0;i<=999;i++)
			{
				myfile << NormalExt(0.8, 0.1 , 0.5, 1.0) << "\n"; // draw from normal distribution, mean=0.8, stdev=0.1, >=0.5, <= 1.0
			
			}
			*/
			
			
		
		double nSumLens = fnSum(pChromosomeSizes , nHaploidChromosomeNum); // get sum of chromosome lengths
		double nLargestChromSize = 0;
		
		for (int i=0;  i < nHaploidChromosomeNum; i++) { // generate markers chromosome by chromosome
			
			double nChrLen = pChromosomeSizes[i];
			int	   nNumMarkers = (int) ceil( ((double) nTotalMarkers) * (nChrLen / nSumLens));
			//fOutFile << nNumMarkers <<endl;
			
			nLargestChromSize = (nLargestChromSize < nChrLen)? nChrLen : nLargestChromSize;
			
			//cout<< "MarkerNum:" << nNumMarkers << endl;
			
			fOutFile << ":chr " << (i+1) << "\tlen = " << nChrLen << "\tcentromere = " << pCentromerePos[i] << endl; // chromosome number
			fOutFile << "\nNum\tPosition Percentage\tPosition Abs.\tFreq. Pop 1\tFreq. Pop 2\n" ;
			
			double * pPos = new double[nNumMarkers];
			for (int j=0;j<nNumMarkers;j++) {
				pPos[j] = UniformGen.Next(); // create random set of uniformly distributed markers
			}
			
			//cout<< "MarkerNum:" << nNumMarkers << endl;
			//Sort markers according to position on chromosome:
			
			qsort (pPos , nNumMarkers , sizeof(double), fnCompare);
			
			//cout<< "MarkerNum:" << pPos[0] << endl;
			
			for (int j=0;j<nNumMarkers;j++) {
				fOutFile << (j + 1);
				fOutFile << '\t';
				fOutFile << pPos[j] ;
				fOutFile << '\t';
				fOutFile << pPos[j] * nChrLen ;
				fOutFile << '\t';
				fOutFile << NormalExt(nPop1AvgAlleleFreq ,nPop1AvgAlleleFreqStdDev , 0.5, 1.0);
				fOutFile << '\t';
				fOutFile <<  NormalExt( 1.0 - nPop2AvgAlleleFreq ,nPop2AvgAlleleFreqStdDev , 0, 0.5);
				fOutFile << endl; 
			}
			
			//cout<< "MarkerNum:" << nNumMarkers << endl;
			
			delete[] pPos;
		}
		
		fOutFile.close();
		
		
		while(true) { //  Get std. dev of the more abundant allele in pop 2
		
			char sRs;
			std::cout << "\nDo you want to generate random recombination probabilities along these chromosomes? y - yes, n - no : ";
			std::cin >> sRs;
			
			if (sRs == 'y') {
				
				UIGenerateRecombinationFreqMap(nHaploidChromosomeNum , nLargestChromSize, pChromosomeSizes, pCentromerePos); 
			}
			
			break;
		}

}

void UIGenerateRecombinationFreqMap(int nChromosomes, double nLargestChromSize, double * pChromosomeSizes, double * pCentromerePos) { //UI generating recombination map
		
		//double nLogMaxAvgRecombRate; // Log(Max(recombination rate));
		char nRate; //1 - uniform, 2 - parabola
		double nTotalExpectedRecombPerMeiosis; // The total expected recombination events per meiosis in the whole genome
		double nParabolaA; // shape of male parabola
		double nParabolaAFemale; //shape of female parabola
		//double nFemCorrA ; // linear function to correct for female rates: y=ax + b; negative values make rates higher near centromere and lower at telomere and vice versa
		//double nFemCorrB ; // This should be larger than one b/c females have higher recombination rate in general. usually 3 
		double nFemFold; //Fold factor compared to max
		double nTotalSamplePoint; //How often to put a point. 
		double nSumLens = fnSum(pChromosomeSizes , nChromosomes); // get sum of chromosome lengths
		
		while(true) { //  Get std. dev of the more abundant allele in pop 2
		
			//std::cout << "\nThe logarithm of the max average recombination rate at the telomere of the longest chromosome per chromosome length unit per gamete (< 0) : ";
			
			std::cout << "\nThe total expected recombination event in the whole genome per meiosis (>= 0), set this to be 4n should be fine : ";
			std::cin >> nTotalExpectedRecombPerMeiosis;
			
			if (nTotalExpectedRecombPerMeiosis < 0 ) {
				std::cout << "Input out of bound, retry." <<endl;
				continue; // 
			}
			
			break;
		}
		/*
		while(true) { //  Get std. dev of the more abundant allele in pop 2
		
			std::cout << "\nDistribution of recombination frequencies 1 - Uniform, 2 - Parabola ... ";
			std::cin >> nRate;
			
			if (nRate < '1' || nRate > '2' ) {
				std::cout << "Input out of bound, retry." <<endl;
				continue; // 
			}
			
			break;
		}
		*/
		nRate = '1';
		
		if (nRate == '2') { 
				while(true) { //  Get std. dev of the more abundant allele in pop 2
		
					std::cout << "\nShape of the male parabola (0, 5] : ";
					std::cin >> nParabolaA;
			
					if (nParabolaA <= 0 || nParabolaA > 5 ) {
						std::cout << "Input out of bound, retry." <<endl;
						continue; // 
					}
			
					break;
				}
		
	
				while(true) { //  Get std. dev of the more abundant allele in pop 2
		
					std::cout << "\nShape of the female parabola (0, 5] : ";
					std::cin >> nParabolaAFemale;
			
					if (nParabolaAFemale <= 0 || nParabolaAFemale > 5 ) {
						std::cout << "Input out of bound, retry." <<endl;
						continue; // 
					}
			
					break;
				}
		}
		
		if (nRate == '1') { // if uniform, just put two since users can specify Use UseUniformRec = yes
			std::cout << "===========================================" <<endl;
			std::cout  <<endl;
			std::cout << "By default Admixem uses the uniform rate" <<endl;
			std::cout << "In addition to files generated here," <<endl;
			std::cout << "Please make sure you ADD THE FOLLOWING OPTION" <<endl;
			std::cout << "in your MAIN CONFIGURATION file:" <<endl;
			std::cout << endl << "        UseUniformRec  =  yes" <<endl;
			std::cout  <<endl;
			std::cout  <<"In case you really want to use the candidate ";
			std::cout  <<"break point approach, refer to the manual. ";
			std::cout << "===========================================" <<endl;

			nTotalSamplePoint = 20 * nChromosomes;
		}
		else {
			while(true) { //  Get std. dev of the more abundant allele in pop 2
		
				std::cout << "\nHow many sample points to put in total?: ";
				std::cin >> nTotalSamplePoint;
			
				if (nTotalSamplePoint <= 0  ) {
					std::cout << "Input out of bound, retry." <<endl;
					continue; // 
				}
			
				break;
			}
		}
		
		while(true) { //  Get std. dev of the more abundant allele in pop 2
		
			std::cout << "\nFold difference between female:male [0.1,10]: ";
			std::cin >> nFemFold;
			
			if (nFemFold < 0.1 || nFemFold > 10 ) {
				std::cout << "Input out of bound, retry." <<endl;
				continue; // 
			}
			
			break;
		}
		
		/*
		while(true) { //  Get std. dev of the more abundant allele in pop 2
		
			std::cout << "\n A linear model y=a*prob(male)+b will be used to generate the female recomb. rate from the male rates. "<<endl;
			std::cout << "A -b < a < 10:"<<endl;
			std::cin >> nFemCorrA;
			std::cout << "B [1,10]:"<<endl;
			std::cin >> nFemCorrB;
			

			
			if (nFemCorrB < 1 || nFemCorrB > 10 ) {
				std::cout << "Parameter B out of bound, retry." <<endl;
				continue; // 
			}
			
				if (nFemCorrA < -nFemCorrB || nFemCorrA > 10 ) {
				std::cout << "Parameter A out of bound, retry." <<endl;
				continue; // 
			}
			
			break;
		}
		*/
		
		
	cout << "Output File name : ";
	string sOutFile;
		
		while(true) { //  Get std. dev of the more abundant allele in pop 2
		
			
		
			
			getline(cin, sOutFile) ;
			
			
			if (sOutFile.size() > 3 ) {
				break;
			}
			
			
			continue;
		}
		
		
		// now prepare file:
		//Figure out how to distribute the markers
		
		ofstream fOutFile;
		fOutFile.open(sOutFile.c_str());
		
		
		for (int i=0;i< nChromosomes;i++) {
			

			
			int nChrIndex =i;
			double nChrLen = pChromosomeSizes[nChrIndex];
			double nCentromerePos =  pCentromerePos[nChrIndex];	
			//double nMaxAvgRecombRate = pow(10, nLogMaxAvgRecombRate);	
			//nLongestArmExpectedRecombPerMeiosis
			
			int	   nNumSamplePoints = (int) ceil( ((double) nTotalSamplePoint) * (nChrLen / nSumLens));

			double nArm1Len =  nCentromerePos;
			double nArm2Len =  nChrLen - nCentromerePos;
				
			fOutFile << ":chr " << (nChrIndex+1) << endl;
			
			double nExpectedMaleRecPerMeiosisArm1 = nTotalExpectedRecombPerMeiosis * ( nArm1Len/nSumLens) / 2;
			double nExpectedMaleRecPerMeiosisArm2 = nTotalExpectedRecombPerMeiosis * ( nArm2Len/nSumLens) / 2;
			double nExpectedFemaleRecPerMeiosisArm1 = nTotalExpectedRecombPerMeiosis * nFemFold * ( nArm1Len/nSumLens) / 2;
			double nExpectedFemaleRecPerMeiosisArm2 = nTotalExpectedRecombPerMeiosis * nFemFold * ( nArm2Len/nSumLens) / 2;
			double nExpectedRecPerMeiosisArm1 = ( nExpectedMaleRecPerMeiosisArm1 + nExpectedFemaleRecPerMeiosisArm1 ) / 2;
			double nExpectedRecPerMeiosisArm2 = ( nExpectedMaleRecPerMeiosisArm2 + nExpectedFemaleRecPerMeiosisArm2 ) / 2;
			
			//fOutFile << "\tmale"<<"\tmale_accumulative\tmale*max\tmale*max_accumulative" << "\tfemale\tfemale_accumulative\tfemale*max\tfemale*max_accumulative" << endl;
			fOutFile << ":ExpectedMaleRecPerMeiosisArm1 = " << nExpectedMaleRecPerMeiosisArm1 << endl;
			fOutFile << ":ExpectedMaleRecPerMeiosisArm2 = " << nExpectedMaleRecPerMeiosisArm2 << endl;
			fOutFile << ":ExpectedFemaleRecPerMeiosisArm1 = " << nExpectedFemaleRecPerMeiosisArm1 << endl;
			fOutFile << ":ExpectedFemaleRecPerMeiosisArm2 = " << nExpectedFemaleRecPerMeiosisArm2 << endl;
			fOutFile << "\tmale"<<"\tmale_accumulative" << "\tfemale\tfemale_accumulative\tMale_recombination_fraction\tFemale_recombination_fraction\tAvg_recombination_fraction\tKosambi_male_interval\tKosambi_male_pos\tKosambi_female_interval\tKosambi_female_pos\tKosambi_both_interval\tKosambi_both_pos" << endl;
		
			std::map<double, double *>  mSamplePoints;

			
			for (int j=0; j < nNumSamplePoints; j++) {
			
				double nPos = nChrLen * UniformGen.Next(); // random position
				
				double nDistanceFromCentromere = ( nPos - nCentromerePos ) / nLargestChromSize; //  distance from centromere standarized by the longest chromosome
				nDistanceFromCentromere = (nDistanceFromCentromere >= 0 )? nDistanceFromCentromere: -nDistanceFromCentromere;
				
				double nArmLen = (nPos > nCentromerePos)? nChrLen - nCentromerePos : nCentromerePos;
				double nDistanceFromCentromereOnArm = (nPos > nCentromerePos)? (nPos - nCentromerePos)/nArmLen : (nCentromerePos - nPos)/nArmLen ;
				
				double nProbAtPosMale =  (nRate=='1')? 0.1 : pow(nParabolaA * nDistanceFromCentromere, 2);
				//double nProbAtPosFemale = nProbAtPosMale * nFemCorrA  + nFemCorrB;
				
				double nProbAtPosFemale = (nRate=='1')? 0.1 * nFemFold : pow(nParabolaAFemale * nDistanceFromCentromere, 2) * nFemFold ; 
				
				double * pProbs = new double[4];
				
				pProbs[0] = NormalExt(nProbAtPosMale , nProbAtPosMale / 2, 0, 1);
				//pProbs[1] = nMaxAvgRecombRate * pProbs[0] ;
				pProbs[2] = NormalExt(nProbAtPosFemale , nProbAtPosFemale / 2, 0, 1);
				//pProbs[3] = nMaxAvgRecombRate * pProbs[2] ;
				
				if (mSamplePoints.find(nPos) == mSamplePoints.end()) {
					mSamplePoints[nPos] = pProbs ;
				}
				else {
					j--;
					continue; // won't work because same position already exists.
				}
				
			}
			
			double nMaleArm1Accum = 0;
			double nFemaleArm1Accum = 0;
			double nMaleArm2Accum = 0;
			double nFemaleArm2Accum = 0;


			
			for (std::map<double,double *>::iterator it = mSamplePoints.begin(); it != mSamplePoints.end(); ++it) {
			
				double nMaleAccum, nFemaleAccum;
				if (it->first < nCentromerePos) { //on arm 1
					/*
					nMaleAccum = nMaleArm1Accum += (it->second)[0];
					nFemaleAccum = nFemaleArm1Accum += (it->second)[2];
					*/
					nMaleArm1Accum += (it->second)[0];
					nFemaleArm1Accum += (it->second)[2];
				}
				else {
					/*
					nMaleAccum = nMaleArm2Accum += (it->second)[0];
					nFemaleAccum = nFemaleArm2Accum += (it->second)[2];
					*/
					 nMaleArm2Accum += (it->second)[0];
					 nFemaleArm2Accum += (it->second)[2];
				}
				
				//fOutFile << setprecision (10) << it->first << '\t' << (it->second)[0] << '\t' << nMaleAccum << '\t' << (it->second)[1] << '\t' << nMaleAccum * nMaxAvgRecombRate  << '\t' << (it->second)[2] << '\t' << nFemaleAccum << '\t' << (it->second)[3] << '\t' << nFemaleAccum * nMaxAvgRecombRate << endl;
			}


			double nScaledMaleArm1Accum = 0;
			double nScaledFemaleArm1Accum = 0;
			double nScaledMaleArm2Accum = 0;
			double nScaledFemaleArm2Accum = 0;

			double nScaledMaleProb = 0;
			double nScaledFemaleProb = 0;
			
			double nPrevScaledMaleArm1Accum =0;
			double nPrevScaledFemaleArm1Accum =0;
			double nPrevScaledMaleArm2Accum =0;
			double nPrevScaledFemaleArm2Accum =0;
			
			double nKosambiMapDistance_male = 0;
			double nKosambiMapDistance_female = 0;
			double nKosambiMapDistance_both = 0;			
			
			double nKosambiMapPos_male = 0;
			double nKosambiMapPos_female = 0;
			double nKosambiMapPos_both = 0;
			
			double nRecombFraction_male = 0;
			double nRecombFraction_female = 0;
			double nRecombFraction_both = 0;
			
	
			for (std::map<double,double *>::iterator it = mSamplePoints.begin(); it != mSamplePoints.end(); ++it) {
			
				double nMaleAccum, nFemaleAccum;
				if (it->first < nCentromerePos) { //on arm 1
					nScaledMaleProb = (it->second)[0] / nMaleArm1Accum;
					nScaledFemaleProb = (it->second)[2] / nFemaleArm1Accum;
					nMaleAccum = nScaledMaleArm1Accum += nScaledMaleProb;
					nFemaleAccum = nScaledFemaleArm1Accum += nScaledFemaleProb;
					
					 nRecombFraction_male = nExpectedMaleRecPerMeiosisArm1 * nScaledMaleProb  ;
					 nKosambiMapPos_male += nKosambiMapDistance_male = 0.25 * log( (1+2 * nRecombFraction_male)/(1-2 * nRecombFraction_male) ) ;
					 nRecombFraction_female = nExpectedFemaleRecPerMeiosisArm1 * nScaledFemaleProb ;
					 nKosambiMapPos_female += nKosambiMapDistance_female = 0.25 * log( (1+2 * nRecombFraction_female)/(1-2 * nRecombFraction_female) ) ;
		
					
					
				}
				else {
					nScaledMaleProb = (it->second)[0] / nMaleArm2Accum;
					nScaledFemaleProb = (it->second)[2] / nFemaleArm2Accum;
					nMaleAccum = nScaledMaleArm2Accum += nScaledMaleProb;
					nFemaleAccum = nScaledFemaleArm2Accum += nScaledFemaleProb;
					
					nRecombFraction_male = nExpectedMaleRecPerMeiosisArm2 * nScaledMaleProb  ;
					 nKosambiMapPos_male += nKosambiMapDistance_male = 0.25 * log( (1+2 * nRecombFraction_male)/(1-2 * nRecombFraction_male) ) ;
					 nRecombFraction_female = nExpectedFemaleRecPerMeiosisArm2 * nScaledFemaleProb  ;
					 nKosambiMapPos_female += nKosambiMapDistance_female = 0.25 * log( (1+2 * nRecombFraction_female)/(1-2 * nRecombFraction_female) ) ;
					 
					
				}
				
				nRecombFraction_both = (nRecombFraction_male + nRecombFraction_female)/2;
				nKosambiMapPos_both += nKosambiMapDistance_both = 0.25 * log( (1+2 * nRecombFraction_both)/(1-2 * nRecombFraction_both) ) ;
					 

				
				fOutFile << setprecision (10) << it->first << '\t' << nScaledMaleProb << '\t' << nMaleAccum << '\t' << nScaledFemaleProb << '\t' << nFemaleAccum;
				fOutFile << setprecision (10) << '\t' << nRecombFraction_male << '\t' << nRecombFraction_female << '\t' << nRecombFraction_both << '\t' << nKosambiMapDistance_male << '\t' << nKosambiMapPos_male << '\t' << nKosambiMapDistance_female << '\t' << nKosambiMapPos_female << '\t' << nKosambiMapDistance_both << '\t' << nKosambiMapPos_both << endl;
			
				double nPrevScaledMaleArm1Accum = nScaledMaleArm1Accum;
				double nPrevScaledFemaleArm1Accum = nScaledFemaleArm1Accum;
				double nPrevScaledMaleArm2Accum = nScaledMaleArm2Accum;
				double nPrevScaledFemaleArm2Accum = nScaledFemaleArm2Accum;
				
			
			}
		
		
		
		}
		
		fOutFile.close();
}


