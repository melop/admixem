#include "config.h"
#include "spline.h"
#include <cstring>
#ifdef _OPENMP
 #include <omp.h>
#endif
SimulationConfigurations SimulConfig; //Global configuration object
extern Normal NormalGen; 
extern Uniform UniformGen;
extern vector< Binomial * > vBinomGen;
extern vector< Normal * > vNomGen;
extern vector< Uniform * > vUniformGen;

using namespace tk;

/*
int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}
*/

#define MAXLINE    250000

void MarkerConfigurations::LoadFromFile(string szConfigFile) {

	FILE *pConfigFile;

	printf("Start loading marker file %s ...\n", szConfigFile.c_str());

	_szConfigFilename = szConfigFile; // save name of config file

	pConfigFile = fopen(szConfigFile.c_str(), "r");

	if (!pConfigFile) {
		printf("%s\n", "Cannot open marker file ...");
		throw "Cannot open marker file!";
	}


	char szBuffer[2048];

	if (fgets(szBuffer, 2048, pConfigFile) == NULL) {
		throw "Error reading marker file!";
	}

	int nHaploidChrNum;

	sscanf(szBuffer, "HaploidChromosomeNum =  %d", &nHaploidChrNum); // Read in chromosome number

	if (nHaploidChrNum == NULL) {
		throw "Haploid Chromosome Number in marker file incorrect!";
	}

	this->_nHaploidChrNum = nHaploidChrNum; //save the haploid chromosome number

	//initialize arrays:
	this->pChrLen = new double[nHaploidChrNum];
	this->pPaternalTransBias = new double[nHaploidChrNum];
	this->pCentromerePos = new double[nHaploidChrNum];

	vector<int> vIndexCounters; //init index counter
	vMarkerIndex.clear();

	for (int i=0;i<_nHaploidChrNum;i++) {
		vIndexCounters.push_back(0);
		map<double, int> oTempMap;
		vMarkerIndex.push_back(oTempMap);
	}
	
	for (int nPop=0;nPop<2;nPop++) 
	{
		vMarkerFreq.insert( vMarkerFreq.end(), new map<double, double>[nHaploidChrNum] );
	}

	int nCurrChr = 0;

	int nNum;
	double nAbsPos, nFreqPop1, nFreqPop2;
	//char sAbsPos[100], sFreqPop1[100], sFreqPop2[100];
	bool bIgnoreMarkerFreq = (SimulConfig.GetConfig("IgnoreMarkerFreq") == "yes");

	while(fgets(szBuffer, 2048, pConfigFile) != NULL) 
	{

		//sscanf(szBuffer, "%d %*f %f %f %f", &nNum, &nAbsPos, &nFreqPop1, &nFreqPop2);

		if (szBuffer[0] == ':') { // it's a chr definition line

			double nChrLen;
			double nCentromerePos;
			double nPaternalTransBias; // when transmitting gametes, the bias towards paternal chr. 0.5 is neutral, 0 is maternal

			int nItemsRead = sscanf(szBuffer, ":chr %d len = %lf centromere = %lf paternaltransbias = %lf", &nCurrChr, &nChrLen , &nCentromerePos, &nPaternalTransBias);

			if (nItemsRead < 4) {
				nPaternalTransBias = 0.5;
			}
			//nCurrChr = (int)strtol( sCurrChr , NULL, 10);
			//nChrLen = strtod( sChrLen, NULL );
			//nCentromerePos = strtod( sCentromerePos , NULL);

			
			nCurrChr--; 

			//printf("%d %f %f", nCurrChr, nChrLen , nCentromerePos );

			if (nCurrChr > nHaploidChrNum) 
			{
				throw "The chromosome number is larger than total chromosome number!";
			}

			// save chromosome information:
			pChrLen[nCurrChr] = nChrLen;
			pPaternalTransBias[nCurrChr] = nPaternalTransBias;
			nTotalGenomeLen += nChrLen;
			pCentromerePos[nCurrChr] = nCentromerePos;
			printf("Reading markers on chromosome %d , length=%f, centromere=%f paternalTransBias=%f ...\n", (nCurrChr + 1), nChrLen, nCentromerePos,nPaternalTransBias);

		}
		else if (sscanf(szBuffer, "%d %*lf %lf %lf %lf", &nNum, &nAbsPos, &nFreqPop1, &nFreqPop2)==4) {

			if (vMarkerFreq[0][nCurrChr].find(nAbsPos) == vMarkerFreq[0][nCurrChr].end()) { // make sure the position is not occupied
				
				vMarkerFreq[0][nCurrChr][nAbsPos] = bIgnoreMarkerFreq? 1: nFreqPop1;
				vMarkerFreq[1][nCurrChr][nAbsPos] = bIgnoreMarkerFreq? 0: nFreqPop2;

				vMarkerIndex[nCurrChr][nAbsPos] = vIndexCounters.at(nCurrChr);
				vIndexCounters.at(nCurrChr) += 1;
			}
			else
			{

				printf("Warning: repetitive definition of position  %f on chromosome %d found. second one will be ignored\n", nAbsPos, nCurrChr );
			}
		}

		else {

		}

		
		/*
		if (sscanf_s(szBuffer, "%d %*f %f %f %f", &nNum, &nAbsPos, &nFreqPop1, &nFreqPop2)!=4) { // if this is not a data line try parsing: :chr 1	len = 3.3e+007	centromere = 1.1e+007
			
			double nChrLen;
			double nCentromerePos;
			
			//sscanf(szBuffer, ":chr %d	len = %f	centromere = %f", &nCurrChr, &nChrLen , &nCentromerePos);

			if (sscanf_s(szBuffer, ":chr %d len = %f centromere = %f", &nCurrChr, &nChrLen , &nCentromerePos) != 3) {
				continue;
			}

			nCurrChr--; 

			//printf("%d %f %f", nCurrChr, nChrLen , nCentromerePos );

			if (nCurrChr > nHaploidChrNum) 
			{
				throw "The chromosome number is larger than total chromosome number!";
			}

			// save chromosome information:
			pChrLen[nCurrChr] = nChrLen;
			pCentromerePos[nCurrChr] = nCentromerePos;

		}
		else { //this is a dataline!
			vMarkerFreq[0][nCurrChr][nAbsPos] = nFreqPop1;
			vMarkerFreq[1][nCurrChr][nAbsPos] = nFreqPop2;
			//printf("%d %f %f %f\n", nNum, nAbsPos, nFreqPop1, nFreqPop2 );
		}

		*/

	}

	fclose(pConfigFile);

	//printf("%f",vMarkerFreq[0][0].begin()->first );

};

void MarkerConfigurations::CalculateMapDistances() { // calculate map distances for each marker
	
	vector<double> * pvMaleMapDistances = ((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->GetBreakPointMapDistances(1);

	//if (!pvMaleMapDistances) return; //The recombination mode is not predetermined mode , not available;

	vector<double> * pvFemaleMapDistances = ((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->GetBreakPointMapDistances(0);
	vector<double> * pvAvgMapDistances = ((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->GetBreakPointMapDistances(2);
	vector<double> * pvBreakpointSamplePositions = ((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->GetBreakPointSamplePositions();
	vAbsToMapDistanceMale.clear();
	vAbsToMapDistanceFemale.clear();
	vAbsToMapDistanceAvg.clear();


	for (int nChr=0 ; nChr < this->_nHaploidChrNum; nChr++ ) {

		double nMaleExpBreakpoints = ((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->HowManyBreakpointsOnArm(true, nChr, 1) + 
										((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->HowManyBreakpointsOnArm(true, nChr, 2);
		double nFemaleExpBreakpoints = ((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->HowManyBreakpointsOnArm(false, nChr, 1) + 
										((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->HowManyBreakpointsOnArm(false, nChr, 2);
		double nAvgExpBreakpoints = (nMaleExpBreakpoints + nFemaleExpBreakpoints) / 2;

		double nChrLen = this->pChrLen[nChr];

		vector<double> * pvBreakpointSamplePositionsOnChr = &pvBreakpointSamplePositions[nChr];
		vector<double> * pvMaleMapDistancesOnChr = &pvMaleMapDistances[nChr];
		vector<double> * pvFemaleMapDistancesOnChr = &pvFemaleMapDistances[nChr];
		vector<double> * pvAvgMapDistancesOnChr = &pvAvgMapDistances[nChr];

		map<double, int> * pmMarkerIndexOnChr = &(vMarkerIndex.at(nChr));
		map<double, double> mAbsToMapDistanceOnChrMale, mAbsToMapDistanceOnChrFemale, mAbsToMapDistanceOnChrAvg ;

		double nPosLastSample = 0;

		for (map<double,int>::iterator itMarker=pmMarkerIndexOnChr->begin(); itMarker!= pmMarkerIndexOnChr->end(); ++itMarker) 
		{

			double nPos = (*itMarker).first;
			double nAbsSampleDistance, nMapMaleSampleDistance, nMapFemaleSampleDistance, nMapAvgSampleDistance, nMapMaleDistance, nMapFemaleDistance, nMapAvgDistance;

			if (((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->IsUseUniform()) { //if user specified using uniform rate, then no need to extrapolate.
				nMapMaleDistance = (nPos / nChrLen ) * nMaleExpBreakpoints;
				nMapFemaleDistance = (nPos / nChrLen ) * nFemaleExpBreakpoints;
				nMapAvgDistance = (nPos / nChrLen ) * nAvgExpBreakpoints;
			}
			else {

				vector<double>::iterator itClosestFollowingBreakpoint = upper_bound(pvBreakpointSamplePositionsOnChr->begin(), pvBreakpointSamplePositionsOnChr->end(), nPos);
				if (itClosestFollowingBreakpoint ==pvBreakpointSamplePositionsOnChr->end() ) {
					int nIndexLastSample = pvBreakpointSamplePositionsOnChr->size()-1;
					int nIndexSecondToLastSample = nIndexLastSample - 1;
					double nPosLastSample = pvBreakpointSamplePositionsOnChr->at(nIndexLastSample);//Last breakpoint sample
					double nPosSecondToLastSample = pvBreakpointSamplePositionsOnChr->at(nIndexSecondToLastSample);
					nAbsSampleDistance = (nPosLastSample - nPosSecondToLastSample);
					nMapMaleSampleDistance = pvMaleMapDistancesOnChr->at(nIndexLastSample);// - pvMaleMapDistancesOnChr->at(nIndexSecondToLastSample);
					nMapFemaleSampleDistance = pvFemaleMapDistancesOnChr->at(nIndexLastSample);// - pvFemaleMapDistancesOnChr->at(nIndexSecondToLastSample);
					nMapAvgSampleDistance = pvAvgMapDistancesOnChr->at(nIndexLastSample);// - pvAvgMapDistancesOnChr->at(nIndexSecondToLastSample);

					//extrapolate
					double nAbsDistanceMarkerToPrevMarker = nPos - nPosLastSample;
					double nPercentFromPrevSample = ( nAbsDistanceMarkerToPrevMarker / nAbsSampleDistance);
					nMapMaleDistance = nMapMaleSampleDistance  * nPercentFromPrevSample;
					nMapFemaleDistance = nMapFemaleSampleDistance * nPercentFromPrevSample;
					nMapAvgDistance = nMapAvgSampleDistance * nPercentFromPrevSample;
				}
				else {
					if (itClosestFollowingBreakpoint != pvBreakpointSamplePositionsOnChr->begin())
					{
						vector<double>::iterator itClosestPrevBreakpoint = itClosestFollowingBreakpoint - 1;
						size_t nIndexClosestFollowingBreakpoint = distance(pvBreakpointSamplePositionsOnChr->begin(),itClosestFollowingBreakpoint );
						size_t nIndexClosestPrevBreakpoint = nIndexClosestFollowingBreakpoint - 1;

						nAbsSampleDistance = ( *itClosestFollowingBreakpoint - *itClosestPrevBreakpoint);
						nMapMaleSampleDistance = pvMaleMapDistancesOnChr->at(nIndexClosestFollowingBreakpoint);// - pvMaleMapDistancesOnChr->at(nIndexClosestPrevBreakpoint);
						nMapFemaleSampleDistance = pvFemaleMapDistancesOnChr->at(nIndexClosestFollowingBreakpoint);// - pvFemaleMapDistancesOnChr->at(nIndexClosestPrevBreakpoint);
						nMapAvgSampleDistance = pvAvgMapDistancesOnChr->at(nIndexClosestFollowingBreakpoint);// - pvAvgMapDistancesOnChr->at(nIndexClosestPrevBreakpoint);


					}
					else {
						nAbsSampleDistance = *itClosestFollowingBreakpoint;
						nMapMaleSampleDistance = pvMaleMapDistancesOnChr->at(0);
						nMapFemaleSampleDistance = pvFemaleMapDistancesOnChr->at(0);
						nMapAvgSampleDistance = pvAvgMapDistancesOnChr->at(0);

					}



					//extrapolate
					double nAbsDistanceMarkerToPrevMarker = nPos - nPosLastSample;;
					double nPercentFromPrevSample = (nAbsDistanceMarkerToPrevMarker  / nAbsSampleDistance);

					nMapMaleDistance = nMapMaleSampleDistance * nPercentFromPrevSample;
					nMapFemaleDistance = nMapFemaleSampleDistance * nPercentFromPrevSample;
					nMapAvgDistance = nMapAvgSampleDistance * nPercentFromPrevSample;
				}

				nPosLastSample = nPos;
				/*
				if (nMapMaleDistance == 0.0) {
					printf("f");
				}
				*/
			}
			mAbsToMapDistanceOnChrMale.insert(pair<double, double>(nPos, nMapMaleDistance * 100));
			mAbsToMapDistanceOnChrFemale.insert(pair<double, double>(nPos, nMapFemaleDistance * 100 ));
			mAbsToMapDistanceOnChrAvg.insert(pair<double, double>(nPos, nMapAvgDistance * 100));
		}

		this->vAbsToMapDistanceAvg.push_back(mAbsToMapDistanceOnChrAvg);
		this->vAbsToMapDistanceMale.push_back(mAbsToMapDistanceOnChrMale);
		this->vAbsToMapDistanceFemale.push_back(mAbsToMapDistanceOnChrFemale);
	}
	//printf("hello");

}

vector<map<double, double> > * MarkerConfigurations::GetMpMarkerMapDistance(int nMode){
	return (nMode==0)? &this->vAbsToMapDistanceFemale : (nMode==1)? &this->vAbsToMapDistanceMale : &this->vAbsToMapDistanceAvg;
}

double MarkerConfigurations::GetChrToGenomeRatio(int nChr) {
	return this->pChrLen[nChr] / this->nTotalGenomeLen;
}

vector< map<double, int> > * MarkerConfigurations::GetMpMarkerIndex() {
	return &this->vMarkerIndex;
}

int MarkerConfigurations::GetHaploidChromosomeNum() {
	return this->_nHaploidChrNum;
}


MarkerConfigurations::MarkerConfigurations(void * pParentConfig) {
	nTotalGenomeLen = 0;
	this->_pParentConfig = pParentConfig;
}

RecombProbConfigurations::RecombProbConfigurations(void * pParentConfig) {
	this->_pParentConfig = pParentConfig;
}

void RecombProbConfigurations::fnEraseNonMonotonic(vector<double> * vMale, vector<double> * vFemale, vector<double> * vPos) {
	
	double nCurrMale , nCurrFemale , nCurrPos;
	nCurrMale = nCurrFemale = nCurrPos = -1;
	vector< size_t > vToDelete;
	for(size_t i=0;i<vPos->size();i++) {

		if (nCurrMale >= vMale->at(i) || nCurrFemale >= vFemale->at(i) || nCurrPos >= vPos->at(i)) {
			vToDelete.push_back(i);
		}
		nCurrMale = vMale->at(i);
		nCurrFemale = vFemale->at(i);
		nCurrPos = vPos->at(i);
	}

	//start deleting from the back;

	for(vector< size_t >::reverse_iterator it=vToDelete.rbegin(); it!=vToDelete.rend();++it) { //iterate from the end
		vMale->erase(vMale->begin() + *it);
		vFemale->erase(vFemale->begin() + *it);
		vPos->erase(vPos->begin() + *it);

	}

}
void RecombProbConfigurations::LoadFromFile(string szConfigFile) {
	#ifdef _OPENMP
		int nTotalCPUCore =  omp_get_max_threads();//omp_get_num_threads();
		//printf("OpenMP enabled\n");
	#else
		int nTotalCPUCore = 1;
		//printf("OpenMP disabled \n");
	#endif

	FILE *pConfigFile;
	char szBuffer[12048];

	this->_nHaploidChrNum = ((SimulationConfigurations *)this->_pParentConfig)->pMarkerConfig->GetHaploidChromosomeNum(); // save a copy of the haploid chromosome number
	//Initialize arrays:
	this->_nExpectedMaleRecPerMeiosisArm1 = new double[_nHaploidChrNum];
	this->_nExpectedMaleRecPerMeiosisArm2 = new double[_nHaploidChrNum];
	this->_nExpectedFemaleRecPerMeiosisArm1 = new double[_nHaploidChrNum];
	this->_nExpectedFemaleRecPerMeiosisArm2 = new double[_nHaploidChrNum];

	this->pvBreakpointSamplePositions = new vector<double>[_nHaploidChrNum];
	this->pvMaleAccuProb = new vector<double>[_nHaploidChrNum];
	this->pvFemaleAccuProb = new vector<double>[_nHaploidChrNum];
	this->pvMaleMapDistance = new vector<double>[_nHaploidChrNum];
	this->pvFemaleMapDistance = new vector<double>[_nHaploidChrNum];
	this->pvAvgMapDistance = new vector<double>[_nHaploidChrNum];
	this->pvMaleRecArm1Gen = new vector< Poisson * >[nTotalCPUCore];
	this->pvMaleRecArm2Gen = new vector< Poisson * >[nTotalCPUCore];
	this->pvFemaleRecArm1Gen = new vector< Poisson * >[nTotalCPUCore];
	this->pvFemaleRecArm2Gen = new vector< Poisson * >[nTotalCPUCore];
	
	
	printf("Start loading marker probability file %s ...\n", szConfigFile.c_str());

	if (((SimulationConfigurations *)this->_pParentConfig)->GetConfig("UseUniformRec") == "yes") {
		printf("Instructed by configuration to use uniform distribution, will ignore candidate breakpoints in recombination file. But expected number of recombinations will still be used.\n");
		this->_bUseUniform = true;
		
	}
	else {
		this->_bUseUniform = false;
	}

	_szConfigFilename = szConfigFile; // save name of config file

	pConfigFile = fopen(szConfigFile.c_str(), "r");

	if (!pConfigFile) {
		printf("%s\n", "Cannot open marker file ...");
		throw "Cannot open marker file!";
	}



	int nCurrChr;
	char szCmd[100];
	double nVal, nPos, nMaleAccu, nFemaleAccu, nMaleMapDistance, nFemaleMapDistance, nAvgMapDistance;

	long nLastBpArm1=0; // the index of the last break point at arm 1
	long nBpIndex = -1;
	double nPrevProb = 0.0; // previous probability. 

	while(fgets(szBuffer, 12048, pConfigFile) != NULL) 
	{

		//sscanf(szBuffer, "%d %*f %f %f %f", &nNum, &nAbsPos, &nFreqPop1, &nFreqPop2);


		if (szBuffer[0] == ':') { // it's metadata line

			if (szBuffer[1] == 'c') { // it's the chr line
				//parse chromosome number:

				sscanf(szBuffer, ":chr %d", &nCurrChr);

				if (nCurrChr > this->_nHaploidChrNum || nCurrChr <=0 ) {
					throw "Chromosome number in marker probability file incorrect!";
				}

				if (nBpIndex > 0) {
					vLastBpOnArm1.push_back(nLastBpArm1==0? nBpIndex : nLastBpArm1); // this is the previous record.
				}
				nLastBpArm1=0;
				nBpIndex = -1;
				nPrevProb = 0.0; // previous probability. 
				nCurrChr--;
				continue; // go to the next line
			}
			else { // it's a command line

				sscanf(szBuffer, "%s = %lf", szCmd, &nVal);
				
				//printf("%s\n",szCmd);

				if (strcmp(szCmd, ":ExpectedMaleRecPerMeiosisArm1")==0) {
					this->_nExpectedMaleRecPerMeiosisArm1[nCurrChr] = nVal;
					for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) {
						Poisson * pPosGen = new Poisson(nVal);
						this->pvMaleRecArm1Gen[nCPU].push_back(pPosGen);
					}
				}
				if (strcmp(szCmd,":ExpectedMaleRecPerMeiosisArm2")==0) {
					this->_nExpectedMaleRecPerMeiosisArm2[nCurrChr] = nVal;
					for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) {
						Poisson * pPosGen = new Poisson(nVal);
						this->pvMaleRecArm2Gen[nCPU].push_back(pPosGen);
					}
				}

				if (strcmp(szCmd,":ExpectedFemaleRecPerMeiosisArm1")==0) {
					this->_nExpectedFemaleRecPerMeiosisArm1[nCurrChr] = nVal;
					for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) {
						Poisson * pPosGen = new Poisson(nVal);
						this->pvFemaleRecArm1Gen[nCPU].push_back(pPosGen);
					}
				}
				if (strcmp(szCmd,":ExpectedFemaleRecPerMeiosisArm2")==0) {
					this->_nExpectedFemaleRecPerMeiosisArm2[nCurrChr] = nVal;
					for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) {
						Poisson * pPosGen = new Poisson(nVal);
						this->pvFemaleRecArm2Gen[nCPU].push_back(pPosGen);
					}
				}

			}

			continue;

		}

		if (szBuffer[0] == '\t') {
			
			continue; //title line, skip
		}
		else { //parse real data

			if (sscanf(szBuffer, "%lf %*lf %lf %*lf %lf %*lf %*lf %*lf %lf %*lf %lf %*lf %lf %*lf", &nPos, &nMaleAccu, &nFemaleAccu, &nMaleMapDistance, &nFemaleMapDistance, &nAvgMapDistance)!=6) {
				continue; // couldn't parse line.	
			};

			nBpIndex++;
			if (nMaleAccu < nPrevProb ) { //switch in prob. crossing the centromere?
				if (nLastBpArm1 != 0) { // this was already set, the sequence is not monotonic.
					throw("Error: accumulative recombination rate must be monotonic within each arm.\n");
				}

				nLastBpArm1 = nBpIndex - 1;//the previous bp is the last one in arm 1.
				
			}
			this->pvBreakpointSamplePositions[nCurrChr].push_back(nPos);
			this->pvMaleAccuProb[nCurrChr].push_back(nMaleAccu);
			this->pvFemaleAccuProb[nCurrChr].push_back(nFemaleAccu);
			this->pvMaleMapDistance[nCurrChr].push_back(nMaleMapDistance);
			this->pvFemaleMapDistance[nCurrChr].push_back(nFemaleMapDistance);
			this->pvAvgMapDistance[nCurrChr].push_back(nAvgMapDistance);

			nPrevProb = nMaleAccu;
		}


	}

	//check that recombination break points on arms are complete:
	if ( pvMaleRecArm1Gen[0].size() != pvMaleRecArm2Gen[0].size() || pvMaleRecArm2Gen[0].size() != pvFemaleRecArm1Gen[0].size() || pvFemaleRecArm1Gen[0].size() != pvFemaleRecArm2Gen[0].size()) {
		throw(new Exception("The number of Recombination break point defined are not equal.\n"));
	}

	vLastBpOnArm1.push_back(nLastBpArm1==0? nBpIndex : nLastBpArm1); // this is the last chromosome.


	//calculate splines
	for(int nCurrChr=0;nCurrChr < this->_nHaploidChrNum; nCurrChr++) {
		double nChrLen = ((SimulationConfigurations*)this->_pParentConfig)->pMarkerConfig->GetChromosomeLength(nCurrChr);
		tk::spline oMaleSplineArm1, oFemaleSplineArm1,oMaleSplineArm2, oFemaleSplineArm2;
		oMaleSplineArm1.set_boundary(spline::first_deriv, 0.0, spline::first_deriv, this->pvBreakpointSamplePositions[nCurrChr].at(this->vLastBpOnArm1[nCurrChr]), false);
		oFemaleSplineArm1.set_boundary(spline::first_deriv, 0.0, spline::first_deriv, this->pvBreakpointSamplePositions[nCurrChr].at(this->vLastBpOnArm1[nCurrChr]), false);
		oMaleSplineArm2.set_boundary(spline::first_deriv, this->pvBreakpointSamplePositions[nCurrChr].at(this->vLastBpOnArm1[nCurrChr]), spline::first_deriv, nChrLen, false);
		oFemaleSplineArm2.set_boundary(spline::first_deriv, this->pvBreakpointSamplePositions[nCurrChr].at(this->vLastBpOnArm1[nCurrChr]), spline::first_deriv, nChrLen, false);
		size_t nTotalBps = this->pvBreakpointSamplePositions[nCurrChr].size();
		vector<double> vPosArm1, vPosArm2, vAccMaleArm1, vAccMaleArm2, vAccFemaleArm1, vAccFemaleArm2;
		copy(this->pvBreakpointSamplePositions[nCurrChr].begin(), this->pvBreakpointSamplePositions[nCurrChr].begin() + this->vLastBpOnArm1[nCurrChr] + 1 ,back_inserter(vPosArm1));
		copy(this->pvBreakpointSamplePositions[nCurrChr].begin() + this->vLastBpOnArm1[nCurrChr] + 1, this->pvBreakpointSamplePositions[nCurrChr].end() ,back_inserter(vPosArm2));
		copy(this->pvMaleAccuProb[nCurrChr].begin(), this->pvMaleAccuProb[nCurrChr].begin() + this->vLastBpOnArm1[nCurrChr] + 1 ,back_inserter(vAccMaleArm1));
		copy(this->pvMaleAccuProb[nCurrChr].begin() + this->vLastBpOnArm1[nCurrChr] + 1, this->pvMaleAccuProb[nCurrChr].end() ,back_inserter(vAccMaleArm2));
		copy(this->pvFemaleAccuProb[nCurrChr].begin(), this->pvFemaleAccuProb[nCurrChr].begin() + this->vLastBpOnArm1[nCurrChr] + 1 ,back_inserter(vAccFemaleArm1));
		copy(this->pvFemaleAccuProb[nCurrChr].begin() + this->vLastBpOnArm1[nCurrChr] + 1, this->pvFemaleAccuProb[nCurrChr].end() ,back_inserter(vAccFemaleArm2));

		// erase elements in case not monotonic:

		if (vPosArm1.size() > 2) { // must have more than 2 pairs of dots to predict a spline.
			this->fnEraseNonMonotonic(&vAccMaleArm1 , &vAccFemaleArm1, &vPosArm1);
			oMaleSplineArm1.set_points(vAccMaleArm1 , vPosArm1);
			oFemaleSplineArm1.set_points(vAccFemaleArm1 , vPosArm1);
		}

		if (vPosArm2.size() > 2) {
			this->fnEraseNonMonotonic(&vAccMaleArm2 , &vAccFemaleArm2, &vPosArm2);

			oMaleSplineArm2.set_points(vAccMaleArm2 , vPosArm2);
			oFemaleSplineArm2.set_points(vAccFemaleArm2 , vPosArm2);
		}


		this->vMaleAccuProbSplineArm1.push_back(oMaleSplineArm1);
		this->vMaleAccuProbSplineArm2.push_back(oMaleSplineArm2);
		this->vFemaleAccuProbSplineArm1.push_back(oFemaleSplineArm1);
		this->vFemaleAccuProbSplineArm2.push_back(oFemaleSplineArm2);
	}

	///printf("%f %f %f", this->pvBreakpointSamplePositions[0].at(0), this->pvMaleAccuProb[0].at(0), this->pvFemaleAccuProb[0].at(0));

};

bool RecombProbConfigurations::IsUseUniform() {
	return this->_bUseUniform;
}

double RecombProbConfigurations::HowManyBreakpointsOnArm(bool bSex, int nChr, int nArm) {

	return  (true==bSex? (nArm==1? _nExpectedMaleRecPerMeiosisArm1[nChr] : _nExpectedMaleRecPerMeiosisArm2[nChr]) : (nArm==1? _nExpectedFemaleRecPerMeiosisArm1[nChr] : _nExpectedFemaleRecPerMeiosisArm2[nChr]));
}

void RecombProbConfigurations::GetBreakPointsByArm(bool bSex, int nChr, int nArm,vector<double> &vRet) { // return a vector of break points for a given arm of a given chromosome for a given sex
	
		#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
			#ifdef DEBUG
			printf("CPU %d: %s : %d\n", nCPU, "SexSelFormulae # ", nCPU);
			#endif
		#else
			int nCPU = 0;
		#endif
	//double nExpectedPoints = (true==bSex? (nArm==1? _nExpectedMaleRecPerMeiosisArm1[nChr] : _nExpectedMaleRecPerMeiosisArm2[nChr]) : (nArm==1? _nExpectedFemaleRecPerMeiosisArm1[nChr] : _nExpectedFemaleRecPerMeiosisArm2[nChr]));
	//int nBreakPointsToPut = (int)nExpectedPoints; // get integral part
	//nBreakPointsToPut += (UniformGen.Next() <= (nExpectedPoints - (double)nBreakPointsToPut)) ? 1:0;
	Poisson * pPoissonGen = (true==bSex? (nArm==1?  pvMaleRecArm1Gen[nCPU][nChr] : pvMaleRecArm2Gen[nCPU][nChr]) : (nArm==1? pvFemaleRecArm1Gen[nCPU][nChr] : pvFemaleRecArm2Gen[nCPU][nChr]));; 
	int nBreakPointsToPut = pPoissonGen->Next();
	
	double nCentromerePos = ((SimulationConfigurations*)this->_pParentConfig)->pMarkerConfig->GetCentromerePosition(nChr);
	double nChrLen = ((SimulationConfigurations*)this->_pParentConfig)->pMarkerConfig->GetChromosomeLength(nChr);

	double nStart = (nArm==1)? 0.0 : nCentromerePos;
	double nEnd   = (nArm==1)? nCentromerePos : nChrLen;

	//vector<double> * pvAccuProb;
	tk::spline * pvAccuProbSpline;
	size_t nCentromereIndex, nLastIndex, nStartIndex, nEndIndex;

	if (!this->_bUseUniform) {
		//pvAccuProb = true==bSex? pvMaleAccuProb : pvFemaleAccuProb;
		if (bSex) {
			pvAccuProbSpline = nArm==1? &vMaleAccuProbSplineArm1[nChr] : &vMaleAccuProbSplineArm2[nChr];
		}
		else {
			pvAccuProbSpline = nArm==1? &vFemaleAccuProbSplineArm1[nChr] : &vFemaleAccuProbSplineArm2[nChr];
		}
		//nCentromereIndex = std::distance( pvBreakpointSamplePositions[nChr].begin(), lower_bound(pvBreakpointSamplePositions[nChr].begin(),pvBreakpointSamplePositions[nChr].end(), nCentromerePos));
		//nLastIndex = std::distance(pvBreakpointSamplePositions[nChr].begin() , pvBreakpointSamplePositions[nChr].end());
		//nStartIndex = (nArm==1)? 0 : nCentromereIndex + 1;
		//nEndIndex = (nArm==1)? nCentromereIndex : nLastIndex;
	}
	
	//printf("Breakpointstoput for sex %d Chr %d Arm %d: %d\n", bSex, nChr, nArm, nBreakPointsToPut);
	//vector<double>::iterator it_start = lower_bound(pvAccuProb[nChr]begin(), pvAccuProb->end(), nStart);

	for (int i=0;i<nBreakPointsToPut;i++) {

		double nRand = vUniformGen[nCPU]->Next();
		if (!this->_bUseUniform) {
			/*
			std::vector<double>::iterator oLowerBound;
			oLowerBound = lower_bound( pvAccuProb[nChr].begin() + nStartIndex, pvAccuProb[nChr].begin() + nEndIndex, nRand);
			size_t nLowerBoundIndex = std::distance(pvAccuProb[nChr].begin(), oLowerBound);
			if (nLowerBoundIndex == 0 || nLowerBoundIndex==nEndIndex) {
				// the first element
				i--;
				continue; //break point cannot be inserted at the tips
			}
			//vector<double>::iterator oPreviousPoint = oLowerBound - 1;
			double nPos1 = pvBreakpointSamplePositions[nChr].at(nLowerBoundIndex - 1);
			double nPos2 = pvBreakpointSamplePositions[nChr].at(nLowerBoundIndex);
			//interpolate
			vRet.push_back( (nPos1 + nPos2) / 2 );
			*/
			//new code that uses spline interpolation
			vRet.push_back((*pvAccuProbSpline)(nRand));
		} else {
			double nArmLen = nEnd - nStart;
			double nBpPos = (nArm==1)? (nArmLen * nRand) : (nCentromerePos + nArmLen * nRand);
			vRet.push_back(nBpPos);
		}
		
	}

	sort(vRet.begin(), vRet.end());


};

vector<double> * RecombProbConfigurations::GetBreakPointMapDistances(int nMode) {
	if (this->_bUseUniform) {
		return NULL;
	}

	return nMode==1? this->pvMaleMapDistance: (nMode==0)? this->pvFemaleMapDistance:this->pvAvgMapDistance;
}

vector<double> * RecombProbConfigurations::GetBreakPointSamplePositions() {
	return pvBreakpointSamplePositions;
}

map<double, double> * MarkerConfigurations::GetMarkersInfoByPop(int nPopId) {
	return this->vMarkerFreq.at(nPopId - 1);
}

double MarkerConfigurations::GetCentromerePosition(int nChr) {
	return this->pCentromerePos[nChr];
}

double MarkerConfigurations::GetChromosomeLength(int nChr) {
	return this->pChrLen[nChr];
}

double MarkerConfigurations::GetPaternalTransBias(int nChr) {
	return this->pPaternalTransBias[nChr];
}

PhenotypeConfigurations::PhenotypeConfigurations(void * pParentConfig) {

	this->_pParentConfig = pParentConfig;

};

void PhenotypeConfigurations::LoadFromFile(string szConfigFile) 
{
	FILE *pConfigFile;
	char szBuffer[2048];
	char szPhenotypeName[100];
	char szFormula[2048];
	
	#ifdef _OPENMP
		int nTotalCPUCore =  omp_get_max_threads();//omp_get_num_threads();
		//printf("OpenMP enabled\n");
	#else
		int nTotalCPUCore = 1;
		//printf("OpenMP disabled \n");
	#endif
	
	printf("Start loading phenotype file %s ...\n", szConfigFile.c_str());

	_szConfigFilename = szConfigFile; // save name of config file

	pConfigFile = fopen(szConfigFile.c_str(), "r");

	if (!pConfigFile) {
		printf("%s\n", "Cannot open phenotype file ...");
		throw "Cannot open phenotype file!";
	}

	if (fgets(szBuffer, 2048, pConfigFile) == NULL) { // skip the header line
		throw "Cannot read phenotype file!";
	}

	while(fgets(szBuffer, 2048, pConfigFile) != NULL) 
	{
		string sPhenotypeName, sFormula;

		sscanf(szBuffer, "%[^\t\n]	%[^\t\n]", szPhenotypeName, szFormula);
		sPhenotypeName = szPhenotypeName;
		sFormula = szFormula;
		if (sPhenotypeName.length() != 0 && sFormula.length() !=0) {
			this->_mpPhenotypes[sPhenotypeName] = sFormula;
			
			Parser * pF = new Parser("");
			try {
				pF->Evaluate(sFormula);
			}
			catch(std::runtime_error e) {
				if ( string(e.what()) == "Divide by zero" ) {
					
					printf("Warning in formula: %s \nErrMsg:%s\nThis could simply be caused by uninitialized values", sFormula.c_str(), e.what());
					
					
				}
				else {
					printf("Error occured when trying to parse formula: %s \nErrMsg:%s\nPlease check formulae grammar!", sFormula.c_str(), e.what());
					exit(100);
				}
			}

			vector< pair<int, double> > vSymbols;
			vector< string > vDadPheSymbols;
			vector< string > vMomPheSymbols;
			vector<string> vSymbolStrings;
			//list all the symbols used in this formula
			for(map<string,double>::iterator it = pF->symbols_.begin(); it != pF->symbols_.end(); ++it) 
			{	
				//Parse chromosome number and position from the symbol
				//cout << (it->first) << "\n"; //list all keys

				int nSymbolChr;
				double nSymbolPos;
				
				if (strlen(it->first.c_str()) < 6) { // no need to add symbol
					continue;
				}

				string sSymbolLabel = it->first;
				size_t nFirstUS = sSymbolLabel.find_first_of("_");
				size_t nLastUS = sSymbolLabel.find_last_of("_");

				if (sSymbolLabel.find_first_of("Phe") == 0) { //if referring to the phenotype of dad or mum
					string sWhichParent = sSymbolLabel.substr(  nFirstUS+1, 3);
					string sPhenotype = sSymbolLabel.substr(  nLastUS+1 );

					if (sWhichParent == "Dad" ) {
						vDadPheSymbols.push_back(sPhenotype);
						//vSymbolStrings.push_back(sSymbolLabel); //don't add this
						this->_mpIsReferToParentPhe[sPhenotypeName] = true;
					}

					if (sWhichParent == "Mom" ) {
						vMomPheSymbols.push_back(sPhenotype);
						//vSymbolStrings.push_back(sSymbolLabel); //don't add this.
						this->_mpIsReferToParentPhe[sPhenotypeName] = true;

					}

					continue;
				}

				string sTemp = it->first.substr(3);	
				nFirstUS = sTemp.find_first_of("_");
				nLastUS = sTemp.find_last_of("_");

				if (nFirstUS != nLastUS) { // there is _a1 _a2 suffixes, delete the suffix
					sTemp = sTemp.substr(0,nLastUS);
					sSymbolLabel = "chr" + sTemp;				
				}

				sTemp = sTemp.replace(sTemp.find_first_of("_"), 1, " ");
				if (sscanf(sTemp.c_str(), "%d %lf", &nSymbolChr, &nSymbolPos)!=2) 
				{ // see if it's another format with _a1, a2 suffixes
					continue;
				};
				pair<int,double> oSymbolPair( nSymbolChr, nSymbolPos);
				vSymbols.push_back(oSymbolPair);
				vSymbolStrings.push_back(sSymbolLabel); // save label too
			}

			for (int nCPU=0;nCPU<nTotalCPUCore;nCPU++) {
				//map<string , Parser *> _mpPhenotypeFormulae;
				this->_mpmpPhenotypeFormulae[nCPU][sPhenotypeName] = new Parser(sFormula);//pF;
			}
			this->_mpPhenotypeFormulaSymbols[sPhenotypeName] = vSymbols;
			this->_mpPhenotypeFormulaDadSymbols[sPhenotypeName] = vDadPheSymbols;
			this->_mpPhenotypeFormulaMomSymbols[sPhenotypeName] = vMomPheSymbols;
			this->_mpPhenotypeFormulaSymbolStrings[sPhenotypeName] = vSymbolStrings;
		}
	}

	for (std::map<string,string>::iterator it = _mpPhenotypes.begin(); it != _mpPhenotypes.end(); ++it) {
		printf("%s %s\n", it->first.c_str(), it->second.c_str());
	}

	printf("%d phenotype(s) loaded. \n", _mpPhenotypes.size());

};

int PhenotypeConfigurations::GetNumPhenotypes() {
	return this->_mpPhenotypes.size();
}

void PhenotypeConfigurations::InitKeys( map<string, double> * pMap) {
	for (std::map<string,string>::iterator it = _mpPhenotypes.begin(); it != _mpPhenotypes.end(); ++it) {
		(*pMap)[it->first] = -1; //initialize values
	}
}

void PhenotypeConfigurations::GetKeys(  vector< string > &vKeys) {
	for (std::map<string,string>::iterator it = _mpPhenotypes.begin(); it != _mpPhenotypes.end(); ++it) {
		
		vKeys.push_back(it->first);
	}
}

vector< pair<int, double> > * PhenotypeConfigurations::GetFormulaSymbols(string sPhenotypeName) {

	return &this->_mpPhenotypeFormulaSymbols[sPhenotypeName];
	
};

vector< string > * PhenotypeConfigurations::GetFormulaParentSymbols(string sPhenotypeName, int nParent) { //nparent = 0 mom, 1 dad
	return (nParent == 0)? &this->_mpPhenotypeFormulaMomSymbols[sPhenotypeName] : &this->_mpPhenotypeFormulaDadSymbols[sPhenotypeName];
};

vector<string> * PhenotypeConfigurations::GetFormulaSymbolStrings(string sPhenotypeName) {

	return &this->_mpPhenotypeFormulaSymbolStrings[sPhenotypeName];
	
};

Parser * PhenotypeConfigurations::GetFormula(string sPhenotypeName) {
	#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
		//printf("%s : %d\n", "SexSelFormulae # ", nCPU);
	#else
		int nCPU = 0;
	#endif

	return this->_mpmpPhenotypeFormulae[nCPU][sPhenotypeName];
};

bool PhenotypeConfigurations::IsReferToParentPhe(string sPhenotypeName) {
	return this->_mpIsReferToParentPhe[sPhenotypeName];
}

NaturalSelectionConfigurations::NaturalSelectionConfigurations(void * pParentConfig) {

	this->_pParentConfig = pParentConfig;

};

void NaturalSelectionConfigurations::LoadFromFile(string szConfigFile)
{
	//FILE *pConfigFile;
	ifstream fsConfigFile;
	string sBuffer;
	char szPopName[MAXLINE];
	char szFormula[MAXLINE];
	
	#ifdef _OPENMP
		int nTotalCPUCore =  omp_get_max_threads();//omp_get_num_threads();
		//printf("OpenMP enabled\n");
	#else
		int nTotalCPUCore = 1;
		//printf("OpenMP disabled \n");
	#endif

	
	
	printf("Start loading natural selection file %s ...\n", szConfigFile.c_str());

	_szConfigFilename = szConfigFile; // save name of config file

	fsConfigFile.open(szConfigFile.c_str() );
	//pConfigFile = fopen(szConfigFile.c_str(), "r");
	

	if (!fsConfigFile.is_open()) {
		printf("%s\n", "Cannot open natural selection file ...");
		throw "Cannot open natural selection file!";
	}

	getline(fsConfigFile,sBuffer);
	if (sBuffer == "") { // skip the header line
		throw "Cannot read natural selection file!";
	}

	string sPopName="";
	string sFormula="";
	string sPrevPopName="";
	string sPrevFormula="";

	while(fsConfigFile.good()) 
	{
		getline(fsConfigFile,sBuffer);
		
		int nGen;
		sBuffer.erase(sBuffer.find_last_not_of(" \n\r\t")+1);//trim
                if (sBuffer == "") {
                        continue;
                }

		sscanf(sBuffer.c_str(), "%[^\t\n]	%d	%[^\t\n]", szPopName, &nGen, szFormula);
		sPrevPopName = sPopName;
		sPrevFormula = sFormula;
		sPopName = szPopName; //make it string so that it saves lots of trouble of fuddling with c strings... i hate c strings.
		sFormula = szFormula;
		if (sPopName==sPrevPopName && sFormula == sPrevFormula) {
			printf("Warning: %s already loaded for pop %s, skipping...\n", sFormula.c_str(), sPopName.c_str());
			continue;
		}

		printf("Parsing natural rule for Pop %s : %s\n", sPopName.c_str(), sFormula.c_str());
		if (sPopName.length() != 0 && sFormula.length() !=0) {
			
			

			if (this->_mpRules.find(sPopName) == this->_mpRules.end()) { // if this pop is not in the map yet
				list< string > vsFormulae;
				list< string > vsFreqDependentFormulae;
				list< pair< Parser * , int> > vpFormulae;
				list< pair< Parser * , int> > vpFreqDependentFormulae;
				list< vector<string> > vvsSymbols;
				list< vector<string> > vvsSymbols2;
				list< vector<string> > vvsSymbols3;
				list< vector<string> > vvsSymbols4;
				list< vector<string> > vvsSymbols5;
				list< vector<string> > vvsSymbols6;
				
				this->_mpRules[sPopName] = vsFormulae;
				this->_mpFreqDependentRules[sPopName] = vsFreqDependentFormulae;
				for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) { // the formula object needs to be duplicated such that it will parallelize
					printf("%d CPU cores found, putting natural selection parser on CPU %d\n",nTotalCPUCore, nCPU );
					if (this->_mpmpRuleFormulae.find(nCPU) == this->_mpmpRuleFormulae.end()) {
						 map<string , list< pair<Parser *, int> > > mpRuleFormulae, _mpFreqDependentRuleFormulae;
						this->_mpmpRuleFormulae[nCPU] =  mpRuleFormulae;
						this->_mpmpFreqDependentRuleFormulae[nCPU] = _mpFreqDependentRuleFormulae;
					}

					this->_mpmpRuleFormulae[nCPU][sPopName] = vpFormulae;
					this->_mpmpFreqDependentRuleFormulae[nCPU][sPopName] = vpFreqDependentFormulae;
				}
				this->_mpRuleFormulaSymbolStrings[sPopName] = vvsSymbols;
				this->_mpRuleFormulaSymbolStringsCourter[sPopName] = vvsSymbols2;
				this->_mpRuleFormulaSymbolStringsSelf[sPopName] = vvsSymbols3;
				this->_mpRuleFormulaSymbolStringsPopWide[sPopName] = vvsSymbols4;
				this->_mpRuleFormulaSymbolStringsPopWideChooser[sPopName] = vvsSymbols5;
				this->_mpRuleFormulaSymbolStringsPopWideCourter[sPopName] = vvsSymbols6;

			}

			
			
			Parser * pF = new Parser("");
			//Parser * pDummyFormula = new Parser("1"); // dummy formula that always return 1.

			try {
				pF->Evaluate(sFormula);
			}
			catch(std::runtime_error e) {
				
				if ( string(e.what()) == "Divide by zero" ) {
					
					printf("Warning in formula: %s \nErrMsg:%s\nThis could simply be caused by uninitialized values", sFormula.c_str(), e.what());
					
					
				}
				else {
					printf("Error occured when trying to parse formula: %s \nErrMsg:%s\nPlease check formulae grammar!", sFormula.c_str(), e.what());
					exit(100);
				}
			}

		
			vector<string> vSymbolStrings; // all symbols
			vector<string> vSymbolStringsCourter; // only symbols referring to an individual courter's trait prefixed with Courter_
			vector<string> vSymbolStringsPopWideCourter; // symbols referring to population-level variables prefixed with PopCourter_
			vector<string> vSymbolStringsSelf; //symbols referring to individual chooser trait prefixed with My_
			vector<string> vSymbolStringsPopWideChooser;//symbols referring to population-level variables prefixed with PopChooser_
			vector<string> vSymbolStringsPopWide; // symbols referring to population-level variables prefixed with Pop_
			//list all the symbols used in this formula

			bool bFreqDependent = false; // does this formula contain frequency dependent symbols?

			for(map<string,double>::iterator it = pF->symbols_.begin(); it != pF->symbols_.end(); ++it) 
			{	
				//Parse chromosome number and position from the symbol
				//cout << (it->first) << "\n"; //list all keys

				
				if (strlen(it->first.c_str()) < 3) { // no need to add symbol
					continue;
				}

				vSymbolStrings.push_back(it->first); // save label too

				//Look at prefixes:
				string sSymbol = it->first;
				
				if (sSymbol.find("Courter_") == 0) {
					throw(new Exception("Courter_ not yet implemented in natural selection."));
					sSymbol.replace(0, 8, "");
					vSymbolStringsCourter.push_back(sSymbol);
				}
				else if (sSymbol.find("PopCourter_") == 0) {
					//throw(new Exception("PopCourter_ not yet implemented in natural selection."));
					sSymbol.replace(0, 11, "");
					vSymbolStringsPopWideCourter.push_back(sSymbol );
					bFreqDependent = true;
				}
				else if (sSymbol.find("My_") == 0) {
					sSymbol.replace(0, 3, "");
					vSymbolStringsSelf.push_back(sSymbol);
				}
				else if (sSymbol.find("PopChooser_") == 0) {
					//throw(new Exception("PopChooser_ not yet implemented in natural selection."));
					sSymbol.replace(0, 11, "");
					vSymbolStringsPopWideChooser.push_back(sSymbol);
					bFreqDependent = true;
				}
				else if (sSymbol.find("Pop_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 4, "");
					vSymbolStringsPopWide.push_back(sSymbol);
					bFreqDependent = true;
				}
				else { //default to be self
					vSymbolStringsSelf.push_back(sSymbol);
				}
			}

			delete pF;
			
			if (bFreqDependent) {
				this->_mpFreqDependentRules[sPopName].push_back(sFormula); // save the formula string for further reference.
				
				
				
				this->_mpRules[sPopName].push_back("1"); // 
				

				for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) { // the formula object needs to be duplicated such that it will parallelize

					this->_mpmpFreqDependentRuleFormulae[nCPU][sPopName].push_back(pair< Parser *, int >( new Parser(sFormula) , nGen )); // literally duplicate the parser object so that it can parallelize
					this->_mpmpRuleFormulae[nCPU][sPopName].push_back(pair< Parser *, int >( new Parser("1") , nGen )); // set dummy formula which always return 1. 
				}
				
			}
			else {

				this->_mpRules[sPopName].push_back(sFormula); // save the formula string for further reference.			

				this->_mpFreqDependentRules[sPopName].push_back("1"); 
				

				for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) { // the formula object needs to be duplicated such that it will parallelize
					this->_mpmpRuleFormulae[nCPU][sPopName].push_back(pair< Parser *, int >( new Parser(sFormula) , nGen ));
					this->_mpmpFreqDependentRuleFormulae[nCPU][sPopName].push_back(pair< Parser *, int >( new Parser("1") , nGen ));// set dummy formula which always return 1.
					
				}

			}
			this->_mpRuleFormulaSymbolStrings[sPopName].push_back(vSymbolStrings);
			this->_mpRuleFormulaSymbolStringsCourter[sPopName].push_back(vSymbolStringsCourter);
			this->_mpRuleFormulaSymbolStringsSelf[sPopName].push_back(vSymbolStringsSelf);
			this->_mpRuleFormulaSymbolStringsPopWide[sPopName].push_back(vSymbolStringsPopWide);
			this->_mpRuleFormulaSymbolStringsPopWideChooser[sPopName].push_back(vSymbolStringsPopWideChooser);
			this->_mpRuleFormulaSymbolStringsPopWideCourter[sPopName].push_back(vSymbolStringsPopWideCourter);

			if (nGen > -1) {
				this->_vSpecialGens.insert(nGen);
			}
		}
	}

	printf("Selection rules for %d populations loaded. \n", this->_mpRules.size());

}

list< vector<string> > * NaturalSelectionConfigurations::GetFormulaSymbolStrings(string sPop) {
	return &this->_mpRuleFormulaSymbolStrings[sPop];
}

list< vector<string> > * NaturalSelectionConfigurations::GetFormulaSymbolStringsCourter(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsCourter[sPop];
}

list< vector<string> > * NaturalSelectionConfigurations::GetFormulaSymbolStringsSelf(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsSelf[sPop];
}

list< vector<string> > * NaturalSelectionConfigurations::GetFormulaSymbolStringsPop(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsPopWide[sPop];
}

list< vector<string> > * NaturalSelectionConfigurations::GetFormulaSymbolStringsPopCourter(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsPopWideCourter[sPop];
}

list< vector<string> > * NaturalSelectionConfigurations::GetFormulaSymbolStringsPopChooser(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsPopWideChooser[sPop];
}

bool NaturalSelectionConfigurations::IgnoreGlobalRules(int nGen) {
	return this->_vSpecialGens.find(nGen) != this->_vSpecialGens.end();
}

list< pair< Parser * , int > > * NaturalSelectionConfigurations::GetFormulae(string sPop) {
	#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
	#else
		int nCPU = 0;
	#endif

	return &this->_mpmpRuleFormulae[nCPU][sPop];
}

list< pair< Parser * , int > > * NaturalSelectionConfigurations::GetFreqDependentFormulae(string sPop) {
	#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
	#else
		int nCPU = 0;
	#endif
	return &this->_mpmpFreqDependentRuleFormulae[nCPU][sPop];
}

map< int, map<string , list< pair<Parser *, int> > > > *  NaturalSelectionConfigurations::GetFormulaeAllCPUs() {
	return &this->_mpmpRuleFormulae;
}

map< int, map<string , list< pair<Parser *, int> > > > *  NaturalSelectionConfigurations::GetFreqDependentFormulaeAllCPUs() {
	return &this->_mpmpFreqDependentRuleFormulae;
}


SexualSelectionConfigurations::SexualSelectionConfigurations(void * pParentConfig) {

	this->_pParentConfig = pParentConfig;

};

void SexualSelectionConfigurations::LoadFromFile(string szConfigFile)
{
	ifstream fsConfigFile;
	string sBuffer;
	char szPopName[MAXLINE];
	char szFormula[MAXLINE];
	
	#ifdef _OPENMP
		int nTotalCPUCore = omp_get_max_threads();//omp_get_num_threads();
	#else
		int nTotalCPUCore = 1;
	#endif
	
	printf("Start loading sexual selection file %s ...\n", szConfigFile.c_str());

	_szConfigFilename = szConfigFile; // save name of config file

	fsConfigFile.open(szConfigFile.c_str());
	//pConfigFile = fopen(szConfigFile.c_str(), "r");
	

	if (!fsConfigFile.is_open()) {
		printf("%s\n", "Cannot open sexual selection file ...");
		throw "Cannot open sexual selection file!";
	}

	getline(fsConfigFile,sBuffer);
	if (sBuffer == "") { // skip the header line
		throw "Cannot read sexual selection file!";
	}


	while(fsConfigFile.good()) 
	{
		getline(fsConfigFile,sBuffer);
		string sPopName, sFormula;
		int nGen = -1;
		sBuffer.erase(sBuffer.find_last_not_of(" \n\r\t")+1);//trim
		if (sBuffer == "") {
			continue;
		}
		sscanf(sBuffer.c_str(), "%[^\t\n]	%d	%[^\t\n]", szPopName, &nGen, szFormula);
		sPopName = szPopName; //make it string so that it saves lots of trouble of fuddling with c strings... i hate c strings.
		sFormula = szFormula;
		printf("Parsing sexual rule for Pop %s : %s for gen %d\n", sPopName.c_str(), sFormula.c_str(), nGen);
		if (sPopName.length() != 0 && sFormula.length() !=0) {
			
			

			if (this->_mpRules.find(sPopName) == this->_mpRules.end()) { // if this pop is not in the map yet
				list< string > vsFormulae;
				list< pair< Parser * , int > > vpFormulae;
				list< vector<string> > vvsSymbols;
				list< vector<string> > vvsSymbols2;
				list< vector<string> > vvsSymbols3;
				list< vector<string> > vvsSymbols4;
				list< vector<string> > vvsSymbols5;
				list< vector<string> > vvsSymbols6;
				list< vector<string> > vvsSymbols7;
				list< vector<string> > vvsSymbols8;
				list< vector<string> > vvsSymbols9;
				list< vector<string> > vvsSymbols10;
				list< vector<string> > vvsSymbols11;
				list< vector<string> > vvsSymbols12;
				list< vector<string> > vvsSymbols13;
				list< vector<string> > vvsSymbols14;

				this->_mpRules[sPopName] = vsFormulae;
				for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) {
					printf("%d CPU cores found, putting sexual selection parser on CPU %d\n",nTotalCPUCore, nCPU );
					if (this->_mpmpRuleFormulae.find(nCPU) == this->_mpmpRuleFormulae.end()) {
						 map<string , list< pair<Parser *, int> > > mpRuleFormulae;
						this->_mpmpRuleFormulae[nCPU] =  mpRuleFormulae;
					}
					this->_mpmpRuleFormulae[nCPU][sPopName] = vpFormulae;
				}
				this->_mpRuleFormulaSymbolStrings[sPopName] = vvsSymbols;
				this->_mpRuleFormulaSymbolStringsCourter[sPopName] = vvsSymbols2;
				this->_mpRuleFormulaSymbolStringsChooser[sPopName] = vvsSymbols3;
				this->_mpRuleFormulaSymbolStringsPopWide[sPopName] = vvsSymbols4;
				this->_mpRuleFormulaSymbolStringsPopWideChooser[sPopName] = vvsSymbols5;
				this->_mpRuleFormulaSymbolStringsPopWideCourter[sPopName] = vvsSymbols6;
				this->_mpRuleFormulaSymbolStringsDad[sPopName] = vvsSymbols7;
				this->_mpRuleFormulaSymbolStringsMom[sPopName] = vvsSymbols8;
				
				this->_mpRuleFormulaSymbolStringsCurrPopWidePrevGen[sPopName] = vvsSymbols9; // Symbols that are population-wide values from prev gen. of current pop
				this->_mpRuleFormulaSymbolStringsCurrPopWidePrevGenCourter[sPopName] = vvsSymbols10; // Symbols that are population-wide courter values from prev gen. of current pop
				this->_mpRuleFormulaSymbolStringsCurrPopWidePrevGenChooser[sPopName] = vvsSymbols11; // Symbols that are population-wide chooser values from prev gen. of current pop
				this->_mpRuleFormulaSymbolStringsPrevPopWidePrevGen[sPopName] = vvsSymbols12; // Symbols that are population-wide values from prev gen. of previous pop, if the individual did not migrate, previous pop = current pop
				this->_mpRuleFormulaSymbolStringsPrevPopWidePrevGenCourter[sPopName] = vvsSymbols13; // Symbols that are population-wide courter values from prev gen. of previous pop, if the individual did not migrate, previous pop = current pop
				this->_mpRuleFormulaSymbolStringsPrevPopWidePrevGenChooser[sPopName] = vvsSymbols14; // Symbols that are population-wide chooser values from prev gen. of previous pop, if the individual did not migrate, previous pop = current pop

			
			}

			
			
			Parser * pF = new Parser("");
			try {
				pF->Evaluate(sFormula);
			}
			catch(std::runtime_error e) {
				if ( string(e.what()) == "Divide by zero" ) {
					
					printf("Warning in formula: %s \nErrMsg:%s\nThis could simply be caused by uninitialized values", sFormula.c_str(), e.what());
					
					
				}
				else {
					printf("Error occured when trying to parse formula: %s \nErrMsg:%s\nPlease check formulae grammar!", sFormula.c_str(), e.what());
					exit(100);
				}
			}

			vector<string> vSymbolStrings; // all symbols
			vector<string> vSymbolStringsCourter; // only symbols referring to an individual courter's trait prefixed with Courter_
			vector<string> vSymbolStringsPopWideCourter; // symbols referring to population-level variables prefixed with PopCourter_
			vector<string> vSymbolStringsChooser; //symbols referring to individual chooser trait prefixed with My_
			vector<string> vSymbolStringsPopWideChooser;//symbols referring to population-level variables prefixed with PopChooser_
			vector<string> vSymbolStringsPopWide; // symbols referring to population-level variables prefixed with Pop_
			vector<string> vSymbolStringsDad; // symbols referring to population-level variables prefixed with Dad_
			vector<string> vSymbolStringsMom; // symbols referring to population-level variables prefixed with Mom_

			vector<string> vSymbolStringsCurrPopWidePrevGen; //PrevGenCurrPop_
			vector<string> vSymbolStringsCurrPopWidePrevGenCourter; //PrevGenCurrPopCourter_
			vector<string> vSymbolStringsCurrPopWidePrevGenChooser; //PrevGenCurrPopChooser_
			vector<string> vSymbolStringsPrevPopWidePrevGen; //PrevGenPrevPop_
			vector<string> vSymbolStringsPrevPopWidePrevGenCourter;  //PrevGenPrevPopCourter_
			vector<string> vSymbolStringsPrevPopWidePrevGenChooser; //PrevGenPrevPopChooser_

			//list all the symbols used in this formula
			for(map<string,double>::iterator it = pF->symbols_.begin(); it != pF->symbols_.end(); ++it) 
			{	
				//Parse chromosome number and position from the symbol
				//cout << (it->first) << "\n"; //list all keys

				
				if (strlen(it->first.c_str()) < 3) { // no need to add symbol
					continue;
				}

				vSymbolStrings.push_back(it->first); // save label too

				//Look at prefixes:
				string sSymbol = it->first;
				
				if (sSymbol.find("Courter_") == 0) {
					//throw(new Exception("Courter_ not yet implemented in natural selection."));
					sSymbol.replace(0, 8, "");
					vSymbolStringsCourter.push_back(sSymbol);
				}
				else if (sSymbol.find("PopCourter_") == 0) {
					//throw(new Exception("PopCourter_ not yet implemented in natural selection."));
					sSymbol.replace(0, 11, "");
					vSymbolStringsPopWideCourter.push_back(sSymbol );
				}
				else if (sSymbol.find("My_") == 0) {
					sSymbol.replace(0, 3, "");
					vSymbolStringsChooser.push_back(sSymbol);
				}
				else if (sSymbol.find("PopChooser_") == 0) {
					//throw(new Exception("PopChooser_ not yet implemented in natural selection."));
					sSymbol.replace(0, 11, "");
					vSymbolStringsPopWideChooser.push_back(sSymbol);
				}
				else if (sSymbol.find("Pop_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 4, "");
					vSymbolStringsPopWide.push_back(sSymbol);
				}
				else if (sSymbol.find("Dad_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 4, "");
					vSymbolStringsDad.push_back(sSymbol);
				}
				else if (sSymbol.find("Mom_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 4, "");
					vSymbolStringsMom.push_back(sSymbol);
				}
				else if (sSymbol.find("PrevGenCurrPop_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 15, "");
					vSymbolStringsCurrPopWidePrevGen.push_back(sSymbol);
				}
				else if (sSymbol.find("PrevGenCurrPopCourter_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 22, "");
					vSymbolStringsCurrPopWidePrevGenCourter.push_back(sSymbol);
				}
				else if (sSymbol.find("PrevGenCurrPopChooser_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 22, "");
					vSymbolStringsCurrPopWidePrevGenChooser.push_back(sSymbol);
				}
				else if (sSymbol.find("PrevGenPrevPop_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 15, "");
					vSymbolStringsPrevPopWidePrevGen.push_back(sSymbol);
				}
				else if (sSymbol.find("PrevGenPrevPopCourter_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 22, "");
					vSymbolStringsPrevPopWidePrevGenCourter.push_back(sSymbol);
				}
				else if (sSymbol.find("PrevGenPrevPopChooser_") == 0) {
					//throw(new Exception("Pop_ not yet implemented in natural selection."));
					sSymbol.replace(0, 22, "");
					vSymbolStringsPrevPopWidePrevGenChooser.push_back(sSymbol);
				}
				else { //default to be self
					vSymbolStringsChooser.push_back(sSymbol);
				}

			}
			
			delete pF;

			this->_mpRules[sPopName].push_back(sFormula); // save the formula string for further reference.
			for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) {
				if (nCPU == 0) {
					printf("Inserting %s to gen %d\n", sFormula.c_str(), nGen);
				}
				this->_mpmpRuleFormulae[nCPU][sPopName].push_back( pair< Parser *, int> ( new Parser(sFormula), nGen) ); //literally copy the pF object, so that there are nTotalCPUCore numbers of copies in the RAM! ;-)
			}
			this->_mpRuleFormulaSymbolStrings[sPopName].push_back(vSymbolStrings);
			this->_mpRuleFormulaSymbolStringsCourter[sPopName].push_back(vSymbolStringsCourter);
			this->_mpRuleFormulaSymbolStringsChooser[sPopName].push_back(vSymbolStringsChooser);
			this->_mpRuleFormulaSymbolStringsPopWide[sPopName].push_back(vSymbolStringsPopWide);
			this->_mpRuleFormulaSymbolStringsPopWideChooser[sPopName].push_back(vSymbolStringsPopWideChooser);
			this->_mpRuleFormulaSymbolStringsPopWideCourter[sPopName].push_back(vSymbolStringsPopWideCourter);
			this->_mpRuleFormulaSymbolStringsDad[sPopName].push_back(vSymbolStringsDad);
			this->_mpRuleFormulaSymbolStringsMom[sPopName].push_back(vSymbolStringsMom);

			this->_mpRuleFormulaSymbolStringsCurrPopWidePrevGen[sPopName].push_back(vSymbolStringsCurrPopWidePrevGen); //PrevGenCurrPop_
			this->_mpRuleFormulaSymbolStringsCurrPopWidePrevGenCourter[sPopName].push_back(vSymbolStringsCurrPopWidePrevGenCourter); //PrevGenCurrPopCourter_
			this->_mpRuleFormulaSymbolStringsCurrPopWidePrevGenChooser[sPopName].push_back(vSymbolStringsCurrPopWidePrevGenChooser); //PrevGenCurrPopChooser_
			this->_mpRuleFormulaSymbolStringsPrevPopWidePrevGen[sPopName].push_back(vSymbolStringsPrevPopWidePrevGen); //PrevGenPrevPop_
			this->_mpRuleFormulaSymbolStringsPrevPopWidePrevGenCourter[sPopName].push_back(vSymbolStringsPrevPopWidePrevGenCourter);  //PrevGenPrevPopCourter_
			this->_mpRuleFormulaSymbolStringsPrevPopWidePrevGenChooser[sPopName].push_back(vSymbolStringsPrevPopWidePrevGenChooser); //PrevGenPrevPopChooser_

			if (nGen > -1) {
				printf("Sexual selection special generation: %d\n", nGen);
				this->_vSpecialGens.insert(nGen);
			}
		}
	}

	printf("Sexual Selection rules for %d populations loaded. \n", this->_mpRules.size());

}

bool SexualSelectionConfigurations::IgnoreGlobalRules(int nGen) {
	return this->_vSpecialGens.find(nGen) != this->_vSpecialGens.end();
}

list< vector<string> > * SexualSelectionConfigurations::GetFormulaSymbolStrings(string sPop) {
	return &this->_mpRuleFormulaSymbolStrings[sPop];
}

list< pair< Parser * , int> > * SexualSelectionConfigurations::GetFormulae(string sPop) {
	#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
		//printf("%s : %d\n", "SexSelFormulae # ", nCPU);
	#else
		int nCPU = 0;
	#endif
	return &this->_mpmpRuleFormulae[nCPU][sPop];
}

map< int, map<string , list< pair<Parser *, int> > > > *  SexualSelectionConfigurations::GetFormulaeAllCPUs() {
	return &this->_mpmpRuleFormulae;
}


list< vector<string> > * SexualSelectionConfigurations::GetFormulaSymbolStringsCourter(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsCourter[sPop];
}

list< vector<string> > * SexualSelectionConfigurations::GetFormulaSymbolStringsSelf(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsChooser[sPop];
}

list< vector<string> > * SexualSelectionConfigurations::GetFormulaSymbolStringsPop(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsPopWide[sPop];
}

list< vector<string> > * SexualSelectionConfigurations::GetFormulaSymbolStringsPopCourter(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsPopWideCourter[sPop];
}

list< vector<string> > * SexualSelectionConfigurations::GetFormulaSymbolStringsPopChooser(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsPopWideChooser[sPop];
}

list< vector<string> > * SexualSelectionConfigurations::GetFormulaSymbolStringsDad(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsDad[sPop];
}

list< vector<string> > * SexualSelectionConfigurations::GetFormulaSymbolStringsMom(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsMom[sPop];
}

list< vector<string> > *  SexualSelectionConfigurations::GetFormulaSymbolStringsPrevGenCurrPop(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsCurrPopWidePrevGen[sPop];
}
list< vector<string> > *  SexualSelectionConfigurations::GetFormulaSymbolStringsPrevGenCurrPopCourter(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsCurrPopWidePrevGenCourter[sPop];
}

list< vector<string> > *  SexualSelectionConfigurations::GetFormulaSymbolStringsPrevGenCurrPopChooser(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsCurrPopWidePrevGenChooser[sPop];	
}
list< vector<string> > *  SexualSelectionConfigurations::GetFormulaSymbolStringsPrevGenPrevPop(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsPrevPopWidePrevGen[sPop];
}

list< vector<string> > *  SexualSelectionConfigurations::GetFormulaSymbolStringsPrevGenPrevPopCourter(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsPrevPopWidePrevGenCourter[sPop];
}

list< vector<string> > *  SexualSelectionConfigurations::GetFormulaSymbolStringsPrevGenPrevPopChooser(string sPop) {
	return &this->_mpRuleFormulaSymbolStringsPrevPopWidePrevGenChooser[sPop];
}


GeneConfigurations::GeneConfigurations(void * pParentConfig) {

	this->_pParentConfig = pParentConfig;

}

void GeneConfigurations::LoadFromFile(string szConfigFile) {
	#ifdef _OPENMP
		int nTotalCPUCore =  omp_get_max_threads();//omp_get_num_threads();
		//printf("OpenMP enabled\n");
	#else
		int nTotalCPUCore = 1;
		//printf("OpenMP disabled \n");
	#endif

	this->_nHaploidChrNum = ((SimulationConfigurations *)this->_pParentConfig)->pMarkerConfig->GetHaploidChromosomeNum();

	//Initialize _mpGenes
	//this->_mpGenes = new map<double, GeneProperties>[_nHaploidChrNum];
	//this->_mpGeneIndex = new map<double, int>[_nHaploidChrNum];
	vector<int> vIndexCounters;
	vector< map<double , GeneProperties> > mpGenes;
	vector< map<double, int> > mpGeneIndex;

	for (int i=0;i<_nHaploidChrNum;i++) {
		map<double , GeneProperties> oPlaceHolderMap1;
		map<double, int> oPlaceHolderMap2;

		mpGenes.push_back(oPlaceHolderMap1);
		mpGeneIndex.push_back(oPlaceHolderMap2);

		vIndexCounters.push_back(0);
	}

	FILE *pConfigFile;
	char szBuffer[12048];
	char szGeneName[100];
	char szMode[20];
	int nChr;
	double  nPos, nDomVal, nRecesVal , nDominantFreqPop1, nDominantFreqPop2;
	char nDominantName, nRecessiveName;

	
	
	printf("Start loading gene file %s ...\n", szConfigFile.c_str());

	_szConfigFilename = szConfigFile; // save name of config file

	pConfigFile = fopen(szConfigFile.c_str(), "r");

	if (!pConfigFile) {
		printf("%s\n", "Cannot open phenotype file ...");
		throw "Cannot open gene file!";
	}

	if (fgets(szBuffer, 2048, pConfigFile) == NULL) { // skip the header line
		throw "Cannot read gene file!";
	}

	

	while(fgets(szBuffer, 12048, pConfigFile) != NULL) 
	{
		
		if (sscanf(szBuffer, "%[^\t\n]	%d	%lf	%c	%lf	%c	%lf	%s %lf %lf", szGeneName, &nChr ,&nPos, &nDominantName, &nDomVal, &nRecessiveName, &nRecesVal, szMode, &nDominantFreqPop1, &nDominantFreqPop2  ) != 10) {
			continue;//empty line
			
		};
		
		if (nChr > this->_nHaploidChrNum || nChr <=0) {
			throw "Chromosome number incorrect in the gene file!";
		}

		GeneProperties oGeneProp;
		oGeneProp.AlleleMode = strcmp(szMode,"Dominant")==0? GeneProperties::Dominant: strcmp(szMode,"Additive")==0? GeneProperties::Additive : GeneProperties::Hemizygous;
		oGeneProp.GeneName = szGeneName;
		oGeneProp.DominantLabel = nDominantName;
		oGeneProp.DominantValue = nDomVal;
		oGeneProp.RecessiveLabel = nRecessiveName;
		oGeneProp.RecessiveValue = nRecesVal;
		oGeneProp.DominantFreqPop1 = nDominantFreqPop1;
		oGeneProp.DominantFreqPop2 = nDominantFreqPop2;

		mpGenes[nChr-1][nPos] = oGeneProp;
		mpGeneIndex[nChr-1][nPos] = vIndexCounters.at(nChr-1);
		vIndexCounters.at(nChr-1) += 1;
		//nIndexCounter++;
		
	}

	
	//Try opening muation file too
	string szMutationFile = ((SimulationConfigurations *)this->_pParentConfig)->GetConfig("GeneMutationFile");
	if (szMutationFile=="") {
		std::printf("No mutation file specified for genes, no mutations will be carried out\n");
	}
	else {
		FILE *pMutationFile;
		char szBuffer[12048];
		char szFormula[12048];
		char szGeneName[100];
		char szPopName[1000];
		int nChr, nStartGen, nEndGen;
		double  nPos, nMuteRate, nLowerBound , nUpperBound;


		pMutationFile = fopen(szMutationFile.c_str(), "r");

		if (!pMutationFile) {
			printf("%s\n", "Cannot open gene mutation file ...");
			throw "Cannot open mutation file!";
		}

		if (fgets(szBuffer, 12048, pMutationFile) == NULL) { // skip the header line
			throw "Cannot read mutation file!";
		}

		while(fgets(szBuffer, 12048, pMutationFile) != NULL) 
		{
			bool bLongLineFormat = true;
			if (sscanf(szBuffer, "%[^\t\n]	%d	%lf	%[^\t\n]	%d	%d	%lf	%lf	%lf %s", szGeneName, &nChr ,&nPos, szPopName, &nStartGen, &nEndGen, &nMuteRate, &nLowerBound, &nUpperBound, szFormula  ) != 10) {
				
				//see if line in short format
				if (sscanf(szBuffer, "%[^\t\n]	%d	%lf	%lf	%lf	%lf %s", szGeneName, &nChr ,&nPos, &nMuteRate, &nLowerBound, &nUpperBound, szFormula  ) != 7) {
					continue;//empty line
				}
				else {
					bLongLineFormat = false;
				}
			}
			

			if (mpGenes.size() < nChr)
			{
				std::printf("Warning: Gene mutation definition applied to an undefined chromosome: %d", nChr);
			}
			else if (mpGenes[nChr-1].find(nPos) == mpGenes[nChr-1].end()) {//
				std::printf("Warning: Gene mutation definition applied to an undefined location: %f on chromosome %d", nPos, nChr);
			}
			else {
				mpGenes[nChr-1][nPos].MutationProb = nMuteRate;
				mpGenes[nChr-1][nPos].LowerBound = nLowerBound;
				mpGenes[nChr-1][nPos].UpperBound = nUpperBound;
				mpGenes[nChr-1][nPos].pFormula = new Parser(szFormula);
				mpGenes[nChr-1][nPos].sFormula = szFormula;
				string sPopName;
				if (bLongLineFormat) {
					mpGenes[nChr-1][nPos].Pops.insert( sPopName=szPopName );
					mpGenes[nChr-1][nPos].StartGen = nStartGen;
					mpGenes[nChr-1][nPos].EndGen = nEndGen;
				}
				else {
					//SHORT FORMAT DEFAULT TO ANY POP, ANY GEN
					mpGenes[nChr-1][nPos].Pops.insert( sPopName = "!ANYPOP" );
					mpGenes[nChr-1][nPos].StartGen = -1;
					mpGenes[nChr-1][nPos].EndGen = -1;
				}
				std::printf("Mutation rule %s (prob=%f) loaded for chr %d : %f, at Generations %d to %d in population %s\n", szFormula, nMuteRate, nChr, nPos, mpGenes[nChr-1][nPos].StartGen, mpGenes[nChr-1][nPos].EndGen, sPopName.c_str() );
			}

			
		};
		
	}

	for( int nCPU=0; nCPU<nTotalCPUCore ; nCPU++) {
		this->_mpmpGenes[nCPU] = mpGenes;
		this->_mpmpGeneIndex[nCPU] = mpGeneIndex;

	}

	//printf("%f \n", this->_mpGenes[23][100.0].DominantFreqPop1);

}

GeneProperties& GeneProperties::operator=(const GeneProperties& oSource) { // overload the copy behavior
	this->AlleleMode = oSource.AlleleMode;
	this->DominantFreqPop1 = oSource.DominantFreqPop1;
	this->DominantFreqPop2 = oSource.DominantFreqPop2;
	this->DominantLabel = oSource.DominantLabel;
	this->DominantValue = oSource.DominantValue;
	this->GeneName = oSource.GeneName;
	this->LowerBound = oSource.LowerBound;
	this->MutationProb = oSource.MutationProb;
	this->RecessiveLabel = oSource.RecessiveLabel;
	this->RecessiveValue = oSource.RecessiveValue;
	this->sFormula = oSource.sFormula;
	this->UpperBound = oSource.UpperBound;
	this->pFormula = new Parser(oSource.sFormula);
	this->StartGen = oSource.StartGen;
	this->EndGen = oSource.EndGen;
	for(set< string >::iterator it=oSource.Pops.begin(); it != oSource.Pops.end(); ++it) {
		this->Pops.insert(*it);
	}
	 return *this ;
}

vector< map<double , GeneProperties> > * GeneConfigurations::GetMpGenes() {
	#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
	#else
		int nCPU = 0;
	#endif

	return &this->_mpmpGenes[nCPU];
}

vector< map<double, int> > *  GeneConfigurations::GetMpGeneIndex() {
	#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
	#else
		int nCPU = 0;
	#endif

	return &this->_mpmpGeneIndex[nCPU];
}

SimulationConfigurations::SimulationConfigurations() {
	this->pMarkerConfig = new MarkerConfigurations(this);
	this->pRecombProbConfig = new RecombProbConfigurations(this);
	this->pPhenotypeConfig = new PhenotypeConfigurations(this);
	this->pGeneConfig = new GeneConfigurations(this);
	this->pNaturalSelConfig = new NaturalSelectionConfigurations(this);
	this->pSexualSelConfig = new SexualSelectionConfigurations(this);
	_szConfigFilename="";
};

int SimulationConfigurations::GetLogVerboseLevel() {
	return this->_nVerboseLogLevel;
}

void SimulationConfigurations::LoadFromFile(string szConfigFile) {


	FILE *pConfigFile;

	printf("%s\n", "Start loading config ...");

	_szConfigFilename = szConfigFile; // save name of config file

	pConfigFile = fopen(szConfigFile.c_str(), "r");

	if (!pConfigFile) {
		printf("%s\n", "Cannot open config file ...");
		throw "Cannot open config file!";
	}

	//Go over each line:
	char szBuffer[1024];


	while(fgets(szBuffer, 1024, pConfigFile) != NULL) 
	{
		char * szKey = new char[1024];
		char * szValue = new char[1024];
		sscanf(szBuffer, "%s = %s %*s", szKey, szValue);
		string sKey(szKey);
		string sVal(szValue);

		this->_mpConfigs[sKey] = sVal;
		printf("%s = %s \n", sKey.c_str(), sVal.c_str());
		//printf("Size of map: %d\n", this->_mpConfigs.size());
	}

	fclose(pConfigFile);

	this->fnParseNumericConfigs();
	//printf("Size of map: %d\n", this->_mpConfigs.size());
	printf("%s = %s \n", "PhenotypeFile", this->GetConfig("PhenotypeFile").c_str());
	printf("%s = %f \n", "pop1_init_size", this->GetNumericConfig("pop1_init_size"));
	//printf("Size of map: %d\n", this->_mpConfigs.size());

	/*
	for (std::map<char *,char *>::iterator it = this->_mpConfigs.begin(); it != this->_mpConfigs.end(); ++it) { // go over the configuration map with string values

		printf("%s%s", it->first, it->second);
	}
	*/
	
	#ifdef _OPENMP

		int nThreads = (int)this->GetNumericConfig("NumThreads");
		if (nThreads <= 0 ) {
			omp_set_num_threads( 2 );
		}
		else {
			omp_set_num_threads( nThreads );
		}
	#endif

	//Load other config files:
	this->pMarkerConfig->LoadFromFile(this->GetConfig("MarkerFile"));
	this->pRecombProbConfig->LoadFromFile(this->GetConfig("MarkerProbFile"));
	this->pPhenotypeConfig->LoadFromFile(this->GetConfig("PhenotypeFile"));
	this->pGeneConfig->LoadFromFile(this->GetConfig("GeneFile"));
	this->pNaturalSelConfig->LoadFromFile(this->GetConfig("NaturalSelection"));
	this->pSexualSelConfig->LoadFromFile(this->GetConfig("SexualSelection"));
	this->pMarkerConfig->CalculateMapDistances();
	this->fnParseSampleGenDef();

	if (this->GetNumericConfig("LogVerbose") == -1) {
		this->_nVerboseLogLevel = 0;//default, not verbose
	}
	else {
		this->_nVerboseLogLevel = this->GetNumericConfig("LogVerbose");
	}

	this->bOutputFirst3Gens = (this->GetConfig("DisableOutputForFirst3Gens")!="yes");

};

bool SimulationConfigurations::IsOutputFirst3Gens() {
	return this->bOutputFirst3Gens;
}

const string SimulationConfigurations::GetConfigFileName() {
	return _szConfigFilename;
}

void SimulationConfigurations::fnParseSampleGenDef() {
	string sDefFile = this->GetConfig("samplegens");
	if (sDefFile == "") {
		return;
	}

	FILE *pDefFile;

	printf("%s : %s\n", "Start loading sampling range file ...", sDefFile.c_str());

	pDefFile = fopen(sDefFile.c_str(), "r");

	if (!pDefFile) {
		printf("%s : %s\n", "Cannot open sampling range file ...", sDefFile.c_str());
		throw "Cannot open sampling range file!";
	}

	//Go over each line:
	char szBuffer[MAXLINE];


	while(fgets(szBuffer, MAXLINE, pDefFile) != NULL) 
	{
		int nStartGen = 0;
		int nEndGen = 0;
		int nReadField = sscanf(szBuffer, "%d %d", &nStartGen, &nEndGen);
		if (nReadField <= 0 ) {
			continue;
		}

		if (nReadField == 1) {
			nEndGen = nStartGen;
		}

		this->_mvSampleGens.push_back( std::pair< int, int>( nStartGen ,  nEndGen ) );

	}
}

bool SimulationConfigurations::IsInUserSpecifiedSamplingRange(int nGen) {
	if ( (this->_mvSampleGens.size()) == 0) {
		return false; // user didn't specify any generations to sample, return false by default.
	}

	for(std::vector< std::pair< int, int > >::iterator it = this->_mvSampleGens.begin(); it != this->_mvSampleGens.end(); ++it) {
		int nBegin = it->first;
		int nEnd = it->second;
		if (nGen >= nBegin && nGen <= nEnd)
		{
			return true;
		}
	}

	return false;
}


void SimulationConfigurations::fnParseNumericConfigs() {

	//printf("Size of map: %d\n", this->_mpConfigs.size());
	printf("%s\n", "Converting numeric config values ...");

	for (std::map<string,string>::iterator it = this->_mpConfigs.begin(); it != this->_mpConfigs.end(); ++it) { // go over the configuration map with string values

		string sKey = it->first;
		string sVal = it->second;

		this->_mpNumericConfigs[sKey] = strtod(sVal.c_str() , NULL);
		//printf("%s = %f \n", pKey, this->_mpNumericConfigs[pKey]);
	}

	printf("%s\n", "Finished ...");



};

const string SimulationConfigurations::GetConfig(string sKey) {

	//printf("Size of map: %d\n", this->_mpConfigs.size());

	

	if (  this->_mpConfigs.find(sKey) == this ->_mpConfigs.end()) { // the requested key is not found in the config

		//printf("The key %s was not found in the configuration file.\n", sKey.c_str());
		return "";

	}
	

	return  this->_mpConfigs[sKey];


};

void SimulationConfigurations::SetConfig(string sKey, string sValue) {
	this->_mpConfigs[sKey] = sValue;
}

const double SimulationConfigurations::GetNumericConfig(string sKey) {
	
	if (  this->_mpNumericConfigs.find(sKey) == this->_mpNumericConfigs.end()) { // the requested key is not found in the config

		//printf("The key %s was not found in the configuration file.\n", sKey.c_str());
		return -1;

	}
	

	return  this->_mpNumericConfigs[sKey];


};

void SimulationConfigurations::SetNumericConfig(string sKey, double nValue) {
	this->_mpNumericConfigs[sKey] = nValue;
}

int SimulationConfigurations::GetCurrGen() {
	return this->nCurrGen;
}

void SimulationConfigurations::SetCurrGen(int nGen) {
	this->nCurrGen = nGen;
}
