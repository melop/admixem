#include "config.h"

SimulationConfigurations SimulConfig; //Global configuration object
extern Normal NormalGen; 
extern Uniform UniformGen;

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

			sscanf(szBuffer, ":chr %d len = %lf centromere = %lf", &nCurrChr, &nChrLen , &nCentromerePos);

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
			nTotalGenomeLen += nChrLen;
			pCentromerePos[nCurrChr] = nCentromerePos;
			printf("Reading markers on chromosome %d , length=%f, centromere=%f ...\n", (nCurrChr + 1), nChrLen, nCentromerePos);

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

				//printf("Warning: repetitive definition of position  %f on chromosome %d found. second one will be ignored\n", nAbsPos, nCurrChr );
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
	vector<double> * pvFemaleMapDistances = ((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->GetBreakPointMapDistances(0);
	vector<double> * pvAvgMapDistances = ((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->GetBreakPointMapDistances(2);
	vector<double> * pvBreakpointSamplePositions = ((SimulationConfigurations *)this->_pParentConfig)->pRecombProbConfig->GetBreakPointSamplePositions();
	vAbsToMapDistanceMale.clear();
	vAbsToMapDistanceFemale.clear();
	vAbsToMapDistanceAvg.clear();


	for (int nChr=0 ; nChr < this->_nHaploidChrNum; nChr++ ) {
		vector<double> * pvBreakpointSamplePositionsOnChr = &pvBreakpointSamplePositions[nChr];
		vector<double> * pvMaleMapDistancesOnChr = &pvMaleMapDistances[nChr];
		vector<double> * pvFemaleMapDistancesOnChr = &pvFemaleMapDistances[nChr];
		vector<double> * pvAvgMapDistancesOnChr = &pvAvgMapDistances[nChr];

		map<double, int> * pmMarkerIndexOnChr = &(vMarkerIndex.at(nChr));
		map<double, double> mAbsToMapDistanceOnChrMale, mAbsToMapDistanceOnChrFemale, mAbsToMapDistanceOnChrAvg ;

		double nPosLastSample = 0;

		for (map<double,int>::iterator itMarker=pmMarkerIndexOnChr->begin(); itMarker!= pmMarkerIndexOnChr->end(); ++itMarker) {
			double nPos = (*itMarker).first;
			double nAbsSampleDistance, nMapMaleSampleDistance, nMapFemaleSampleDistance, nMapAvgSampleDistance, nMapMaleDistance, nMapFemaleDistance, nMapAvgDistance;

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

void RecombProbConfigurations::LoadFromFile(string szConfigFile) {

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
	
	
	printf("Start loading marker probability file %s ...\n", szConfigFile.c_str());

	_szConfigFilename = szConfigFile; // save name of config file

	pConfigFile = fopen(szConfigFile.c_str(), "r");

	if (!pConfigFile) {
		printf("%s\n", "Cannot open marker file ...");
		throw "Cannot open marker file!";
	}



	int nCurrChr;
	char szCmd[100];
	double nVal, nPos, nMaleAccu, nFemaleAccu, nMaleMapDistance, nFemaleMapDistance, nAvgMapDistance;

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

				nCurrChr--;
				continue; // go to the next line
			}
			else { // it's a command line

				sscanf(szBuffer, "%s = %lf", szCmd, &nVal);
				
				//printf("%s\n",szCmd);

				if (strcmp(szCmd, ":ExpectedMaleRecPerMeiosisArm1")==0) {
					this->_nExpectedMaleRecPerMeiosisArm1[nCurrChr] = nVal;
				}
				if (strcmp(szCmd,":ExpectedMaleRecPerMeiosisArm2")==0) {
					this->_nExpectedMaleRecPerMeiosisArm2[nCurrChr] = nVal;
				}

				if (strcmp(szCmd,":ExpectedFemaleRecPerMeiosisArm1")==0) {
					this->_nExpectedFemaleRecPerMeiosisArm1[nCurrChr] = nVal;
				}
				if (strcmp(szCmd,":ExpectedFemaleRecPerMeiosisArm2")==0) {
					this->_nExpectedFemaleRecPerMeiosisArm2[nCurrChr] = nVal;
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
			this->pvBreakpointSamplePositions[nCurrChr].push_back(nPos);
			this->pvMaleAccuProb[nCurrChr].push_back(nMaleAccu);
			this->pvFemaleAccuProb[nCurrChr].push_back(nFemaleAccu);
			this->pvMaleMapDistance[nCurrChr].push_back(nMaleMapDistance);
			this->pvFemaleMapDistance[nCurrChr].push_back(nFemaleMapDistance);
			this->pvAvgMapDistance[nCurrChr].push_back(nAvgMapDistance);

		}


	}

	///printf("%f %f %f", this->pvBreakpointSamplePositions[0].at(0), this->pvMaleAccuProb[0].at(0), this->pvFemaleAccuProb[0].at(0));

};

void RecombProbConfigurations::GetBreakPointsByArm(bool bSex, int nChr, int nArm,vector<double> &vRet) { // return a vector of break points for a given arm of a given chromosome for a given sex
	
	double nExpectedPoints = (true==bSex? (nArm==1? _nExpectedMaleRecPerMeiosisArm1[nChr] : _nExpectedMaleRecPerMeiosisArm2[nChr]) : (nArm==1? _nExpectedFemaleRecPerMeiosisArm1[nChr] : _nExpectedFemaleRecPerMeiosisArm2[nChr]));
	int nBreakPointsToPut = (int)nExpectedPoints; // get integral part
	nBreakPointsToPut += (UniformGen.Next() <= (nExpectedPoints - (double)nBreakPointsToPut)) ? 1:0;
	double nCentromerePos = ((SimulationConfigurations*)this->_pParentConfig)->pMarkerConfig->GetCentromerePosition(nChr);
	double nChrLen = ((SimulationConfigurations*)this->_pParentConfig)->pMarkerConfig->GetChromosomeLength(nChr);

	double nStart = (nArm==1)? 0.0 : nCentromerePos;
	double nEnd   = (nArm==1)? nCentromerePos : nChrLen;

	
	vector<double> * pvAccuProb = true==bSex? pvMaleAccuProb : pvFemaleAccuProb;
	
	size_t nCentromereIndex = std::distance( pvBreakpointSamplePositions[nChr].begin(), lower_bound(pvBreakpointSamplePositions[nChr].begin(),pvBreakpointSamplePositions[nChr].end(), nCentromerePos));
	size_t nLastIndex = std::distance(pvBreakpointSamplePositions[nChr].begin() , pvBreakpointSamplePositions[nChr].end());
	size_t nStartIndex = (nArm==1)? 0 : nCentromereIndex + 1;
	size_t nEndIndex = (nArm==1)? nCentromereIndex : nLastIndex;
	
	//printf("Breakpointstoput for sex %d Chr %d Arm %d: %d\n", bSex, nChr, nArm, nBreakPointsToPut);
	//vector<double>::iterator it_start = lower_bound(pvAccuProb[nChr]begin(), pvAccuProb->end(), nStart);

	for (int i=0;i<nBreakPointsToPut;i++) {

		double nRand = UniformGen.Next();
		std::vector<double>::iterator oLowerBound;
		oLowerBound = lower_bound( pvAccuProb[nChr].begin() + nStartIndex, pvAccuProb[nChr].begin() + nEndIndex, nRand);
		size_t nLowerBoundIndex = std::distance(pvAccuProb[nChr].begin(), oLowerBound);
		if (nLowerBoundIndex == 0 || nLowerBoundIndex==nEndIndex) {
			// the first element
			continue; //break point cannot be inserted at the tips
		}
		//vector<double>::iterator oPreviousPoint = oLowerBound - 1;
		double nPos1 = pvBreakpointSamplePositions[nChr].at(nLowerBoundIndex - 1);
		double nPos2 = pvBreakpointSamplePositions[nChr].at(nLowerBoundIndex);
		//interpolate
		vRet.push_back( (nPos1 + nPos2) / 2 );
		
	}

	sort(vRet.begin(), vRet.end());


};

vector<double> * RecombProbConfigurations::GetBreakPointMapDistances(int nMode) {
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

PhenotypeConfigurations::PhenotypeConfigurations(void * pParentConfig) {

	this->_pParentConfig = pParentConfig;

};

void PhenotypeConfigurations::LoadFromFile(string szConfigFile) 
{
	FILE *pConfigFile;
	char szBuffer[2048];
	char szPhenotypeName[100];
	char szFormula[2048];
	
	
	
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
			pF->Evaluate(sFormula);
			vector< pair<int, double> > vSymbols;
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

				string sTemp = it->first.substr(3);				
				sTemp = sTemp.replace(sTemp.find_first_of("_"), 1, " ");
				if (sscanf(sTemp.c_str(), "%d %lf", &nSymbolChr, &nSymbolPos)!=2) 
				{
					continue;
				};
				pair<int,double> oSymbolPair( nSymbolChr, nSymbolPos);
				vSymbols.push_back(oSymbolPair);
				vSymbolStrings.push_back(it->first); // save label too
			}
			this->_mpPhenotypeFormulae[sPhenotypeName] = pF;
			this->_mpPhenotypeFormulaSymbols[sPhenotypeName] = vSymbols;
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

vector< pair<int, double> > * PhenotypeConfigurations::GetFormulaSymbols(string sPhenotypeName) {

	return &this->_mpPhenotypeFormulaSymbols[sPhenotypeName];
	
};

vector<string> * PhenotypeConfigurations::GetFormulaSymbolStrings(string sPhenotypeName) {

	return &this->_mpPhenotypeFormulaSymbolStrings[sPhenotypeName];
	
};

Parser * PhenotypeConfigurations::GetFormula(string sPhenotypeName) {
	return this->_mpPhenotypeFormulae[sPhenotypeName];
};


GeneConfigurations::GeneConfigurations(void * pParentConfig) {

	this->_pParentConfig = pParentConfig;

}

void GeneConfigurations::LoadFromFile(string szConfigFile) {

	this->_nHaploidChrNum = ((SimulationConfigurations *)this->_pParentConfig)->pMarkerConfig->GetHaploidChromosomeNum();

	//Initialize _mpGenes
	//this->_mpGenes = new map<double, GeneProperties>[_nHaploidChrNum];
	//this->_mpGeneIndex = new map<double, int>[_nHaploidChrNum];
	vector<int> vIndexCounters;

	for (int i=0;i<_nHaploidChrNum;i++) {
		map<double , GeneProperties> oPlaceHolderMap1;
		map<double, int> oPlaceHolderMap2;

		this->_mpGenes.push_back(oPlaceHolderMap1);
		this->_mpGeneIndex.push_back(oPlaceHolderMap2);

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

		this->_mpGenes[nChr-1][nPos] = oGeneProp;
		this->_mpGeneIndex[nChr-1][nPos] = vIndexCounters.at(nChr-1);
		vIndexCounters.at(nChr-1) += 1;
		//nIndexCounter++;
		
	}

	//printf("%f \n", this->_mpGenes[23][100.0].DominantFreqPop1);

}

vector< map<double , GeneProperties> > * GeneConfigurations::GetMpGenes() {
	return &this->_mpGenes;
}

vector< map<double, int> > *  GeneConfigurations::GetMpGeneIndex() {
	return &this->_mpGeneIndex;
}

SimulationConfigurations::SimulationConfigurations() {
	this->pMarkerConfig = new MarkerConfigurations(this);
	this->pRecombProbConfig = new RecombProbConfigurations(this);
	this->pPhenotypeConfig = new PhenotypeConfigurations(this);
	this->pGeneConfig = new GeneConfigurations(this);
	_szConfigFilename="";
};

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

	//Load other config files:
	this->pMarkerConfig->LoadFromFile(this->GetConfig("MarkerFile"));
	this->pRecombProbConfig->LoadFromFile(this->GetConfig("MarkerProbFile"));
	this->pPhenotypeConfig->LoadFromFile(this->GetConfig("PhenotypeFile"));
	this->pGeneConfig->LoadFromFile(this->GetConfig("GeneFile"));
	this->pMarkerConfig->CalculateMapDistances();
};

const string SimulationConfigurations::GetConfigFileName() {
	return _szConfigFilename;
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

		printf("The key %s was not found in the configuration file.\n", sKey.c_str());
		return NULL;

	}
	

	return  this->_mpConfigs[sKey];


};

const double SimulationConfigurations::GetNumericConfig(string sKey) {
	
	if (  this->_mpNumericConfigs.find(sKey) == this->_mpNumericConfigs.end()) { // the requested key is not found in the config

		printf("The key %s was not found in the configuration file.\n", sKey.c_str());
		return NULL;

	}
	

	return  this->_mpNumericConfigs[sKey];


};