#ifdef _OPENMP
 #include <omp.h>
#endif

#include "Population.h"

extern Normal NormalGen; 
extern Uniform UniformGen;
extern SimulationConfigurations SimulConfig;

Population::Population(void)
{
	this->_nCurrMales=0;
	this->_nCurrFemales = 0;
	this->_bBred = false;
}


Population::~Population(void)
{
}

void Population::Init(string sPopName,int nPopId,  char nAncestryLabel, int nPopInitSize, int nPopMaxSize, double nMaleRatio) {

	this->_nPopId = nPopId;
	this->_sPopName = sPopName;
	this->_nAncestryLabel = nAncestryLabel;
	this->_nPopInitSize = nPopInitSize;
	this->_nMaleRatio = nMaleRatio;
	this->_nPopMaxSize = nPopMaxSize;

	int nMalesToCreate = (int)((double)nPopInitSize * nMaleRatio);
	int nFemalesToCreate = nPopInitSize - nMalesToCreate;

	printf("To create... Male: %d  Female %d\n", nMalesToCreate, nFemalesToCreate );

	while(true)

	{
		
		printf("Created Male: %d  Female %d\n", this->_nCurrMales, this->_nCurrFemales );
		if (nMalesToCreate <= this->_nCurrMales && nFemalesToCreate <= this->_nCurrFemales) {
			break;
		}

		Individual * pNewInd = new Individual(this, nAncestryLabel);

		if (pNewInd->GetSex() == Individual::Male) { // If the created individual is male

			if (nMalesToCreate > this->_nCurrMales) { // if still have room for new males
				this->_mpMales.insert(this->_mpMales.begin(), pNewInd); // put the new male in the population;
				this->_nCurrMales++;
			}
			else {
				delete pNewInd;
				continue;
			}
		}

		if (pNewInd->GetSex() == Individual::Female) { // If the created individual is female

			if (nFemalesToCreate > this->_nCurrFemales) { // if still have room for new females
				this->_mpFemales.insert(this->_mpFemales.begin(), pNewInd); // put the new female in the population;
				this->_nCurrFemales++;
			}
			else {
				delete pNewInd;
				continue;
			}
		}
	}
};

void Population::SummarizePhenotype() {
	vector< string > vPhenotypes;
	SimulConfig.pPhenotypeConfig->GetKeys(vPhenotypes);

	for(vector< string >::iterator it=vPhenotypes.begin(); it != vPhenotypes.end(); ++it) { // go over each key
		string sKey = *it;
		int nPop=0, nMale=0, nFemale=0;
		double nSum=0.0, nSumMale=0.0, nSumFemale=0.0;
		double nMean=0.0, nMeanMale=0.0, nMeanFemale=0.0;
		double nSqrDiffSum=0.0, nSqrDiffSumMale=0.0, nSqrDiffSumFemale=0.0;

		//go over males
		for(vector< Individual *>::iterator itInd = _mpMales.begin(); itInd !=_mpMales.end(); ++itInd ) {
			nPop++;
			nMale++;
			double nVal = (*itInd)->GetPhenotype(sKey);
			double nValSqr = nVal * nVal;
			nSum += nVal;
			nSumMale += nVal;
		}

				//go over females
		for(vector< Individual *>::iterator itInd = _mpFemales.begin(); itInd !=_mpFemales.end(); ++itInd ) {
			nPop++;
			nFemale++;
			double nVal = (*itInd)->GetPhenotype(sKey);
			double nValSqr = nVal * nVal;
			nSum += nVal;
			nSumFemale += nVal;
		}

		nMean = (nPop == 0)? 0.0:  nSum / nPop;
		nMeanMale = (nMale == 0)? 0.0: nSumMale / nMale;
		nMeanFemale = (nFemale == 0) ? 0.0:nSumFemale / nFemale;

		//go over males
		for(vector< Individual *>::iterator itInd = _mpMales.begin(); itInd !=_mpMales.end(); ++itInd ) {
			double nVal = (*itInd)->GetPhenotype(sKey);
			double nDiff = nVal - nMean;
			double nDiffMale = nVal - nMeanMale;
			nSqrDiffSum += nDiff * nDiff;
			nSqrDiffSumMale += nDiffMale * nDiffMale;
		}

				//go over females
		for(vector< Individual *>::iterator itInd = _mpFemales.begin(); itInd !=_mpFemales.end(); ++itInd ) {
			double nVal = (*itInd)->GetPhenotype(sKey);
			double nDiff = nVal - nMean;
			double nDiffFemale = nVal - nMeanFemale;
			nSqrDiffSum += nDiff * nDiff;
			nSqrDiffSumFemale += nDiffFemale * nDiffFemale;
		}


		
		this->_mpSumPhenotype[sKey] = pair< double , double >( nMean  , sqrt((nPop==0)? 0.0:(nSqrDiffSum / nPop))); // mean, standard deviation
		this->_mpSumPhenotypeMale[sKey] = pair< double , double >( nMeanMale , sqrt((nMale==0)? 0.0:(nSqrDiffSumMale / nMale)) ); // mean, standard deviation
		this->_mpSumPhenotypeFemale[sKey] = pair< double , double >( nMeanFemale , sqrt((nFemale==0)? 0.0: (nSqrDiffSumFemale / nFemale)) ); // mean, standard deviation

	}
}

bool Population::Breed() {
	
	
	#ifdef _OPENMP
		int nThreads = (int)SimulConfig.GetNumericConfig("NumThreads");
		if (nThreads <= 0 ) {
			omp_set_num_threads( 2 );
		}
		else {
			omp_set_num_threads( nThreads );
		}
	#endif
	
	int nCurrGen = SimulConfig.GetCurrGen(); // tell what is the current generation
	bool bIgnoreGlobalRules = SimulConfig.pSexualSelConfig->IgnoreGlobalRules(nCurrGen); // is there special sexual selection rules for this generation?
	bool bIgnoreGlobalRulesNa = SimulConfig.pNaturalSelConfig->IgnoreGlobalRules(nCurrGen);
	double nSampleMate = SimulConfig.GetNumericConfig("SampleMate");

	if (_mpMales.size()==0 || _mpFemales.size()==0) {
		printf("Pop %s cannot breed because one of the sexes is 0...\n", _sPopName.c_str());
		return false;
	}

	printf("Pop %s breeding...\n", _sPopName.c_str());
	int nNewOffSpringCount = 0;
	int nNumFemales = _mpFemales.size();
	double nAvgKidPerFemale = (double)_nPopMaxSize / (double)nNumFemales;

	//for (vector<Individual *>::iterator itFemale = _mpFemales.begin(); itFemale!= _mpFemales.end(); ++itFemale) {
	//enable omp. because this is random access, ideal for parallel
	#pragma omp parallel shared(nNumFemales, nSampleMate, bIgnoreGlobalRules, nAvgKidPerFemale, bIgnoreGlobalRulesNa, nNewOffSpringCount)
	{

	#pragma omp for 
	for(int i=0;i<nNumFemales;i++) {
		//Go over each female so that they can mate.
		Individual * pFemale = this->_mpFemales[fnGetRandIndex(nNumFemales)]; //get random female

		int nCourters =  (int)NormalExt(nSampleMate , 1, 0, 100);

		
			for (int j=0;j<nCourters;j++) {
				//printf("Breed::beforerandom\n");
				//printf("_mpMale.size() %d", _mpMales.size());
				Individual * pCourter = this->_mpMales.at(fnGetRandIndex(this->_mpMales.size() ));
				//printf("Breed::afterrandom\n");
				#pragma omp critical
				{
					pFemale->HandleCourter(pCourter , bIgnoreGlobalRules);
				}
				//printf("Breed::aftercourterhandler\n");
			}
		

		//After mating, get the offspring out!

		vector<Individual *> vOffSprings;

		#pragma omp critical
		{
			pFemale->GiveBirth(vOffSprings, round(NormalExt(nAvgKidPerFemale,nAvgKidPerFemale/4, 0,100)), bIgnoreGlobalRulesNa); // to save memory, natural selection that isn't frequency dependent is carried out in the GiveBirth Function!
		}

		for (vector<Individual *>::iterator itOffSpring = vOffSprings.begin(); itOffSpring!=vOffSprings.end(); ++itOffSpring) {

			if (nNewOffSpringCount <= this->_nPopMaxSize)
			{
				if ( (*itOffSpring)->GetSex() == Individual::Male) {

					#pragma omp critical
					{
						this->_mpNewGenMales.push_back(*itOffSpring);
					}
				}
				else {
					#pragma omp critical
					{
						this->_mpNewGenFemales.push_back(*itOffSpring);
					}
				}
				#pragma omp critical
				{
					nNewOffSpringCount++;
				}
			}
			else {
				break; // no need to do anything more, already exceeded limit;
			}
		}

	}

	} //end parallel block

	this->_bBred = true; //set flag
	return true;

}

void Population::KillOldGen() { //

	printf("Pop %s's old generation was killed off...\n", _sPopName.c_str());

	if (!this->_bBred) { // the pop hasn't bred yet. don't kill the old ones!!
		throw "There is no new generation, so cannot kill the olds";
	}

	for (vector<Individual *>::iterator itMale = _mpMales.begin(); itMale!=_mpMales.end(); ++itMale) {
		delete (*itMale);
	}

	for (vector<Individual *>::iterator itFemale = _mpFemales.begin(); itFemale!=_mpFemales.end(); ++itFemale) {
		delete (*itFemale);
	}

	_mpMales.clear();
	_mpFemales.clear();

	//printf("_mpNewGenMales.size() %d\n", _mpNewGenMales.size());
	//printf("_mpNewGenFemales.size() %d\n", _mpNewGenFemales.size());
	copy(_mpNewGenMales.begin(), _mpNewGenMales.end(), back_inserter(_mpMales));
	copy(_mpNewGenFemales.begin(), _mpNewGenFemales.end(), back_inserter(_mpFemales));
	//printf("_mpMales.size() %d\n", _mpMales.size());
	//printf("_mpFemales.size() %d\n", _mpFemales.size());

	_mpNewGenMales.clear();
	_mpNewGenFemales.clear();

	this->_bBred = false;
}

Individual * Population::Emigrate() {
	int nMales = _mpMales.size();
	int nFemales = _mpFemales.size();
	int nTotal = nMales + nFemales;
	bool bWhichSex = (UniformGen.Next() <= (double)nMales / (double)nTotal)? true:false;
	vector< Individual *> * pIndividuals = bWhichSex? &_mpMales:&_mpFemales;

	if (pIndividuals->size() == 0) { // no more individuals to return
		return NULL;
	}

	int nEmigrant = fnGetRandIndex(pIndividuals->size());
	//printf("1\n");
	Individual * pEmigrant = pIndividuals->at(nEmigrant);
	//printf("2\n");
	pIndividuals->erase(pIndividuals->begin() + nEmigrant);
	//printf("3\n");
	return pEmigrant;
}

bool Population::Immigrate(Individual * pImmigrant, bool bForceExistingDie) { // Force current members to die or not if exceeding population limits
	
	if (!pImmigrant) {
		printf("Immigrant cannot be NULL\n");
		return false;
	}
	//printf("1\n");
	pImmigrant->ChangePopulation(this);
	Individual::Sex bWhichSex = pImmigrant->GetSex();
	//printf("2\n");
	vector< Individual *> * pIndividuals = (bWhichSex==Individual::Male)? &_mpMales:&_mpFemales;
	//printf("3\n");
	if (pIndividuals->size() < _nPopMaxSize) { // Still within limit of max size
		//printf("4\n");
		pIndividuals->push_back(pImmigrant); //add this individual to the population.
		//printf("5\n");
		return true;
	}
	//printf("MaxPopSize: %d \n", _nPopMaxSize );
	if ( (pIndividuals->size() >= _nPopMaxSize) && bForceExistingDie) {
		// Randomly delete one individual from the current pop.
		//printf("Immigrate::6\n");
		pIndividuals->erase(pIndividuals->begin() + fnGetRandIndex( pIndividuals->size() ));
		//printf("Immigrate::7\n");
		pIndividuals->push_back(pImmigrant);
		//printf("Immigrate::8\n");
		return true;
	}

	return false; // already full, cannot accept any more.

}

int Population::GetPopSize(int nMode) {
	switch(nMode) {
		case 0: return _mpFemales.size();
		case 1: return _mpMales.size();
		case 2: return _mpNewGenFemales.size();
		case 3: return _mpNewGenMales.size();
		case 4: return _mpFemales.size()+_mpMales.size();
		case 5: return _mpNewGenFemales.size()+_mpNewGenMales.size();
	}
}

void Population::Sample(ofstream &fMarkerOutFile, ofstream &fGeneOutFile,  ofstream &fPhenotypeOutFile,  ofstream &fPhenoSumOutFile ) {
	// now dump all the males:
	this->fnSample(fMarkerOutFile, fGeneOutFile,  fPhenotypeOutFile , Individual::Male);
	// dump all females:
	this->fnSample(fMarkerOutFile, fGeneOutFile,  fPhenotypeOutFile , Individual::Female);

	this->fnSamplePhenotypeStats(fPhenoSumOutFile);
}

void Population::fnSamplePhenotypeStats(ofstream &fPhenoSumOutFile) {

	map< string , pair< double, double> >::iterator itMale = this->_mpSumPhenotypeMale.begin();
	map< string , pair< double, double> >::iterator itFemale = this->_mpSumPhenotypeFemale.begin();

	for ( map< string , pair< double, double> >::iterator itAll = this->_mpSumPhenotype.begin(); itAll != _mpSumPhenotype.end(); ++itAll) {

		fPhenoSumOutFile << this->GetPopId() << '\t' << itAll->first << '\t';
		
		fPhenoSumOutFile << itAll->second.first << '[' << itAll->second.second << "]\t"
			<< itMale->second.first << '[' << itMale->second.second << "]\t"
			<< itFemale->second.first << '[' << itFemale->second.second << ']' << endl;

		++itMale;
		++itFemale;
	}
}

void Population::fnSample(ofstream &fMarkerOutFile, ofstream &fGeneOutFile,  ofstream &fPhenotypeOutFile, Individual::Sex bSex) {

	bool bWriteMarkers = !(SimulConfig.GetConfig("MarkerOutput") == "Off");
	
	vector< Individual *> * pMpIndividuals = (bSex == Individual::Male)? &this->_mpMales : &this->_mpFemales;
	
	

	int i=0;

	for (vector< Individual *>::iterator itInd=pMpIndividuals->begin(); itInd!=pMpIndividuals->end(); ++ itInd) {

		if (bWriteMarkers) {
			this->fnWriteIndividualMarkers(fMarkerOutFile, *itInd);
		}
		this->fnWriteIndividualGenes(fGeneOutFile, *itInd);

		if (i==0) {
				// phenotype file needs header.
				fPhenotypeOutFile << "id\tpop\tsex\tFatherId\tMotherId\t";
				(*itInd)->WritePhenotypeHeader(fPhenotypeOutFile);
				fPhenotypeOutFile << endl;
		}

		this->fnWriteIndividualPhenotypes(fPhenotypeOutFile, *itInd);

		i++;
	}

}

void Population::fnWriteIndividualMarkers(ofstream &fOutFile, Individual * pInd) {

	this->fnWriteIndividualID(fOutFile , pInd);
	pInd->DumpMarkers(fOutFile, 0);
	fOutFile << endl;

	this->fnWriteIndividualID(fOutFile , pInd);
	pInd->DumpMarkers(fOutFile, 1);
	fOutFile << endl;
}

void Population::fnWriteIndividualGenes(ofstream &fOutFile, Individual * pInd) {

	this->fnWriteIndividualID(fOutFile , pInd);
	pInd->DumpGenes(fOutFile, 0);
	fOutFile << endl;

	this->fnWriteIndividualID(fOutFile , pInd);
	pInd->DumpGenes(fOutFile, 1);
	fOutFile << endl;
}

void Population::fnWriteIndividualPhenotypes(ofstream &fOutFile, Individual * pInd) {

	this->fnWriteIndividualID(fOutFile , pInd);
	pInd->DumpPhenotypes(fOutFile);
	fOutFile << endl;

}

void Population::fnWriteIndividualID(ofstream &fOutFile, Individual * pInd) {
	fOutFile << pInd->GetID() << '\t';
	fOutFile << this->_nPopId << '\t';
	fOutFile << ((pInd->GetSex() == Individual::Male)? "1":"0") << '\t';
	fOutFile << pInd->GetFatherId() << '\t';
	fOutFile << pInd->GetMotherId() << '\t';
}

int Population::GetPopId() {
	return this->_nPopId;
};

string  Population::GetPopName() {
	return this->_sPopName;
};