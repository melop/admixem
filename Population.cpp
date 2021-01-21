#ifdef _OPENMP
 #include <omp.h>
#endif

#include "Population.h"
extern int nTotalCPUCore;
extern Normal NormalGen; 
extern Uniform UniformGen;
extern SimulationConfigurations SimulConfig;

Population::Population(void)
{
	this->_nCurrMales=0;
	this->_nCurrFemales = 0;
	this->_bBred = false;
	this->_bRecordNatSelProb = (SimulConfig.GetConfig("DumpNatSelProb") == "On");

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

	this->_mpvNewGenMales = new vector< Individual *>[nTotalCPUCore];
	this->_mpvNewGenFemales = new vector< Individual *>[nTotalCPUCore];


	printf("To create... Male: %d  Female %d\n", nMalesToCreate, nFemalesToCreate );
	this->_mpMales.reserve(nMalesToCreate);
	this->_mpFemales.reserve(nFemalesToCreate);
	#pragma omp parallel shared(nMalesToCreate, nFemalesToCreate) 
	//private(pFemale, nCourters, pCourter, vOffSprings, itOffSpring)
	{
	


	#pragma omp for 
	for (int nCPU=0; nCPU<nTotalCPUCore;nCPU++)
	{
		
		//printf("Created Male: %d  Female %d\n", this->_nCurrMales, this->_nCurrFemales );
		while(true) {
			if (nMalesToCreate <= this->_nCurrMales && nFemalesToCreate <= this->_nCurrFemales) {
				break;
			}

			Individual * pNewInd = new Individual(this, nAncestryLabel);

			if (pNewInd->GetSex() == Individual::Male) { // If the created individual is male

				if (nMalesToCreate > this->_nCurrMales) { // if still have room for new males
					#pragma omp critical 
					{
					this->_mpMales.push_back(pNewInd); // put the new male in the population;
					this->_nCurrMales++;
					}
				}
				else {
					delete pNewInd;
					continue;
				}
			}

			if (pNewInd->GetSex() == Individual::Female) { // If the created individual is female

				if (nFemalesToCreate > this->_nCurrFemales) { // if still have room for new females
					#pragma omp critical 
					{
						this->_mpFemales.push_back(pNewInd); // put the new female in the population;
						this->_nCurrFemales++;
					}
				}
				else {
					delete pNewInd;
					continue;
				}
			}
		}//while end.
	} //for end

	} //end parallel
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

		/* Save the phenotype for the previous generation , used for oblique imprinting*/
		this->_mpPrevGenSumPhenotype[sKey] = (_mpSumPhenotype.find(sKey) == _mpSumPhenotype.end())? pair< double, double> ( 0.0 , 0.0 ) : this->_mpSumPhenotype[sKey];
		this->_mpPrevGenSumPhenotypeMale[sKey] = (_mpSumPhenotypeMale.find(sKey) == _mpSumPhenotypeMale.end())? pair< double, double> ( 0.0 , 0.0 ) : this->_mpSumPhenotypeMale[sKey]; 
		this->_mpPrevGenSumPhenotypeFemale[sKey] = (_mpSumPhenotypeFemale.find(sKey) == _mpSumPhenotypeFemale.end())? pair< double, double> ( 0.0 , 0.0 ) : this->_mpSumPhenotypeFemale[sKey];  ; 
		
		this->_mpSumPhenotype[sKey] = pair< double , double >( nMean  , sqrt((nPop==0)? 0.0:(nSqrDiffSum / nPop))); // mean, standard deviation
		this->_mpSumPhenotypeMale[sKey] = pair< double , double >( nMeanMale , sqrt((nMale==0)? 0.0:(nSqrDiffSumMale / nMale)) ); // mean, standard deviation
		this->_mpSumPhenotypeFemale[sKey] = pair< double , double >( nMeanFemale , sqrt((nFemale==0)? 0.0: (nSqrDiffSumFemale / nFemale)) ); // mean, standard deviation

	}
}

bool Population::Breed() {
	
	
	int nCurrGen = SimulConfig.GetCurrGen(); // tell what is the current generation
	bool bIgnoreGlobalRules = SimulConfig.pSexualSelConfig->IgnoreGlobalRules(nCurrGen); // is there special sexual selection rules for this generation?
	bool bIgnoreGlobalRulesNa = SimulConfig.pNaturalSelConfig->IgnoreGlobalRules(nCurrGen);
	double nSampleMate = SimulConfig.GetNumericConfig("SampleMate");
	string sOffSpringCountFunc = SimulConfig.GetConfig("kids_per_female_func");
	int nOffSpringCountFunc = (sOffSpringCountFunc=="Poisson")? 1:0; //1 - Poisson, 0- normal distribution with variance being 1/4 of mean

	if (_mpMales.size()==0 || _mpFemales.size()==0) {
		printf("Pop %s cannot breed because one of the sexes is 0...\n", _sPopName.c_str());
		this->_bBred = true;
		return false;
	}

	//Check if natural selection prob or sexual selection prob is zero, if so, no need to continue.
	list< pair< Parser *, int> > * pqParsers = SimulConfig.pNaturalSelConfig->GetFormulae(this->GetPopName());

	for (list< pair< Parser *, int> >::iterator itParser= pqParsers->begin(); itParser != pqParsers->end() ; ++itParser) {
		int nGen = itParser->second;
		if ((nGen ==-1 && !bIgnoreGlobalRules) || (bIgnoreGlobalRulesNa && nGen == nCurrGen) ) {
			if (itParser->first->IsZero) {
				this->_bBred = true;
				return true; //no need to continue since surviving prob is zero.
			}
		}
	}

	printf("Pop %s breeding...\n", _sPopName.c_str());
	int nNewOffSpringCount = 0;
	int nNumFemales = _mpFemales.size();
	
	vector< Poisson* > vOffSpringPoissonGen; 
	

	//set population level variables for freq dependent sexual selection:

	
	map< int, map<string , list< pair<Parser *, int> > > > * mppqParsers = SimulConfig.pSexualSelConfig->GetFormulaeAllCPUs();

	list< vector<string> > * pqvPopCourterSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsPopCourter( this->_sPopName);
	list< vector<string> > * pqvPopChooserSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsPopChooser( this->_sPopName);
	list< vector<string> > * pqvPopSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsPop( this->_sPopName);



	for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) { //set global current population parameters for all the CPUs

		list< pair< Parser *, int > > * pqParsers = &(*mppqParsers)[nCPU][this->_sPopName];
		list< vector<string> >::iterator itPopCourterSymbols= pqvPopCourterSymbols->begin();
		list< vector<string> >::iterator itPopChooserSymbols= pqvPopChooserSymbols->begin();
		list< vector<string> >::iterator itPopSymbols= pqvPopSymbols->begin();

		for (list< pair< Parser *, int > > ::iterator itParser= pqParsers->begin(); itParser != pqParsers->end() ; ++itParser) {
				vector<string> vSymbolsPopCourter = *itPopCourterSymbols;
				vector<string> vSymbolsPopChooser = *itPopChooserSymbols;
				vector<string> vSymbolsPop = *itPopSymbols;

				Parser * pParser = itParser->first;

				//Set courter symbol values
				for(vector<string>::iterator itSymbol=vSymbolsPop.begin();itSymbol!=vSymbolsPop.end();++itSymbol)
				{
					string sSymbol = (*itSymbol).substr(4);
					string sType = (*itSymbol).substr(0, 3); //either Avg or Std
					if (sType != "Avg" && sType != "Std") {
						throw new Exception("Population parameter type unknown!");
					}
					if (this->_mpSumPhenotype.find(sSymbol) == _mpSumPhenotype.end()) {
						throw new Exception("Unable to find population symbol");
					}
			
					pParser->symbols_[string("Pop_Avg_"+sSymbol)] = this->_mpSumPhenotype[sSymbol].first;
					pParser->symbols_[string("Pop_Std_"+sSymbol)] = this->_mpSumPhenotype[sSymbol].second;
					cout << "CPU:" << nCPU << "Pop_Avg_" << sSymbol.c_str() << "=" << this->_mpSumPhenotype[sSymbol].first << "\n";
					cout << "CPU:" << nCPU << "Pop_Std_" << sSymbol.c_str() << "=" << this->_mpSumPhenotype[sSymbol].second << "\n";
				}

				for(vector<string>::iterator itSymbol=vSymbolsPopCourter.begin();itSymbol!=vSymbolsPopCourter.end();++itSymbol)
				{
					string sSymbol = (*itSymbol).substr(4);
					string sType = (*itSymbol).substr(0, 3); //either Avg or Std
					if (sType != "Avg" && sType != "Std") {
						throw new Exception("Population parameter type unknown!");
					}
					if (this->_mpSumPhenotypeMale.find(sSymbol) == _mpSumPhenotypeMale.end()) {
						throw new Exception("Unable to find population symbol: male");
					}
			
					pParser->symbols_[string("PopCourter_Avg_"+sSymbol)] = this->_mpSumPhenotypeMale[sSymbol].first;
					pParser->symbols_[string("PopCourter_Std_"+sSymbol)] = this->_mpSumPhenotypeMale[sSymbol].second;
					cout << "CPU:" << nCPU << "PopCourter_Avg_" << sSymbol.c_str() << "=" << this->_mpSumPhenotypeMale[sSymbol].first << "\n";
					cout << "CPU:" << nCPU << "PopCourter_Std_" << sSymbol.c_str() << "=" << this->_mpSumPhenotypeMale[sSymbol].second << "\n";
				}

				for(vector<string>::iterator itSymbol=vSymbolsPopChooser.begin();itSymbol!=vSymbolsPopChooser.end();++itSymbol)
				{
					string sSymbol = (*itSymbol).substr(4);
					string sType = (*itSymbol).substr(0, 3); //either Avg or Std
					if (sType != "Avg" && sType != "Std") {
						throw new Exception("Population parameter type unknown!");
					}
					if (this->_mpSumPhenotypeFemale.find(sSymbol) == _mpSumPhenotypeFemale.end()) {
						throw new Exception("Unable to find population symbol: female");
					}
			
					pParser->symbols_[string("PopChooser_Avg_"+sSymbol)] = this->_mpSumPhenotypeFemale[sSymbol].first;
					pParser->symbols_[string("PopChooser_Std_"+sSymbol)] = this->_mpSumPhenotypeFemale[sSymbol].second;
					cout << "CPU:" << nCPU << "PopChooser_Avg_" << sSymbol.c_str() << "=" << this->_mpSumPhenotypeFemale[sSymbol].first << "\n";
					cout << "CPU:" << nCPU << "PopChooser_Std_" << sSymbol.c_str() << "=" << this->_mpSumPhenotypeFemale[sSymbol].second << "\n";


				}
			
				++itPopCourterSymbols;
				++itPopChooserSymbols;
				++itPopSymbols;
		}
	}

	//for (vector<Individual *>::iterator itFemale = _mpFemales.begin(); itFemale!= _mpFemales.end(); ++itFemale) {
	//enable omp. because this is random access, ideal for parallel
	std::set<int> stSampledFemales;//keep a list of the females already sampled so that they're not sampled twice.
	std::set<int> stExhaustedFemales; // a list of females without further gametes
	//int nExhaustedFemales=0;
	bool bCourterHandeled = false;
	int nMaxLoop = SimulConfig.GetNumericConfig("FemaleGiveBirthMaxIterations");
	nMaxLoop = (nMaxLoop==-1)? 1000:nMaxLoop; //do at most 20 loops.
	int nLoopCount = 0;
	int nSameAvgKidIterCount = 0;
	double nPrevKidPerFemale = -1.0;
	do 
	{
		nLoopCount++;
		if (nNewOffSpringCount >= this->_nPopMaxSize) {
			printf("Max Pop Size met.\n");
			break; //no need to continue;
		}
		//#pragma omp critical
		//{
			if (stSampledFemales.size() != 0) {
				stSampledFemales.clear();
			}
		//}
	double nAvgKidPerFemale = (double)(_nPopMaxSize - nNewOffSpringCount) / (double)nNumFemales; // average kid per female, given the number offsprings left to fill the pop
	vOffSpringPoissonGen.clear();
		
	if (nPrevKidPerFemale == nAvgKidPerFemale ) {
		//no addition;
		if (nSameAvgKidIterCount > nMaxLoop) {
			printf("No more breeding possible\n");
			break;
		} else {
			nSameAvgKidIterCount++;	
		}
	}
	nPrevKidPerFemale = nAvgKidPerFemale;
		
	for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) { //set global current population parameters for all the CPUs

		//initialize poisson generators for offspring numbers
		if (nOffSpringCountFunc == 1 ) {
			Poisson * OffSpPois = new Poisson(nAvgKidPerFemale);
			vOffSpringPoissonGen.push_back(OffSpPois);
		}
	}

	printf("AvgKidPerFemale: %f\n", nAvgKidPerFemale);

	#pragma omp parallel shared(bCourterHandeled, stExhaustedFemales, stSampledFemales, nNumFemales, nSampleMate, bIgnoreGlobalRules, nAvgKidPerFemale, bIgnoreGlobalRulesNa, nNewOffSpringCount) 
	//private(pFemale, nCourters, pCourter, vOffSprings, itOffSpring)
	{
	


	#pragma omp for 
		//reduction(+: nNewOffSpringCount) 
	for(int i=0;i<nNumFemales;i++) 
	{
		if (nNewOffSpringCount >= this->_nPopMaxSize || stExhaustedFemales.size() >= this->_mpFemales.size()) {
			continue; // first see if new pop already filled up by other threads. if so then do nothing.
		}
		
		#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
			#ifdef DEBUG
			printf("CPU %d: %s : %d\n", nCPU, "SexSelFormulae # ", nCPU);
			#endif
		#else
			int nCPU = 0;
		#endif

		//Go over each female so that they can mate.
		int nRandFemaleInd=fnGetRandIndex(nNumFemales);
		bool bGetRandSuccess = true;
		
			while((stSampledFemales.find(nRandFemaleInd)!= stSampledFemales.end()) && (stExhaustedFemales.find(nRandFemaleInd) == stExhaustedFemales.end())) {
				nRandFemaleInd=fnGetRandIndex(nNumFemales);
				if (stExhaustedFemales.size() >= this->_mpFemales.size()) {
					bGetRandSuccess = false;
					break;
				}
			}
			if (bGetRandSuccess) {
				#pragma omp critical
				{
					stSampledFemales.insert(nRandFemaleInd);
				}
			}
		#ifdef DEBUG
		printf("CPU %d: Female Id: %d\n", nCPU, nRandFemaleInd);
		#endif
		if (!bGetRandSuccess) {
			continue;
		}

		Individual * pFemale = this->_mpFemales.at(nRandFemaleInd); //get random female

		if (pFemale->GetGameteNum() == 0) 
		{
				#pragma omp critical
				{
			
			
						stExhaustedFemales.insert(nRandFemaleInd); // this female cannot produce more offsprings.

				}
				continue; //if already no gamete then don't bother.
				//nExhaustedFemales++;
		}

		if (!bCourterHandeled) {
		int nCourters = (int)nSampleMate; //(int)NormalExt(nSampleMate , 1, 0, 100);

		
			for (int j=0;j<nCourters;j++) {
				#ifdef DEBUG
				printf("CPU %d: Breed::beforerandom\n", nCPU);
				printf("CPU %d: _mpMale.size() %d\n", nCPU, _mpMales.size());
				#endif
				
				Individual * pCourter = this->_mpMales.at(fnGetRandIndex(this->_mpMales.size() ));

				#ifdef DEBUG
				printf("CPU %d: Breed::afterrandom\n", nCPU);
				#endif
				//#pragma omp critical
				//{
					pFemale->HandleCourter(pCourter , bIgnoreGlobalRules);
				//}
				#ifdef DEBUG
				printf("CPU %d : Breed::aftercourterhandler\n", nCPU);
				#endif
			}
		}
		

		//After mating, get the offspring out!

		vector<Individual *> vOffSprings;

		int nTargetOffSpringCount =  round(NormalExt(nAvgKidPerFemale,nAvgKidPerFemale/4, 0,10000));

		if ( nOffSpringCountFunc == 1 ) { //use poisson distribution for offspring number
			#ifdef DEBUG
			printf("CPU %d : Getting Target offspring...\n", nCPU );
			#endif
			nTargetOffSpringCount = vOffSpringPoissonGen[nCPU]->Next();
			#ifdef DEBUG
			printf("CPU %d : Target offspring: %d\n", nCPU, nTargetOffSpringCount );
			#endif
		}
		//#pragma omp critical
		//{
		if (nTargetOffSpringCount == 0) {
			continue;
		}

		#ifdef DEBUG

		printf("CPU %d: Breed::GiveBirth before\n", nCPU);
		#endif
		
		pFemale->GiveBirth(vOffSprings, nTargetOffSpringCount, bIgnoreGlobalRulesNa); // to save memory, natural selection that isn't frequency dependent is carried out in the GiveBirth Function!
		
		#ifdef DEBUG

		printf("CPU %d: Breed::GiveBirth after\n", nCPU);
		#endif
		//}

		if (pFemale->GetGameteNum() == 0) 
		{
			#ifdef DEBUG

			printf("CPU %d: Exhausted female before\n", nCPU);
			#endif
			#pragma omp critical
			{
			
			
				stExhaustedFemales.insert(nRandFemaleInd); // this female cannot produce more offsprings.
				//nExhaustedFemales++;
			}

			#ifdef DEBUG
			printf("CPU %d: Exhausted female after\n", nCPU);
			#endif
		}

		#ifdef DEBUG
		printf("CPU %d: Start insert offspring\n", nCPU);
		#endif
		//#pragma omp critical
		//{
		for (vector<Individual *>::iterator itOffSpring = vOffSprings.begin(); itOffSpring!=vOffSprings.end(); ++itOffSpring) 
		{

			if (nNewOffSpringCount < this->_nPopMaxSize)
			{
				if ( (*itOffSpring)->GetSex() == Individual::Male) {

					//#pragma omp critical
						this->_mpvNewGenMales[nCPU].push_back(*itOffSpring);
					//}
				}
				else {
					//#pragma omp critical
					//{
						this->_mpvNewGenFemales[nCPU].push_back(*itOffSpring);
					//}
				}
				#pragma omp critical
				{
					nNewOffSpringCount++;
				}
			}
			else {

				printf("No need to insert more offsprings, break;\n");
				break; // no need to do anything more, already exceeded limit;
			}
		}

		//} // end omp
		#ifdef DEBUG
		printf("CPU %d: end insert offspring\n", nCPU);
		#endif

	}

		
	} //end parallel block

	bCourterHandeled = true;
		
		for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) { //set global current population parameters for all the CPUs
			//delete poisson generators for offspring numbers
			if (nOffSpringCountFunc == 1 ) {
				#ifdef DEBUG
				printf("CPU %d: Realease poisson generator\n", nCPU);
				#endif
				delete vOffSpringPoissonGen[nCPU];
				vOffSpringPoissonGen[nCPU] = NULL;
			}
		}


	} // end do block
	while( (nLoopCount <= nMaxLoop) && (nNewOffSpringCount < this->_nPopMaxSize) && (stExhaustedFemales.size() < this->_mpFemales.size()) ) ; //end do while block.
	/* needs more work
	int nExcessOffsprings = nNewOffSpringCount - this->_nPopMaxSize;
	int nExcessMales = nExcessOffsprings / 2;
	int nExcessFemales = nExcessOffsprings - nExcessMales;
	int nDeletedMales = 0;
	int nDeletedFemales = 0;
	
	if (nExcessOffsprings > 0) {
		printf("Delete %d extra males and %d females (total %d)\n", nExcessMales,  nExcessFemales, nExcessOffsprings);	
	}
	
		for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) { //set global current population parameters for all the CPUs

			if (nDeletedMales < nExcessMales && (!this->_mpvNewGenMales[nCPU].empty())) {
				Individual *pToDelete = (this->_mpvNewGenMales[nCPU].back());
				this->_mpvNewGenMales[nCPU].pop_back();
				delete pToDelete;
				nDeletedMales++;
			}
                        if (nDeletedFemales < nExcessFemales && (!this->_mpvNewGenFemales[nCPU].empty())) {
                                Individual *pToDelete = (this->_mpvNewGenFemales[nCPU].back());
                                this->_mpvNewGenFemales[nCPU].pop_back();
				delete pToDelete;
                                nDeletedFemales++;
                        }
		}
	*/
	
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
	_mpImCacheMales.clear();
	_mpImCacheFemales.clear();

	//printf("_mpNewGenMales.size() %d\n", _mpNewGenMales.size());
	//printf("_mpNewGenFemales.size() %d\n", _mpNewGenFemales.size());
	_mpMales.reserve(this->GetPopSize(2));
	_mpFemales.reserve(this->GetPopSize(3));

	for(int nCpu=0; nCpu<nTotalCPUCore;nCpu++) {
		copy(_mpvNewGenMales[nCpu].begin(), _mpvNewGenMales[nCpu].end(), back_inserter(_mpMales));
		copy(_mpvNewGenFemales[nCpu].begin(), _mpvNewGenFemales[nCpu].end(), back_inserter(_mpFemales));

		_mpvNewGenMales[nCpu].clear();
		_mpvNewGenFemales[nCpu].clear();
	}
	//printf("_mpMales.size() %d\n", _mpMales.size());
	//printf("_mpFemales.size() %d\n", _mpFemales.size());
	
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

bool Population::Immigrate(Individual * pImmigrant) { // Force current members to die or not if exceeding population limits
	
	if (!pImmigrant) {
		printf("Immigrant cannot be NULL\n");
		return false;
	}
	//printf("1\n");
	pImmigrant->ChangePopulation(this);
	Individual::Sex bWhichSex = pImmigrant->GetSex();
	//printf("2\n");
	vector< Individual *> * pIndividuals = (bWhichSex==Individual::Male)?  &_mpImCacheMales : &_mpImCacheFemales; //&_mpMales:&_mpFemales;
	//printf("3\n");

		pIndividuals->push_back(pImmigrant); //add this individual to the population.


	return true; // already full, cannot accept any more.

}

bool Population::ImmigrateConfirm(bool bForceExistingDie) {
	int nCacheMaleSize = _mpImCacheMales.size();
	int nCacheFemaleSize = _mpImCacheFemales.size();
	int nTotalCached = nCacheMaleSize + nCacheFemaleSize;
	int nTotalToAdd;
	int nCurrPopSize = this->GetPopSize(4);
	int nRoomLeft = _nPopMaxSize-nCurrPopSize >=0? _nPopMaxSize-nCurrPopSize : 0;
	double nProbToAdd;

	if (bForceExistingDie) {
		nTotalToAdd = nTotalCached < _nPopMaxSize? nTotalCached : _nPopMaxSize;
		nProbToAdd = (double)nTotalToAdd / (double)nTotalCached;
	} else {
		if (nRoomLeft == 0) {
			return false; //no room to add any immigrants
		}
		nTotalToAdd = nTotalCached < nRoomLeft? nTotalCached : nRoomLeft;
		nProbToAdd = (double)nTotalToAdd / (double)nTotalCached;
	}
		//printf("MaxPopSize: %d \n", _nPopMaxSize );
	printf("Candidate immigrating males %d\n", nCacheMaleSize);
	printf("Candidate immigrating females %d\n", nCacheFemaleSize);

	int nAddedMales=0;
	int nAddedFemales=0;

	if (bForceExistingDie && nRoomLeft < nTotalToAdd) {
		int nMakeMoreRoom = nTotalToAdd - nRoomLeft;
		for(int i=0;i< nMakeMoreRoom;i++) {
			Individual * pVictim = this->Emigrate();
			if (!pVictim) {
				
			} else {
				delete pVictim;
			}
		}
	}

	for(vector< Individual * >::iterator it=this->_mpImCacheMales.begin(); it!=this->_mpImCacheMales.end();++it) {
		if (UniformGen.Next() < nProbToAdd) {
			//add

			_mpMales.push_back( *it );
			nAddedMales++;
		}
	}

	for(vector< Individual * >::iterator it=this->_mpImCacheFemales.begin(); it!=this->_mpImCacheFemales.end();++it) {
		if (UniformGen.Next() < nProbToAdd) {
			//add
			_mpFemales.push_back( *it );
			nAddedFemales++;
		}
	}

	printf("Added immigrating males %d\n", nAddedMales);
	printf("Added immigrating females %d\n", nAddedFemales);
	//do final check on pop size:
	nCurrPopSize = this->GetPopSize(4);
	if (nCurrPopSize > _nPopMaxSize) {
		//randomly remove some individuals
		while(this->GetPopSize(4) > _nPopMaxSize) {
			Individual * pVictim = this->Emigrate();
			if (!pVictim) {
				
			} else {
				delete pVictim;
			}
		}
	}
}

int Population::GetPopSize(int nMode) {
	switch(nMode) {
		case 0: return _mpFemales.size();
		case 1: return _mpMales.size();
		case 2: {int nCount=0; for( int nCpu=0;nCpu<nTotalCPUCore;nCpu++) {nCount+= _mpvNewGenMales[nCpu].size();} return nCount;};
		case 3: {int nCount=0; for( int nCpu=0;nCpu<nTotalCPUCore;nCpu++) {nCount+= _mpvNewGenFemales[nCpu].size();} return nCount;};
		case 4: return _mpFemales.size()+_mpMales.size();
		case 5: {int nCount=0; for( int nCpu=0;nCpu<nTotalCPUCore;nCpu++) {nCount+=_mpvNewGenMales[nCpu].size()+ _mpvNewGenFemales[nCpu].size();} return nCount;};
	}
}

void Population::Sample(ofstream &fMarkerOutFile, ofstream &fGeneOutFile,  ofstream &fPhenotypeOutFile,  ofstream &fPhenoSumOutFile, ofstream &fOffSpringNatSelProb ) {
	// now dump all the males:
	this->fnSample(fMarkerOutFile, fGeneOutFile,  fPhenotypeOutFile , Individual::Male);
	// dump all females:
	this->fnSample(fMarkerOutFile, fGeneOutFile,  fPhenotypeOutFile , Individual::Female);

	this->fnSamplePhenotypeStats(fPhenoSumOutFile);

	this->fnDumpNaturalProb(fOffSpringNatSelProb);

	this->_mpOffSpringNaturalProb.clear();
}

void Population::FreqDependentNaturalSelection() {

	int nCurrGen = SimulConfig.GetCurrGen();
	bool bIgnoreGlobalRules = SimulConfig.pNaturalSelConfig->IgnoreGlobalRules(nCurrGen);
	string sPop = this->_sPopName;

	map< int, map<string , list< pair<Parser *, int> > > > * mppqParsers = SimulConfig.pNaturalSelConfig->GetFreqDependentFormulaeAllCPUs();


	//list< vector<string> > * pqvCourterSymbols = SimulConfig.pNaturalSelConfig->GetFormulaSymbolStringsCourter( sPop);
	list< vector<string> > * pqvSelfSymbols = SimulConfig.pNaturalSelConfig->GetFormulaSymbolStringsSelf( sPop);
	list< vector<string> > * pqvPopSymbols = SimulConfig.pNaturalSelConfig->GetFormulaSymbolStringsPop( sPop);
	list< vector<string> > * pqvPopCourterSymbols = SimulConfig.pNaturalSelConfig->GetFormulaSymbolStringsPopCourter( sPop);
	list< vector<string> > * pqvPopChooserSymbols = SimulConfig.pNaturalSelConfig->GetFormulaSymbolStringsPopChooser( sPop);

	//set population level parameters


	for(int nCPU=0;nCPU<nTotalCPUCore;nCPU++) { //set global current population parameters for all the CPUs

		list< pair< Parser *, int > > * pqParsers = &(*mppqParsers)[nCPU][this->_sPopName];
		list< vector<string> >::iterator itPopSymbols = pqvPopSymbols->begin();
		list< vector<string> >::iterator itPopCourterSymbols = pqvPopCourterSymbols->begin();
		list< vector<string> >::iterator itPopChooserSymbols = pqvPopChooserSymbols->begin();

		for (list< pair< Parser *, int> >::iterator itParser= pqParsers->begin(); itParser != pqParsers->end() ; ++itParser) {
			vector<string> vSymbolsPop = *itPopSymbols;
			vector<string> vSymbolsPopCourter = *itPopCourterSymbols;
			vector<string> vSymbolsPopChooser = *itPopChooserSymbols;

			Parser * pParser = itParser->first;
			//set population level symbols

			for(vector<string>::iterator itSymbol=vSymbolsPop.begin();itSymbol!=vSymbolsPop.end();++itSymbol)
			{
				string sSymbol = (*itSymbol).substr(4);
				string sType = (*itSymbol).substr(0, 3); //either Avg or Std
				if (sType != "Avg" && sType != "Std") {
					throw new Exception("Population parameter type unknown!");
				}
				if (this->_mpSumPhenotype.find(sSymbol) == _mpSumPhenotype.end()) {
					throw new Exception("Unable to find population symbol");
				}
			
				pParser->symbols_[string("Pop_Avg_"+sSymbol)] = this->_mpSumPhenotype[sSymbol].first;
				pParser->symbols_[string("Pop_Std_"+sSymbol)] = this->_mpSumPhenotype[sSymbol].second;

			}

			for(vector<string>::iterator itSymbol=vSymbolsPopCourter.begin();itSymbol!=vSymbolsPopCourter.end();++itSymbol)
			{
				string sSymbol = (*itSymbol).substr(4);
				string sType = (*itSymbol).substr(0, 3); //either Avg or Std
				if (sType != "Avg" && sType != "Std") {
					throw new Exception("Population parameter type unknown!");
				}
				if (this->_mpSumPhenotypeMale.find(sSymbol) == _mpSumPhenotypeMale.end()) {
					throw new Exception("Unable to find population symbol: male");
				}
			
				pParser->symbols_[string("PopCourter_Avg_"+sSymbol)] = this->_mpSumPhenotypeMale[sSymbol].first;
				pParser->symbols_[string("PopCourter_Std_"+sSymbol)] = this->_mpSumPhenotypeMale[sSymbol].second;

			}

			for(vector<string>::iterator itSymbol=vSymbolsPopChooser.begin();itSymbol!=vSymbolsPopChooser.end();++itSymbol)
			{
				string sSymbol = (*itSymbol).substr(4);
				string sType = (*itSymbol).substr(0, 3); //either Avg or Std
				if (sType != "Avg" && sType != "Std") {
					throw new Exception("Population parameter type unknown!");
				}
				if (this->_mpSumPhenotypeFemale.find(sSymbol) == _mpSumPhenotypeFemale.end()) {
					throw new Exception("Unable to find population symbol: female");
				}
			
				pParser->symbols_[string("PopChooser_Avg_"+sSymbol)] = this->_mpSumPhenotypeFemale[sSymbol].first;
				pParser->symbols_[string("PopChooser_Std_"+sSymbol)] = this->_mpSumPhenotypeFemale[sSymbol].second;

			}

		
			++itPopSymbols;
			++itPopCourterSymbols;
			++itPopChooserSymbols;
		}

	}
	//now go through each individual 
	vector<int>::size_type nMale = this->_mpMales.size();
	vector<int>::size_type nFemale = this->_mpFemales.size();

	#pragma omp parallel 
	{

	#pragma omp for
	//for (vector< Individual * >::size_type i=0; i< nMale; i++ ) {
	for (signed long i=0; i< nMale; i++ ) {
		Individual * pInd = _mpMales[i];

		#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
		#else
		int nCPU = 0;
		#endif

		list< pair< Parser *, int > > * pqParsers = &(*mppqParsers)[nCPU][this->_sPopName];
		list< vector<string> >::iterator itSelfSymbols = pqvSelfSymbols->begin();
		for (list< pair< Parser *, int> >::iterator itParser= pqParsers->begin(); itParser != pqParsers->end() ; ++itParser) {
			vector<string> vSymbolsSelf = *itSelfSymbols;
			Parser * pParser = itParser->first;
			int nGen = itParser->second;

			if ((nGen ==-1 && !bIgnoreGlobalRules) || (bIgnoreGlobalRules && nGen == nCurrGen) ) {

				for(vector<string>::iterator itSymbol=vSymbolsSelf.begin();itSymbol!=vSymbolsSelf.end();++itSymbol)
				{
					//#pragma omp critical
					//{
						pParser->symbols_[string("My_"+(*itSymbol))] = pInd->GetPhenotype(*itSymbol);
						pParser->symbols_[*itSymbol] = pInd->GetPhenotype(*itSymbol); //set both variables
					//}
				}

				bool bLive = true;
				//#pragma omp critical 
				//{
					bLive = (UniformGen.Next() <= pParser->Evaluate())? true : false;
				//}

				if (!bLive) {

					#pragma omp critical 
					{
						delete pInd; // uhoh, dead!!
						_mpMales[i] = NULL; //mark it
						//_mpMales.erase(_mpMales.begin() + i);
					}
					break;
				}
			}
		
			++itSelfSymbols;
		}


	}

	#pragma omp for
	//for (vector< Individual * >::size_type i=0; i< nFemale; i++ ) {
	for (signed long i=0; i< nFemale; i++ ) {
		Individual * pInd = _mpFemales[i];

		#ifdef _OPENMP
		int nCPU = omp_get_thread_num();
		#else
		int nCPU = 0;
		#endif

		list< pair< Parser *, int > > * pqParsers = &(*mppqParsers)[nCPU][this->_sPopName];

		list< vector<string> >::iterator itSelfSymbols = pqvSelfSymbols->begin();
		for (list< pair< Parser *, int> >::iterator itParser= pqParsers->begin(); itParser != pqParsers->end() ; ++itParser) {
			vector<string> vSymbolsSelf = *itSelfSymbols;
			Parser * pParser = itParser->first;
			int nGen = itParser->second;

			if ((nGen ==-1 && !bIgnoreGlobalRules) || (bIgnoreGlobalRules && nGen == nCurrGen) ) {

				for(vector<string>::iterator itSymbol=vSymbolsSelf.begin();itSymbol!=vSymbolsSelf.end();++itSymbol)
				{
					//#pragma omp critical 
					//{
						pParser->symbols_[string("My_"+(*itSymbol))] = pInd->GetPhenotype(*itSymbol);
						pParser->symbols_[*itSymbol] = pInd->GetPhenotype(*itSymbol); //set both variables
					//}
				}

				bool bLive = true;
				//#pragma omp critical 
				//{
					bLive = (UniformGen.Next() <= pParser->Evaluate())? true : false;
				//}

				if (!bLive) {
					#pragma omp critical 
					{
						delete pInd; // uhoh, dead!!
						_mpFemales[i] = NULL;
						//_mpFemales.erase(_mpFemales.begin() + i);
					}
					break;
				}
			}
		
			++itSelfSymbols;
		}


	}

	} //end parallel block

	//clear NULL pointers from the arrays:

	for (vector< Individual * >::size_type i=0; i< nMale; i++ ) {
		if (! _mpMales[i]) {
			_mpMales.erase(_mpMales.begin() + i);
			nMale--; // readjust upper bound
			i--;
		}
	}
	
	for (vector< Individual * >::size_type i=0; i< nFemale; i++ ) {
		if (! _mpFemales[i]) {
			_mpFemales.erase(_mpFemales.begin() + i);
			nFemale--; // readjust upper bound
			i--;
		}
	}

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

void Population::fnDumpNaturalProb(ofstream &fNaturalProbOutFile) {
	if (!_bRecordNatSelProb) {
		return;
	}
	fNaturalProbOutFile << "Rule\tIndID\tDadID\tDadPrevPop\tMomID\tMomPrevPop\tNatSelProb_Pop" << this->GetPopId() << "\tSurvived" << endl;
	for (vector< vector< double > >::iterator itInd=this->_mpOffSpringNaturalProb.begin(); itInd!=this->_mpOffSpringNaturalProb.end(); ++ itInd) {
		for(vector< double >::iterator itFields=itInd->begin(); itFields!=itInd->end(); ++ itFields) {
			fNaturalProbOutFile << *itFields << "\t";
		}
		fNaturalProbOutFile << endl;
	}
}

void Population::fnAddNatSelProb(double nParserCount, double nIndID, double nDadID, double nDadPop, double nMomID, double nMomPop, double nProb, bool bSurvived) {
	if (!_bRecordNatSelProb) {
		return;
	}
	#pragma omp critical 
	{
		vector< double > vRecord;
		vRecord.push_back(nParserCount);
		vRecord.push_back(nIndID);
		vRecord.push_back(nDadID);
		vRecord.push_back(nDadPop);
		vRecord.push_back(nMomID);
		vRecord.push_back(nMomPop);
		vRecord.push_back(nProb);
		vRecord.push_back(bSurvived? 1.0:0.0);
		_mpOffSpringNaturalProb.push_back(vRecord);
	}
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
