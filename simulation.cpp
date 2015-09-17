#include "simulation.h"
#include "maths.h"
#include <time.h>
#include <random>

#ifdef _OPENMP
 #include <omp.h>
#endif

int nTotalCPUCore = 1;


extern Normal NormalGen; 
extern Uniform UniformGen;
extern double nRandSeed;
std::default_random_engine RandEngine;

extern SimulationConfigurations SimulConfig; // defined in config.cpp

//random generators that are parallelized.
vector< Binomial * > vBinomGen;
vector< Normal * > vNomGen;
vector< Uniform * > vUniformGen;

void fnInitBinom(vector< Binomial * > &vBinom) {
	for(int nCpu=0; nCpu < nTotalCPUCore; nCpu++) {
		Binomial * pBinom = new Binomial(1,0.5);
		vBinom.push_back(pBinom );
	}
}

void fnInitNom(vector< Normal * > &vNom) {
	for(int nCpu=0; nCpu < nTotalCPUCore; nCpu++) {
		Normal * pNom = new Normal();
		vNom.push_back(pNom );
	}
}

void fnInitUniform(vector< Uniform * > &vNom) {
	for(int nCpu=0; nCpu < nTotalCPUCore; nCpu++) {
		Uniform * pNom = new Uniform();
		vNom.push_back(pNom );
	}
}



void UISimulation() {


		UILoadConfig();


		clock_t t;
		t = clock();
		PerformSimulation();
		t = clock() - t;
		printf ("Execution time %f seconds.\n",((float)t)/CLOCKS_PER_SEC);


};

string PopId2PopString(int nId) {
	ostringstream sRet;
	if (nId == 1 || nId == 2) {
		sRet << "pop" << nId;
	}
	else if (nId == 3) {
		sRet << "hybrid";
	}
	else {
		sRet << "hybrid" << (nId - 3);
	}

	return sRet.str();
}

void PerformSimulation() {
	int nGenerations = (int)SimulConfig.GetNumericConfig("generations");
	vector< Population * > vPops;

	double nRandSeed = SimulConfig.GetNumericConfig("RandomSeed");
	if (nRandSeed != -1 && nRandSeed>0.0 && nRandSeed<1.0) {
		Random::Set(nRandSeed); //initialize random number generator
		#ifndef SRANDCALLED		
		srand((int) (nRandSeed * (double)RAND_MAX));
		#define SRANDCALLED
		#endif

		printf("User-specified random seed is %f \n" , nRandSeed);
	}
	else {
		#ifndef SRANDCALLED
		srand(time(NULL)); // set seed with time.
		#define SRANDCALLED
		#endif
		nRandSeed = (double)rand()/(double)RAND_MAX;
		if (nRandSeed == 0.0) {
			nRandSeed = 0.0001;
		}

		if (nRandSeed == 1.0) {
			nRandSeed = 0.9999;
		}

		Random::Set(nRandSeed); //initialize random number generator
		//srand((int) (nRandSeed * (double)RAND_MAX));
		printf("System chosen random seed is %f \n" , nRandSeed);
		SimulConfig.SetNumericConfig("RandomSeed" , nRandSeed);
	}

	#ifdef _OPENMP
		nTotalCPUCore =  omp_get_max_threads();//omp_get_num_threads();
		//printf("OpenMP enabled\n");
	#else
		nTotalCPUCore = 1;
		//printf("OpenMP disabled \n");
	#endif

	fnInitBinom(vBinomGen);
	fnInitNom(vNomGen);
	fnInitUniform(vUniformGen);

	string sOutFolder = SimulConfig.GetConfig("OutputFolder");
	
	if (!fnFileExists(sOutFolder.c_str())) {
		fnMakeDir(sOutFolder.c_str()); // make output folder
	}
	else { //output folder already exists, try another name
		int nTryFolder=1;
		while(true) {
			string sCurrFolderNum = fnIntToString(nTryFolder);
			sOutFolder = SimulConfig.GetConfig("OutputFolder") + "_" + sCurrFolderNum;
				if (!fnFileExists(sOutFolder.c_str())) {
					fnMakeDir(sOutFolder.c_str()); // make output folder
					SimulConfig.SetConfig("OutputFolder" , sOutFolder);
					break;
				}

			nTryFolder++;
		}

	}

	string sOutFileBase = sOutFolder + "/" + "Gen";

	int nPopId = 1;//start from 1
	while(true) { // now try to parse in populations.

		string sPopString = PopId2PopString(nPopId);
		Population * pNewPop = new Population();
		printf((sPopString+"_size_limit\n").c_str());
		if (SimulConfig.GetNumericConfig(sPopString+"_size_limit")==-1) {
			printf("%d populations defined in configuration.\n", (nPopId-1));
			break;
		}

		if (nPopId==1 || nPopId==2) { // parental pops
				printf(string("Initializing Parental " + sPopString + "...\n").c_str());
				pNewPop->Init(
					SimulConfig.GetConfig(sPopString+"_name"),nPopId, 
					SimulConfig.GetConfig(sPopString+"_ancestry_label").c_str()[0],
					(int)SimulConfig.GetNumericConfig(sPopString+"_init_size"),
					(int)SimulConfig.GetNumericConfig(sPopString+"_size_limit"),
					SimulConfig.GetNumericConfig(sPopString+"_male_ratio")
				);
		}
		else {
				printf(string("Initializing Hybrid pop " + sPopString + "...\n").c_str());
				pNewPop->Init(
					SimulConfig.GetConfig(sPopString+"_name"),nPopId, 
					NULL,
					0,
					(int)SimulConfig.GetNumericConfig(sPopString+"_size_limit"),
					0
				);
		}

		vPops.push_back(pNewPop);

		nPopId++;
	}




	// Now start simulation:


	// Do migration:
	int nSampleFreq = (int)SimulConfig.GetNumericConfig("samplefreq");

	for (int nCurrGen=0; nCurrGen <= nGenerations; nCurrGen++) {

		clock_t t;
		t = clock();
		
		
		printf("\n----------------------------------------\nGeneration: %d\n", nCurrGen);
		SimulConfig.SetCurrGen(nCurrGen);

		//Summarize phenotypes
		printf("Calculating phenotype distributions...\n");
		for (vector< Population * >::iterator itpPop = vPops.begin() ; itpPop != vPops.end(); ++itpPop)
		{
			(*itpPop)->SummarizePhenotype();
		}
	
		bool bAllOutputOff = (SimulConfig.GetConfig("AllOutput")=="Off");

		if (bAllOutputOff) {
			printf("Warning: AllOutput is set to Off, no output will be written. \n");
		}

		if ( (!bAllOutputOff) && ((nSampleFreq > 0  && (nCurrGen % nSampleFreq == 0)) || nCurrGen==0 || nCurrGen==1 || nCurrGen==2 || SimulConfig.IsInUserSpecifiedSamplingRange(nCurrGen))) {

		
			printf("Writing to disk...\n");
			/*
			char szCurrGen[100];
			_itoa(nCurrGen , szCurrGen, 10);
			*/
			string szCurrGen = fnIntToString(nCurrGen);
			//FILE * pOutFile = fopen( (sOutFileBase + szCurrGen + ".txt").c_str(), "w" );
			ofstream fMarkerOutFile, fGeneOutFile, fPhenotypeOutFile, fPhenoSumFile, fOffSpringNatSelProb;
			fMarkerOutFile.open((sOutFileBase + szCurrGen + "_markers.txt").c_str());// output file stream
			fGeneOutFile.open((sOutFileBase + szCurrGen + "_genes.txt").c_str());// output file stream
			fPhenotypeOutFile.open((sOutFileBase + szCurrGen + "_phenotypes.txt").c_str());// output file stream
			fPhenoSumFile.open((sOutFileBase + szCurrGen + "_phenostats.txt").c_str());// output file stream
			fOffSpringNatSelProb.open((sOutFileBase + szCurrGen + "_natselprobdump.txt").c_str());// output file stream

			fMarkerOutFile.setf(std::ios_base::fixed);
			fGeneOutFile.setf(std::ios_base::fixed);
			fPhenotypeOutFile.setf(std::ios_base::fixed);
			fOffSpringNatSelProb.setf(std::ios_base::fixed);

			//dump everything to disk:
			for (vector< Population * >::iterator itpPop = vPops.begin() ; itpPop != vPops.end(); ++itpPop)
			{
				(*itpPop)->Sample(fMarkerOutFile, fGeneOutFile, fPhenotypeOutFile, fPhenoSumFile, fOffSpringNatSelProb);
			}
		
			fMarkerOutFile.close();
			fGeneOutFile.close();
			fPhenotypeOutFile.close();
		}
		// do frequency-dependent selection

		for (vector< Population * >::iterator itpPop = vPops.begin() ; itpPop != vPops.end(); ++itpPop)
		{
				printf(string("FreqDependent Selection " + (*itpPop)->GetPopName() + "...\n").c_str());
				(*itpPop)->FreqDependentNaturalSelection();
				
		}



		//do migration:
		//try to do all pairwise migration, if the parameter is not specified, then it's deemed 0
		bool bMigrateOnlyFirstGen = SimulConfig.GetConfig("migration_only_first_gen") == "yes";

		if (bMigrateOnlyFirstGen && nCurrGen!=0) {//no migration then.

		}
		else {
			string szCurrGen = fnIntToString(nCurrGen);
			string sMigrationParamLabelPrefix = "gen" + szCurrGen + "_"; //generation specific rule prefix
			for (vector< Population * >::iterator itpPop1 = vPops.begin() ; itpPop1 != vPops.end(); ++itpPop1)
			{
				for (vector< Population * >::iterator itpPop2 = vPops.begin() ; itpPop2 != vPops.end(); ++itpPop2) {

					if ((*itpPop1)->GetPopId() == (*itpPop2)->GetPopId()) {
						continue; // same pops, no migration needed
					}
					string sMigrationParamLabel = sMigrationParamLabelPrefix + PopId2PopString((*itpPop1)->GetPopId()) + "_to_" + PopId2PopString((*itpPop2)->GetPopId()); 
					string sMigrationParamLabelGeneral = PopId2PopString((*itpPop1)->GetPopId()) + "_to_" + PopId2PopString((*itpPop2)->GetPopId()); 
					if (SimulConfig.GetNumericConfig(sMigrationParamLabel)==-1 && SimulConfig.GetNumericConfig(sMigrationParamLabelGeneral)==-1) {
						continue; // didn't find the configuration for this migration, move on
					}
					else {
						double nNumberMigrants = SimulConfig.GetNumericConfig(sMigrationParamLabel);
						double nNumberMigrantsGeneral = SimulConfig.GetNumericConfig(sMigrationParamLabelGeneral);
						if (nNumberMigrants == -1 && nNumberMigrantsGeneral!=-1) { // IF generation specific rule not available, use general rule
							nNumberMigrants = nNumberMigrantsGeneral;
						}
						else if (nNumberMigrants == -1 && nNumberMigrantsGeneral==-1) {
							continue;
						}

						if (SimulConfig.GetConfig("UseBinomMigrateRate") == "yes" && nNumberMigrants!=0 && (*itpPop1)->GetPopSize(4)>0) {
							if (nNumberMigrants <=0 || nNumberMigrants > 1) {
								printf("ERROR: When UseBinomMigrateRate is turned on, migrational rates must be 0-1.\n");
								Throw(Out_of_range("ERROR: When UseBinomMigrateRate is turned on, migrational rates must be 0-1."));
								
							}
							//Binomial oBinom((*itpPop1)->GetPopSize(4) ,  nNumberMigrants);
							std::binomial_distribution<int> oBinomDist((*itpPop1)->GetPopSize(4) ,  nNumberMigrants);
							printf("Migrant number drawn from Binom(n=%d , p=%f)\n", (*itpPop1)->GetPopSize(4) , nNumberMigrants);
							nNumberMigrants = round( oBinomDist(RandEngine) );
							
						}
						printf("Migration %s : %d migrants...\n", sMigrationParamLabel.c_str(), (int)nNumberMigrants);
						for (int i=0; i<nNumberMigrants; i++) {
							if ((*itpPop1)->GetPopSize(4) == 0) {
								printf("Population size of %s is now 0, only migrated %d individuals.\n", (*itpPop1)->GetPopName().c_str(), i);
								break;							
							}
							(*itpPop2)->Immigrate( (*itpPop1)->Emigrate() );
						}

					}
				}
				
			}

			//confirm immigrants for every pop
			for (vector< Population * >::iterator itpPop = vPops.begin() ; itpPop != vPops.end(); ++itpPop)
			{
				(*itpPop)->ImmigrateConfirm(true);
			}
		}

		//breed and kill
		
		bool bAtLeastOneBred = false;

		for (vector< Population * >::iterator itpPop = vPops.begin() ; itpPop != vPops.end(); ++itpPop)
		{
				
				if ((*itpPop)->Breed())
				{
					
					bAtLeastOneBred = true;
				}

				(*itpPop)->KillOldGen();
				
		}


		if (!bAtLeastOneBred ) {
			printf("All populations went extinct because all of them failed to breed.");
			exit(200);
		}

		ostringstream sPopSizeReport;
		sPopSizeReport<< "Population size: ";
		for (vector< Population * >::iterator itpPop = vPops.begin() ; itpPop != vPops.end(); ++itpPop)
		{
			sPopSizeReport <<  (*itpPop)->GetPopName() << "  ";
			sPopSizeReport << ((*itpPop)->GetPopSize(4));
			sPopSizeReport << " (m:f=" ;
			sPopSizeReport << ((double)(*itpPop)->GetPopSize(1) / (double)(*itpPop)->GetPopSize(0));
			sPopSizeReport << ") ";
		}
		sPopSizeReport << "\n";
		 
		printf(sPopSizeReport.str().c_str());

		t = clock() - t;
		printf ("Execution time for generation : %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

	}
	
};

void UIExportResults() {
	UILoadConfig();
		
	int nGen = 0; // generation to sample:
	int nRep = 0; //which replicate to sample from?
	string sExportFolder;
	int nPop1Size, nPop2Size, nPop3Size , nMarkerNum; // how many individuals to sample
	double nRandSeed;

	while(true) { //
		
			std::cout << "\n Which generation to sample? :" ;
			std::cin >> nGen;
			
			if (nGen < 0 ) {
				std::cout << "\n Make sure the generation is correct. \n";
				continue; // 
			}
			else {
				break;
			}
			
			
			
		}

	while(true) { //
		
			std::cout << "\n Which replicate to sample? Starting from 0 :" ;
			std::cin >> nRep;
			
			if (nRep < 0 ) {
				std::cout << "\n Make sure the replicate number is correct. \n";
				continue; // 
			}
			else {
				break;
			}
			
			
			
		}

	while(true) { //
		
			std::cout << "\n Export to folder :" ;
			std::cin >> sExportFolder;
			
			if (sExportFolder.size() == 0 ) {
				std::cout << "\n Make sure the folder name is correct. \n";
				continue; // 
			}
			else {
				fnMakeDir(sExportFolder.c_str()); // make dir.
				break;
			}
			
			
			
		}

	/*
	while(true) { //
		
			std::cout << "\n How many individuals to sample for each pop (-1 to include all) in the order of pop1 pop2 pophybrid :" ;
			if(scanf("%d %d %d", &nPop1Size, &nPop2Size, &nPop3Size)!=3) {
				continue;
			}
			else {
				break;
			}
			
		}
		*/
	nPop1Size = -1;
	nPop2Size = -1;
	nPop3Size = -1;

	while(true) { //
		
			std::cout << "\n How many markers to sample? max: 10000 :" ;
			if(scanf("%d", &nMarkerNum)!=1) {
				continue;
			}
			else {
				break;
			}
			
		}

	while(true) { //
		
			std::cout << "\n Random Seed? 0-1 :" ;
			if(scanf("%lf", &nRandSeed)!=1) {
				continue;
			}
			else {
				break;
			}
			
		}


	PerformExport(nGen, nRep, sExportFolder, nPop1Size, nPop2Size, nPop3Size, nMarkerNum, nRandSeed);
}

void PerformExport(int nGen, int nRep, string sExportFolder, int nPop1Size, int nPop2Size, int nPop3Size, int nMarkerNum, double nRandSeed) {

	Random::Set(nRandSeed);
	string sOutFolder = SimulConfig.GetConfig("OutputFolder");
	sExportFolder = sExportFolder + "/Gen" + convertInt(nGen);
	fnMakeDir(sExportFolder.c_str());

	ofstream fPriorAlleleFreq, fOutcomeVars, fGenotypes, fEtaPriors, fLoci, fSampledMarkers;// Create priorallelefreqs.txt
	//fPriorAlleleFreq.open(sExportFolder + "/priorallelefreqs.txt");
	//fOutcomeVars.open(sExportFolder + "/outcomevars.txt");
	//fGenotypes.open(sExportFolder + "/genotypes.txt");
	//fEtaPriors.open(sExportFolder + "/etapriors.txt");
	fLoci.open((sExportFolder + "/lociChr.txt").c_str());
	fSampledMarkers.open("sampledmarkers.txt"); // this is an instruction to the php script regarding what markers to draw

	ifstream fMarkers, fPhenotypes;
	/*char szCurrGen[100];
	_itoa(nGen , szCurrGen, 10);
	*/
	string sGenHeader = sOutFolder;
	sGenHeader = (nRep==0)? sGenHeader : (sGenHeader + "_") + fnIntToString(nRep);
	sGenHeader = sGenHeader + "/Gen";
	sGenHeader = sGenHeader + fnIntToString(nGen);//szCurrGen;

	fMarkers.open((sGenHeader + "_markers.txt").c_str() , ios_base::in );

	fPhenotypes.open((sGenHeader + "_phenotypes.txt").c_str() , ios_base::in );
	//convert the phenotype file to outcomevars.txt
	string sCmd = "php export_phenotype.php -- -s \"" + sGenHeader + "_phenotypes.txt\" -o \"" + sExportFolder + "/outcomevars.txt\"";
	printf("Calling: %s ...\n", sCmd.c_str());
	system(sCmd.c_str());

	//write loci.txt file
	printf("Writing locus information...\n");
	fLoci << "\"LocusName\"\t\"NumberOfAlleles\"\t\"cM\"\t\"Chr\"" <<endl; // Header
	
	fLoci << fixed << showpoint; 
	fLoci << setprecision(6);
	

	//vector<vector<int>> vRandomMarkerSubset;
	vector<map<double, double> > * pvMarkerMapDistance = SimulConfig.pMarkerConfig->GetMpMarkerMapDistance(2); // use average map distance

	size_t nTotalChrSize = pvMarkerMapDistance->size();
	for( size_t nChr=0; nChr < nTotalChrSize ; nChr++) {
		map<double, double> * pmMarkerMapDistanceOnChr = & (pvMarkerMapDistance->at(nChr));
		size_t nCurrMarker = 1;
		size_t nSampledMarkerCount = 1;
		double nPrevMarkerPos=0.0;

		double nExpectedMarkerNumOnChr = (double)nMarkerNum * SimulConfig.pMarkerConfig->GetChrToGenomeRatio( nChr );
		double nProbToKeep = nExpectedMarkerNumOnChr / (double)pmMarkerMapDistanceOnChr->size();

		for (map<double,double>::iterator itMarker = pmMarkerMapDistanceOnChr->begin(); itMarker!=pmMarkerMapDistanceOnChr->end(); ++itMarker) {

			if (UniformGen.Next() <= nProbToKeep) {

				double nCurrPos = itMarker->second; // itMarker->first  * 10.0e-6;

				fLoci << "chr" << (nChr+1) << "_" << nCurrMarker << "\t2\t";
				fSampledMarkers << nChr << "\t" << ( nCurrMarker -1 ) << "\t" << nCurrPos << endl;
			
				if ((nSampledMarkerCount==1)) {
				
					fLoci << "NA\t"; 
				}
				else {
				
					fLoci <<  (nCurrPos - nPrevMarkerPos) << "\t";
					
				
				}
				nPrevMarkerPos = nCurrPos ;
				fLoci << (nChr+1) << endl;
				nSampledMarkerCount++;
			}
			else {
				//nPrevMarkerPos += itMarker->first; //accumulate distance
			}

			nCurrMarker++;
		}
	}

	fSampledMarkers.close();

	sCmd = "php export_genotype.php -- -s \"" + sGenHeader + "_markers.txt\" -o \"" + sExportFolder + "\" -i sampledmarkers.txt -pop1 100000 -pop2 100000 -pop3 999999";
	printf("Calling: %s ...\n", sCmd.c_str());
	system(sCmd.c_str());



	/*
	string sLn;
	int iLn = 0;
	if (std::getline(fPhenotypes, sLn)) {// Header

	}

	else {
		printf("Phenotype output error.\n");
		return;
	}

	while(std::getline(fPhenotypes, sLn)) {
		
		int nId, nPop, nSex, nFatherId, nMotherId;
		sscanf(sLn.c_str(), "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\n]", &nId, &nPop, &nSex, &nFatherId, &nMotherId);
		string sPrelude
		iLn++;
	}
	*/
}

void UILoadConfig() {
		string szFilename;
		
		if (SimulConfig.GetConfigFileName() == "") {

	
		while(true) { //Get random seed from user
		
			std::cout << "\n Enter simulation configuration file name: ";
			std::cin >> szFilename;
			
			if (szFilename.size() == 0 ) {
				std::cout << "\n Please enter file name! ";
				continue; // 
			}
			else {
				break;
			}
			
			
			
		}

		SimulConfig.LoadFromFile(szFilename);
		//printf("%s\n", szFilename.c_str());
		}
}

string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

void fnMakeDir(const char * sPath)
{
	#ifdef _WIN32
	_mkdir(sPath);
	#elif __linux__
	string szCmd("mkdir -m 777 -p ");
	system((szCmd + "\"" + sPath + "\"").c_str());
	#else
	string szCmd("mkdir -m 777 -p ");
	system((szCmd + "\"" + sPath + "\"").c_str());
	#endif
}

bool fnFileExists(const char *filename)
{
  ifstream ifile(filename);
  return ifile.is_open();
}