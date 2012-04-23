#include "simulation.h"
#include "maths.h"

extern Normal NormalGen; 
extern Uniform UniformGen;
extern double nRandSeed;

extern SimulationConfigurations SimulConfig; // defined in config.cpp

void UISimulation() {


		UILoadConfig();


		PerformSimulation();

};

void PerformSimulation() {
	int nGenerations = (int)SimulConfig.GetNumericConfig("generations");
	Population * pPop1 = new Population();
	Population * pPop2 = new Population();
	Population * pPop_hybrid = new Population();

	Random::Set(SimulConfig.GetNumericConfig("RandomSeed")); //initialize random number generator
	srand((int)SimulConfig.GetNumericConfig("RandomSeed") * 10);

	string sOutFolder = SimulConfig.GetConfig("OutputFolder") + "/";
	string sOutFileBase = sOutFolder + "Gen";
	fnMakeDir(sOutFolder.c_str()); // make output folder

	printf("Initializing Pop 1...\n");
	pPop1->Init(
					SimulConfig.GetConfig("pop1_name"),1, 
					SimulConfig.GetConfig("pop1_ancestry_label").c_str()[0],
					(int)SimulConfig.GetNumericConfig("pop1_init_size"),
					(int)SimulConfig.GetNumericConfig("pop1_size_limit"),
					SimulConfig.GetNumericConfig("pop1_male_ratio")
				);

	printf("Initializing Pop 2...\n");
	pPop2->Init(
					SimulConfig.GetConfig("pop2_name"),2, 
					SimulConfig.GetConfig("pop2_ancestry_label").c_str()[0],
					(int)SimulConfig.GetNumericConfig("pop2_init_size"),
					(int)SimulConfig.GetNumericConfig("pop2_size_limit"),
					SimulConfig.GetNumericConfig("pop2_male_ratio")
				);
	printf("Initializing Pop hybrid...\n");
	pPop_hybrid->Init(
					SimulConfig.GetConfig("hybrid_name"),3, 
					NULL,
					0,
					(int)SimulConfig.GetNumericConfig("hybrid_size_limit"),
					0
				);

	// Now start simulation:
	// Do migration:
	int nSampleFreq = (int)SimulConfig.GetNumericConfig("samplefreq");

	for (int nCurrGen=0; nCurrGen <= nGenerations; nCurrGen++) {

		printf("\n----------------------------------------\nGeneration: %d\n", nCurrGen);
		SimulConfig.SetCurrGen(nCurrGen);

		if ( (nSampleFreq > 0  && (nCurrGen % nSampleFreq == 0)) || nCurrGen==0 || nCurrGen==1 || nCurrGen==2) {

		
			printf("Writing to disk...\n");
			/*
			char szCurrGen[100];
			_itoa(nCurrGen , szCurrGen, 10);
			*/
			string szCurrGen = fnIntToString(nCurrGen);
			//FILE * pOutFile = fopen( (sOutFileBase + szCurrGen + ".txt").c_str(), "w" );
			ofstream fMarkerOutFile, fGeneOutFile, fPhenotypeOutFile;
			fMarkerOutFile.open((sOutFileBase + szCurrGen + "_markers.txt").c_str());// output file stream
			fGeneOutFile.open((sOutFileBase + szCurrGen + "_genes.txt").c_str());// output file stream
			fPhenotypeOutFile.open((sOutFileBase + szCurrGen + "_phenotypes.txt").c_str());// output file stream

			//dump everything to disk:
			pPop1->Sample(fMarkerOutFile, fGeneOutFile, fPhenotypeOutFile);
			pPop2->Sample(fMarkerOutFile, fGeneOutFile, fPhenotypeOutFile);
			pPop_hybrid->Sample(fMarkerOutFile, fGeneOutFile, fPhenotypeOutFile);
		
			fMarkerOutFile.close();
			fGeneOutFile.close();
			fPhenotypeOutFile.close();
		}

		//do migration:
		int nPop1_to_hybrid = (nCurrGen==0)? (int)SimulConfig.GetNumericConfig("gen1_pop1_to_hybrid"):(int)SimulConfig.GetNumericConfig("pop1_to_hybrid") ;
		int nPop2_to_hybrid = (nCurrGen==0)? (int)SimulConfig.GetNumericConfig("gen1_pop2_to_hybrid") : (int)SimulConfig.GetNumericConfig("pop2_to_hybrid");
		bool bMigrateOnlyFirstGen = SimulConfig.GetConfig("migration_only_first_gen") == "yes";

		
		if ( (bMigrateOnlyFirstGen && nCurrGen==0) || !bMigrateOnlyFirstGen ) {
			printf("Start Migration...\n");
			printf("Pop1_to_hybrid : %d...\n", nPop1_to_hybrid);
			for (int i=0; i<nPop1_to_hybrid; i++) {
				pPop_hybrid->Immigrate( pPop1->Emigrate(), true );
			}

			printf("Pop2_to_hybrid : %d...\n", nPop2_to_hybrid);
			for (int i=0; i<nPop2_to_hybrid; i++) {
				pPop_hybrid->Immigrate( pPop2->Emigrate(), true );
			}
			printf("End Migration...\n");
		}
		
		bool bAtLeastOneBred = false;
		if (pPop1->Breed())
		{
			pPop1->KillOldGen();
			bAtLeastOneBred = true;
		}
		if (pPop2->Breed())
		{
			pPop2->KillOldGen();
			bAtLeastOneBred = true;
		}
		if (pPop_hybrid->Breed()){

			pPop_hybrid->KillOldGen();
			bAtLeastOneBred = true;
		}

		if (!bAtLeastOneBred ) {
			printf("All populations went extinct because all of them failed to breed.");
			exit(200);
		}

		printf("Population size: pop1 %d (m:f=%3f) pop2 %d (m:f=%3f) pop3 %d (m:f=%3f) \n", 
			pPop1->GetPopSize(4),
			(double)pPop1->GetPopSize(1) / (double)pPop1->GetPopSize(0),
			pPop2->GetPopSize(4),
			(double)pPop2->GetPopSize(1) / (double)pPop2->GetPopSize(0),
			pPop_hybrid->GetPopSize(4),
			(double)pPop_hybrid->GetPopSize(1) / (double)pPop_hybrid->GetPopSize(0)
			);

	}
	
};

void UIExportResults() {
	UILoadConfig();
		
	int nGen = 0; // generation to sample:
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

	while(true) { //
		
			std::cout << "\n How many individuals to sample for each pop (-1 to include all) in the order of pop1 pop2 pophybrid :" ;
			if(scanf("%d %d %d", &nPop1Size, &nPop2Size, &nPop3Size)!=3) {
				continue;
			}
			else {
				break;
			}
			
		}

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


	PerformExport(nGen, sExportFolder, nPop1Size, nPop2Size, nPop3Size, nMarkerNum, nRandSeed);
}

void PerformExport(int nGen, string sExportFolder, int nPop1Size, int nPop2Size, int nPop3Size, int nMarkerNum, double nRandSeed) {

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
	string sGenHeader = sOutFolder + "/Gen";
	sGenHeader = sGenHeader + fnIntToString(nGen);//szCurrGen;

	fMarkers.open((sGenHeader + "_markers.txt").c_str() , ios_base::in );

	fPhenotypes.open((sGenHeader + "_phenotypes.txt").c_str() , ios_base::in );
	//convert the phenotype file to outcomevars.txt
	string sCmd = "php export_phenotype.php -- -s \"" + sGenHeader + "_phenotypes.txt\" -o \"" + sExportFolder + "/outcomevars.txt\"";
	printf("Calling: %s ...\n", sCmd.c_str());
	system(sCmd.c_str());

	//write loci.txt file
	printf("Writing locus information...\n");
	fLoci << "\"LocusName\"\t\"NumberOfAlleles\"\t\"Mb\"\t\"Chr\"" <<endl; // Header
	
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

				fLoci << "chr" << (nChr+1) << "_" << nCurrMarker << "\t2\t";
				fSampledMarkers << nChr << "\t" << ( nCurrMarker -1 ) << "\t" << itMarker->first << endl;
			
				if ((nSampledMarkerCount==1)) {
				
					fLoci << "NA\t"; 
				}
				else {
				
					fLoci <<  (itMarker->first - nPrevMarkerPos) * 10.0e-6 << "\t";
					
				
				}
				nPrevMarkerPos = itMarker->first ;
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
	#endif
}