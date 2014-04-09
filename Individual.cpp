/*
 *  Individual.cpp
 *  admixsimul
 *
 *  Created by Rosenthal Lab on 11/21/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "Individual.h"
#include "Population.h"

extern SimulationConfigurations SimulConfig;
extern Normal NormalGen; 
extern Uniform UniformGen;
extern Binomial BinomGen(1, 0.5);

int nCurrIndividualId = 1; //Global counter
/*
void Genome::BuildGenome(Chromosome * pPaternalGamete, Chromosome * pMaternalGamete) {

};
*/

Individual::Individual(void) {
	this->_bDead = false;
}

Individual::~Individual(void) {
	//delete[] _arrMarkers;
	//delete[] _pmGenes;
	//delete[] _arrGenes;
	//delete &_mpPhenotypes;
	//delete &_mpEnvPhenotypes;
}

Individual::Individual(void * pPop, char nAncestryLabel) { //Initializing a founder
	this->_bDead = false;
	int nPopId = ((Population *)pPop)->GetPopId();
	this->_pPrevPop = NULL;
	this->_pPop = pPop;
	if (nPopId == 3) {
		throw "You cannot directly create hybrid founders.";
	}

	this->_nFatherId = 0; //founder, father unknown.
	this->_nMotherId = 0; //founder, mother unknown
	this -> _nId = nCurrIndividualId;
	nCurrIndividualId++;

	this->_nHaploidChrNum = SimulConfig.pMarkerConfig->GetHaploidChromosomeNum();
	this->_nTotalChrNum = 2* this->_nHaploidChrNum ;

	SimulConfig.pPhenotypeConfig->InitKeys(&this->_mpPhenotypes);
	SimulConfig.pPhenotypeConfig->InitKeys(&this->_mpEnvPhenotypes);

	map<double, double> * pMarkerInfo = SimulConfig.pMarkerConfig->GetMarkersInfoByPop(nPopId);
	vector< map<double , GeneProperties> > * pGenes = SimulConfig.pGeneConfig->GetMpGenes();

	//this->_pmMarkers = new map<double, Marker>[_nTotalChrNum];
	//this->_arrMarkers = new vector<Marker>* [_nTotalChrNum];
	//this->_pmGenes   = new map<double, Gene>[_nTotalChrNum];
	//this->_arrGenes = new vector<Gene>* [_nTotalChrNum];

	//fill individual genome with markers
	for (int i=0; i<this->_nTotalChrNum; i++) {

		int nCurrChrNum;
		if (i % 2 ==0) {
			nCurrChrNum = i / 2;
		}
		else {
			nCurrChrNum = (i-1) / 2;
		}

		//map<double, Marker> oMarkerOnChr;
		vector<Marker> oMarkerOnChr;
	
		for (map<double, double>::iterator it=pMarkerInfo[nCurrChrNum].begin(); it!=pMarkerInfo[nCurrChrNum].end(); ++it)
		{
			

			double nFrequency = it->second;

			Marker oMarker;
			oMarker.Ancestry = nAncestryLabel;

			/*
			if (UniformGen.Next() <= nFrequency) {
				oMarker.Allele = 'A';
			}
			else {
				oMarker.Allele = 'a';
			}
			*/
			oMarker.Allele = (UniformGen.Next() <= nFrequency)? 'A':'a';
			/*
			char nAbundantAllele = (nPopId==1)? 'A':'a';
			char nRareAllele = (nPopId==1)? 'a':'A'; // make sure these are ancestry informative.

			oMarker.Allele = (UniformGen.Next() <= nFrequency )? nAbundantAllele : nRareAllele;
			*/
			//oMarkerOnChr.insert(oMarkerOnChr.end(),pair<double, Marker>(it->first , oMarker));
			oMarkerOnChr.push_back(oMarker);
				
		}

		//this->_pmMarkers[i] = oMarkerOnChr;
		this->_arrMarkers.push_back(oMarkerOnChr);

		//map<double, Gene> oGeneOnChr;
		vector<Gene> oGeneOnChr;

		for (map<double, GeneProperties>::iterator it2=pGenes->at(nCurrChrNum).begin(); it2!=pGenes->at(nCurrChrNum).end(); ++it2) // genes
		{
			double nFrequency = (nPopId==1)? it2->second.DominantFreqPop1 : it2->second.DominantFreqPop2;
			Gene oGene;
			oGene.Ancestry = nAncestryLabel;
			/*
				if (i==19) {
						printf("%s\n", (it2->second.AlleleMode == GeneProperties::Hemizygous && (i % 2)==1)? "true":"false" );
					}
					*/
			if (it2->second.AlleleMode == GeneProperties::Hemizygous && (i % 2)==1 ) { // if this gene is hemizygous, then only need to randomize for the first chromosome.
				//if ( this->_pmGenes[i-1][it->first].Allele == it->second.DominantLabel ) { // if the homologous chrom. already has the dominant allele
					oGene.Allele = it2->second.RecessiveLabel; // make the current one recessive.
					oGene.Value = it2->second.RecessiveValue;

				//}
			}
			/*
			else if (UniformGen.Next() <= nFrequency) {
				oGene.Allele = it2->second.DominantLabel;
			}
			else 
			{
				oGene.Allele = it2->second.RecessiveLabel;
			}
			*/
			/*
			else {
				oGene.Allele = (rand() % 100 <= nFrequency * 100)?  it2->second.DominantLabel:it2->second.RecessiveLabel;
			}
			*/
			else {
				oGene.Allele = (UniformGen.Next()<= nFrequency)?  it2->second.DominantLabel:it2->second.RecessiveLabel;
				oGene.Value = (oGene.Allele == it2->second.DominantLabel)?  it2->second.DominantValue:it2->second.RecessiveValue;
				/*if ( nCurrChrNum == 9 ) {
					printf("something wrong.\n");
				}*/
			}

			//oGeneOnChr[it2->first] = oGene;
			//oGeneOnChr.insert(oGeneOnChr.end(),pair<double, Gene>(it2->first , oGene));
			oGeneOnChr.push_back(oGene);
		}

		//this->_pmGenes[i] = oGeneOnChr;
		this->_arrGenes.push_back(oGeneOnChr);

	}

	this->fnDeterminePhenotypes();
	this->fnDetermineSex();
	this->fnDetermineNumGametes();
	/*
	vector<vector<Marker> > vTempMarkers;
	vector<vector<Gene> > vTempGenes;
	this->GetGamete(vTempMarkers, vTempGenes);
	//printf("hello");
	*/

}

void Individual::fnDetermineSex() 
{
	//this->_bSex = UniformGen.Next()>=0.5? Sex::Male:Sex::Female; // temp function.
	this->_bSex = this->_mpPhenotypes["Sex"]==0.0? Female:Male;
}

void Individual::fnDetermineNumGametes() {
	if (this->_bSex == Male) {
		this->_nAvailableGametes = INT_MAX;
		return;
	}

	//Deal with females
	int nExpected = SimulConfig.GetNumericConfig("avg_female_gamete");
	int nStdDev	= SimulConfig.GetNumericConfig("std_female_gamete");
	this->_nAvailableGametes = (double)nStdDev * NormalGen.Next() + (double)nExpected;
}

void Individual::fnDeterminePhenotypes() { // Calculate the phenotypic values from genotypes

	for(std::map<string, double >::iterator it=this->_mpPhenotypes.begin(); it!=this->_mpPhenotypes.end(); ++it) {

		//Go over all phenotypes
		string sPhenoName = it->first;
		Parser *pFormula = SimulConfig.pPhenotypeConfig->GetFormula(sPhenoName);
		vector< pair<int, double> > * pvSymbols = SimulConfig.pPhenotypeConfig->GetFormulaSymbols(sPhenoName);
		vector<string> * pvSymbolStrings = SimulConfig.pPhenotypeConfig->GetFormulaSymbolStrings(sPhenoName);

		vector<string>::iterator oSymbolString = pvSymbolStrings->begin();

		for (vector< pair<int, double> >::iterator oSymbol = pvSymbols->begin(); oSymbol != pvSymbols->end(); ++oSymbol) 
		{
			
			int nChr = oSymbol->first; // Chromosome num
			double nPos = oSymbol->second; // Position of gene
			map<double, GeneProperties>  oGeneList = SimulConfig.pGeneConfig->GetMpGenes()->at(nChr-1);
			map<double, int> oGeneIndexList = SimulConfig.pGeneConfig->GetMpGeneIndex()->at(nChr-1);
			//look up information of this gene
			int nIndex2 = nChr  * 2 - 1;
			int nIndex1 = nIndex2 - 1;
			Gene oAllele1 = this->_arrGenes.at(nIndex1).at(oGeneIndexList[nPos]); 
			Gene oAllele2 = this->_arrGenes.at(nIndex2).at(oGeneIndexList[nPos]);
			GeneProperties oGeneProp = oGeneList[nPos];

			double nVal=0.0;
			/*
			if (sPhenoName == "ComplexEpi2") {
				printf("break");
			}
			*/

			if (oGeneProp.AlleleMode==GeneProperties::Hemizygous || oGeneProp.AlleleMode==GeneProperties::Dominant) { // if this gene is hemizygous or dominant
				if (oAllele1.Allele  != oAllele2.Allele) { // if one of the alleles is dominant
					nVal = (oAllele1.Allele == oGeneProp.DominantLabel )? oAllele1.Value : oAllele2.Value;//oGeneProp.DominantValue;
				}
				/*else if (oAllele1.Allele == oGeneProp.DominantLabel && oAllele2.Allele == oGeneProp.DominantLabel){ // both dominant
					nVal = ( oAllele1.Value + oAllele2.Value)/2;//oGeneProp.RecessiveValue; take the mean. 
				}*/
				else {//both recessive or dominant
					nVal = ( oAllele1.Value + oAllele2.Value)/2;
				}


			}
			else { //this gene is additive
				/*nVal  = oAllele1.Allele == oGeneProp.DominantLabel? oGeneProp.DominantValue:oGeneProp.RecessiveValue;
				  nVal += oAllele2.Allele == oGeneProp.DominantLabel? oGeneProp.DominantValue:oGeneProp.RecessiveValue;*/
				nVal = ( oAllele1.Value + oAllele2.Value);
			}
			
			
			pFormula->symbols_[*oSymbolString] = nVal; // put value to symbol;
			pFormula->symbols_[(*oSymbolString)+"_a1"] = oAllele1.Value;
			pFormula->symbols_[(*oSymbolString)+"_a2"] = oAllele2.Value;
			oSymbolString++;
		}

		this->_mpPhenotypes[sPhenoName] = pFormula->Evaluate(); // re-evaluate.
	}

}


Individual::Individual(Individual * pFather, Individual * pMother, bool &bSuccess) {
	this->_bDead = false;
	this->_pPrevPop = NULL;
	this->_pPop = pMother->GetPop();

	if (pFather->GetSex() != Male || pMother -> GetSex() != Female) { // gay sex not implemented.
		bSuccess = false;
		//throw "Cannot create a new individual from same-sex parents!";
		return;
	}
	
	this -> _nId = nCurrIndividualId;
	nCurrIndividualId++;
	
	this->_nFatherId = pFather->GetID();
	this->_nMotherId = pMother->GetID();

	
	this->_nHaploidChrNum = SimulConfig.pMarkerConfig->GetHaploidChromosomeNum();
	this->_nTotalChrNum = 2* this->_nHaploidChrNum ;

	SimulConfig.pPhenotypeConfig->InitKeys(&this->_mpPhenotypes);
	SimulConfig.pPhenotypeConfig->InitKeys(&this->_mpEnvPhenotypes);
	
	this->_mpDadPhenotypes.insert(pFather->_mpPhenotypes.begin(), pFather->_mpPhenotypes.end()); // save parent's phenotypes
	this->_mpMomPhenotypes.insert(pMother->_mpPhenotypes.begin(), pMother->_mpPhenotypes.end());
	
	vector< vector<Marker> > vFatherMarkers, vMotherMarkers;
	vector< vector<Gene> > vFatherGenes, vMotherGenes;
	#pragma omp critical 
	{
		pFather->GetGamete(vFatherMarkers, vFatherGenes);
		pMother->GetGamete(vMotherMarkers, vMotherGenes);
	}

	//printf("vFatherMarkers.size() %d\n", vFatherMarkers.size());
	//printf("vFatherGenes.size() %d\n", vFatherGenes.size());
	this->_arrMarkers.reserve(vFatherMarkers.size() * 2);
	this->_arrGenes.reserve(vFatherGenes.size() * 2);

	//put chromosomes in order:
	for (int i=0;i<vFatherMarkers.size();i++) {
		this->_arrMarkers.push_back(vFatherMarkers.at(i));
		this->_arrMarkers.push_back(vMotherMarkers.at(i));

	}

	for (int i=0;i<vFatherGenes.size();i++) {
		this->_arrGenes.push_back(vFatherGenes.at(i));
		this->_arrGenes.push_back(vMotherGenes.at(i));
	}

	//printf("_arrMarkers.size() %d\n", _arrMarkers.size());
	//printf("_arrGenes.size() %d\n", _arrGenes.size());
	this->fnDeterminePhenotypes();
	this->fnDetermineSex();
	this->fnDetermineNumGametes();
	bSuccess = true;
};

bool Individual::Court(Individual * pChooser) {
	if (this->_bSex != Male) {
		throw "Females courting is not implemented yet!";
	}

	int nCurrGen = SimulConfig.GetCurrGen(); // tell what is the current generation
	bool bIgnoreGlobalRules = SimulConfig.pSexualSelConfig->IgnoreGlobalRules(nCurrGen); // is there special sexual selection rules for this generation?

	return pChooser->HandleCourter(this , bIgnoreGlobalRules);
}

int Individual::HandleCourter(Individual * pCourter , bool bIgnoreGlobalRules) {

	int nCurrGen = SimulConfig.GetCurrGen();

	if (pCourter->IsDead()) {
		throw(new Exception("Cannot mate with dead courter. Something went wrong."));
	}

	if (this->_nAvailableGametes == 0) {
		//no more gametes, byebye.
		return  0;
	}
	//entry point for sexual selection.

	string sPop = ((Population*)(this->_pPop))->GetPopName();
	
	list< vector<string> > * pqvCourterSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsCourter( sPop);
	list< vector<string> > * pqvSelfSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsSelf( sPop);
	list< vector<string> > * pqvDadSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsDad( sPop);
	list< vector<string> > * pqvMomSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsMom( sPop);

	//Traits for oblique imprinting; note that things for peer imprinting has already been set at the pop level. see Population::Breed()
	list< vector<string> > * pqvPrevGenCurrPopSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsPrevGenCurrPop(sPop);
	list< vector<string> > * pqvPrevGenCurrPopCourterSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsPrevGenCurrPopCourter(sPop);
	list< vector<string> > * pqvPrevGenCurrPopChooserSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsPrevGenCurrPopChooser(sPop);

	list< vector<string> > * pqvPrevGenPrevPopSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsPrevGenPrevPop(sPop);
	list< vector<string> > * pqvPrevGenPrevPopCourterSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsPrevGenPrevPopCourter(sPop);
	list< vector<string> > * pqvPrevGenPrevPopChooserSymbols = SimulConfig.pSexualSelConfig->GetFormulaSymbolStringsPrevGenPrevPopChooser(sPop);

	list< vector<string> >::iterator itSelfSymbols= pqvSelfSymbols->begin();
	list< vector<string> >::iterator itCourterSymbols= pqvCourterSymbols->begin();
	list< vector<string> >::iterator itDadSymbols= pqvDadSymbols->begin();
	list< vector<string> >::iterator itMomSymbols= pqvMomSymbols->begin();

	list< vector<string> >::iterator itPrevGenCurrPopSymbols		=        pqvPrevGenCurrPopSymbols->begin();
	list< vector<string> >::iterator itPrevGenCurrPopCourterSymbols = 		 pqvPrevGenCurrPopCourterSymbols->begin();
	list< vector<string> >::iterator itPrevGenCurrPopChooserSymbols =		 pqvPrevGenCurrPopChooserSymbols->begin();
	
	list< vector<string> >::iterator itPrevGenPrevPopSymbols		= 				 pqvPrevGenPrevPopSymbols->begin();
	list< vector<string> >::iterator itPrevGenPrevPopCourterSymbols =  		 pqvPrevGenPrevPopCourterSymbols->begin();
	list< vector<string> >::iterator itPrevGenPrevPopChooserSymbols =		 pqvPrevGenPrevPopChooserSymbols->begin();

	bool bAccept = true;

	//#pragma omp critical 
	//{
		list< pair< Parser *, int > > * pqParsers = SimulConfig.pSexualSelConfig->GetFormulae(sPop);
		for (list< pair< Parser *, int > > ::iterator itParser= pqParsers->begin(); itParser != pqParsers->end() ; ++itParser) {
			vector<string> vSymbolsSelf = *itSelfSymbols;
			vector<string> vSymbolsCourter = *itCourterSymbols;
			vector<string> vSymbolsDad = *itDadSymbols;
			vector<string> vSymbolsMom = *itMomSymbols;

			vector<string> vSymbolsPrevGenCurrPop		=	*itPrevGenCurrPopSymbols	;	
			vector<string> vSymbolsPrevGenCurrPopCourter = *itPrevGenCurrPopCourterSymbols;
			vector<string> vSymbolsPrevGenCurrPopChooser = *itPrevGenCurrPopChooserSymbols;
			 
			vector<string> vSymbolsPrevGenPrevPop		 =	*itPrevGenPrevPopSymbols		;
			vector<string> vSymbolsPrevGenPrevPopCourter = *itPrevGenPrevPopCourterSymbols;
			vector<string> vSymbolsPrevGenPrevPopChooser = *itPrevGenPrevPopChooserSymbols;

			Parser * pParser = itParser->first;
			int nGen = itParser->second;

			if ((nGen ==-1 && !bIgnoreGlobalRules) || (bIgnoreGlobalRules && nGen == nCurrGen) ) { // two scenarios: apply global rule or apply generation-specific rule
				

				//Set self symbol values
				for(vector<string>::iterator itSymbol=vSymbolsSelf.begin();itSymbol!=vSymbolsSelf.end();++itSymbol)
				{
					pParser->symbols_[string("My_"+(*itSymbol))] = this->GetPhenotype(*itSymbol);
					pParser->symbols_[*itSymbol] = this->GetPhenotype(*itSymbol); //set both variables
				}

				//Set courter symbol values
				for(vector<string>::iterator itSymbol=vSymbolsCourter.begin();itSymbol!=vSymbolsCourter.end();++itSymbol)
				{
					pParser->symbols_[string("Courter_"+(*itSymbol))] = pCourter->GetPhenotype(*itSymbol);				
				}

				//Set dad symbol values
				bool bSkipRule = false; // rules containing references to parental phenotypes will be skipped in generation 0
				for(vector<string>::iterator itSymbol=vSymbolsDad.begin();itSymbol!=vSymbolsDad.end();++itSymbol)
				{
					if (nCurrGen == 0) {
						bSkipRule = true;
						break;
					}

					pParser->symbols_[string("Dad_"+(*itSymbol))] = this->_mpDadPhenotypes[*itSymbol];				
				}

				//Set mom symbol values
				for(vector<string>::iterator itSymbol=vSymbolsMom.begin();itSymbol!=vSymbolsMom.end();++itSymbol)
				{
					if (nCurrGen == 0) {
						bSkipRule = true;
						break;
					}

					pParser->symbols_[string("Mom_"+(*itSymbol))] = this->_mpMomPhenotypes[*itSymbol];				
				}

				//Set PrevGenCurrPop symbol values
				bSkipRule = bSkipRule || this->fnSetParserPopWide(pParser , vSymbolsPrevGenCurrPop, string("PrevGenCurrPop"), ((Population*)this->GetPop())->_mpPrevGenSumPhenotype);
				//Set PrevGenCurrPopCourter symbol values
				bSkipRule = bSkipRule || this->fnSetParserPopWide(pParser , vSymbolsPrevGenCurrPopCourter , string("PrevGenCurrPopCourter"), ((Population*)this->GetPop())->_mpPrevGenSumPhenotypeMale);
				//Set PrevGenCurrPopChooser
				bSkipRule = bSkipRule || this->fnSetParserPopWide(pParser , vSymbolsPrevGenCurrPopChooser, string("PrevGenCurrPopChooser"), ((Population*)this->GetPop())->_mpPrevGenSumPhenotypeFemale);
				
				//set PrevGenPrevPop
				bSkipRule = bSkipRule || this->fnSetParserPopWide(pParser , vSymbolsPrevGenPrevPop, string("PrevGenPrevPop"), ((Population*)this->GetPrevPop())->_mpPrevGenSumPhenotype);
				// set PrevGenPrevPopCourter
				bSkipRule = bSkipRule || this->fnSetParserPopWide(pParser , vSymbolsPrevGenPrevPopCourter, string("PrevGenPrevPopCourter"), ((Population*)this->GetPrevPop())->_mpPrevGenSumPhenotypeMale);
				//set PrevGenPrevPopChooser
				bSkipRule = bSkipRule || this->fnSetParserPopWide(pParser , vSymbolsPrevGenPrevPopChooser, string("PrevGenPrevPopChooser"), ((Population*)this->GetPrevPop())->_mpPrevGenSumPhenotypeFemale);
				


				if (bSkipRule) continue;
				

				bAccept = (UniformGen.Next() <= pParser->Evaluate())? true : false;

				if (!bAccept) {
					return 0; //reject mate
				}

			}

			/*
			if (!bAccept) {
				break;
			}
			*/
			++itSelfSymbols;
			++itCourterSymbols;
			++itDadSymbols;
			++itMomSymbols;

			++itPrevGenCurrPopSymbols;	
			++itPrevGenCurrPopCourterSymbols;
			++itPrevGenCurrPopChooserSymbols;
			
			++itPrevGenPrevPopSymbols;
			++itPrevGenPrevPopCourterSymbols;
			++itPrevGenPrevPopChooserSymbols;

		}
	
	//}//end critical

	if (!bAccept) {
		return 0;
	}

	int nSuccess;
	#pragma omp critical 
	{
		if (pCourter) {
			this->_arrOtherParentsForOffsprings.push_back(pCourter);
			nSuccess = 1;
		}
		else {
			nSuccess = 0;
		}
	}//end critical

	return nSuccess;
	
}

bool Individual::fnSetParserPopWide(Parser * pParser,  vector<string> &vSymbols, const string &sPrefix,  map< string , pair< double, double> > &mpSumPhenotypes) { // returns skipped.
	int nCurrGen = SimulConfig.GetCurrGen();
	for(vector<string>::iterator itSymbol=vSymbols.begin();itSymbol!=vSymbols.end();++itSymbol)
				{
					if (nCurrGen == 0) {
						return true; // skipped this rule
					}

					string sSymbol = (*itSymbol).substr(4);
					string sType = (*itSymbol).substr(0, 3); //either Avg or Std

					if (mpSumPhenotypes.find(sSymbol) == mpSumPhenotypes.end()) {
						throw new Exception("Unable to find population symbol");
					}

					if (sType == "Avg") {
						pParser->symbols_[string(sPrefix + "_Avg_"+sSymbol)] = mpSumPhenotypes[sSymbol].first;
					}
					else if (sType == "Std") {
						pParser->symbols_[string(sPrefix + "_Std_"+sSymbol)] = mpSumPhenotypes[sSymbol].second;
					}

	}

	return false;
}

int Individual::GetMateNumber() {
	return this->_arrOtherParentsForOffsprings.size();
}

int Individual::GetGameteNum() {
	return this->_nAvailableGametes;
}

void Individual::GiveBirth(vector<Individual *> &vOffSprings, int nNum, bool bIgnoreGlobalRules) {

	if (this->_nAvailableGametes <= 0) {
		//printf("Not enough gametes left.\n");
		return; //no more offsprings.
	}

	int nCurrGen = SimulConfig.GetCurrGen();
	//Go through the parenthood of each inseminated gamete:
	if (this->_bSex == Male) {
		throw "Males cannot give birth, r u nuts??";
	}

	nNum = nNum==-1? this->_arrOtherParentsForOffsprings.size():nNum;
	nNum = (nNum <= this->_nAvailableGametes)? nNum : this->_nAvailableGametes;
	//printf("Female %d: %d offsprings\n", this->GetID(), nNum);
	// do frequency-independent natural selection here!

	string sPop = ((Population*)(this->_pPop))->GetPopName();

	list< pair< Parser *, int> > * pqParsers;
	//#pragma omp critical 
	//{
		pqParsers = SimulConfig.pNaturalSelConfig->GetFormulae(sPop);
	

	//list< vector<string> > * pqvCourterSymbols = SimulConfig.pNaturalSelConfig->GetFormulaSymbolStringsCourter( sPop);
	list< vector<string> > * pqvSelfSymbols = SimulConfig.pNaturalSelConfig->GetFormulaSymbolStringsSelf( sPop);

	
	for(int i=0; i<nNum;i++) { 
		if (this->_arrOtherParentsForOffsprings.size()==0) return;

		
		Individual * pOffSpring = NULL;
		Individual * pDad = NULL;

		do {
			#pragma omp critical 
			{
				int nRandDad = fnGetRandIndex(this->_arrOtherParentsForOffsprings.size() );
				pDad = this->_arrOtherParentsForOffsprings.at(nRandDad);
			}
		}
		while(pDad==NULL);

		bool bSuccess = false;
		pOffSpring = new Individual( pDad, this , bSuccess); // create a new kid
		if (!bSuccess) {
			printf("Something wrong when creating new individual. return false\n");
			i--;
			continue;
		}

		if (!pOffSpring) {
			printf("Something wrong when creating new individual. return NULL\n");
			i--;
			continue;
		}

		this->_nAvailableGametes--; // one less gamete!

		// See if it's lucky enough to survive the cruel nature!!
		// Go through each selection rule:
		
		//list< vector<string> >::iterator itCourterSymbols= pqvCourterSymbols->begin();
		list< vector<string> >::iterator itSelfSymbols= pqvSelfSymbols->begin();
		int nParserCount = -1;
		for (list< pair< Parser *, int> >::iterator itParser= pqParsers->begin(); itParser != pqParsers->end() ; ++itParser) {
			nParserCount++;
			vector<string> vSymbolsSelf = *itSelfSymbols;
			Parser * pParser = itParser->first;
			int nGen = itParser->second;

			if ((nGen ==-1 && !bIgnoreGlobalRules) || (bIgnoreGlobalRules && nGen == nCurrGen) ) {
				//Set self symbol values
				bool bLive = false;

				for(vector<string>::iterator itSymbol=vSymbolsSelf.begin();itSymbol!=vSymbolsSelf.end();++itSymbol)
				{
					pParser->symbols_[string("My_"+(*itSymbol))] = pOffSpring->GetPhenotype(*itSymbol);
					pParser->symbols_[*itSymbol] = pOffSpring->GetPhenotype(*itSymbol); //set both variables
				}

				double nSurvivalProb = pParser->Evaluate();
				bLive = (UniformGen.Next() <= nSurvivalProb)? true : false;

				((Population*)this->_pPop)->fnAddNatSelProb(nParserCount, pOffSpring->GetID(), pDad->GetID(), ((Population*)pDad->GetPrevPop())->GetPopId(),
															this->GetID(), ((Population*)this->GetPrevPop())->GetPopId(),
															nSurvivalProb , bLive);

				/*if (SimulConfig.GetConfig("DumpNaturalSelProb") == "On") {
					printf("%f\n",nSurvivalProb);
				}*/

				if (!bLive) {
					delete pOffSpring; // uhoh, dead!!
					pOffSpring = NULL;
					break;
				}
			}

			++itSelfSymbols;
		}
		
		

		if (!pOffSpring) {
			
		}
		else {

				vOffSprings.push_back(pOffSpring);
				// just for now, inseminate one egg:			

		}

		
		//_arrOtherParentsForOffsprings.erase(_arrOtherParentsForOffsprings.begin() + nRandDad);
	}

	//} //end critical

	// now clear all gametes
	//this->_arrOtherParentsForOffsprings.clear();

}

double Individual::GetPhenotype(string sPhenotype) {
	return this->_mpPhenotypes[sPhenotype];
}

/*
void Individual::GetGamete(Chromosome * pGamete) {

	
}
*/
void Individual::ChangePopulation(void * pPop) {
	this->_pPrevPop = this->_pPop;
	this->_pPop = pPop;
}

void * Individual::GetPop() {
	return this->_pPop;
}

void * Individual::GetPrevPop() {
	return (this->_pPrevPop == NULL)? this->_pPop:this->_pPrevPop;
}

Individual::Sex Individual::GetSex() {
	return this-> _bSex;
};

unsigned int Individual::GetID() {
	return this-> _nId;
}

unsigned int Individual::GetFatherId() {
	return this->_nFatherId;
}

unsigned int Individual::GetMotherId() {
	return this->_nMotherId;
}

void Individual::GetGamete(vector< vector<Marker> > &vMarkers, vector< vector<Gene> > &vGenes ) {

	//do recombination. 

	bool bSex = (this->_bSex == Male);
	vector< map<double, int> > * pMarkerIndex = SimulConfig.pMarkerConfig->GetMpMarkerIndex();
	vector< map<double, int> > * pGeneIndex = SimulConfig.pGeneConfig->GetMpGeneIndex();

	for (int nChr=0;nChr<this->_nHaploidChrNum;nChr++) { // go over each chromosome
		

		int nMainChr, nOtherChr;
		//int nWhichToPick = (UniformGen.Next() <=0.5)? 0:1;
		int nWhichToPick = (int)BinomGen.Next();
		nMainChr = nChr*2 + nWhichToPick;
		nOtherChr = nWhichToPick==0? nMainChr+1:nMainChr-1;
		map<double, int> * pMarkerIndexOnChr = & (pMarkerIndex->at(nChr));
		map<double, int> * pGeneIndexOnChr = & (pGeneIndex->at(nChr));
		double nCentromerePos = SimulConfig.pMarkerConfig->GetCentromerePosition(nChr);
		double nChrLen = SimulConfig.pMarkerConfig->GetChromosomeLength(nChr);

		vector<Marker> vNewMarkers;
		vNewMarkers.reserve(pMarkerIndexOnChr->size());
		vector<Gene> vNewGenes;
		vNewGenes.reserve(pGeneIndexOnChr->size());

					
		size_t nLastIndex = pMarkerIndexOnChr->size()-1;
		map<double,int>::iterator oMarkerCentromereIt = pMarkerIndexOnChr->upper_bound(nCentromerePos);
		//nLastIndex =  pMarkerIndexOnChr->size()-1;

		size_t nCentromereIndex = (oMarkerCentromereIt == pMarkerIndexOnChr->end())? nLastIndex + 1 :oMarkerCentromereIt->second ;
			
		size_t nCentromereIndex_gene=-1, nLastIndex_gene=-1;

		//if (pGeneIndexOnChr->size() != 0)
		//{
			map<double,int>::iterator oGeneCentromereIt = pGeneIndexOnChr->upper_bound(nCentromerePos);
			nLastIndex_gene =  pGeneIndexOnChr->size()-1;
			nCentromereIndex_gene = (oGeneCentromereIt == pGeneIndexOnChr->end())? nLastIndex_gene + 1 :oGeneCentromereIt->second ;
				
		//}

		for (int nArm=1;nArm<=2;nArm++) 
		{


			vector<double> oBreakPoints;
			SimulConfig.pRecombProbConfig->GetBreakPointsByArm(bSex, nChr, nArm, oBreakPoints);
			
			/*
			if (pBreakPoints->size() > 0){
				printf("%d break points, first: %f\n",pBreakPoints->size(), *(pBreakPoints->begin()));
			}
			else {
				printf("No break point on Chr %d arm %d\n", nChr, nArm);
			}
			*/
			int nBreakPoints = oBreakPoints.size();

			bool bSelf = (nArm==1)? (nBreakPoints % 2 ==0) : true; //if it's arm 1, determine which chromosome to start
			int nBpCount = 0;




			int nStartIndex = (nArm == 1)? ((nCentromereIndex!=0)?  0 : INT_MAX ) : nCentromereIndex; // start copying from index 0
			int nEndIndex = 0;

			int nStartIndex_gene = (nArm == 1)? ((nCentromereIndex_gene!=0)?  0 : INT_MAX ): nCentromereIndex_gene; // start copying from index 0
			int nEndIndex_gene = 0;
			/*
			if (nChr == 5 ) {
				printf("break");
			}
			*/
			if (nBreakPoints == 0) { // if no break point, copy the whole damn arm.
				nEndIndex = (nArm == 1)? nCentromereIndex - 1 : nLastIndex;
				nEndIndex_gene = (nArm == 1)? nCentromereIndex_gene - 1 : nLastIndex_gene ;

				vector<Marker> * pvSourceMarkers = &this->_arrMarkers[nMainChr];
				vector<Gene> * pvSourceGenes = &this->_arrGenes[nMainChr];
				if (nStartIndex <= nEndIndex) {
					copy(pvSourceMarkers->begin()+nStartIndex,pvSourceMarkers->begin()+nEndIndex+1,back_inserter(vNewMarkers));
				}
				if (nStartIndex_gene <= nEndIndex_gene) {
					copy(pvSourceGenes->begin()+nStartIndex_gene,pvSourceGenes->begin()+nEndIndex_gene+1,back_inserter(vNewGenes));
				}
				continue; // go to next arm
			}

			bool bNoNeedContinueMarker = false;
			bool bNoNeedContinueGene = false;
			

			for (vector<double>::iterator nBreakPos=oBreakPoints.begin(); nBpCount <= nBreakPoints;) 
			{//Go over each break point
				//find index using position:
				int nChrIndex = bSelf? nMainChr:nOtherChr;
				vector<Marker> * pvSourceMarkers = &this->_arrMarkers[nChrIndex];
				vector<Gene> * pvSourceGenes = &this->_arrGenes[nChrIndex];

				//printf("max size %d\n",  pvSourceMarkers->max_size());

				if (nBpCount == nBreakPoints) // if it's beyond the last break point
				{
					nEndIndex = (nArm == 1)? nCentromereIndex - 1 : nLastIndex;
					nEndIndex_gene = (nArm == 1)? nCentromereIndex_gene - 1 : nLastIndex_gene;
				}
				else {

					map<double, int>::iterator oMarkerBreakPosIt = pMarkerIndexOnChr->upper_bound(*nBreakPos);
					map<double, int>::iterator oGeneBreakPosIt = pGeneIndexOnChr->upper_bound( *nBreakPos);

					if (oMarkerBreakPosIt == pMarkerIndexOnChr->end()) { // there is no more markers beyond the break point
						
						nEndIndex = (nArm == 1)? nCentromereIndex - 1 : nLastIndex;
						if (nStartIndex <= nEndIndex && !bNoNeedContinueMarker) {
							copy(pvSourceMarkers->begin()+nStartIndex,pvSourceMarkers->begin()+nEndIndex+1,back_inserter(vNewMarkers));
						}
						bNoNeedContinueMarker = true;
					}
					else {
						if (!bNoNeedContinueMarker) {
							nEndIndex = (oMarkerBreakPosIt->second < nCentromereIndex)? oMarkerBreakPosIt->second - 1: nCentromereIndex - 1; // need to minus 1 so that all markers on left but not right of Break point are copied
						}
					}

					if (oGeneBreakPosIt == pGeneIndexOnChr->end()) {
						nEndIndex_gene = (nArm == 1)? nCentromereIndex_gene - 1 : nLastIndex_gene;
						if (nStartIndex_gene <= nEndIndex_gene && !bNoNeedContinueGene ) {
							copy(pvSourceGenes->begin()+nStartIndex_gene,pvSourceGenes->begin()+nEndIndex_gene+1,back_inserter(vNewGenes));
						}
						bNoNeedContinueGene = true;
						
					}
					else {
						if (!bNoNeedContinueGene) {
							nEndIndex_gene = ( oGeneBreakPosIt->second < nCentromereIndex_gene)? oGeneBreakPosIt->second - 1: nCentromereIndex_gene - 1 ;
						}
					}

				}



				if (nStartIndex <= nEndIndex && !bNoNeedContinueMarker) {
					copy(pvSourceMarkers->begin()+nStartIndex,pvSourceMarkers->begin()+nEndIndex+1,back_inserter(vNewMarkers));
				}
				if (nStartIndex_gene <= nEndIndex_gene  && !bNoNeedContinueGene) {
					copy(pvSourceGenes->begin()+nStartIndex_gene,pvSourceGenes->begin()+nEndIndex_gene+1,back_inserter(vNewGenes));
				}

				bSelf = !bSelf; // flip
				nStartIndex = nEndIndex + 1;
				nStartIndex_gene = nEndIndex_gene + 1;
				nBpCount++;
				if (nBreakPos != oBreakPoints.end()) {
					++nBreakPos;
				}

	
			}

			
		}

		vMarkers.push_back(vNewMarkers);
		vGenes.push_back(vNewGenes);



	}

	//Perform mutation of genes:
	vector< map<double , GeneProperties> > * pGenes = SimulConfig.pGeneConfig->GetMpGenes();

	for (int nChr=0;nChr<this->_nHaploidChrNum;nChr++) {
		vector< Gene > * pChrmToMutate = &vGenes.at(nChr);
		int nInx = 0;
		for (map<double, GeneProperties>::iterator it2=pGenes->at(nChr).begin(); it2!=pGenes->at(nChr).end(); ++it2) {// genes
			//see if apply mutation rule:
			if ( (
					(*it2).second.Pop == "!ANYPOP"
						||
					(*it2).second.Pop == ((Population*)this->GetPop())->GetPopName()
				  )
				   &&
				  (
						(*it2).second.StartGen==-1
							||
						(
							SimulConfig.GetCurrGen() >= (*it2).second.StartGen
							&&
							SimulConfig.GetCurrGen() <= (*it2).second.EndGen
						)

				  )
				
				)
			{
			
				double nMuteProb = (*it2).second.MutationProb;
				if (nMuteProb==0) {//no mutation needed for this locus

				}
				else if( UniformGen.Next() <= nMuteProb){ // if need to mutate this locus
					#pragma omp critical
					{
						double nCurrVal = pChrmToMutate->at(nInx).Value;
						(*it2).second.pFormula->symbols_["CurrVal"] = nCurrVal; //give value to symbol;
						double nNewVal = (*it2).second.pFormula->Evaluate();
						pChrmToMutate->at(nInx).Value = (nNewVal > (*it2).second.UpperBound)? (*it2).second.UpperBound : ((nNewVal < (*it2).second.LowerBound)? (*it2).second.LowerBound : nNewVal);
					}
				}
			}
			nInx++;
		}
	}
	

	
}

void Individual::DumpMarkers(ofstream &fOutFile, int nChromosomeSide) { // nChromosomeSide is either 0 or 1
	for (int i=nChromosomeSide; i< this->_arrMarkers.size(); i+=2) {
		vector<Marker> * pMarkers = &(this->_arrMarkers.at(i));
		fOutFile << -1 << '\t'; // output chromosome spacer

		for (vector<Marker>::iterator itMarker=pMarkers->begin(); itMarker != pMarkers->end() ; ++itMarker) {
			fOutFile << (*itMarker).Allele <<   (*itMarker).Ancestry << '\t';
		}

	}
}

void Individual::DumpGenes(ofstream &fOutFile, int nChromosomeSide) { // nChromosomeSide is either 0 or 1
	for (int i=nChromosomeSide; i< this->_arrGenes.size(); i+=2) {
		vector<Gene> * pGenes = &(this->_arrGenes.at(i));
		fOutFile << -1 << '\t'; // output chromosome spacer

		for (vector<Gene>::iterator itGene=pGenes->begin(); itGene != pGenes->end() ; ++itGene) {
			fOutFile << (*itGene).Allele <<  (*itGene).Ancestry << '\t' << (*itGene).Value << '\t';
		}

	}
}

void Individual::DumpPhenotypes(ofstream &fOutFile) { // nChromosomeSide is either 0 or 1


		for (map<string, double >::iterator itPheno=this->_mpPhenotypes.begin(); itPheno != this->_mpPhenotypes.end() ; ++itPheno) {
			fOutFile << (*itPheno).second << '\t';
		}


}

void Individual::WritePhenotypeHeader(ofstream &fOutFile) { // nChromosomeSide is either 0 or 1


		for (map<string, double >::iterator itPheno=this->_mpPhenotypes.begin(); itPheno != this->_mpPhenotypes.end() ; ++itPheno) {
			fOutFile << (*itPheno).first << '\t';
		}


}

bool Individual::IsDead() {
	return this->_bDead;
}

void Individual::Die() {
	if (this->_bDead) {
		throw(new Exception("Individual is already dead, cannot die again!"));
	}
	else {
		this->_bDead = true;
	}
}