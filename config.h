#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include <set>
#include "parser/parser.h"
#include "newran02/newran.h"
#include "maths.h"
#include "spline.h"


#define MAXTEXT 1000000


#pragma once
#ifndef _CONFIG_FILE
#define _CONFIG_FILE 0

using namespace std;

class MarkerConfigurations {

	public:
		MarkerConfigurations(void * pParentConfig);
		void LoadFromFile(string szConfigFile);
		int GetHaploidChromosomeNum();
		map<double, double> * GetMarkersInfoByPop(int nPopId);
		double GetCentromerePosition(int nChr);
		double GetChromosomeLength(int nChr);
		double GetPaternalTransBias(int nChr);
		vector<map<double, int> > * GetMpMarkerIndex();
		void CalculateMapDistances();
		vector<map<double, double> > * GetMpMarkerMapDistance(int nMode);
		double GetChrToGenomeRatio(int nChr);

		
	private:
		void * _pParentConfig; // convert to SimulationConfigurations!
		int	_nHaploidChrNum;
		string _szConfigFilename;
		double * pChrLen; // array containing chromosome lengths.
		double * pPaternalTransBias;
		double nTotalGenomeLen;
		double * pCentromerePos; //array containing centromere positions
		vector< map<double, double> * > vMarkerFreq; // array containing array. first dimension is population, second dimension is chromosome, and the key for the map is marker position, value of map is frequency of the dominant allele
		vector< map<double, int> > vMarkerIndex; //divided by chromosomes, key are absolute position, values are index
		vector< map<double, double> > vAbsToMapDistanceMale; // neucleotide position to map distance 
		vector< map<double, double> > vAbsToMapDistanceFemale;
		vector< map<double, double> > vAbsToMapDistanceAvg;
};

class RecombProbConfigurations {

	public:
		RecombProbConfigurations(void * pParentConfig);
		void LoadFromFile(string szConfigFile);
		void GetBreakPointsByArm(bool bSex, int nChr, int nArm, vector<double> &vRet); // sex: true male, false female, nChr starts from 0 arm is either 1 or 2.
		vector<double> * GetBreakPointMapDistances(int nMode);
		vector<double> * GetBreakPointSamplePositions();
		bool IsUseUniform();
		double HowManyBreakpointsOnArm(bool bSex, int nChr, int nArm);

	private:
		void * _pParentConfig;// convert to SimulationConfigurations!
		int	_nHaploidChrNum;
		string _szConfigFilename;
		double * _nExpectedMaleRecPerMeiosisArm1; //array storing expected recombinations per meiosis on each chromosome's arm 1 for males
		double * _nExpectedMaleRecPerMeiosisArm2;
		double * _nExpectedFemaleRecPerMeiosisArm1;
		double * _nExpectedFemaleRecPerMeiosisArm2;
		bool _bUseUniform;

		vector<double> * pvBreakpointSamplePositions; // array of vectors containing absolute positions of the sample break points
		vector<double> * pvMaleAccuProb; //array of vectors containing male accumulative rec. prob. restart at centromere
		vector<double> * pvFemaleAccuProb;
		vector<tk::spline> vMaleAccuProbSplineArm1;
		vector<tk::spline> vFemaleAccuProbSplineArm1;
		vector<tk::spline> vMaleAccuProbSplineArm2;
		vector<tk::spline> vFemaleAccuProbSplineArm2;
		vector<double> vLastBpOnArm1;
		vector<double> * pvMaleMapDistance; // array of map distance of current sample point to previous one. in the unit of Morgan, not centiMorgan.
		vector<double> * pvFemaleMapDistance;
		vector<double> * pvAvgMapDistance;

};

class PhenotypeConfigurations {

	public:
		PhenotypeConfigurations(void * pParentConfig);
		void LoadFromFile(string szConfigFile);
		int GetNumPhenotypes();
		void InitKeys( map<string, double> * pMap);
		void GetKeys( vector<string> &vKeys);
		vector<pair<int, double> > * GetFormulaSymbols(string sPhenotypeName);
		vector< string > * GetFormulaParentSymbols(string sPhenotypeName, int nParent);
		vector< string > * GetFormulaSymbolStrings(string sPhenotypeName);
		Parser * GetFormula(string sPhenotypeName);
		bool IsReferToParentPhe(string sPhenotypeName); //does this phenotype depend on parent's phenotype?

	private:
		string _szConfigFilename;
		map<string, bool > _mpIsReferToParentPhe;
		void * _pParentConfig;// convert to SimulationConfigurations!
		map<string , string> _mpPhenotypes; //definition of phenotypes as functions of genotypes
		map< int, map<string , Parser *> > _mpmpPhenotypeFormulae; // similar as above, but the expression is parsed. first level key is CPU core, second level key is name of phenotype.
		map<string , vector< pair<int, double> > > _mpPhenotypeFormulaSymbols; // Symbols used by each formula
		map<string , vector< string > > _mpPhenotypeFormulaDadSymbols; // Symbols used by each formula, referring to dad's phenotypes
		map<string , vector< string > > _mpPhenotypeFormulaMomSymbols; // Symbols used by each formula, referring to dad's phenotypes
		map<string , vector<string> > _mpPhenotypeFormulaSymbolStrings; // Symbols as strings
};

class NaturalSelectionConfigurations {

	public:
		NaturalSelectionConfigurations(void * pParentConfig);
		void LoadFromFile(string szConfigFile);
		list< vector<string> > * GetFormulaSymbolStrings(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsCourter(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsSelf(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPop(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPopCourter(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPopChooser(string sPop);
		list< pair< Parser *, int> > * GetFormulae(string sPop);
		list< pair< Parser *, int> > * GetFreqDependentFormulae(string sPop);
		map< int, map<string , list< pair<Parser *, int> > > > * GetFormulaeAllCPUs(); // 
		map< int, map<string , list< pair<Parser *, int> > > > * GetFreqDependentFormulaeAllCPUs(); // 
		bool IgnoreGlobalRules(int nGen);

	private:
		string _szConfigFilename;
		void * _pParentConfig;// convert to SimulationConfigurations!
		map<string , list< string > > _mpRules; //definition of selection rules (freq independent rules). key is population name
		map<string , list< string > >  _mpFreqDependentRules; //definition of selection rules (freq dependent rules).  key is population name
		map< int, map<string , list< pair< Parser *, int> > > > _mpmpRuleFormulae; // similar as above, but the expression is parsed  first level key is CPU core, second level key is population name
		map< int, map<string , list< pair< Parser *, int> > > > _mpmpFreqDependentRuleFormulae; // similar as above, but the expression is parsed
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStrings; // Symbols as strings
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsCourter; // Symbols that are courter values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsSelf; // Symbols that are self (female) values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsPopWide; // Symbols that are population-wide values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsPopWideCourter; // Symbols that are population-wide courter values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsPopWideChooser; // Symbols that are population-wide chooser values
		set< int > _vSpecialGens; //generations with special rules. if the generation is in here, all rules with -1 will be ignored.
};

class SexualSelectionConfigurations {

	public:
		SexualSelectionConfigurations(void * pParentConfig);
		void LoadFromFile(string szConfigFile);
		map< int, map<string , list< pair<Parser *, int> > > > * GetFormulaeAllCPUs(); // 
		list< vector<string> > * GetFormulaSymbolStrings(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsCourter(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsSelf(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPop(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPopCourter(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPopChooser(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsDad(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsMom(string sPop);

		list< vector<string> > * GetFormulaSymbolStringsPrevGenCurrPop(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPrevGenCurrPopCourter(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPrevGenCurrPopChooser(string sPop);

		list< vector<string> > * GetFormulaSymbolStringsPrevGenPrevPop(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPrevGenPrevPopCourter(string sPop);
		list< vector<string> > * GetFormulaSymbolStringsPrevGenPrevPopChooser(string sPop);

		list< pair< Parser * , int> > * GetFormulae(string sPop);
		bool IgnoreGlobalRules(int nGen);

	private:
		string _szConfigFilename;
		void * _pParentConfig;// convert to SimulationConfigurations!
		map<string , list< string > > _mpRules; //definition of selection rules. key is population name
		map< int, map<string , list< pair<Parser *, int> > > > _mpmpRuleFormulae; // similar as above, but the expression is parsed, level one key is CPU core, level two key is pop name
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStrings; // Symbols as strings
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsCourter; // Symbols that are courter values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsChooser; // Symbols that are self (female) values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsPopWide; // Symbols that are population-wide values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsPopWideCourter; // Symbols that are population-wide courter values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsPopWideChooser; // Symbols that are population-wide chooser values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsDad; // Symbols that are father's values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsMom; // Symbols that are mother's values
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsCurrPopWidePrevGen; // Symbols that are population-wide values from prev gen. of current pop
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsCurrPopWidePrevGenCourter; // Symbols that are population-wide courter values from prev gen. of current pop
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsCurrPopWidePrevGenChooser; // Symbols that are population-wide chooser values from prev gen. of current pop
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsPrevPopWidePrevGen; // Symbols that are population-wide values from prev gen. of previous pop, if the individual did not migrate, previous pop = current pop
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsPrevPopWidePrevGenCourter; // Symbols that are population-wide courter values from prev gen. of previous pop, if the individual did not migrate, previous pop = current pop
		map<string , list< vector<string> > > _mpRuleFormulaSymbolStringsPrevPopWidePrevGenChooser; // Symbols that are population-wide chooser values from prev gen. of previous pop, if the individual did not migrate, previous pop = current pop


		set< int > _vSpecialGens; //generations with special rules. if the generation is in here, all rules with -1 will be ignored.
};

struct GeneProperties {
	
	GeneProperties() {
		this->MutationProb = 0;
	}

	enum Mode {Hemizygous, Additive, Dominant};
	
	string GeneName;
	char DominantLabel;
	char RecessiveLabel;
	double DominantValue;
	double RecessiveValue;
	double DominantFreqPop1;
	double DominantFreqPop2;
	Mode AlleleMode;

	//Mutation properties, optional
	double MutationProb;
	double LowerBound;
	double UpperBound;
	Parser * pFormula;
	string sFormula;
	set<string> Pops;
	int StartGen;
	int EndGen;

	GeneProperties& operator=(const GeneProperties& oSource);
};


class GeneConfigurations {
	public:
		GeneConfigurations(void * pParentConfig);
		void LoadFromFile(string szConfigFile);
		vector< map<double , GeneProperties> > * GetMpGenes();
		vector< map<double, int> > * GetMpGeneIndex();

	private:
		string _szConfigFilename;
		int	_nHaploidChrNum;
		void * _pParentConfig;// convert to SimulationConfigurations!
		map< int, vector< map<double , GeneProperties> > >  _mpmpGenes; //array by chromosome. each element is a map of genes, key is the position of genes
		map< int, vector< map<double, int> > > _mpmpGeneIndex; // same structured as above, this records the consecutive index of each gene on the chromosome.
};


class SimulationConfigurations {

	public:
		SimulationConfigurations();
		void LoadFromFile(string szConfigFile);
		const string GetConfig(string sKey);
		void SetConfig(string sKey, string sValue);
		const double GetNumericConfig(string sKey);
		void SetNumericConfig(string sKey, double nValue);
		MarkerConfigurations * pMarkerConfig;
		RecombProbConfigurations * pRecombProbConfig;
		PhenotypeConfigurations * pPhenotypeConfig;
		GeneConfigurations * pGeneConfig;
		NaturalSelectionConfigurations * pNaturalSelConfig; 
		SexualSelectionConfigurations * pSexualSelConfig; 
		const string GetConfigFileName();
		int GetCurrGen();
		void SetCurrGen( int nGen);
		bool IsInUserSpecifiedSamplingRange(int nGen);
		int GetLogVerboseLevel();

	private:
		int _nVerboseLogLevel;
		string _szConfigFilename;
		std::map<string, string> _mpConfigs; 
		std::map<string, double> _mpNumericConfigs;
		std::vector< std::pair< int, int > > _mvSampleGens;
		int nCurrGen;

		void fnParseNumericConfigs();
		void fnParseSampleGenDef() ;


};



#endif