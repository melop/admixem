#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include "parser/parser.h"
#include "newran02/newran.h"
#include "maths.h"



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
		vector<map<double, int> > * GetMpMarkerIndex();
		void CalculateMapDistances();
		vector<map<double, double> > * GetMpMarkerMapDistance(int nMode);
		double GetChrToGenomeRatio(int nChr);
		
	private:
		void * _pParentConfig; // convert to SimulationConfigurations!
		int	_nHaploidChrNum;
		string _szConfigFilename;
		double * pChrLen; // array containing chromosome lengths.
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

	private:
		void * _pParentConfig;// convert to SimulationConfigurations!
		int	_nHaploidChrNum;
		string _szConfigFilename;
		double * _nExpectedMaleRecPerMeiosisArm1; //array storing expected recombinations per meiosis on each chromosome's arm 1 for males
		double * _nExpectedMaleRecPerMeiosisArm2;
		double * _nExpectedFemaleRecPerMeiosisArm1;
		double * _nExpectedFemaleRecPerMeiosisArm2;

		vector<double> * pvBreakpointSamplePositions; // array of vectors containing absolute positions of the sample break points
		vector<double> * pvMaleAccuProb; //array of vectors containing male accumulative rec. prob. restart at centromere
		vector<double> * pvFemaleAccuProb;
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
		vector<pair<int, double> > * GetFormulaSymbols(string sPhenotypeName);
		vector<string> * GetFormulaSymbolStrings(string sPhenotypeName);
		Parser * GetFormula(string sPhenotypeName);

	private:
		string _szConfigFilename;
		void * _pParentConfig;// convert to SimulationConfigurations!
		map<string , string> _mpPhenotypes; //definition of phenotypes as functions of genotypes
		map<string , Parser *> _mpPhenotypeFormulae; // similar as above, but the expression is parsed
		map<string , vector< pair<int, double> > > _mpPhenotypeFormulaSymbols; // Symbols used by each formula
		map<string , vector<string> > _mpPhenotypeFormulaSymbolStrings; // Symbols as strings
};

struct GeneProperties {
	
	enum Mode {Hemizygous, Additive, Dominant};
	
	string GeneName;
	char DominantLabel;
	char RecessiveLabel;
	double DominantValue;
	double RecessiveValue;
	double DominantFreqPop1;
	double DominantFreqPop2;
	Mode AlleleMode;
	
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
		vector< map<double , GeneProperties> >  _mpGenes; //array by chromosome. each element is a map of genes, key is the position of genes
		vector< map<double, int> > _mpGeneIndex; // same structured as above, this records the consecutive index of each gene on the chromosome.
};


class SimulationConfigurations {

	public:
		SimulationConfigurations();
		void LoadFromFile(string szConfigFile);
		const string GetConfig(string sKey);
		const double GetNumericConfig(string sKey);
		MarkerConfigurations * pMarkerConfig;
		RecombProbConfigurations * pRecombProbConfig;
		PhenotypeConfigurations * pPhenotypeConfig;
		GeneConfigurations * pGeneConfig;
		const string GetConfigFileName();

	private:
		string _szConfigFilename;
		std::map<string, string> _mpConfigs; 
		std::map<string, double> _mpNumericConfigs; 

		void fnParseNumericConfigs();


};



#endif