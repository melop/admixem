#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include "Individual.h"
#include "maths.h"

#pragma once
#ifndef _POPULATION_FILE
#define _POPULATION_FILE 0

using namespace std;
//extern int nCurrIndividualId; // defined in individual.cpp

class Population
{
public:
	Population(void);
	~Population(void);
	void Init(string sPopName, int nPopId, char nAncestryLabel, int nPopInitSize, int nPopMaxSize, double nMaleRatio);
	bool Immigrate(Individual * pImmigrant, bool bForceExistingDie);
	Individual * Emigrate(); // delete individual from the population and return it 
	bool Breed();
	void KillOldGen();
	void Sample(ofstream &fMarkerOutFile, ofstream &fGeneOutFile,  ofstream &fPhenotypeOutFile);
	void NaturalSelection();
	int GetPopSize(int nMode);
	int GetPopId();
	string  GetPopName();

private:
	string _sPopName;
	int		_nPopId;
	char   _nAncestryLabel;
	int	   _nPopInitSize;
	int    _nPopMaxSize;
	int    _nCurrMales;
	int    _nCurrFemales;
	double _nMaleRatio;
	bool _bBred; //already bred.
	vector< Individual *> _mpMales;  // array storing males
	vector< Individual *> _mpFemales; //array storing females
	vector< Individual *> _mpNewGenMales;  // array storing males
	vector< Individual *> _mpNewGenFemales; //array storing females

	void fnWriteIndividualMarkers(ofstream &fOutFile, Individual * pInd);
	void fnWriteIndividualGenes(ofstream &fOutFile, Individual * pInd);
	void fnWriteIndividualPhenotypes(ofstream &fOutFile, Individual * pInd);
	void fnWriteIndividualID(ofstream &fOutFile, Individual * pInd);
	void fnSample(ofstream &fMarkerOutFile, ofstream &fGeneOutFile,  ofstream &fPhenotypeOutFile, Individual::Sex bSex);
};

#endif