/*
 *  Individual.h
 *  admixsimul
 *
 *  Created by Rosenthal Lab on 11/21/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <iterator>
#include <map>
//#include "def.h"
#include "makemarkerfile.h"
#include "newran02/newran.h"
#include "config.h"
#include "maths.h"
//#include "Population.h"


#pragma once
#ifndef _INDIVIDUAL_FILE
#define _INDIVIDUAL_FILE 0

using namespace std;

struct GeneticLocus {
		
	//int Chromosome;
	//double Position;
	char Ancestry;
	char Allele;
};

struct Marker: public GeneticLocus {

};

struct Gene: public GeneticLocus {
	double Value;
};
/*
class Chromosome {

	private:
	double _nTotalLen; 
	double _nCentromerePos;  // position of the centromere.
};

class Genome {

	public:
	void BuildGenome( Chromosome * pPaternalGamete,  Chromosome * pMaternalGamete);

	private:
	int _n; // n = ? ploidy is always 2.
	Chromosome* Chromosomes; // size will be 2n. Odd number and the next even number are homologs.
	
	
};
*/

class Individual {

    public:
	Individual(void);
	Individual(void * pPop, char nAncestryLabel); //constructor for founders
	Individual(Individual * pFather, Individual * pMother, bool &bSuccess); //Constructor . Fertilize a new individual given parents
	
	enum Sex {Male, Female};
	//Genome CurrentGenome;
	
	
	~Individual(void); 
	
	//Public methods
	bool Court(Individual * pChooser); // Only available to males
	double EvaluateCourter(Individual * pCourter); //Assign a probability for the female to mate with the given courter. This is the preference function of the female.
	int HandleCourter(Individual * pCourter, bool bIgnoreGlobalRules); // Handling the courtship. Only available to females. EvaluateCourter() will be called, then self condition will also be considered. If accepted, will inseminate some eggs. return value are the eggs inseminated
	Sex	GetSex();
	unsigned int GetID();
	unsigned int GetFatherId();
	unsigned int GetMotherId();
	int GetGameteNum();
	int GetMateNumber();
	double GetPhenotype(string sPhenotype);

	void GetGamete(vector< vector<Marker> > &vMarkers, vector< vector<Gene> > &vGenes );
	void GiveBirth(vector<Individual *> &vOffSprings, int nNum, bool bIgnoreGlobalRules); // give birth to any number inseminated eggs, pass -1 to nNum to get all. this takes a lot of memory so be careful.
	
	void DumpMarkers(ofstream &fOutFile, int nChromosomeSide); // write all markers to file stream. each chromosome separated by -1.
	void DumpGenes(ofstream &fOutFile, int nChromosomeSide); // write all markers to file stream. each chromosome separated by -1.
	void DumpPhenotypes(ofstream &fOutFile);
	void WritePhenotypeHeader(ofstream &fOutFile);
	bool IsDead();
	void Die(); // sets the dead flag
	void ChangePopulation(void * pPop);
	void * GetPop();
	void * GetPrevPop(); //What is the prev pop? If the individual migrated, then prev pop will differ from the current.
	//friend class Individual;

	//void GetGamete( Chromosome * pGamete);

	private:
	 void fnDetermineSex();
	 void fnDeterminePhenotypes();
	 void fnDetermineNumGametes(); //determines how many gametes does this female have.
	 bool fnSetParserPopWide(Parser * pParser,  vector<string> &vSymbols, const string &sPrefix,  map< string , pair< double, double> > &mpSumPhenotypes);

	 unsigned int _nFatherId; //id of father
	 unsigned int _nMotherId; //id of mother
	 unsigned int _nId; // fish id
	 void * _pPop; //pointer to the current population it's living in.
	 void * _pPrevPop; // pointer to the previous population, before migration. this is relavant in oblique imprinting, where the individual learned before migration.
	 Sex _bSex;
	 bool _bMatured;
	 bool _bDead;
	 double _nCondition; //Condition of fish
	 double _nAge; //Age of the individual
	 unsigned int _nAvailableGametes; // how many more gametes available? Only meaningful for females. males have unlimited num of gametes (-1)
	 
	 int _nHaploidChrNum;
	 int _nTotalChrNum;

	 std::map<string, double >  _mpPhenotypes; // map recording phenotypic values directly calculated from genotypic values;
	 std::map<string, double >  _mpEnvPhenotypes; // Phenotypes with added enviromental variations NOT USED
	 std::map<string, double >  _mpDadPhenotypes; //father's phenotypes, useful for imprinting
	 std::map<string, double >  _mpMomPhenotypes; //mother's phenotypes, useful for imprinting

	 //map<double, Marker> * _pmMarkers; // 2n chromosomes, each element in vector is one chromosome, the map are indexed by absolute positions on that chromosome
	 vector< vector<Marker> >  _arrMarkers; // 2n chromosomes, each
	 //map<double, Gene> * _pmGenes; // same structured as above
	 vector< vector<Gene> >  _arrGenes; // 2n chromosomes, each
	 vector<Individual *> _arrOtherParentsForOffsprings;// the fatherhoods of each inseminated egg.
	
};

#endif