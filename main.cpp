#include <iostream>
#include <fstream>
#include <string.h>
//#include "def.h"
#include "makemarkerfile.h"
#include "testexpression.h"
#include "simulation.h"
#include <math.h>
#include <string>


#pragma once

using namespace std;
extern SimulationConfigurations SimulConfig; // defined in config.cpp

int main (int argc, char * const argv[]) {

	if (argc == 2) {
		string sConfig = argv[1];
		SimulConfig.LoadFromFile(sConfig);
		UISimulation();
		return 0;
	}

    // insert code here...
    //std::cout << "Hello, World!\n";
	
	int nProgramMode;
	
	/*
	string sTest = "";
	
	getline(cin , sTest);
	
	std::cout << "Length: " << sTest << endl;
	
	if (sTest.length() > 0 )  {
		cout << "lalala" << endl;
	
	}
	*/

	cout << floor(-0.1) <<endl;
	
	
	while(true)
	{
		std::cout << "\nWhat do you want to do? 0 - Cancel,  1 - Generate markers, 2 - Export simulated results, 3 - Simulation: ";
	
		std::cin >>  nProgramMode;
		
		switch(nProgramMode) {
		
			case 0: return 0;
			
			case 1: UIMakeNewMarkerFile(); continue; //if user wants to create a new marker file

			case 4: UITestExpression(); continue;

			case 3: UISimulation(); continue;

			case 2: UIExportResults(); continue;
		
		}
		
		break;
	
	}

	

	
    return 0;
}


