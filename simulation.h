#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#ifdef _WIN32
#include <direct.h>
#endif

/*
#include <sys/stat.h>
*/
#include "config.h"
#include "Population.h"
#include "maths.h"

#pragma once

using namespace std;
void UISimulation();
void PerformSimulation();
void UILoadConfig() ;
void UIExportResults();
void PerformExport(int nGen, int nRep, string sExportFolder, int nPop1Size, int nPop2Size, int nPop3Size, int nMarkerNum, double nRandSeed);

string convertInt(int number);

void fnMakeDir(const char * sPath);
bool fnFileExists(const char *filename);