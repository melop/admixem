/*
 *  testexpression.cpp
 *  admixsimul
 *
 *  Created by Rosenthal Lab on 9/21/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include "testexpression.h"
#include "parser/parser.h"

using namespace std;

void UITestExpression() {

	Parser p("");
	
	
	p["chr1_30"] = 1;
	p["chr1_20"] = 0;
	
	
	double result = p.Evaluate ("if(chr1_30==chr1_20, 10000 , chr1_30 + chr1_20 * 2)");   // use in expression 
    
	cout << ( p.symbols_.find("chr1_31") == p.symbols_.end()? "true": "false") << endl;
	
	cout << result << endl;
	
	result = p.Evaluate ("if(chr1_30==0, sqrt(20000) , chr1_30 + chr1_10 * 2)");   // use in expression 
    
	cout << result << endl;

	result = p.Evaluate ("(chr1_30 == chr1_20)*3");   // use in expression 
    
	cout << result << endl;
	
	
	for(std::map<std::string,double>::iterator it = p.symbols_.begin(); it != p.symbols_.end(); ++it) 
	{
		
		cout << (it->first) << "\n"; //list all keys
	}


}