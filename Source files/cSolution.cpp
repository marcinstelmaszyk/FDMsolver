#include "cSolution.h"

cSolution::cSolution() : cGrid() {	
	mCaseSolved = false;
}

void cSolution::addValue(double inValue, eSolutionData valueType){
	switch(valueType){
		case DISTANCE:
			mSolution.d.push_back(inValue);
			break;
		case PHI:
			mSolution.phi.push_back(inValue);						
	}
}

double cSolution::getValue(int k, eSolutionData valueType){
	switch(valueType){
		case DISTANCE:
			return mSolution.d[k];
			break;
		case PHI:
			return mSolution.phi[k];						
	}
}

double cSolution::getValue(int i, int j, eSolutionData valueType){
	int k = i%nodeCount(PROFILE) + j*nodeCount(PROFILE);	//Modulo used to make grid continuous e.g. i=N == i=0
	
	return getValue(k, valueType);
}

void cSolution::SolutionDone(){
	mCaseSolved = true;
}