#include <iostream>
#include "cSolver.h"


int main(int argc, char * argv[]){
	//Grid parametres
	cPoint	S				= cPoint(0.5, 0.0);	//Center of the circle
	double	circleRadius	= 0.7; 
	int		gridDensity		= 100;

	cSolver CASE;

	//Read profile from a file
	CASE.addProfile("../data/NACA_0012.dat");

	//Generate grid based on profile
	CASE.generateGrid(S, circleRadius, gridDensity);

	//Solve problem
	std::cout << CASE.solve(POISSON) << std::endl;
	
	//Save solution to a file
	CASE.saveSolution(DISTANCE,"../data/NACA_0012_poisson.dat");

	system("PAUSE");
	return 0;
}