#ifndef _functions_h_
#define _functions_h_

#include "cGrid.h"
#include "cPoint.h"
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>

struct sCircle{
	cPoint mS;
	double mR;
public:
	sCircle(cPoint & inP = cPoint(), double inR = 0.0){ mS = inP; mR = inR; } 
};


//Read points from profile
cGrid & readPoints(std::string fileName, cGrid & inGrid){
	std::fstream file;
	file.open(fileName, std::ios::in);

	if(file.good() == false)
		std::cout << "Error";

	double x,y;
	cPoint p;

	while(!file.eof()){
		file >> x >> y;
		p.setXY(x, y);
		inGrid.addNode(p);
	}

	file.close();

	return inGrid;
}

//Save generated points to file
void writePoints(std::string fileName, cGrid & inGrid){
	std::ofstream file;
	cPoint p;

	file.open(fileName, std::ios::out);

	if(file.good() == false)
		std::cout << "Error";

	int n = inGrid.nodeCount();

	for(int i = 1; i <= n; ++i)
		file << std::setw(15) << inGrid.getNode(i).getX() << std::setw(15) << inGrid.getNode(i).getY() << std::endl;
	
	file.close();
}

/*--------- Grid generation function---------------
 1) Create point on circle
 2) Create vector from point on profile and on circle, and normalize it
 3) Create points on line using unit vector 
 4) Repeat for full circle */
cGrid & gridGeneration(cGrid & inGrid, sCircle & inCircle){
	int N = inGrid.nodeCount();	//Number of airfoil's nodes
	double divisions = 10.0;	//Number of divisions along line
	double lineLength, x, y;	//
	cPoint p;
	cVector u;
	
	const double PI = 4.0*atan(1.0);	//Define pi
	
	for(int i = 1; i <= N; ++i){
		//Create circle
		p.setXY(inCircle.mS.getX() + inCircle.mR*cos((i-1)*2*PI/N), inCircle.mS.getY() + inCircle.mR*sin((i-1)*2*PI/N));
		inGrid.addNode(p);
		
		u.setByPoints(inGrid.getNode(i), p);	//Create vector using point p
		lineLength = u.getLength();				//Length of line connecting node on circle and airfoil
		u = u/lineLength;						//Unit vector
		
		for(int j = 1; j < divisions; ++j){
			x = inGrid.getNode(i).getX() + j*lineLength/divisions*u.getX();
			y = inGrid.getNode(i).getY() + j*lineLength/divisions*u.getY();
			p.setXY(x, y);
			inGrid.addNode(p);	//Add node to grid
		}
	}

	return inGrid;
}

#endif