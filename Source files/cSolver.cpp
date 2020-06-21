#include "cSolver.h"
#include <string>
#include <fstream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctime>
#include <eigen3/Eigen/LU>


double cSolver::distanceBetween(cPoint & p_1, cPoint & p_2){
	double x = p_2.getX() - p_1.getX();
	double y = p_2.getY() - p_1.getY();
	
	return sqrt(x*x + y*y);
}

void cSolver::dialog(std::string message){
	std::cout << message << std::endl;
}

void cSolver::addProfile(std::string adress){
	std::fstream file;
	file.open(adress, std::ios_base::in);

	if(!file.good())
		std::cout << "Error opening a file";

	double x,y;
	cPoint point;

	//Read pairs of coordinates until end of file
	while(!file.eof()){
		file >> x >> y;
		mSolution.addNode(point.setXY(x,y), PROFILE);
	}

	file.close();
}

void cSolver::generateGrid(cPoint & circleC, double circleR, int linesDensity){
	double x, y;
	cPoint profilePoint, circlePoint, point;
	
	mSolution.setRows(linesDensity);
	dialog("Generating grid");

	int N = mSolution.nodeCount(PROFILE);	//Number of nodes on the profile
	int M = mSolution.rowCount();			//Number of rows
	for(int j = 1; j <= M; ++j){
		for(int i = 0; i < N; ++i){
			circlePoint.setX(circleC.getX() + circleR*cos(i*2.0*M_PI/N));
			circlePoint.setY(circleC.getY() + circleR*sin(i*2.0*M_PI/N));
			profilePoint = mSolution(i);
			
			x = profilePoint.getX() + j*(circlePoint.getX() - profilePoint.getX())/linesDensity;
			y = profilePoint.getY() + j*(circlePoint.getY() - profilePoint.getY())/linesDensity;

			mSolution.addNode(circlePoint.setXY(x,y), GRID);
		}
	}
	dialog("Grid generated!");
}

void cSolver::showGrid(eGridMode mode) const {
	for(int i = 0; i < mSolution.nodeCount(mode); i++)
		std::cout << "[" << i << "] " << mSolution(i) << std::endl;
}

void cSolver::saveGrid(std::string filePath, eFileExtension extension){
	std::fstream file;

	dialog("Saving grid...");

	file.open(filePath, std::ios_base::out);
	if(!file.good()){
		std::cout << "Error saving grid!";
		exit(1);
	}
		
	if(extension == TECPLOT){
		file << "VARIABLES = \"X\", \"Y\"," << std::endl << "ZONE I=";
		file << mSolution.nodeCount(PROFILE) << ", J=";
		file << mSolution.rowCount() + 1 << ", DATAPACKING=POINT" << std::endl << std::endl;
	}

	int N = mSolution.nodeCount();	//Number of nodes to save

	for(int i = 0; i < N; ++i)
		file << std::setw(15) << mSolution(i).getX() << std::setw(15) << mSolution(i).getY() << std::endl;

	dialog("Grid saved!");
}

void cSolver::solveBruteForce(){
	double d_min	= 0.0;	//Minimal distance computed for particular node (by far)
	double d		= 0.0;	//Current computed distance between nodes

	//Distances of nodes on profile
	for(int i = 0; i < mSolution.nodeCount(PROFILE); ++i)
			mSolution.addValue(0.0);

	//Distances of nodes 
	for(int j = 1; j <= mSolution.rowCount(); ++j)					//Iterate by rows
		for(int i = 0; i < mSolution.nodeCount(PROFILE); ++i){		//Iterate by nodes in j-row
			d_min = distanceBetween(mSolution(i,j), mSolution(0));	//Initial value of d_min

			for(int k = 0; k < mSolution.nodeCount(PROFILE); ++k){  //Iterate by all profile's nodes
				d =	distanceBetween(mSolution(i,j), mSolution(k));	//Distance computed for (i,j) node
				d_min = (d < d_min) ? d : d_min;					
			}
			mSolution.addValue(d_min);	//Set computed minimal distance for (i,j) node
		}
}

void cSolver::solvePoisson(){
	//Count nodes
	int n_prof	= mSolution.nodeCount(PROFILE);		//Nodes on the profile
	int n_grid	= mSolution.nodeCount();			//Nodes in grid
	int m		= mSolution.rowCount();				//Number of rows	

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(6*n_grid);
	f.resize(n_grid);

	//Set Dirichlet boundary condition k=i, j=0
	for(int k = 0; k < n_prof; ++k){
		tripletList.push_back(T(k,k,1.0));	//Element on diagonal		
		f(k) = 0.0;							//Right hand side value
	}

	//Set Neumann boundary conditions 
	int k, j = mSolution.rowCount();	
	for(int i = 0; i < n_prof; ++i){
		k = i + j*n_prof;								//Nodes on boundary
		tripletList.push_back(T(k, k, 1.0));
		tripletList.push_back(T(k, k - n_prof, -1.0)); 		
		f(k) = 0.0;										//Right hand side value
	}

	//Compute matrix coefficients	
	double A,B,C,D,a,b,c,d,e;
	for(int j = 1; j < m; ++j){
		for(int i = 0; i < n_prof; ++i){
			k = i + j*n_prof;
			
			//Poisson equation coefficients
			A = 1.0/partial(k, h_1);
			B = 1.0/partial(k, h_2);
			C = 1.0/partial(k, h_1)/partial(k, h_2)*(partial(k, x_v)*partial(k, x_uv) + partial(k, y_v)*partial(k, y_uv));
			C -= 1.0/partial(k, h_1)/partial(k, h_1)*(partial(k, x_u)*partial(k, x_uu) + partial(k, y_u)*partial(k, y_uu));
			D = 1.0/partial(k, h_1)/partial(k, h_2)*(partial(k, x_u)*partial(k, x_uv) + partial(k, y_u)*partial(k, y_uv));
			D -= 1.0/partial(k, h_2)/partial(k, h_2)*(partial(k, x_v)*partial(k, x_vv) + partial(k, y_v)*partial(k, y_vv));
			
			//Coefficient of PHI(i,j-1)
			a = B/h_v/h_v - D/2.0/h_v;
			tripletList.push_back(T(k, k - n_prof, a));
			
			//Coefficient of PHI(i-1,j)
			b = A/h_u/h_u - C/2.0/h_u;
			if(i == 0)	//Left PHI doesn't exist, use PHI from the end of the row
				tripletList.push_back(T(k, k + n_prof - 1, b));
			else
				tripletList.push_back(T(k, k-1, b));
			
			//Coefficient of PHI(i,j)
			c = -2.0*(A/h_u/h_u + B/h_v/h_v);
			tripletList.push_back(T(k, k, c));
				
			//Coefficient of PHI(i+1,j)
			d = A/h_u/h_u + C/2.0/h_u;
			if(i == n_prof - 1)	//Right PHI doesn't exist, use PHI from the begining of the row
				tripletList.push_back(T(k, k - n_prof + 1, d));
			else
				tripletList.push_back(T(k, k+1, d));
			
			//Coefficient of PHI(i,j+1)
			e = B/h_v/h_v + D/2.0/h_v;
			tripletList.push_back(T(k, k + n_prof, e));
			
			f(k) = -1.0;	//Right hand side value			
		}
	}
	
	//Insert coefficients into sparse matrix
	K.resize(n_grid, n_grid);
	K.setFromTriplets(tripletList.begin(), tripletList.end());
	
	//Initial guess
	x_guess.resize(n_grid);
	double step = 1./m;
	for(int i = 0; i < n_grid; ++i)
		x_guess(i) = i%n_prof*step;

	//Solve matrix equation
	solver.setTolerance(0.1);	//Tolerance
	K.makeCompressed();
	x = solver.compute(K).solveWithGuess(f, x_guess);	
	std::cout << "    #iterations: " << solver.iterations() << "   estimated error: " << solver.error() << std::endl;		
	
	//Save PHI in mSolution
	for(int j = 0; j <= m; ++j)
		for(int i = 0; i < n_prof; ++i)
			mSolution.addValue(x(i + j*n_prof), PHI);
	
	computeDistance();
}

void cSolver::solveEikonal(){

}

double cSolver::partial(int k, ePartial partial){
	int j = k/mSolution.nodeCount(PROFILE);
	int i = k%mSolution.nodeCount(PROFILE);
	
	double X_u, X_v, Y_u, Y_v;
	
	switch(partial){
	case h_1:
		X_u = (mSolution(i+1,j).getX() - mSolution(i-1,j).getX())/2.0/h_u;		
		Y_u = (mSolution(i+1,j).getY() - mSolution(i-1,j).getY())/2.0/h_u;		
		return X_u*X_u + Y_u*Y_u;	//Returns h_1^2
	case h_2:
		X_v = (mSolution(i,j+1).getX() - mSolution(i,j-1).getX())/2.0/h_v;		
		Y_v = (mSolution(i,j+1).getY() - mSolution(i,j-1).getY())/2.0/h_v;		
		return X_v*X_v + Y_v*Y_v;	//Returns h_2^2

	case x_u:
		return (mSolution(i+1,j).getX() - mSolution(i-1,j).getX())/2.0/h_u;		
	case x_v:
		return (mSolution(i,j+1).getX() - mSolution(i,j-1).getX())/2.0/h_v;
	case x_uu:
		return (mSolution(i+1,j).getX() - 2.0*mSolution(i,j).getX() + mSolution(i-1,j).getX())/h_u/h_u;
	case x_vv:
		return (mSolution(i,j+1).getX() - 2.0*mSolution(i,j).getX() + mSolution(i,j-1).getX())/h_v/h_v;
	case x_uv:
		return (mSolution(i+1,j+1).getX() - mSolution(i+1,j-1).getX() - mSolution(i-1,j+1).getX() + mSolution(i-1,j-1).getX())/4.0/h_u/h_v;

	case y_u:
		return (mSolution(i+1,j).getY() - mSolution(i-1,j).getY())/2.0/h_u;		
	case y_v:
		return (mSolution(i,j+1).getY() - mSolution(i,j-1).getY())/2.0/h_v;
	case y_uu:
		return (mSolution(i+1,j).getY() - 2.0*mSolution(i,j).getY() + mSolution(i-1,j).getY())/h_u/h_u;
	case y_vv:
		return (mSolution(i,j+1).getY() - 2.0*mSolution(i,j).getY() + mSolution(i,j-1).getY())/h_v/h_v;
	case y_uv:
		return (mSolution(i+1,j+1).getY() - mSolution(i+1,j-1).getY() - mSolution(i-1,j+1).getY() + mSolution(i-1,j-1).getY())/4.0/h_u/h_v;	
	}
}

double cSolver::solve(eSolutionType type){
	double elapsedTime = 0.0;			//Time of complete distance evaluation
	clock_t start_time, stop_time;
	
	if(type == BRUTE_FORCE){
		dialog("Solving brute force...");

		start_time = clock();
		solveBruteForce();
		stop_time = clock();	

	}else if(type == POISSON){	
		dialog("Solving Poisson...");

		start_time = clock();
		solvePoisson();		
		stop_time = clock();
		
	}else{
		dialog("Solving Eikonal...");

		start_time = clock();
		solveEikonal();
		stop_time = clock();
	}

	elapsedTime = static_cast<double>(stop_time - start_time) / CLOCKS_PER_SEC;
	mSolution.SolutionDone();
	dialog("Solution solved!");

	return elapsedTime;
}


void cSolver::saveSolution(eSolutionData dataType, std::string filePath, eFileExtension extension){
	std::fstream file;

	dialog("Saving solution...");

	file.open(filePath, std::ios_base::out);
	if(!file.good())
		std::cout << "Error saving grid!";

	int N = mSolution.nodeCount();

	//Add header (format needed by TecPlot360)
	if(extension == TECPLOT){
		file << "VARIABLES = \"X\", \"Y\", \"U\"," << std::endl << "ZONE I=";
		file << mSolution.nodeCount(PROFILE) + 1  << ", J=";
		file << mSolution.rowCount() + 1 << ", DATAPACKING=POINT" << std::endl << std::endl;
	}

	int n = mSolution.nodeCount(PROFILE);
	
	for(int i = 0; i < N; ++i){
		file << std::setw(15) << mSolution(i).getX() << std::setw(15) << mSolution(i).getY() << std::setw(15) << mSolution.getValue(i, dataType) << std::endl;

		//Add additional point at the end of row (closes grid in TecPlot360)
		if(i != 0 && (i+1)%n == 0)
			file << std::setw(15) << mSolution(i+1-n).getX() << std::setw(15) << mSolution(i+1-n).getY() << std::setw(15) << mSolution.getValue(i+1-n, dataType) << std::endl;
	}
	dialog("Solution saved!");
	file.close();
}

void cSolver::computeDistance(){
	//Transform coordinates (computational_space -> physical_space)
	Eigen::Matrix2d J;			//Transformation matrix x->u
	Eigen::Matrix2d I;			//Inverse transformation matrix u->x
	Eigen::Vector2d phi_xy;
	Eigen::Vector2d phi_uv;

	double phi_u, phi_v, phi_x, phi_y, d;
	bool invertible;
	int k;	//Node number

	int n_prof = mSolution.nodeCount(PROFILE);	//Nodes on profile
	int n_rows = mSolution.rowCount();			//Number of rows

	//Insert zero distance to solution vector (nodes on the profile)
	for(int i = 0; i < n_prof; ++i)
		mSolution.addValue(0.0, DISTANCE);
	
	for(int j = 1; j <= n_rows; ++j){
		for(int i = 0; i < n_prof; ++i){
			k = i + j*n_prof;

			J(0,0) = partial(k, x_u);
			J(0,1) = partial(k, y_u);
			J(1,0) = (mSolution(i,j).getX() - mSolution(i,j-1).getX())/h_v;
			J(1,1) = (mSolution(i,j).getY() - mSolution(i,j-1).getY())/h_v;
			
			J.computeInverseWithCheck(I, invertible);
			
			phi_uv(0) = (mSolution.getValue(i+1,j,PHI) - mSolution.getValue(i-1,j,PHI))/2.0/h_u;
			phi_uv(1) = (mSolution.getValue(i,j,PHI) - mSolution.getValue(i,j-1,PHI))/h_v;

			phi_xy = I*phi_uv;			//Transform!
			d =	phi_xy.dot(phi_xy);		//Dot product == phi_x*phi_x + phi_y*phi_y
			d = sqrt(d + 2*mSolution.getValue(i, j, PHI)) - sqrt(d);	//Distance formula

			mSolution.addValue(d, DISTANCE);	//Save computed distance		
		}
	}
}