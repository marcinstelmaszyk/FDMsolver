#ifndef _CSOLVER_H_
#define _CSOLVER_H_

#include <eigen3/Eigen/Sparse>
#include "cSolution.h"
#include "cGrid.h"

enum eFileExtension {DAT, TECPLOT};	//Output data format
enum ePartial {h_1, h_2, x_u, x_v, y_u, y_v, x_uu, x_uv, x_vv, y_uu, y_uv, y_vv};

//Shorten names
typedef Eigen::SparseMatrix<double>						Eigen_SparseMatrix;
typedef Eigen::MatrixXd									Eigen_DenseMatrix;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double>>	Eigen_BiCGSTAB;

class cSolver{
	//Results
	cSolution	mSolution;						//Cointainer for solutions
	void		dialog(std::string message);	//Show message

	//Brute force method
	double		distanceBetween(cPoint & p_1, cPoint & p_2);	//Distance between two points
	void		solveBruteForce();								//Use bruteforce method

	//Differential methods objects & functions
	Eigen_SparseMatrix	K;				//Coefficient matrix
	Eigen::VectorXd		f;				//Right hand side
	Eigen::VectorXd		x;				//Phi solution
	Eigen::VectorXd		x_guess;
	Eigen_BiCGSTAB		solver;			//Biconjugative iterative method
	
	//Poisson method
	void		solvePoisson();	

	//Eikonal method
	void		solveEikonal();			//Not implemented

	double		h_u, h_v;
	double		partial(int k, ePartial partial);
	void		computeDistance();
public:
	cSolver(double H_U = 1.0, double H_V = 1.0) : h_u(H_U), h_v(H_V) {};

	//Profile
	void		addProfile(std::string filePath);

	//Grid
	void		generateGrid(cPoint & circleCenter, double circleRadius, int linesDensity);
	void		showGrid(eGridMode mode = GRID) const;
	void		saveGrid(std::string filePath = "grid.dat", eFileExtension extension = TECPLOT);

	//Solving
	double		solve(eSolutionType type = POISSON);

	//Saving
	void		saveSolution(eSolutionData dataType = DISTANCE, std::string filePath = "solution.dat", eFileExtension extension = TECPLOT);	
};

#endif