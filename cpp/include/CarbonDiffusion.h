
#ifndef CARBONDIFFUSION_CARBONDIFFUSION2D_H
#define CARBONDIFFUSION_CARBONDIFFUSION2D_H

#include <iostream>
#include "../lib/Eigen/Sparse"
#include "../lib/Eigen/IterativeLinearSolvers"
#include "InfoMatrix.h"
#include "AssembleMatrix.h"

using std::cout;
using std::endl;
using std::vector;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::ConjugateGradient;

class CarbonDiffusion {
private:
    //return results
    vector< VectorXd > carbon_concentration;
    vector<double> max_error;
   //info matrix
    vector< vector<double> > P;
    vector< vector<double> > T;
    vector< vector<double> > Pb;
    vector< vector<double> > Tb;
    vector< vector<double> > boundary_nodes;

    AssembleMatrix assembleMatrix;
    ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;

    SparseMatrix<double> A;
    VectorXd b, x;

    double factor_D(double temperature);
    double initial_function(double x, double y);
    double boundary_function(vector<double> omega, double x, double y, double t);
    double exact_function(double x, double y, double t);

    SparseMatrix<double> treat_boundary_A(vector< vector<int> > boundary_nodes, SparseMatrix<double> A);
    VectorXd treat_boundary_b(vector<double> omega, vector< vector<int> > boundary_nodes, VectorXd b, vector< vector<double> > Pb, double t);

public:
    void solve(vector<double> omega, vector<double> h, int Nt, int dt, string basis_type, int temperature);
    vector< VectorXd > getCarbonConcentration();
    vector<double> getMaxError();
};



#endif //CARBONDIFFUSION_CARBONDIFFUSION2D_H
