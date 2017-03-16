
#ifndef CARBONDIFFUSION_ASSEMBLEMATRIX_H
#define CARBONDIFFUSION_ASSEMBLEMATRIX_H

#include <vector>
#include "../lib/Eigen/Sparse"
using std::vector;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

class AssembleMatrix {
private:
    vector< std::vector<double> > vertices;

    double basis_function(double x, double y, vector< vector<double> > vertices, int basis_index, int derivative_degree_x, int derivative_degree_y);
    double func_f(double x, double y, double t, double D);
    double gauss2d_integral_trial_test(vector< vector<double> > vertices, int basis_index_trial, int basis_index_test,
                                        int derivative_degree_x, int derivative_degree_y);
    double gauss2d_integral_test(vector< vector<double> > vertices, int basis_index,
                                 int derivative_degree_x, int derivative_degree_y, double time, double D);
public:
    AssembleMatrix();
    SparseMatrix<double> solve_A(vector< vector<double> > P, vector< vector<int> > T,
            vector< vector<double> > Pb, vector< vector<int> > Tb,
            int N, int Nb, int Nlb, int derivative_degree_x, int derivative_degree_y);
    VectorXd solve_b(vector< vector<double> > P, vector< vector<int> > T,
            vector< vector<double> > Pb, vector< vector<int> > Tb,
            int N, int Nb, int Nlb, double time, double D, int derivative_degree_x, int derivative_degree_y);
};

#endif //CARBONDIFFUSION_ASSEMBLEMATRIX_H
