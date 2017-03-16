#include <iostream>
#include <vector>
#include "../lib/Eigen/Sparse"
#include "../lib/Eigen/IterativeLinearSolvers"
#include "../include/CarbonDiffusion.h"

const double PI = atan(1.0)*4;

using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::ConjugateGradient;
using std::cout;
using std::vector;
using std::endl;

double factor_D(double temperature);
SparseMatrix<double> treat_boundary_A(vector< vector<int> > boundary_nodes, SparseMatrix<double> A);
VectorXd treat_boundary_b(vector<double> omega, vector< vector<int> > boundary_nodes, VectorXd b, vector< vector<double> > Pb, double t);
double boundary_function(vector<double> omega, double x, double y, double t);

int main() {
    vector<double> omega = {-1, 1, -1, 1};
    vector<double> h = {1.0/2, 1.0/2};
    string basis_type = "linear";
    double t_min = 0;
    double t_max = 3600;
    double dt = 1800;
    int Nt = (t_max - t_min) / dt;
    double temperature = 1183;
    double D = factor_D(temperature);
    InfoMatrix generator = InfoMatrix(omega, h, basis_type);
    vector< vector<double> > P = generator.getCoordinateMatrix_mesh();
    vector< vector<int> > T = generator.getNodeIndicesMatrix_mesh();
    vector< vector<double> > Pb = generator.getCoordinateMatrix_fe();
    vector< vector<int> > Tb = generator.getNodeIndicesMatrix_fe();
    vector< vector<int> > boundary_nodes = generator.getBoundaryNode();

    cout << "Matrix P:" << endl;
    for (int i = 0; i < P.size(); i++) {
        for (int j = 0; j < 2; j++) {
            cout << P[i][j] << " ";
        }
        cout << endl;
    }
    cout << "Matrix T:" << endl;
    for (int i = 0; i < T.size(); i++) {
        for (int j = 0; j < 3; j++) {
            cout << T[i][j] << " ";
        }
        cout << endl;
    }
    cout << "Boundary Nodes:" << endl;
    for (int i = 0; i < boundary_nodes.size(); i++) {
        for (int j = 0; j < 2; j++) {
            cout << boundary_nodes[i][j] << " ";
        }
        cout << endl;
    }
    AssembleMatrix assembleMatrix;
    int N = T.size();
    int Nb = Pb.size();
    int Nlb = T[0].size();
    cout << "N:" << N << " Nb: " << Nb << " Nlb: " << Nlb << endl;
    SparseMatrix<double> A1 = assembleMatrix.solve_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 0);
    SparseMatrix<double> A2 = assembleMatrix.solve_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 1);
    SparseMatrix<double> A3 = assembleMatrix.solve_A(P, T, Pb, Tb, N, Nb, Nlb, 1, 0);
    cout << "A1:" << endl;
    cout << A1 << endl;
    cout << "A2:" << endl;
    cout << A2 << endl;
    cout << "A3:" << endl;
    cout << A3 << endl;
    cout << "A1 + dt * D * (A2 + A3)" << endl;
    SparseMatrix<double> A = A1 + dt * D * (A2 + A3);
    cout << A << endl;
    cout << "b:" << endl;
    VectorXd b = assembleMatrix.solve_b(P, T, Pb, Tb, N, Nb, Nlb, 0, D, 0, 0);
    cout << b << endl;
    cout << "treat bounday node:" << endl;
    A = treat_boundary_A(boundary_nodes, A);
    b = treat_boundary_b(omega, boundary_nodes, b, Pb, 0);
    cout << "A:" << endl;
    cout << A << endl;
    cout << b << endl;
    cout << "solve:" << endl;

    ConjugateGradient<SparseMatrix<double> > solver;
    VectorXd x = solver.compute(A).solve(b);
    cout << x << endl;

    return 0;
}

double factor_D(double temperature) {
    double D0 = 0.162 * 100;
    double Q = 137800;
    double R = 8.314;
    double T = temperature;
    double D = D0 * exp(-Q / (R * T));

    return D;
}

SparseMatrix<double> treat_boundary_A(vector< vector<int> > boundary_nodes, SparseMatrix<double> A) {
    long nbn = boundary_nodes.size();
    for (int i = 0; i < nbn; i++) {
        if (boundary_nodes[i][0] == -1) {
            int k = boundary_nodes[i][1];
            for (int i = 0; i < A.innerSize(); i++) {
                A.coeffRef(k, i) = 0;
            }
            A.coeffRef(k, k) = 1;
        }
    }
    return A;
}

VectorXd treat_boundary_b(vector<double> omega, vector< vector<int> > boundary_nodes, VectorXd b, vector< vector<double> > Pb, double t) {
    long nbn = boundary_nodes.size();
    for (int i = 0; i < nbn; i++) {
        if (boundary_nodes[i][0] == -1) {
            int k = boundary_nodes[i][1];
            b[k] = boundary_function(omega, Pb[i][0], Pb[i][1], t);
        }
    }
    return b;
}

double boundary_function(vector<double> omega, double x, double y, double t) {
    double left = omega[0];
    double right = omega[1];
    double bottom = omega[2];
    double top = omega[3];
    double result;

    if (x == left) {
        result = 0;
    } else if (x == right) {
        result = 0;
    } else if (y == bottom) {
        result = -t * sin(x);
    } else if (y == top) {
        result = -t * sin(x);
    }
    return result;
}

double exact_function(double x, double y, double t) {
    double result;
    result = t * sin(x) * cos(y);

    return result;
}