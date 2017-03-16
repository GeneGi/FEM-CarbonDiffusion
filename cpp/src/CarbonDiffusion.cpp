
#include "../include/CarbonDiffusion.h"


void CarbonDiffusion::solve(vector<double> omega, vector<double> h, int Nt, int dt, string basis_type, int temperature) {
    //generate info matrix and boundary nodes
    cout << "generate info matrix" << endl;
    InfoMatrix generator = InfoMatrix(omega, h, basis_type);
    vector< vector<double> > P = generator.getCoordinateMatrix_mesh();
    vector< vector<int> > T = generator.getNodeIndicesMatrix_mesh();
    vector< vector<double> > Pb = generator.getCoordinateMatrix_fe();
    vector< vector<int> > Tb = generator.getNodeIndicesMatrix_fe();
    vector< vector<int> > boundary_nodes = generator.getBoundaryNode();
    int N = T.size();
    int Nb = Pb.size();
    int Nlb = T[0].size();
    double D = factor_D(temperature);

    // assemble stiffness matrix A
    cout << "assemble stiffness matrix" << endl;
    SparseMatrix<double> A1 = assembleMatrix.solve_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 0);
    SparseMatrix<double> A2 = assembleMatrix.solve_A(P, T, Pb, Tb, N, Nb, Nlb, 0, 1);
    SparseMatrix<double> A3 = assembleMatrix.solve_A(P, T, Pb, Tb, N, Nb, Nlb, 1, 0);
    A = A1 + dt * D * (A2 + A3);
    //cout << "A: " << A << endl;
    A = treat_boundary_A(boundary_nodes, A);
//    cout << "A: " << A << endl;
    solver.compute(A);
    // init result C
    cout << "init Carbon" << endl;
    vector< VectorXd > C(Nt, VectorXd(Nb));
    for (int i = 0; i < Nb; i++) {
        C[0][i] = initial_function(Pb[i][0], Pb[i][1]);
    }
    // time iteration
    for (int k = 0; k < Nt-1; k++) {
        cout << "time iteration " << k << endl;
        VectorXd b1 = assembleMatrix.solve_b(P, T, Pb, Tb, N, Nb, Nlb, (k+1)*dt, D, 0, 0);
//        cout << b1 << endl;
        VectorXd b2 = A1 * C[k];
//        cout << "b2: " << b2 << endl;
        b = dt * b1 + b2;
//        cout << "before b:" << b << endl;
        b = treat_boundary_b(omega, boundary_nodes, b, Pb, (k+1)*dt);
//        cout << " after b:" << b << endl;
        C[k+1] = solver.solve(b);
//        cout << "C[k+1]: " << endl;
//        cout << C[k+1] << endl;
    }
    vector<double> E(Nt);
    for (int i = 0; i < Nt; i++) {
        double max = 0;
        for (int j = 0; j < Nb; j++) {
            double exact = exact_function(Pb[j][0], Pb[j][1], dt*i);
            double error = fabs(exact - C[i][j]);
            if (error > max) {
               max = error;
            }
        }
        E[i] = max;
    }
    carbon_concentration = C;
    max_error = E;
}

vector< VectorXd > CarbonDiffusion::getCarbonConcentration() {
    return carbon_concentration;
}

vector<double> CarbonDiffusion::getMaxError() {
    return max_error;
}

SparseMatrix<double> CarbonDiffusion::treat_boundary_A(vector<vector<int> > boundary_nodes,
                                                              SparseMatrix<double> A) {
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

VectorXd CarbonDiffusion::treat_boundary_b(vector<double> omega, vector<vector<int> > boundary_nodes, VectorXd b,
                                                  vector<vector<double> > Pb, double t) {
    long nbn = boundary_nodes.size();
    for (int i = 0; i < nbn; i++) {
        if (boundary_nodes[i][0] == -1) {
            int k = boundary_nodes[i][1];
            b[k] = boundary_function(omega, Pb[k][0], Pb[k][1], t);
        }
    }
    return b;
}

double CarbonDiffusion::factor_D(double temperature) {
    double D0 = 0.162 * 100;
    double Q = 137800;
    double R = 8.314;
    double T = temperature;
    double D = D0 * exp(-Q / (R * T));

    return D;
}

double CarbonDiffusion::initial_function(double x, double y) {
    return 0;
}

double CarbonDiffusion::boundary_function(vector<double> omega, double x, double y, double t) {
    double left = omega[0];
    double right = omega[1];
    double bottom = omega[2];
    double top = omega[3];
    double result = 0;

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

double CarbonDiffusion::exact_function(double x, double y, double t) {
    double result;
    result = t * sin(x) * cos(y);

    return result;
}