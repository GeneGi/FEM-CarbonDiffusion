
#include "../include/AssembleMatrix.h"
AssembleMatrix::AssembleMatrix() {
    vertices.resize(3);
    for (int i = 0; i < 3; i++) {
        vertices[i].resize(2);
    }
}

Eigen::SparseMatrix<double> AssembleMatrix::solve_A(vector <vector<double> > P, vector <vector<int> > T, vector <vector<double> > Pb,
                           vector <vector<int> > Tb, int N, int Nb, int Nlb, int derivative_degree_x,
                           int derivative_degree_y) {
    Eigen::SparseMatrix<double> StiffnessMatrix(Nb, Nb);
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vertices[k][0] = P[T[i][k]][0];
            vertices[k][1] = P[T[i][k]][1];
        }
        for (int alpha = 0; alpha < Nlb; alpha++) {
            for (int beta = 0; beta < Nlb; beta++) {
                double temp = gauss2d_integral_trial_test(vertices, alpha, beta, derivative_degree_x, derivative_degree_y);
                StiffnessMatrix.coeffRef((long) Tb[i][beta], (long) Tb[i][alpha]) += temp;
            }
        }
   }
    return StiffnessMatrix;
}

Eigen::VectorXd AssembleMatrix::solve_b(vector<vector<double> > P, vector<vector<int> > T,
                                                    vector<vector<double> > Pb, vector<vector<int> > Tb, int N,
                                                    int Nb, int Nlb, double time, double D, int derivative_degree_x,
                                                    int derivative_degree_y) {
    Eigen::VectorXd right_side(Nb);
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vertices[k][0] = P[T[i][k]][0];
            vertices[k][1] = P[T[i][k]][1];
        }
        for (int beta = 0; beta < Nlb; beta++) {
            double temp = gauss2d_integral_test(vertices, beta, derivative_degree_x, derivative_degree_y, time, D);
            right_side(Tb[i][beta]) += temp;
        }
    }
    return right_side;
}
double AssembleMatrix::gauss2d_integral_trial_test(vector <vector<double> > vertices, int basis_index_trial,
                                                   int basis_index_test, int derivative_degree_x,
                                                   int derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];

    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
        quad = quad + gauss_weight[i] * basis_function(x, y, vertices, basis_index_trial, derivative_degree_x, derivative_degree_y)
                                      * basis_function(x, y, vertices, basis_index_test, derivative_degree_x, derivative_degree_y);
    }
    quad = area * quad;
    return quad;
}

double AssembleMatrix::gauss2d_integral_test(vector<vector<double> > vertices, int basis_index,
                                             int derivative_degree_x, int derivative_degree_y,
                                             double time, double D) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];
    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
        quad = quad + gauss_weight[i] * func_f(x, y, time, D) * basis_function(x, y, vertices, basis_index, derivative_degree_x, derivative_degree_y);
    }
    quad = quad * area;
    return quad;
}

double AssembleMatrix::basis_function(double x, double y, vector <vector<double> > vertices, int basis_index,
                                      int derivative_degree_x, int derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];
    double J = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    double result = 0;
    if (derivative_degree_x == 0 && derivative_degree_y == 0) {
        if (basis_index == 0) {
            result = - ((y3 - y1) * (x - x1) - (x3 - x1) * (y - y1)) / J
                     - (-(y2 - y1) * (x - x1) + (x2 - x1) * (y - y1)) / J + 1;
        } else if (basis_index == 1) {
            result = ((y3 - y1) * (x - x1) - (x3 - x1) * (y - y1)) / J;
        } else if (basis_index == 2) {
            result = (-(y2 - y1) * (x - x1) + (x2 - x1) * (y - y1)) / J;
        } else {
        }
    } else if (derivative_degree_x == 1 && derivative_degree_y == 0) {
        if (basis_index == 0) {
            result = (y2 - y3) / J;
        } else if (basis_index == 1) {
            result = (y3 - y1) / J;
        } else if (basis_index == 2) {
            result = -(y2 - y1) / J;
        } else {
        }
    } else if (derivative_degree_x == 0 && derivative_degree_y == 1) {
        if (basis_index == 0) {
            result = (x3 - x2) / J;
        } else if (basis_index == 1) {
            result = - (x3 - x1) / J;
        } else if (basis_index == 2) {
            result = (x2 - x1) / J;
        }
    }
    return result;
}

double AssembleMatrix::func_f(double x, double y, double t, double D) {
    return (1 + 2 * D * t) * sin(x) * cos(y);
}