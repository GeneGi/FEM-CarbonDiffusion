#include <iostream>

#include "../include/CarbonDiffusion.h"
using std::cout;
using std::endl;
const double PI = atan(1.0)*4;

int main()
{
    vector<double> omega = {-PI, PI, -PI, PI};
    vector< vector<double> > H = {{1.0/4 * PI, 1.0/4 * PI},
                                 {1.0/8 * PI, 1.0/8 * PI},
                                 {1.0/16 * PI, 1.0/16 * PI},
                                 {1.0/32 * PI, 1.0/32 * PI}};
    int H_num = H.size();
    string basis_type = "linear";
    double t_min = 0;
    double t_max = 3600;
    double dt = 900;
    double Nt = (t_max - t_min) / dt + 1;
    int temperature = 1183;

    CarbonDiffusion solver;
    vector< vector<double> > E(H_num, vector<double>(Nt));
    vector<double> Order(H_num-1);
    for (int i = 0; i < H_num; i++) {
       vector<double> h = H[i];
       solver.solve(omega, h, Nt, dt, basis_type, temperature);
       vector< VectorXd > C = solver.getCarbonConcentration();
       E[i] = solver.getMaxError();
        for (int i = 0; i < E.size(); i++) {
            for (int j = 0; j < E[0].size(); j++) {
                cout << E[i][j] << " ";
            }
            cout << endl;
        }
    }
    for (int k = 0; k < Nt; k++) {
        for (int i = 0; i < H_num-1; i++) {
            Order[i] = log(E[i][k] / E[i+1][k]) / log(H[0][i] / H[0][i+1]);
        }
    }
    for (int i = 0; i < Order.size(); i++) {
        cout << Order[i] << endl;
    }

    return 0;
}