
#ifndef CARBONDIFFUSION_INFOMATRIX_H
#define CARBONDIFFUSION_INFOMATRIX_H

#include <vector>
#include <string>
using std::vector;
using std::string;

class InfoMatrix {
private:
    double left, right, bottom, top;
    double hx, hy;
    string basis_type;

    int Nx, Ny, N;

    vector< vector<double> > CoordinateMatrix_mesh;
    vector< vector<int> > NodeIndicesMatrix_mesh;
    vector< vector<double> > CoordinateMatrix_fe;
    vector< vector<int> > NodeIndicesMatrix_fe;
    vector< vector<int> > boundary_nodes;

    void generateCoordinateMatrix_mesh();
    void generateNodeIndicesMatrix_mesh();
    void generateCoordinateMatrix_fe();
    void generateNodeIndicesMatrix_fe();
    void generateBoundaryNodes();
public:
    InfoMatrix(vector<double> omega, vector<double> h, string basis_type);
    ~InfoMatrix();
    vector< vector<double> > getCoordinateMatrix_mesh();
    vector< vector<int> > getNodeIndicesMatrix_mesh();
    vector< vector<double> > getCoordinateMatrix_fe();
    vector< vector<int> > getNodeIndicesMatrix_fe();
    vector< vector<int> > getBoundaryNode();
};


#endif //CARBONDIFFUSION_INFOMATRIX_H
