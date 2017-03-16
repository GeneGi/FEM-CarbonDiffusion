#include "../include/InfoMatrix.h"

InfoMatrix::InfoMatrix(vector<double> omega, vector<double> h, string basis_type) {
    this->left = omega[0];
    this->right = omega[1];
    this->bottom = omega[2];
    this->top = omega[3];
    this->hx = h[0];
    this->hy = h[1];
    this->basis_type = basis_type;

    Nx = (int) ((right - left) / hx);
    Ny = (int) ((top - bottom) / hy);
    N = (Nx+1) * (Ny+1);

    generateCoordinateMatrix_mesh();
    generateNodeIndicesMatrix_mesh();
    generateCoordinateMatrix_fe();
    generateNodeIndicesMatrix_fe();
    generateBoundaryNodes();
}
InfoMatrix::~InfoMatrix() {}

vector< vector<double> > InfoMatrix::getCoordinateMatrix_mesh() {
    return CoordinateMatrix_mesh;
}

vector< vector<int> > InfoMatrix::getNodeIndicesMatrix_mesh() {
    return NodeIndicesMatrix_mesh;
}

vector< vector<double> > InfoMatrix::getCoordinateMatrix_fe() {
    return CoordinateMatrix_fe;
}

vector< vector<int> > InfoMatrix::getNodeIndicesMatrix_fe() {
    return NodeIndicesMatrix_fe;
}

vector< vector<int> > InfoMatrix::getBoundaryNode(){
    return boundary_nodes;
}

void InfoMatrix::generateCoordinateMatrix_mesh() {
    CoordinateMatrix_mesh.resize(N);
    for (int i = 0; i < N; i++) {
        CoordinateMatrix_mesh[i].resize(2);
    }

    for (int i = 0; i < Nx+1; i++) {
        for (int j = 0; j < Ny+1; j++) {
            CoordinateMatrix_mesh[i*(Ny+1) + j][0] = left + i * hx;
            CoordinateMatrix_mesh[i*(Ny+1) + j][1] = bottom + j * hy;
        }
    }
}

void InfoMatrix::generateNodeIndicesMatrix_mesh() {
    NodeIndicesMatrix_mesh.resize((unsigned long) (2 * Nx * Ny));
    for (int i = 0; i < 2*Nx*Ny; i++) {
        NodeIndicesMatrix_mesh[i].resize(3);
    }

    vector< vector<int> > temp_NodeIndices;
    temp_NodeIndices.resize((unsigned long) (Nx + 1));
    for (int i = 0; i < Nx+1; i++) {
        temp_NodeIndices[i].resize((unsigned long) (Ny + 1));
    }
    for (int i = 0; i < Nx+1; i++) {
        for (int j = 0; j < Ny+1; j++) {
            temp_NodeIndices[i][j] = i * (Ny+1) + j;
        }
    }
    int row, column;
    for (int i = 0; i < Nx*Ny; i++) {
        if ((i+1) % Ny == 0) {
            row = Ny - 1;
            column = i / Ny;
        } else {
            row = i % Ny;
            column = i / Ny;
        }
        NodeIndicesMatrix_mesh[2*i][0]= temp_NodeIndices[column][row];
        NodeIndicesMatrix_mesh[2*i][1] = temp_NodeIndices[column+1][row];
        NodeIndicesMatrix_mesh[2*i][2] = temp_NodeIndices[column][row+1];
        NodeIndicesMatrix_mesh[2*i+1][0] = temp_NodeIndices[column][row+1];
        NodeIndicesMatrix_mesh[2*i+1][1] = temp_NodeIndices[column+1][row];
        NodeIndicesMatrix_mesh[2*i+1][2] = temp_NodeIndices[column+1][row+1];
    }
}

void InfoMatrix::generateCoordinateMatrix_fe() {
    if (basis_type == "linear") {
        CoordinateMatrix_fe = CoordinateMatrix_mesh;
    } else if (basis_type == "quadratic") {
        //TODO basis type quadratic
    }
}

void InfoMatrix::generateNodeIndicesMatrix_fe() {
    if (basis_type == "linear") {
        NodeIndicesMatrix_fe = NodeIndicesMatrix_mesh;
    } else if (basis_type == "quadratic") {
        //TODO basis type quadratic
    }
}

void InfoMatrix::generateBoundaryNodes() {
    if (basis_type == "quadratic") {
        Nx = 2 * Nx;
        Ny = 2 * Ny;
    }
    int nbn = 2 * (Nx+Ny);
    boundary_nodes.resize(nbn);
    for (int i = 0; i < nbn; i++) {
        boundary_nodes[i].resize(2);
    }
    for (int i = 0; i < nbn; i++) {
        boundary_nodes[i][0] = -1;
    }
    // left side
    for (int i = 0; i < Nx; i++) {
        boundary_nodes[i][1] = i * (Ny+1);
    }
    // top side
    for (int i = Nx; i < Nx+Ny; i++) {
        boundary_nodes[i][1] = Nx * Ny + i;
    }
    // right side
    for (int i = Nx+Ny; i < 2*Nx+Ny; i++) {
        boundary_nodes[i][1] = (2*Nx + Ny + 1 - i) * (Ny + 1) - 1;
    }
    // bottom side
    for (int i = 2*Nx+Ny; i < nbn; i++) {
        boundary_nodes[i][1] = 2*Nx + 2*Ny - i;
    }
}