#include "Lagrange.h"
#include <vector>

Lagrange::Lagrange(double **ptrMatrix, int dimension) {

    this->dimension = dimension;
    
    this->copyMatrix(ptrMatrix);

}

void Lagrange::solve() {
    
}

void Lagrange::copyMatrix(double **ptrMatrix) {
    this->distanceMatrix = vvi(this->dimension);
    this->modifiedMatrix = vvi(this->dimension);

    std::vector<double> aux(dimension);
    for(int i = 0; i < this->dimension; i++) {
        for(int j = 0; j < this->dimension; j++) {
            aux[j] = ptrMatrix[i][j];
        }
        distanceMatrix[i] = aux;
    }
    this->modifiedMatrix = this->distanceMatrix;
}