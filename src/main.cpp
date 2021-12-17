#include <stdio.h>
#include "data.h"
#include "Lagrange.h"

int main(int argc, char *argv[]) {

    Data data(argc, argv[1]);
    data.readData();

    double **matrix = data.getMatrixCost();
    int dimension = data.getDimension();

    vvi distanceMatrix(dimension);
    std::vector<double> aux(dimension);
    for(int i = 0; i < dimension; i++) {
        for(int j = 0; j < dimension; j++) {
            aux[j] = matrix[i][j];
        }
        distanceMatrix[i] = aux;
    }

    double upperbound;
    printf("Informe o limitante primal: ");
    int n = scanf("%lf", &upperbound);
    getchar();

    Lagrange lagrange(&distanceMatrix, dimension, upperbound);
    lagrange.solve();
    
    printf("Custo: %.2lf\n", lagrange.getCost());
    printf("Upperbound: %.2lf\n", lagrange.getUpperbound());

    return 0;
}