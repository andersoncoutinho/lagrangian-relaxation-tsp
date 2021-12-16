#include <stdio.h>
#include "data.h"
#include "Lagrange.h"

int main(int argc, char *argv[]) {

    Data data(argc, argv[1]);
    data.readData();

    double **matrix = data.getMatrixCost();
    int dimension = data.getDimension();

    double upperbound;
    printf("Informe o limitante primal: ");
    int n = scanf("%lf", &upperbound);
    getchar();

    Lagrange lagrange(matrix, dimension, upperbound);
    lagrange.solve();
    
    printf("Custo: %.2lf\n", lagrange.getCost());
    printf("Upperbound: %.2lf\n", lagrange.getUpperbound());

    return 0;
}