#include "data.h"

int main(int argc, char *argv[]) {

    Data data(argc, argv[1]);
    data.readData();

    double **matrix = data.getMatrixCost();
    int dimension = data.getDimension();;



}