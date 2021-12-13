#include "Kruskal.h"

class Lagrange {
    public:
        Lagrange(double **ptrMatrix, int dimension);
        void solve();
    private:
        vvi distanceMatrix;
        vvi modifiedMatrix;
        std::vector<int> subgradients;
        std::vector<double> u;
        int dimension;
        void copyMatrix(double **ptrMatrix);
        void calculateNodeDegrees();
        void calculateSubgradients();

};