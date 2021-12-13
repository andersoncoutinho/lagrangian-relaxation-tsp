#include "Kruskal.h"

class Lagrange {
    public:
        Lagrange(double **ptrMatrix, int dimension, double upperbound);
        void solve();
    private:
        vvi distanceMatrix;
        vvi modifiedMatrix;
        std::vector<int> subgradients;
        std::vector<double> u;
        double cost;
        double upperbound;
        int dimension;
        void copyMatrix(double **ptrMatrix);
        vector<int> calculateNodeDegrees();
        void calculateSubgradients(vector<int> nodesDegrees);
        void calculateU();

};