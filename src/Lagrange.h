#include "Kruskal.h"

#define CORRECTION_FACTOR 0.000001
#define MIN_EPSILON 0.0005
#define QTD_ITERATIONS 30

class Lagrange {
    public:
        Lagrange(double **ptrMatrix, int dimension, double upperbound);
        void solve();
        double getCost();
        double getUpperbound();
    private:
        vvi distanceMatrix;
        vvi modifiedMatrix;
        std::vector<int> subgradients;
        std::vector<double> u;
        std::vector<ii> forbiddenEdges;
        std::vector<int> nodesDegree;
        double cost;
        double upperbound;
        int dimension;
        int iterations;
        double L;
        double EPSILON;
        bool feasible;
        void copyMatrix(double **ptrMatrix);
        void calculateNodesDegree();
        void calculateSubgradients();
        bool isFeasible();
        void calculateU();
        void modifyMatrix();
        void generateForbiddenEdges();
        vii edges;

};