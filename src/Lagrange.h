#include "Kruskal.h"

#define CORRECTION_FACTOR 0.00001
#define MIN_EPSILON 0.0005
#define QTD_ITERATIONS 30

extern double upperbound;
extern int dimension;

class Lagrange {
    public:
        Lagrange(const vvi &matrix, vector<double> u, ii forb);
        void solve();
        double getCost();
        vii getForbiddenEdges();
        vector<double> getU();
        bool cut();
        vvi* getMatrix();
    private:
        vvi distanceMatrix;
        vvi modifiedMatrix;
        std::vector<int> subgradients;
        std::vector<double> u;
        std::vector<ii> forbiddenEdges;
        double cost;
        int iterations;
        double L;
        double EPSILON;
        bool feasible;
        void modifyMatrix();
        void calculateSubgradients();
        bool isFeasible();
        void calculateU();
        void generateForbiddenEdges();
        vii edges;

        vector<double> best_U;

};