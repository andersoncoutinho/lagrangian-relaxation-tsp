#include "Kruskal.h"

class Lagrange {
    public:
        Lagrange(double **ptrMatrix, int dimension);
        void solve();
    private:
        vvi distanceMatrix;
        vvi modifiedMatrix;
        void copyMatrix(double **ptrMatrix);

};