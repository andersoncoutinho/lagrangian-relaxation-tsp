#include <vector>
#include "Lagrange.h"
#include "data.h"

Lagrange::Lagrange(double **ptrMatrix, int dimension) {

    this->subgradients = vector<int>(dimension);
    this->u = vector<double>(dimension);
    this->dimension = dimension;
    
    this->copyMatrix(ptrMatrix);
    

}

void Lagrange::solve() {
    this->calculateSubgradients(this->calculateNodeDegrees());
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

vector<int> Lagrange::calculateNodeDegrees() {
    
    Kruskal kruskal(this->modifiedMatrix);
    double cost = kruskal.MST(this->dimension-1);
    vii edges = kruskal.getEdges();

    printf("Custo: %.0lf\n", cost);
    ii bestCost, bestNodes;
    bestCost.first = bestCost.second = INFINITE;
    for(int j = 0; j < dimension; j++) {
        if(modifiedMatrix[0][j] < bestCost.first) {
            bestCost.second = bestCost.first;
            bestNodes.second = bestNodes.first;
            bestCost.first = modifiedMatrix[0][j];
            bestNodes.first = j;
        } else if (modifiedMatrix[0][j] < bestNodes.second) {
            bestCost.second = modifiedMatrix[0][j]; 
            bestNodes.second = j;
        }
    }
    edges.push_back({0, bestNodes.first});
    edges.push_back({0, bestNodes.second});
    cost += (modifiedMatrix[0][bestNodes.first] + modifiedMatrix[0][bestNodes.second]);
    printf("Custo: %.0lf\n", cost);
    
    vector<int> nodeDegrees(dimension);
    for(int i = 0; i < edges.size(); i++) {
        nodeDegrees[edges[i].first]++;
        nodeDegrees[edges[i].second]++;
    }

    for(int i = 0; i < dimension; i++) {
        printf("Grau do nó %d: %d\n", i, nodeDegrees[i]);
    }

    return nodeDegrees;

}

void Lagrange::calculateSubgradients(vector<int> nodesDegrees) {

    for(int i = 0; i < nodesDegrees.size(); i++) {
        subgradients[i] = 2 - nodesDegrees[i];
        printf("Subgradiente do nó %d: %d\n", i, subgradients[i]);
    }
    
}