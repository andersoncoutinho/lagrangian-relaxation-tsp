#include <vector>
#include "Lagrange.h"
#include "data.h"

Lagrange::Lagrange(double **ptrMatrix, int dimension, double upperbound) {
    /*
    this->subgradients = vector<int>(dimension);
    this->u = vector<double>(dimension);
    this->dimension = dimension;
    this->upperbound = upperbound;
    
    this->copyMatrix(ptrMatrix);
    */
   this->subgradients = vector<int>(5);
    this->u = vector<double>(5);
    this->dimension = 5;
    this->upperbound = upperbound;
    this->distanceMatrix = vvi(5);
    distanceMatrix[0] = {INFINITE, 30, 26, 50, 40};
    distanceMatrix[1] = {30, INFINITE, 24, 40, 50};
    distanceMatrix[2] = {26, 24, INFINITE, 24, 26};
    distanceMatrix[3] = {50, 40, 24, INFINITE, 30};
    distanceMatrix[4] = {40, 50, 26, 30, INFINITE};
    modifiedMatrix = distanceMatrix;

}

void Lagrange::solve() {

    this->calculateSubgradients(this->calculateNodeDegrees());
    this->calculateU();
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
    cost = kruskal.MST(this->dimension-1);
    vii edges = kruskal.getEdges();

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
    
    vector<int> nodeDegrees(dimension);
    for(int i = 0; i < edges.size(); i++) {
        nodeDegrees[edges[i].first]++;
        nodeDegrees[edges[i].second]++;
    }

    return nodeDegrees;

}

void Lagrange::calculateSubgradients(vector<int> nodesDegrees) {

    for(int i = 0; i < nodesDegrees.size(); i++) {
        subgradients[i] = 2 - nodesDegrees[i];
        printf("Subgradiente [%d]: %d\n", i, subgradients[i]);
    }
}

void Lagrange::calculateU() {

    double powSubgrad = 0;
    for(int i = 0; i < this->subgradients.size(); i++) {
        powSubgrad += (subgradients[i]*subgradients[i]);
    }
    
    double step = ((this->upperbound - this->cost) / powSubgrad);
    printf("PAsso: %.2lf\n", step);
    for(int i = 0; i < u.size(); i++) {
        u[i] += (step * subgradients[i]);
        printf("Multiplicador de lagrangre [%d]: %.3lf\n", i, u[i]);
    }
}