#include <vector>
#include "Lagrange.h"
#include "data.h"

Lagrange::Lagrange(vvi *matrix, int dimension, double upperbound) {

    this->subgradients = vector<int>(dimension);
    this->u = vector<double>(dimension);
    this->dimension = dimension;
    this->upperbound = upperbound;
    this->iterations = 0;
    this->L = 0;
    this->EPSILON = 2;
    this->feasible = false;
    this->nodesDegree = vector<int>(dimension);
   
    this->distanceMatrix = *matrix;
    for(int i = 0; i < this->forbiddenEdges.size(); i++) {
        this->distanceMatrix[forbiddenEdges[i].first][forbiddenEdges[i].second] = 
        this->distanceMatrix[forbiddenEdges[i].second][forbiddenEdges[i].first] = INFINITE;
    }

    this->modifiedMatrix = this->distanceMatrix;
}

void Lagrange::solve() {
    
    while(true) {

        this->calculateNodesDegree();
        this->calculateSubgradients();

        if(this->isFeasible()) {
            this->feasible = true;
            this->upperbound = this->cost;
            break;
        }

        if(this->cost >= this->L) {
            this->L = this->cost;
            iterations = -1;
        } else if(iterations == QTD_ITERATIONS) {
            this->EPSILON /= 2;
            if(this->EPSILON < MIN_EPSILON) {
                this->generateForbiddenEdges();
                break;
            }
            iterations = -1;
        }
        
        this->calculateU();
        this->modifyMatrix();

        iterations++;
    }
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

void Lagrange::calculateNodesDegree() {
    
    Kruskal kruskal(modifiedMatrix);
    cost = kruskal.MST(dimension-1);
    edges = kruskal.getEdges();
    
    ii bestCost, bestNodes;
    bestCost.first = bestCost.second = INFINITE;
    for(int j = 1; j < dimension; j++) {
        if(modifiedMatrix[0][j] < bestCost.first) {
            bestCost.second = bestCost.first;
            bestNodes.second = bestNodes.first;
            bestCost.first = modifiedMatrix[0][j];
            bestNodes.first = j;
        } else if (modifiedMatrix[0][j] < bestCost.second) {
            bestCost.second = modifiedMatrix[0][j]; 
            bestNodes.second = j;
        }
    }
    edges.push_back({0, bestNodes.first});
    edges.push_back({0, bestNodes.second});
    cost += (modifiedMatrix[0][bestNodes.first] + modifiedMatrix[0][bestNodes.second]);
    double sum = 0;
    for(int i = 0; i < u.size(); i++) {
        sum += u[i];
    } 
    cost += (2 * sum);

    fill(nodesDegree.begin(), nodesDegree.end(), 0);
    for(int i = 0; i < edges.size(); i++) {
        nodesDegree[edges[i].first]++;
        nodesDegree[edges[i].second]++;
    }
}

void Lagrange::calculateSubgradients() {
    for(int i = 0; i < nodesDegree.size(); i++) {
        subgradients[i] = 2 - nodesDegree[i];
    }
}

void Lagrange::calculateU() {

    double powSubgrad = 0;
    for(int i = 0; i < this->subgradients.size(); i++) {
        powSubgrad += (subgradients[i]*subgradients[i]);
    }
    
    double step = EPSILON * ((this->upperbound - this->cost) / powSubgrad);

    for(int i = 0; i < u.size(); i++) {
        u[i] += (step * subgradients[i]);
    }
}

void Lagrange::modifyMatrix() {   
    for(int i = 0; i < this->dimension; i++) {
        for(int j = 0; j < this->dimension; j++) {
            if(i!=j) {
                this->modifiedMatrix[i][j] = this->distanceMatrix[i][j] - (u[i] + u[j]);
            }
        }
    }
    
}

bool Lagrange::isFeasible() {
    for(int i = 0; i < subgradients.size(); i++) {
        if(subgradients[i] != 0) {
            return false;
        }
    }
    return true;
}

double Lagrange::getCost() {
    return this->cost;
}

double Lagrange::getUpperbound() {
    return this->upperbound;
}

void Lagrange::generateForbiddenEdges() {
    int index = -1;
    int degree = 0;
    for(int i = 0; i < this->nodesDegree.size(); i++) {
        if(this->nodesDegree[i] > degree) {
            degree = this->nodesDegree[i];
            index = i;
        }
    }

    for(int i = 0; i < edges.size(); i++) {
        if(edges[i].first == index || edges[i].second == index) {
            this->forbiddenEdges.push_back(edges[i]);
        }
    }
}