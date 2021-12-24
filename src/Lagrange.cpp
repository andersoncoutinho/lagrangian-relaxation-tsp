#include <vector>
#include "Lagrange.h"
#include "data.h"

double upperbound = 0;
int dimension = 0;

Lagrange::Lagrange(const vvi &matrix, vector<double> u, ii forb) {

    this->subgradients = vector<int>(dimension);
    this->u = u;
    this->iterations = 0;
    this->L = 0;
    this->EPSILON = 2;
    this->feasible = false;
   
    this->distanceMatrix = matrix;
    distanceMatrix[forb.first][forb.second] = distanceMatrix[forb.second][forb.first] = INFINITE;
    this->modifiedMatrix = this->distanceMatrix;
}

void Lagrange::solve() {
    
    while(true) {
        
        this->modifyMatrix();
        this->calculateSubgradients();
        
        if(this->L < this->cost) {
            this->L = this->cost;
            this->best_U = this->u;
        }

        if(this->isFeasible()) {
            this->feasible = true;
            upperbound = this->cost;
            break;
        }

        if(this->cost <= this->L) {
            iterations++;

            if(iterations >= QTD_ITERATIONS) {
                iterations = 0;
                this->EPSILON /= 2;
                if(this->EPSILON <= MIN_EPSILON) {
                    break;
                }
            }
        } else {
            iterations = 0;
        }
        if(upperbound - this->L <= CORRECTION_FACTOR) {
            break;
        }
        
        this->calculateU();
         
    }
    this->u = best_U;
    this->modifyMatrix();
    this->calculateSubgradients();
    this->generateForbiddenEdges();
        
}

void Lagrange::calculateSubgradients() {
    
    Kruskal kruskal(modifiedMatrix);
    cost = kruskal.MST(dimension);
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

    fill(subgradients.begin(), subgradients.end(), 2);
    for(int i = 0; i < edges.size(); i++) {
        subgradients[edges[i].first]--;
        subgradients[edges[i].second]--;
    }
}

void Lagrange::calculateU() {

    int powSubgrad = 0;
    for(int i = 0; i < this->subgradients.size(); i++) {
        powSubgrad += (subgradients[i]*subgradients[i]);
    }
    
    double step = EPSILON * ((upperbound - this->cost) / powSubgrad);

    for(int i = 0; i < u.size(); i++) {
        u[i] += (step * subgradients[i]);
    }
}

void Lagrange::modifyMatrix() {   
    for(int i = 0; i < dimension; i++) {
        double u_i = u[i];
        for(int j = i+1; j < dimension; j++) {
            if(i!=j) {
                this->modifiedMatrix[i][j] = this->distanceMatrix[i][j] - (u_i + u[j]);
                this->modifiedMatrix[j][i] = this->distanceMatrix[i][j] - (u_i + u[j]);
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

void Lagrange::generateForbiddenEdges() {
    int index = -1;
    int degree = 0;
    for(int i = 0; i < this->subgradients.size(); i++) {
        if(this->subgradients[i] < degree) {
            degree = this->subgradients[i];
            index = i;
        }
    }

    for(int i = 0; i < edges.size(); i++) {
        if(edges[i].first == index || edges[i].second == index) {
            this->forbiddenEdges.push_back(edges[i]);
        }
    }
}

bool Lagrange::cut() {
    return this->feasible;
}

vii Lagrange::getForbiddenEdges() {
    return this->forbiddenEdges;
}

vector<double> Lagrange::getU() {
    return this->u;
}

vvi* Lagrange::getMatrix() {
    return &this->distanceMatrix;
}