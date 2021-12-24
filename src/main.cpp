#include <stdio.h>
#include <list>
#include <chrono>
#include "data.h"
#include "Lagrange.h"

double primal_bound_calc(vvi cost);
int node_degree(const vii & edges, int node);
bool cycle_search(vii edges, int dimension);
vii tree_node_children(vii edges, int parent);
bool cycle_search_path(vii edges, int v, bool visited[], int parent);

int main(int argc, char *argv[]) {

    Data data(argc, argv[1]);
    data.readData();

    double **matrix = data.getMatrixCost();
    dimension = data.getDimension();

    vvi distanceMatrix(dimension, std::vector<double>(dimension));  
    for(int i = 0; i < dimension; i++) {
        for(int j = 0; j < dimension; j++) {
            distanceMatrix[i][j] = matrix[i][j];
        }
    }
    
    double bestCost = INFINITE;
    
    std::cout << "Informe o limitante primal: ";
    std::cin >> upperbound;
    
    printf("Limitante primal: %.2lf\n", upperbound);

    auto t1 = std::chrono::high_resolution_clock::now();
    Lagrange root(distanceMatrix, vector<double>(dimension), ii());

    std::list<Lagrange> tree;
    tree.push_front(root);

    while(!tree.empty()) {
        Lagrange node = tree.front();
        tree.pop_front();
    
        node.solve();
                    
        if(!node.cut() && node.getCost() < upperbound) {
            vii forbidden = node.getForbiddenEdges();
            for(int i = 0; i < forbidden.size(); i++) {
                
                Lagrange newNode(*node.getMatrix(), node.getU(), forbidden[i]);
                tree.push_back(newNode);
            }
        } else if(node.getCost() < bestCost) {
            bestCost = node.getCost();
        }

    }
    
    auto t2 = std::chrono::high_resolution_clock::now();
    printf("Cost: %.3lf\n", bestCost);
    auto exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    printf("TIME: %.10lf\n", (double)exec_time / 10e2);
    
    return 0;
}