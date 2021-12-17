#include <stdio.h>
#include <list>
#include <chrono>
#include "data.h"
#include "Lagrange.h"

int main(int argc, char *argv[]) {

    Data data(argc, argv[1]);
    data.readData();

    double **matrix = data.getMatrixCost();
    int dimension = data.getDimension();

    vvi distanceMatrix(dimension);
    std::vector<double> aux(dimension);
    for(int i = 0; i < dimension; i++) {
        for(int j = 0; j < dimension; j++) {
            aux[j] = matrix[i][j];
        }
        distanceMatrix[i] = aux;
    }
    double bestCost = INFINITE;
    double upperbound;
    printf("Informe o limitante primal: ");
    int n = scanf("%lf", &upperbound);
    getchar();

    auto t1 = std::chrono::high_resolution_clock::now();

    Lagrange root(&distanceMatrix, dimension, upperbound, vector<double>(dimension), ii());
    root.solve();

    std::list<Lagrange> tree;
    tree.push_front(root);

    while(!tree.empty()) {
        Lagrange node = tree.front();
        tree.pop_front();

                    
        if(!node.cut() && node.getCost() < upperbound) {
            vii forbidden = node.getForbiddenEdges();
            for(int i = 0; i < forbidden.size(); i++) {
                
                Lagrange newNode(node.getMatrix(), dimension, node.getUpperbound(), node.getU(), forbidden[i]);
                newNode.solve();
                if(newNode.cut()) {
                    if(newNode.getCost() < bestCost) {
                        upperbound = newNode.getCost();
                        bestCost = newNode.getCost();
                    }
                } else {
                    if(newNode.getCost() < upperbound) {
                        tree.push_front(newNode);
                    }
                }
            }
        } else if(node.getCost() < bestCost) {
            bestCost = node.getCost();
        }

    }
    auto t2 = std::chrono::high_resolution_clock::now();
    printf("Cost: %.3lf\n", bestCost);
    auto exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "TIME: " << (double)exec_time / 10e2  << std::endl;
    return 0;
}