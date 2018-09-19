#ifndef PERFORMANCEMODELS_H
#define PERFORMANCEMODELS_H

#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include "Traversals.h"

typedef typename Traversals::traversalNames TRAVERSAL;

// Function of (cutoff, density) -> time estimate
typedef const std::function<double(double, double)> Model;

class PerformanceModels {

    // TODO: Maybe read models from a file instead of hard coding them.
    std::map<TRAVERSAL, Model> models = {
            {TRAVERSAL::C04,    [](double cutoff, double density) {
                return 2 * cutoff + 4 * density; //TODO Replace example with real model
            }},
            {TRAVERSAL::C08,    [](double cutoff, double density) {
                return 3 * cutoff + 2 * density; //TODO Replace example with real model
            }},
            {TRAVERSAL::SLICED, [](double cutoff, double density) {
                return 5 * cutoff + 1 * density; //TODO Replace example with real model
            }}
    };

public:

    PerformanceModels() = default;
    ~PerformanceModels() = default;

    /**
     * Predicts the best traversal based on the stored performance models.
     * @param cutoff Current cutoff
     * @param density Current molecule density of the domain/rank
     * @param applicableTraversals List of traversals that should be considered (Defaults to all)
     * @return Prediction for the best traversal
     */
    TRAVERSAL predictBest(double cutoff, double density,
            std::vector<TRAVERSAL> applicableTraversals = Traversals::AllTraversals){

        TRAVERSAL bestTraversal = applicableTraversals.front();
        double bestTime = std::numeric_limits<double>::infinity();

        // Search through applicable traversals and find the one with the shortest predicted time
        for(auto traversal : applicableTraversals){

            auto model = models.at(traversal);
            double time = model(cutoff, density);

            if (time < bestTime){
                bestTraversal = traversal;
                bestTime = time;
            }
        }

        return bestTraversal;

    }

};


#endif //PERFORMANCEMODELS_H
