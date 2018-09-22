#ifndef PERFORMANCEMODELS_H
#define PERFORMANCEMODELS_H

#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <vector>
#include <math.h>

#include "Traversals.h"

typedef typename Traversals::traversalNames TRAVERSAL;

// Function of (cutoff, density) -> time estimate
typedef const std::function<double(double, double)> Model;

class PerformanceModels {

    // TODO: Maybe read models from a file instead of hard coding them. e.g. parse them with ExprTk
    std::map<TRAVERSAL, Model> models = {
            {TRAVERSAL::C04,    [](double cutoff, double density) {
                return 80.8166 - 262.153 * log2(cutoff) + 1812.23 * (sqrt(density));
            }},
            {TRAVERSAL::C08,    [](double cutoff, double density) {
                return 514.694 - 690.874 * sqrt(log2(cutoff)) + 1883.07 * (sqrt(density));
            }},
            {TRAVERSAL::SLICED, [](double cutoff, double density) {
                return 57.0112 - 250.564 * log2(cutoff) + 1823.25 * (sqrt(density));
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
                          std::vector<TRAVERSAL> applicableTraversals = Traversals::AllTraversals) {

        // Currently: Search minimum as we model time. Change this to maximum if e.g. molsteps/s are modeled

        TRAVERSAL bestTraversal = applicableTraversals.front();
        double bestTime = std::numeric_limits<double>::infinity();

        // Search through applicable traversals and find the one with the shortest predicted time
        for (auto traversal : applicableTraversals) {

            auto model = models.at(traversal);
            double time = model(cutoff, density);

            if (time < bestTime) {
                bestTraversal = traversal;
                bestTime = time;
            }
        }

        return bestTraversal;

    }

};


#endif //PERFORMANCEMODELS_H
