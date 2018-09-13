#ifndef TRAVERSALS_H
#define TRAVERSALS_H

#include <vector>

struct Traversals {
    enum traversalNames {
        ORIGINAL = 0,
        C08 = 1,
        C04 = 2,
        SLICED = 3,
        HS = 4,
        MP = 5,
        QSCHED = 6,
    };

    static const std::vector<traversalNames> AllTraversals;
};

#endif //TRAVERSALS_H
