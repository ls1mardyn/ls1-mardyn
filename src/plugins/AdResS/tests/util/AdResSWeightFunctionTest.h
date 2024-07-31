//
// Created by alex on 7/26/24.
//

#ifndef MARDYN_ADRESSWEIGHTFUNCTIONTEST_H
#define MARDYN_ADRESSWEIGHTFUNCTIONTEST_H

#include "utils/Testing.h"
#include "plugins/AdResS/util/WeightFunction.h"

class AdResSWeightFunctionTest : public utils::Test {
TEST_SUITE(AdResSWeightFunctionTest);
    TEST_METHOD(testBounds);
    TEST_METHOD(testCenter);
    TEST_METHOD(testOther);
TEST_SUITE_END;
public:
    AdResSWeightFunctionTest();
    virtual ~AdResSWeightFunctionTest() = default;

    void testBounds();
    void testCenter();
    void testOther();
    // TODO add test for vectorized weight
private:
    using d3 = std::array<double, 3>;
    using test_points = std::array<d3, 5>;
    void checkBounds(Weight::function_t fun);
    void checkCenter(Weight::function_t fun);
    void checkOther(Weight::function_t fun);
    Resolution::FPRegion _region;
    // FP inter C inter CG
    const std::array<test_points , 26> _test_sets {
        // single dirs
        test_points {d3 {2.1, 2.5, 2.5}, d3 {1.75, 2.5, 2.5}, d3 {1.5, 2.5, 2.5}, d3 {1.25, 2.5, 2.5}, d3 {0.5, 2.5, 2.5}}, // Left
        test_points {d3 {2.9, 2.5, 2.5}, d3 {3.25, 2.5, 2.5}, d3 {3.5, 2.5, 2.5}, d3 {3.75, 2.5, 2.5}, d3 {4.5, 2.5, 2.5}}, // Right
        test_points {d3 {2.5, 2.1, 2.5}, d3 {2.5, 1.75, 2.5}, d3 {2.5, 1.5, 2.5}, d3 {2.5, 1.25, 2.5}, d3 {2.5, 0.5, 2.5}}, // Down
        test_points {d3 {2.5, 2.9, 2.5}, d3 {2.5, 3.25, 2.5}, d3 {2.5, 3.5, 2.5}, d3 {2.5, 3.75, 2.5}, d3 {2.5, 4.5, 2.5}}, // Up
        test_points {d3 {2.5, 2.5, 2.1}, d3 {2.5, 2.5, 1.75}, d3 {2.5, 2.5, 1.5}, d3 {2.5, 2.5, 1.25}, d3 {2.5, 2.5, 0.5}}, // Front
        test_points {d3 {2.5, 2.5, 2.9}, d3 {2.5, 2.5, 3.25}, d3 {2.5, 2.5, 3.5}, d3 {2.5, 2.5, 3.75}, d3 {2.5, 2.5, 4.5}}, // Back
        // x and y
        test_points {d3 {2.1, 2.1, 2.5}, d3 {1.75, 1.75, 2.5}, d3 {1.5, 1.5, 2.5}, d3 {1.25, 1.25, 2.5}, d3 {0.5, 0.5, 2.5}}, // Left-Down
        test_points {d3 {2.9, 2.1, 2.5}, d3 {3.25, 1.75, 2.5}, d3 {3.5, 1.5, 2.5}, d3 {3.75, 1.25, 2.5}, d3 {4.5, 0.5, 2.5}}, // Right-Down
        test_points {d3 {2.1, 2.9, 2.5}, d3 {1.75, 3.25, 2.5}, d3 {1.5, 3.5, 2.5}, d3 {1.25, 3.75, 2.5}, d3 {0.5, 4.5, 2.5}}, // Left-Up
        test_points {d3 {2.9, 2.9, 2.5}, d3 {3.25, 3.25, 2.5}, d3 {3.5, 3.5, 2.5}, d3 {3.75, 3.75, 2.5}, d3 {4.5, 4.5, 2.5}}, // Right-Up
        // x and z
        test_points {d3 {2.1, 2.5, 2.1}, d3 {1.75, 2.5, 1.75}, d3 {1.5, 2.5, 1.5}, d3 {1.25, 2.5, 1.25}, d3 {0.5, 2.5, 0.5}}, // Left-Front
        test_points {d3 {2.9, 2.5, 2.1}, d3 {3.25, 2.5, 1.75}, d3 {3.5, 2.5, 1.5}, d3 {3.75, 2.5, 1.25}, d3 {4.5, 2.5, 0.5}}, // Right-Front
        test_points {d3 {2.1, 2.5, 2.9}, d3 {1.75, 2.5, 3.25}, d3 {1.5, 2.5, 3.5}, d3 {1.25, 2.5, 3.75}, d3 {0.5, 2.5, 4.5}}, // Left-Back
        test_points {d3 {2.9, 2.5, 2.9}, d3 {3.25, 2.5, 3.25}, d3 {3.5, 2.5, 3.5}, d3 {3.75, 2.5, 3.75}, d3 {4.5, 2.5, 4.5}}, // Right-Back
        // y and z
        test_points {d3 {2.5, 2.1, 2.1}, d3 {2.5, 1.75, 1.75}, d3 {2.5, 1.5, 1.5}, d3 {2.5, 1.25, 1.25}, d3 {2.5, 0.5, 0.5}}, // Down-Front
        test_points {d3 {2.5, 2.9, 2.1}, d3 {2.5, 3.25, 1.75}, d3 {2.5, 3.5, 1.5}, d3 {2.5, 3.75, 1.25}, d3 {2.5, 4.5, 0.5}}, // Up-Front
        test_points {d3 {2.5, 2.1, 2.9}, d3 {2.5, 1.75, 3.25}, d3 {2.5, 1.5, 3.5}, d3 {2.5, 1.25, 3.75}, d3 {2.5, 0.5, 4.5}}, // Down-Back
        test_points {d3 {2.5, 2.9, 2.9}, d3 {2.5, 3.25, 3.25}, d3 {2.5, 3.5, 3.5}, d3 {2.5, 3.75, 3.75}, d3 {2.5, 4.5, 4.5}}, // Up-Back
        // all dirs
        test_points {d3 {2.1, 2.1, 2.1}, d3 {1.75, 1.75, 1.75}, d3 {1.5, 1.5, 1.5}, d3 {1.25, 1.25, 1.25}, d3 {0.5, 0.5, 0.5}}, // Left-Down-Front
        test_points {d3 {2.9, 2.1, 2.1}, d3 {3.25, 1.75, 1.75}, d3 {3.5, 1.5, 1.5}, d3 {3.75, 1.25, 1.25}, d3 {4.5, 0.5, 0.5}}, // Right-Down-Front
        test_points {d3 {2.1, 2.9, 2.1}, d3 {1.75, 3.25, 1.75}, d3 {1.5, 3.5, 1.5}, d3 {1.25, 3.75, 1.25}, d3 {0.5, 4.5, 0.5}}, // Left-Up-Front
        test_points {d3 {2.9, 2.9, 2.1}, d3 {3.25, 3.25, 1.75}, d3 {3.5, 3.5, 1.5}, d3 {3.75, 3.75, 1.25}, d3 {4.5, 4.5, 0.5}}, // Right-Up-Front
        test_points {d3 {2.1, 2.1, 2.9}, d3 {1.75, 1.75, 3.25}, d3 {1.5, 1.5, 3.5}, d3 {1.25, 1.25, 3.75}, d3 {0.5, 0.5, 4.5}}, // Left-Down-Back
        test_points {d3 {2.9, 2.1, 2.9}, d3 {3.25, 1.75, 3.25}, d3 {3.5, 1.5, 3.5}, d3 {3.75, 1.25, 3.75}, d3 {4.5, 0.5, 4.5}}, // Right-Down-Back
        test_points {d3 {2.1, 2.9, 2.9}, d3 {1.75, 3.25, 3.25}, d3 {1.5, 3.5, 3.5}, d3 {1.25, 3.75, 3.75}, d3 {0.5, 4.5, 4.5}}, // Left-Up-Back
        test_points {d3 {2.9, 2.9, 2.9}, d3 {3.25, 3.25, 3.25}, d3 {3.5, 3.5, 3.5}, d3 {3.75, 3.75, 3.75}, d3 {4.5, 4.5, 4.5}}  // Right-Up-Back
    };
};
#endif //MARDYN_ADRESSWEIGHTFUNCTIONTEST_H
