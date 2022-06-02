/*
 * ExtendedProfileSamplingTest.cpp
 *
 *  Created on: Jun 2022
 *      Author: homes
 */

#include "ExtendedProfileSamplingTest.h"

TEST_SUITE_REGISTRATION(ExtendedProfileSamplingTest);

ExtendedProfileSamplingTest::ExtendedProfileSamplingTest() {}

ExtendedProfileSamplingTest::~ExtendedProfileSamplingTest() {}

void ExtendedProfileSamplingTest::testEPSampling() {

    double delta = 1e-6;  // Tolerate deviation between expected and actual value

    // Expected temperature
    std::array<double, 20> temp_expected = {1.7848196, 1.8002271, 1.8008320, 1.76791660, 1.7873646,
                                            1.7978649, 1.7972279, 1.7947228, 1.79220220, 1.7932661,
                                            1.7981219, 1.7815903, 1.7938292, 1.79564790, 1.7968447,
                                            1.7974076, 1.7866776, 1.7894757, 1.78767060, 1.7948688};
    // Expected density
    std::array<double, 20> rho_expected  = {0.49, 0.5200, 0.4875, 0.5075, 0.4650,
                                            0.50, 0.5075, 0.4650, 0.4850, 0.4875,
                                            0.53, 0.4900, 0.5075, 0.5375, 0.4975,
                                            0.51, 0.4925, 0.5200, 0.5075, 0.4925};
    // Expected kinetic energy
    std::array<double, 20> ekin_expected = {2.6885133, 2.7107971, 2.7045695, 2.7001697, 2.7007511,
                                            2.6991446, 2.6999478, 2.6984594, 2.7011729, 2.6954951,
                                            2.7040797, 2.6960042, 2.7014050, 2.7001036, 2.7017692,
                                            2.7000050, 2.6986238, 2.6989247, 2.7000430, 2.7013787};


    const char* filename = "LJTS2-5_equilibrated.inp";
    double cutoff = 2.5;
    std::unique_ptr<ParticleContainer> container{
        initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff)};

    std::unique_ptr<ExtendedProfileSampling> plugin {new ExtendedProfileSampling()};

    plugin->init(container.get(), _domainDecomposition, _domain);
    plugin->afterForces(container.get(), _domainDecomposition, 1);

    // Only root knows correct values; see plugin for explanation
    if (_domainDecomposition->getRank() == 0) {
        // Default binwidth is 1 and global domain length of test system is 20; therefore 20 bins in total
        for (unsigned int i = 0; i < 20; i++) {

            const double temp_actual_all = plugin->getQuantity(_domainDecomposition, "T", i);
            const double rho_actual_all  = plugin->getQuantity(_domainDecomposition, "rho", i);
            const double ekin_actual_all = plugin->getQuantity(_domainDecomposition, "ekin", i);
            ASSERT_DOUBLES_EQUAL_MSG("Temperature at index "+std::to_string(i)+" not as expected", temp_expected[i], temp_actual_all, delta);
            ASSERT_DOUBLES_EQUAL_MSG("Temperature at index "+std::to_string(i)+" not as expected", temp_expected[i], temp_actual_all, delta);
            ASSERT_DOUBLES_EQUAL_MSG("Temperature at index "+std::to_string(i)+" not as expected", temp_expected[i], temp_actual_all, delta);

            // As there is only one component in the test system, the componentwise and overall values are supposed to be equal
            const double temp_actual_cid1 = plugin->getQuantity(_domainDecomposition, "T", i+20u);
            const double rho_actual_cid1  = plugin->getQuantity(_domainDecomposition, "rho", i+20u);
            const double ekin_actual_cid1 = plugin->getQuantity(_domainDecomposition, "ekin", i+20u);
            ASSERT_DOUBLES_EQUAL(temp_actual_all, temp_actual_cid1, delta);
            ASSERT_DOUBLES_EQUAL(rho_actual_all, rho_actual_cid1, delta);
            ASSERT_DOUBLES_EQUAL(ekin_actual_all, ekin_actual_cid1, delta);
        }
        cout << endl;
    }

}
