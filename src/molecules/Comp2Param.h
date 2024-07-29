#ifndef COMP2PARAM_H_
#define COMP2PARAM_H_

#include <vector>

#include "molecules/Component.h"
#include "molecules/Array2D.h"
#include "molecules/ParaStrm.h"
#include "molecules/mixingrules/MixingRuleBase.h"


/** Comp2Param provides a flexible hash table interface to interaction parameter streams.
 *
 * @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
 */
class Comp2Param {

    using MixRulesType = std::map<int,std::map<int,MixingRuleBase*>>;

    public:
        /** Create a new empty parameter stream. */
        Comp2Param() : m_numcomp(0), m_ssparatbl(0,0) {}

        /** Create a new interaction parameter hash table initialized with
         * the given components and parameters.
         */
        Comp2Param(const std::vector<Component>& components,
				   const MixRulesType & mixcoeff,
				   double epsRF, double rc, double rcLJ) :
			m_numcomp(components.size()), m_ssparatbl(m_numcomp,m_numcomp)
		{
            initialize(components, mixcoeff, epsRF, rc, rcLJ);
        }

        /** Get parameter stream for the interaction between component i and j. */
        ParaStrm& operator()(unsigned int i, unsigned int j) {
            return m_ssparatbl(i,j);
        }

        /** @brief Initialize the parameter streams for each component-component table entry.
         *
         *   The order of the entries must correspond to the
         *   PotForce() function found in potforce.h reading the stream
         */
        void initialize(const std::vector<Component>& components, const MixRulesType & mixcoeff, double epsRF, double rc, double rcLJ);

    private:
        unsigned int m_numcomp;  /**< number of components */

        /** table for parameter streams
         * for each component-component combination (vector<vector<double>* >, row order)
         * and here each site-site combination (vector<double>, row order) store
         * 24*epsilon11, sigma11^2,24*epsilon12, sigma12^2,...
         * absMy1*absMy1,absMy1*absMy2,...
         * 1.5*absMy1*absQ1,1.5*absMy1*absQ2,...
         * 1.5*absQ1*absMy1,1.5*absQ1*absMy2,...
         * 0.75*absQ1*absQ1,0.75*absQ1*absQ2,...
         */
        Array2D<ParaStrm> m_ssparatbl;  /**< table for parameter streams */
};
#endif /* COMP2PARAM_H_ */
