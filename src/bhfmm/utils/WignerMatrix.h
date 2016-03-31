/*
 * WignerMatrix.h
 *
 *  Created on: Dec 11, 2014
 *      Author: uwe
 */

#ifndef WIGNERMATRIX_H_
#define WIGNERMATRIX_H_

#include "CuboidPyramidalMatrix.h"

namespace bhfmm {

class WignerMatrix {
public:
	/* constructor */
	WignerMatrix() : order(0), mat()  {}

	WignerMatrix(unsigned ord, bool clean = false) :
		order(ord), mat(ord + 1, clean)  {}

	/* destructor */
	~WignerMatrix() {}


	void clear() {mat.clear();}
	unsigned get_num_entries() const { return mat.get_num_entries();}
	unsigned get_order() const {return order;}

	void evaluate(double theta);

	/* prints all matrix elements to order maxl */
    void print(int maxl);

    //const accesssor
    inline double acc_c(unsigned l, unsigned m, int k) const {return mat.access_const(l,m,k);}

    //direct accessor
    inline double & acc(unsigned l, unsigned m, int k) {return mat.access(l,m,k);}

private:

    //sequential accessor
    inline double & acc_seq(unsigned l) {return mat.access_seq(l);}

    //sequential const accessor
    inline double acc_c_seq(unsigned l) const {return mat.access_seq_const(l);}

    void evaluate_0_Pip2(double theta);

    void mirror_k();

    template <class PowPair>
    void apply_minus_one_pow();


	unsigned order;
	CuboidPyramidalMatrix mat;

	// helper classes to avoid code duplication
	class l_m {
	public:
		inline static int sum(int l, int m, int k) {
			return l + m;
		}
	};

	class l_k {
	public:
		inline static int sum(int l, int m, int k) {
			return l + k;
		}
	};

	class m_k {
	public:
		inline static int sum(int l, int m, int k) {
			return m + k;
		}
	};

};

} // namespace bhfmm


#endif /* WIGNERMATRIX_H_ */
