//
// Created by alex on 04.01.24.
//

#ifndef MARDYN_INTERPOLATION_H
#define MARDYN_INTERPOLATION_H

#include "utils/mardyn_assert.h"

#include <vector>
#include <string>
#include <array>
#include <complex>
#include <algorithm>

namespace Interpolation {
	/**
	* Represents one interpolated function using Cubic Hermite Splines
	* */
	struct Function {
		//! @brief number of knots
		unsigned long n;
		//! @brief distance between each knot
		std::vector<double> step_width;
		//! @brief starting point of samples, x_0
		double begin;
		//! @brief samples of f(x)
		std::vector<double> function_values;
		//! @brief samples of f'(x)
		std::vector<double> gradients;

		//! @brief write this function to file in XML format
		void writeXML(const std::string &filename) const;

		//! @brief write this function to file in text format
		void writeTXT(const std::string &filename) const;

		//! @brief load function from file
		void loadTXT(const std::string &filename);
	};

	// TODO
	struct Function3D {
	};

	/**
	 * Matrix NxM N rows, M columns
	 * Row-First Memory alignment
	 * */
	struct Matrix {
		Matrix() : _dim0(0), _dim1(0) {}

		Matrix(unsigned long dim0, unsigned long dim1) : _dim0(dim0), _dim1(dim1) {
			_data.resize(dim0 * dim1, 0.0);
		}

		/**
		 * @param vec single column
		 * @param repeat how often vec should be repeated horizontally
		 * */
		Matrix(const std::vector<double> &vec, unsigned long repeat) : Matrix(vec.size(), repeat) {
#if defined(_OPENMP)
#pragma omp parallel for collapse(2)
#endif
			for (unsigned long n = 0; n < _dim0; n++) {
				for (unsigned long m = 0; m < _dim1; m++) {
					setAt(n, m, vec[n]);
				}
			}
		}

		/**
		 * @param vec single row
		 * @param repeat how often vec should be repeated vertically
		 * */
		Matrix(unsigned long repeat, const std::vector<double> &vec) : Matrix(repeat, vec.size()) {
#if defined(_OPENMP)
#pragma omp parallel for collapse(2)
#endif
			for (unsigned long n = 0; n < _dim0; n++) {
				for (unsigned long m = 0; m < _dim1; m++) {
					setAt(n, m, vec[m]);
				}
			}
		}

		std::vector<double> operator*(const std::vector<double> &vec) const {
			mardyn_assert((vec.size() == _dim1));
			std::vector<double> result;
			result.resize(_dim0, 0.0);
			auto *raw_result = std::data(result);

#if defined(_OPENMP)
#pragma omp parallel for simd reduction(+:raw_result[:_dim0]) collapse(2)
#endif
			for (unsigned long n = 0; n < _dim0; n++) {
				for (unsigned long m = 0; m < _dim1; m++) {
					raw_result[n] += vec[m] * _data[n * _dim1 + m];
				}
			}
			return result;
		}

		void setAt(unsigned long n, unsigned long m, double value) {
			_data[n * _dim1 + m] = value;
		}

		double getAt(unsigned long n, unsigned long m) {
			return _data[n * _dim1 + m];
		}

		unsigned long _dim0;
		unsigned long _dim1;
		std::vector<double> _data;
	};

	// TODO
	struct Tensor {
	};

	/**
	 * Numerically computes the gradient of the input vector. Uses finite difference coefficients (based on Lagrange Polynomials).
	 * Input and output have equal size. The output on the borders use forward or backward differences respectively, the rest central.
	 * @param input sample points of f(x)
	 * @param output sample points of f'(x)
	 * */
	[[maybe_unused]] void computeGradient(std::vector<double> &input, std::vector<double> &output);

	/**
	 * Solves the equation system Ax=d, where A is a tridiagonal matrix composed of the diagonals a, b, and c.
	 * The lengths of a,b, and c should be equal. Pad with 0.
	 * Solves it in O(n).
	 * Implementation from: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	 * @param a lower diagonal (0, a_2 to a_n) padded with 0 in front
	 * @param b middle diagonal (b_1 to b_n)
	 * @param c upper diagonal (c_1 to c_(n-1), 0) padded with 0 in the end
	 * @param x output vector
	 * @param d right vector
	 * */
	[[maybe_unused]] void solveTriDiagonalMatrix(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c,
												 std::vector<double> &x, std::vector<double> &d);

	/**
	 * Evaluates the Cubic Hermite interpolation spline at position x.
	 * @param x evaluation point
	 * @param fun function representation
	 * */
	[[maybe_unused]] double computeHermiteAt(double x, Function &fun);

	/**
	 * Creates the Cubic Hermite interpolation spline based on the specified knots and stores them in fun.
	 * fun will take ownership of knots and steps.
	 * Assumes that the boundary gradients are zero.
	 * @param begin starting point of knots, x_0
	 * @param fVals sample points of f(x)
	 * @param steps distance between each knot, i.e. x_i - x_(i+1), must be of size: samples - 1.
	 * @param samples total number of samples
	 * @param fun output buffer
	 * */
	[[maybe_unused]] void
	computeHermite(double begin, std::vector<double> &fVals, std::vector<double> &steps, int samples, Function &fun);

	/**
	 * Integrates f by integrating each spline piece of f symbolically.
	 * Constant offset is assumed to be zero, but can be added manually later.
	 * @param f function f(x)
	 * @param F integral function F(x), with F'(x) = f(x)
	 * */
	[[maybe_unused]] void computeIntegral(Function &f, Function &F);

	/**
	 * Computes the derivative of F piece-wise.
	 * @param F function F(x)
	 * @param f derivative of F, with f(x) = F'(x)
	 * */
	[[maybe_unused]] void computeGradient(Function &F, Function &f);

	/**
	 * Gaussian Kernel for kernel smoothing. The term 2b**2 was replaced with sigma**2
	 * */
	[[maybe_unused]] inline double gaussian_kernel(double x, double x_i, double sigma);

	/**
	 * Generates as matrix to be used for smoothing sampled data, if sampled data is y=f(x)
	 * @param begin begin of x
	 * @param end end of x
	 * @param step_width sampling step width
	 * @param sigma filter strength
	 * @param output output
	 * */
	[[maybe_unused]] void
	createGaussianMatrix(double begin, double end, double step_width, double sigma, Matrix &output);

	/**
	 * Evaluates the provided function in the range begin:step_width:end to in- or decrease resolution.
	 * Results is stored back into same function.
	 * @param begin begin inclusive
	 * @param end end exclusive
	 * @param step_width distance between each step
	 * @param function in/out function
	 * */
	[[maybe_unused]] void resampleFunction(double begin, double end, double step_width, Function &function);

	/**
	 * Generates a Gaussian Mixture Model directly by using the provided vector of points a roots.
	 * Since every point is part of a gaussian function evaluated at each position i * step_width,
	 * the result is approximated using a cubic hermite spline to avoid a complexity of O(N*M) for each GMM evaluation.
	 * See: https://arxiv.org/pdf/1504.07351.pdf page 3 (right)
	 * @param begin begin of x
     * @param end end of x
     * @param samples number of samples for hermite spline
     * @param xi filter strength
     * @param centers gaussian centers
     * @param function output function
	 * */
	[[maybe_unused]] void
	createGMM(double begin, double end, int samples, double xi, const std::vector<double> &centers, Function &function);

	/**
	 * Transforms a real valued and periodic function f into its spectral space, whilst using limited frequencies.
	 * @param R container of vectors r_i
	 * @param k_max highest frequency to analyze for
	 * @param output buffer to store results
	 * @param T cycle length
	 * */
	[[maybe_unused]] void
	realFT(const std::vector<double> &R, unsigned int k_max, double T, std::vector<std::complex<double>> &output);

	/**
	 * Filters out high frequencies from F
	 * TODO: dummy function at the moment, applies gaussian to low frequencies
	 * */
	[[maybe_unused]] void filterFT(std::vector<std::complex<double>> &F);

	/**
	 * Computes the inverse Real FT for the coefficients F.
	 * This should transform F back into real space f, which is real, periodic and continuous.
	 * To represent f, it is sampled and then stored as a Hermite Spline
	 * @param F fourier coefficients
	 * @param begin start of sampling
	 * @param end end of sampling
	 * @param samples number of samples
	 * @param function buffer to store f
	 * @param highest frequency of F
	 * @param T cycle length
	 * */
	[[maybe_unused]] void
	ift(const std::vector<std::complex<double>> &F, unsigned int k_max, double T, double begin, double end, int samples,
		Function &function);

	/**
	 * Finite Elements
	 * */
	namespace FE {
		using i3 = std::array<int, 3>;
		using d3 = std::array<double, 3>;
		using idx_t = unsigned long;

		/**
		 * Singular node of FE, stores its location in domain space and property to be measured.
		 * */
		template<typename T>
		class Node {
		public:
			Node() = default;

			/**
			 * Get immutable reference to node position
			 * */
			[[nodiscard]] const d3 &getPos() const { return _pos; }

			/**
			 * Sets position to @param p
			 * */
			void setPos(const d3 &p) { _pos = p; }

			/**
			 * Sets position to @params x, y, z
			 * */
			void setPos(double x, double y, double z) { _pos = {x, y, z}; }

			/**
			 * Gets a mutable reference to the stored property
			 * */
			T &data() { return _data; }

		private:
			/// position in domain space
			d3 _pos = {0, 0, 0};
			/// measured property
			T _data;
		};

		/**
		 * Singular element of FE, stores its lower and upper bounds inside, as well as its respective nodes as 1D coordinates
		 * */
		class Element {
		public:
			Element() = default;

			Element(const std::array<idx_t, 8> &indices, const d3 &lower, const d3 &upper) : _node_indices(indices),
																							 _lower_bound(lower),
																							 _upper_bound(upper) {}

			/**
			 * Get immutable reference to node indices
			 * */
			[[nodiscard]] const std::array<idx_t, 8> &getNodes() const { return _node_indices; }

			/**
			 * Sets the index of node @param pos to index @param idx
			 * pos must be a value in range of 0 (inclusive) to 8 (exclusive)
			 * */
			void setNodeIdx(idx_t idx, int pos) { _node_indices[pos] = idx; }

			/**
			 * @brief Transforms the provided position into a percentage scale from 0 to 1 for each dimension, based on where it
			 * it inside of this element.
			 * If the point is outside the element, behaviour is undefined.
			 * @param position query position
			 * */
			d3 getRelativeLocalPosition(const d3& position) {
				d3 local_point;
				local_point[0] = 2.0 * (position[0]-_lower_bound[0]) / (_upper_bound[0] - _lower_bound[0]) - 1.0;
				local_point[1] = 2.0 * (position[1]-_lower_bound[1]) / (_upper_bound[1] - _lower_bound[1]) - 1.0;
				local_point[2] = 2.0 * (position[2]-_lower_bound[2]) / (_upper_bound[2] - _lower_bound[2]) - 1.0;
				return local_point;
			}

			/**
			 * Sets all indices
			 * */
			void setNodeIndices(idx_t i0, idx_t i1, idx_t i2, idx_t i3, idx_t i4, idx_t i5, idx_t i6,
								idx_t i7) { _node_indices = {i0, i1, i2, i3, i4, i5, i6, i7}; }

			/**
			 * Get immutable reference to lower bound
			 * */
			[[nodiscard]] const d3 &getLow() const { return _lower_bound; }

			/**
			 * Get immutable reference to upper bound
			 * */
			[[nodiscard]] const d3 &getUpper() const { return _upper_bound; }

			/**
			 * Sets lower bound
			 * */
			void setLow(const d3 &low) { _lower_bound = low; }

			/**
			 * Sets lower bound
			 * */
			void setLow(double x, double y, double z) { _lower_bound = {x, y, z}; }

			/**
			 * Sets upper bound
			 * */
			void setUpper(const d3 &upper) { _upper_bound = upper; }

			/**
			 * Sets upper bound
			 * */
			void setUpper(double x, double y, double z) { _upper_bound = {x, y, z}; }

		private:
			/// nodes indices, idx into Node Vector
			std::array<idx_t, 8> _node_indices;
			/// lower bound in domain space
			d3 _lower_bound;
			/// upper bound in domain space
			d3 _upper_bound;
		};

		/**
		 * Domain of FE is discretized into N elements. Each element consists of M nodes and can be an arbitrary 3D object.
		 * Neighbouring elements, share the same nodes on touching faces.
		 * For now we assume, that elements are cubes.
		 * This container, stores all nodes of the FE domain.
		 * T is the type of the measured property.
		 * */
		template<typename T>
		class Nodes {
		public:
			Nodes() = default;

			Node<T> *data() { return _data.data(); }

			void resize(idx_t size, const Node<T> &data = Node<T>()) {
				_data.resize(size, data);
				_node_count = size;
			}

			Node<T> &operator[](idx_t idx) { return _data[idx]; }

			Node<T> &at(idx_t idx) { return _data.at(idx); }

			Node<T> &at(idx_t idx_x, idx_t idx_y, idx_t idx_z) { return _data[get1DIndex(idx_x, idx_y, idx_z)]; }

			[[nodiscard]] idx_t get1DIndex(idx_t idx_x, idx_t idx_y, idx_t idx_z) const {
				return idx_x + _nodes_per_dim[0] * idx_y + _nodes_per_dim[0] * _nodes_per_dim[1] * idx_z;
			}

			idx_t size() { return _node_count; }

			[[nodiscard]] const i3 &shape() const { return _nodes_per_dim; }

			void setShape(const i3 &shape) { _nodes_per_dim = shape; }

			typename std::vector<Node<T>>::iterator begin() { return _data.begin(); }

			typename std::vector<Node<T>>::const_iterator cbegin() { return _data.cbegin(); }

			typename std::vector<Node<T>>::iterator end() { return _data.end(); }

			typename std::vector<Node<T>>::const_iterator cend() { return _data.cend(); }

		private:
			/// storage of all nodes
			std::vector<Node<T>> _data;
			/// count of nodes for each dimension
			i3 _nodes_per_dim = {0, 0, 0};
			/// same as _data.size()
			idx_t _node_count = 0;
		};

		/**
		 * Domain of FE is discretized into N elements. Each element consists of M nodes and can be an arbitrary 3D object.
		 * Neighbouring elements, share the same nodes on touching faces.
		 * For now we assume, that elements are cubes.
		 * This container, stores all elements of the FE domain.
		 * */
		class Elements {
		public:
			Elements() = default;

			Element *data() { return _data.data(); }

			void resize(idx_t size, const Element &data = Element()) {
				_data.resize(size, data);
				_elem_count = size;
			}

			Element &operator[](idx_t idx) { return _data[idx]; }

			Element &at(idx_t idx) { return _data.at(idx); }

			Element &at(idx_t idx_x, idx_t idx_y, idx_t idx_z) { return _data[get1DIndex(idx_x, idx_y, idx_z)]; }

			[[nodiscard]] idx_t get1DIndex(idx_t idx_x, idx_t idx_y, idx_t idx_z) const {
				return idx_x + _elems_per_dim[0] * idx_y + _elems_per_dim[0] * _elems_per_dim[1] * idx_z;
			}

			[[nodiscard]] idx_t size() const { return _elem_count; }

			[[nodiscard]] const i3 &shape() const { return _elems_per_dim; }

			void setShape(const i3 &shape) { _elems_per_dim = shape; }

			[[nodiscard]] const d3 &getElementSize() const { return _elem_size; }

			void setElementSize(const d3 &size) {
				_elem_size = size;
				_elem_volume = _elem_size[0] * _elem_size[1] * _elem_size[2];
			}

			[[maybe_unused]] [[nodiscard]] double getElementVolume() const { return _elem_volume; }

			std::vector<Element>::iterator begin() { return _data.begin(); }

			std::vector<Element>::const_iterator cbegin() { return _data.cbegin(); }

			std::vector<Element>::iterator end() { return _data.end(); }

			std::vector<Element>::const_iterator cend() { return _data.cend(); }

		private:
			/// storage of all elements
			std::vector<Element> _data;
			/// count of elements for each dimension
			i3 _elems_per_dim = {0, 0, 0};
			/// width of each element in each dimension
			d3 _elem_size = {0, 0, 0};
			/// same as _data.size()
			idx_t _elem_count = 0;
			/// same as _elem_size[0]*_elem_size[1]*_elem_size[2]
			double _elem_volume = 0;
		};

		/**
		 * Spans a cuboid subregion within a grid by defining the 3D coordinates of the lower and upper nodes
		 * */
		class SubGrid {
		public:
			SubGrid(const i3 &lower, const i3 &upper) : _lower(lower), _upper(upper) {}

		private:
			/// lower bound of cuboid
			i3 _lower;
			/// upper bound of cuboid
			i3 _upper;
		};

		/**
		  * The grid class creates a highly structured grid of box elements with nodes on the corners. It is composed
		  * of nodes and elements (no edges nor faces). It must be defined using the number of elements per dimension
		  * as well as the lower and upper coordinates of the region to be meshed.
		  * T is the type of the measured property.
		 * */
		template<typename T>
		class Grid {
		public:
			Grid() = default;

			void init(const d3 &lower, const d3 &upper, int num_elem_x, int num_elem_y, int num_elem_z);

			Nodes<T> &getNodes() { return _nodes; }

			Elements &getElements() { return _elements; }

			/// lower bound of entire grid
			[[nodiscard]] const d3 &getLower() const { return _lower_bound; }

			/// upper bound of entire grid
			[[nodiscard]] const d3 &getUpper() const { return _upper_bound; }

			/**
		 	 * Will find the element in which the point resides in. If the point is outside of the range spanned by the grid,
		 	 * behaviour is undefined.
		 	 * */
			Element& getElementOf(const d3& point);

			template<typename R>
			friend std::ostream &operator<<(std::ostream &out, Grid<R> &grid);

		private:
			/// container for all nodes
			Nodes<T> _nodes;
			/// container for all elements
			Elements _elements;
			/// lower bound of the entire grid
			d3 _lower_bound = {0, 0, 0};
			/// upper bound of the entire grid
			d3 _upper_bound = {0, 0, 0};
		};

		template<typename T>
		std::ostream &operator<<(std::ostream &out, Grid<T> &grid) {
			static const std::string prefix = "//[Grid]: ";
			std::stringstream ss;
			ss << prefix << "Meshed region is: ("
			   << grid._lower_bound[0] << "," << grid._lower_bound[1] << "," << grid._lower_bound[2] << ")x("
			   << grid._upper_bound[0] << "," << grid._upper_bound[1] << "," << grid._upper_bound[2] << ")\n";
			ss << prefix << "Total elements: " << grid._elements.size() << "\n";
			ss << prefix << "Length per direction: ["
			   << grid._elements.getElementSize()[0] << ","
			   << grid._elements.getElementSize()[1] << ","
			   << grid._elements.getElementSize()[2] << "]\n";
			ss << prefix << "Elements per direction: ["
			   << grid._elements.shape()[0] << ","
			   << grid._elements.shape()[1] << ","
			   << grid._elements.shape()[2] << "]\n";
			ss << prefix << "Volume per element: " << grid._elements.getElementVolume() << "\n";
			ss << prefix << "Largest Element Index: " << grid._elements.size() - 1 << "\n";
			ss << prefix << "Total nodes: " << grid._nodes.size() << "\n";
			ss << prefix << "Nodes per direction: ("
			   << grid._nodes.shape()[0] << ","
			   << grid._nodes.shape()[1] << ","
			   << grid._nodes.shape()[2] << ")\n";

			ss << " \n\n";
			ss << "idx\t x\t y\t z \n";
			int total_nodes = grid._nodes.size();
			for (int i = 0; i < total_nodes; i++) {
				ss << i << "\t"
				   << grid._nodes[i].getPos()[0] << "\t"
				   << grid._nodes[i].getPos()[1] << "\t"
				   << grid._nodes[i].getPos()[2] << "\n";
			}

			ss << "\n\n";
			ss << "idx\t 0\t 1\t 2\t 3\t 4\t 5\t 6\t 7\n";
			int total_elements = grid._elements.size();
			for (int i = 0; i < total_elements; i++) {
				const std::array<int, 8> &nodes = grid._elements[i].getNodes();
				ss << i << "\t";
				for (int j = 0; j < 8; j++) ss << nodes[j] << "\t";
				ss << "\n";

			}

			out << ss.str();
			return out;
		}

		template<typename T>
		void Grid<T>::init(const d3 &lower, const d3 &upper, int num_elem_x, int num_elem_y, int num_elem_z) {
			_lower_bound = lower;
			_upper_bound = upper;

			// Element construction
			_elements.setShape({num_elem_x, num_elem_y, num_elem_z});
			_elements.resize(num_elem_x * num_elem_y * num_elem_z);
			d3 element_size = {
					(upper[0] - lower[0]) / static_cast<double>(num_elem_x),
					(upper[1] - lower[1]) / static_cast<double>(num_elem_y),
					(upper[2] - lower[2]) / static_cast<double>(num_elem_z),
			};
			_elements.setElementSize(element_size);
			for (int z = 0; z < num_elem_z; z++) {
				for (int y = 0; y < num_elem_y; y++) {
					for (int x = 0; x < num_elem_x; x++) {
						auto &e = _elements.at(x, y, z);
						e.setLow(x * element_size[0], y * element_size[1], z * element_size[2]);
						e.setUpper((x + 1) * element_size[0], (y + 1) * element_size[1], (z + 1) * element_size[2]);
					}
				}
			}

			// Node construction
			_nodes.setShape({num_elem_x + 1, num_elem_y + 1, num_elem_z + 1});
			_nodes.resize((num_elem_x + 1) * (num_elem_y + 1) * (num_elem_z + 1));
			for (int z = 0; z < num_elem_z + 1; z++) {
				for (int y = 0; y < num_elem_y + 1; y++) {
					for (int x = 0; x < num_elem_x + 1; x++) {
						auto &n = _nodes.at(x, y, z);
						n.setPos(x * element_size[0], y * element_size[1], z * element_size[2]);
					}
				}
			}

			// Create links from nodes to elements
			for (int z = 0; z < num_elem_z; z++) {
				for (int y = 0; y < num_elem_y; y++) {
					for (int x = 0; x < num_elem_x; x++) {
						auto &e = _elements.at(x, y, z);

						int counter = 0;
						for (int dz = 0; dz < 2; dz++) {
							for (int dy = 0; dy < 2; dy++) {
								for (int dx = 0; dx < 2; dx++) {
									e.setNodeIdx(_nodes.get1DIndex(x + dx, y + dy, z + dz), counter);
									counter++;
								}
							}
						}
					}
				}
			}
		}

		/**
		 * Will compute gradient over provided grid.
		 * EX_SRC loads the property of each cell over which the gradient should be computed.
		 * Its signature must be equivalent to: D (*)(const n_data_t&)
		 * EX_DST stores the three gradient values back into the grid.
		 * Its signature must be equivalent to: void (*)(n_data_t&, int dim, D grad)
		 * dim will have values 0, 1, 2
		 *
		 * Gradients on the border will use forward/backward differences, inner ones central differences.
		 * */
		template<typename n_data_t, typename EX_SRC, typename EX_DST, typename D=double>
		void computeGradient(Grid<n_data_t>& grid, EX_SRC src_extractor, EX_DST dst_extractor) {
			auto& nodes = grid.getNodes();
			const auto& shape = nodes.shape();
			using i2 = std::array<int, 2>;
			const std::array<i2, 3> dim_map = {
					i2{2, 1}, // dim 0
					i2{2, 0}, // dim 1
					i2{1, 0}  // dim 2
			};
			auto& mesh_width = grid.getElements().getElementSize();

			// do forward and backward FD for border
			for (int dim = 0; dim < 3; dim++) {
				for (idx_t j = 0; j < shape[dim_map[dim][0]]; j++) {
					for (idx_t i = 0; i < shape[dim_map[dim][1]]; i++) {
						// forward FD
						std::array<idx_t, 3> coords0 { 0 };
						coords0[dim_map[dim][0]] = j;
						coords0[dim_map[dim][1]] = i;
						std::array<idx_t, 3> coords1 = coords0;
						coords1[dim] += 1;

						auto& node0_fw = nodes.at(coords0[0], coords0[1], coords0[2]);
						auto& node1_fw = nodes.at(coords1[0], coords1[1], coords1[2]);
						D fd_fw = (src_extractor(node1_fw) - src_extractor(node0_fw)) / mesh_width[dim];
						dst_extractor(node0_fw, dim, fd_fw);

						// backward FD
						coords0[dim] = shape[dim] - 1;
						coords1[dim] = coords0[dim] - 1;

						auto& node0_bw = nodes.at(coords0[0], coords0[1], coords0[2]);
						auto& node1_bw = nodes.at(coords1[0], coords1[1], coords1[2]);
						D fd_bw = (src_extractor(node0_bw) - src_extractor(node1_bw)) / mesh_width[dim];
						dst_extractor(node0_bw, dim, fd_bw);
					}
				}
			}

			// do central FD for inner nodes
			for (idx_t z = 0; z < shape[2]; z++) {
				for (idx_t y = 0; y < shape[1]; y++) {
					for (idx_t x = 0; x < shape[0]; x++) {
						for (int dim = 0; dim < 3; dim++) {
							std::array<idx_t, 3> coords { x, y, z };
							if (coords[dim] == 0 || coords[dim] == shape[dim] - 1) continue;

							std::array<idx_t, 3> coords0 = coords;
							std::array<idx_t, 3> coords1 = coords;
							coords0[dim] -= 1;
							coords1[dim] += 1;

							auto& node_c = nodes.at(coords[0], coords[1], coords[2]);
							auto& node0 = nodes.at(coords0[0], coords0[1], coords0[2]);
							auto& node1 = nodes.at(coords1[0], coords1[1], coords1[2]);
							D fd_c = (src_extractor(node1) - src_extractor(node0)) / (2 * mesh_width[dim]);
							dst_extractor(node_c, dim, fd_c);
						}
					}
				}
			}
		}

		/**
		 * Will find the element in which the point resides in. If the point is outside of the range spanned by the grid,
		 * behaviour is undefined.
		 * */
		template<typename T>
		Element &Grid<T>::getElementOf(const d3 &point) {
			const auto& elem_size = _elements.getElementSize();
			std::array<idx_t, 3> coord { 0 };
			for (int dim = 0; dim < 3; dim++) coord[dim] = point[dim] / elem_size[dim];

			const auto& elem_shape = _elements.shape();
			for (int dim = 0; dim < 3; dim++) std::clamp<idx_t>(coord[dim], 0, elem_shape[dim] - 1);

			return _elements.at(coord[0], coord[1], coord[2]);
		}
	};
}
#endif //MARDYN_INTERPOLATION_H
