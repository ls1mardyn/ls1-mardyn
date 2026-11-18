/**
 * @file DomainBoxTest.cpp
 * @author Georg von Bismarck
 * @date 17.11.2026
 */

#include "DomainBoxTest.h"

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "parallel/DomainBox.h"

TEST_SUITE_REGISTRATION(DomainBoxTest);

void DomainBoxTest::assertBoxEquals(const DomainBox& expected, const DomainBox& actual, double delta) {
	for (int dim = 0; dim < DIMgeom; ++dim) {
		ASSERT_DOUBLES_EQUAL(expected.min[dim], actual.min[dim], delta);
		ASSERT_DOUBLES_EQUAL(expected.max[dim], actual.max[dim], delta);
	}
}

void DomainBoxTest::testConstructors() {
	const DomainBox defaultBox;
	for (int dim = 0; dim < DIMgeom; ++dim) {
		ASSERT_EQUAL(0.0, defaultBox.min[dim]);
		ASSERT_EQUAL(0.0, defaultBox.max[dim]);
	}

	const DomainPoint lower{1.0, 2.0, 3.0};
	const DomainPoint upper{4.0, 5.0, 6.0};
	const DomainBox fromPoints(lower, upper);
	for (int dim = 0; dim < DIMgeom; ++dim) {
		ASSERT_EQUAL(dim + 1.0, fromPoints.min[dim]);
		ASSERT_EQUAL(dim + 1.0 + 3.0, fromPoints.max[dim]);
	}

	const DomainBox fromValues(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
	for (int dim = 0; dim < DIMgeom; ++dim) {
		ASSERT_EQUAL(dim + 1.0, fromValues.min[dim]);
		ASSERT_EQUAL(dim + 1.0 + 3.0, fromValues.max[dim]);
	}
}

void DomainBoxTest::testFromArrayToArray() {
	const std::array<double, 6> data{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

	const DomainBox box = DomainBox::fromArray(data);
	for (int dim = 0; dim < DIMgeom; ++dim) {
		ASSERT_EQUAL(dim + 1.0, box.min[dim]);
		ASSERT_EQUAL(dim + 1.0 + 3.0, box.max[dim]);
	}

	//array -> box -> array
	const std::array<double, 6> roundTrip = box.toArray();
	for (std::size_t i = 0; i < data.size(); ++i) {
		ASSERT_EQUAL(data[i], roundTrip[i]);
	}

	//box -> array -> box
	const DomainBox original(-1.5, 0.0, 2.25, 3.5, 7.0, 9.75);
	assertBoxEquals(original, DomainBox::fromArray(original.toArray()), 0.0);
}

void DomainBoxTest::testFromFlatArray() {
	ASSERT_EQUAL(static_cast<std::size_t>(0), DomainBox::fromFlatArray({}).size());

	const std::vector<DomainBox> one = DomainBox::fromFlatArray({1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
	ASSERT_EQUAL(static_cast<std::size_t>(1), one.size());
	assertBoxEquals(DomainBox(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), one[0], 0.0);

	const std::vector<DomainBox> two = DomainBox::fromFlatArray(
		{0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
		 1.0, 0.0, 0.0, 2.0, 1.0, 1.0});
	ASSERT_EQUAL(static_cast<std::size_t>(2), two.size());
	assertBoxEquals(DomainBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0), two[0], 0.0);
	assertBoxEquals(DomainBox(1.0, 0.0, 0.0, 2.0, 1.0, 1.0), two[1], 0.0);
}

void DomainBoxTest::testCreateAsRegularGrid() {
	{
		const DomainPoint domainLength{10.0, 20.0, 30.0};
		const DomainGridPoint gridSize{2, 4, 3};

		const DomainBox first = DomainBox::createAsRegularGrid(domainLength, gridSize, {0, 0, 0});
		assertBoxEquals(DomainBox(0.0, 0.0, 0.0, 5.0, 5.0, 10.0), first);

		const DomainBox inner = DomainBox::createAsRegularGrid(domainLength, gridSize, {0, 1, 1});
		assertBoxEquals(DomainBox(0.0, 5.0, 10.0, 5.0, 10.0, 20.0), inner);

		const DomainBox last = DomainBox::createAsRegularGrid(domainLength, gridSize, {1, 3, 2});
		assertBoxEquals(DomainBox(5.0, 15.0, 20.0, 10.0, 20.0, 30.0), last);
	}

	// non exactly representable cell sizes (1/3, 2/7, ...)
	const DomainPoint domainLength{1.0, 2.0, 3.0};
	const DomainGridPoint gridSize{3, 7, 5};

	double totalVolume = 0.0;
	for (int x = 0; x < gridSize[0]; ++x) {
		for (int y = 0; y < gridSize[1]; ++y) {
			for (int z = 0; z < gridSize[2]; ++z) {
				const DomainGridPoint coords{x, y, z};
				const DomainBox box = DomainBox::createAsRegularGrid(domainLength, gridSize, coords);
				ASSERT_TRUE(box.isValid());
				totalVolume += box.volume();

				for (int dim = 0; dim < DIMgeom; ++dim) {
					// the lower domain boundary is matched exactly
					if (coords[dim] == 0) {
						ASSERT_EQUAL(0.0, box.min[dim]);
					}

					// the upper domain boundary is matched exactly
					if (coords[dim] == gridSize[dim] - 1) {
						ASSERT_EQUAL(domainLength[dim], box.max[dim]);
						continue;
					}

					// neighbouring cells share their boundary exactly (no gaps, no overlap)
					DomainGridPoint neighbourCoords = coords;
					++neighbourCoords[dim];
					const DomainBox neighbour =
						DomainBox::createAsRegularGrid(domainLength, gridSize, neighbourCoords);
					ASSERT_EQUAL(box.max[dim], neighbour.min[dim]);
				}
			}
		}
	}

	// all cells together cover the whole domain
	ASSERT_DOUBLES_EQUAL(domainLength[0] * domainLength[1] * domainLength[2], totalVolume, 1e-12);
}

void DomainBoxTest::testDimensionsAndVolume() {
	const DomainBox box(1.0, 2.0, 3.0, 4.0, 6.0, 8.0);
	ASSERT_DOUBLES_EQUAL(3.0, box.width(), 0.0);
	ASSERT_DOUBLES_EQUAL(4.0, box.height(), 0.0);
	ASSERT_DOUBLES_EQUAL(5.0, box.depth(), 0.0);
	ASSERT_DOUBLES_EQUAL(60.0, box.volume(), 0.0);

	const DomainBox negative(-2.0, -3.0, -4.0, 0.0, 0.0, 0.0);
	ASSERT_DOUBLES_EQUAL(24.0, negative.volume(), 0.0);

	// flat box: valid, but without volume
	const DomainBox flat(0.0, 0.0, 0.0, 1.0, 1.0, 0.0);
	ASSERT_TRUE(flat.isValid());
	ASSERT_DOUBLES_EQUAL(0.0, flat.volume(), 0.0);

	// invalid box: volume is 0 instead of a negative value
	const DomainBox invalid(2.0, 0.0, 0.0, 1.0, 1.0, 1.0);
	ASSERT_DOUBLES_EQUAL(-1.0, invalid.width(), 0.0);
	ASSERT_DOUBLES_EQUAL(0.0, invalid.volume(), 0.0);
}

void DomainBoxTest::testIsValid() {
	ASSERT_TRUE(DomainBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0).isValid());

	// min == max is still valid (degenerate box)
	ASSERT_TRUE(DomainBox(1.0, 1.0, 1.0, 1.0, 1.0, 1.0).isValid());
	ASSERT_TRUE(DomainBox().isValid());

	// min > max in exactly one dimension is invalid
	for (int dim = 0; dim < DIMgeom; ++dim) {
		DomainBox box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
		box.min[dim] = 2.0;
		ASSERT_TRUE_MSG("min > max in dimension " + std::to_string(dim), !box.isValid());
	}
}

void DomainBoxTest::testIntersection() {
	struct IntersectionCase {
		std::string name;
		DomainBox other;
		std::optional<DomainBox> expectedIntersection;  // std::nullopt if no intersection is expected
	};

	const DomainBox baseBox(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

	std::vector<IntersectionCase> cases{
		{"partial overlap", DomainBox(1.0, 1.0, 1.0, 3.0, 3.0, 3.0), DomainBox(1.0, 1.0, 1.0, 2.0, 2.0, 2.0)},
		{"other inside a", DomainBox(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), DomainBox(0.5, 0.5, 0.5, 1.5, 1.5, 1.5)},
		{"a inside other", DomainBox(-1.0, -1.0, -1.0, 3.0, 3.0, 3.0), baseBox},
		{"identical", baseBox, baseBox},
		// touching boxes: intersection() returns a box with zero thickness, not std::nullopt
		{"touching", DomainBox(2.0, 0.0, 0.0, 4.0, 2.0, 2.0), DomainBox(2.0, 0.0, 0.0, 2.0, 2.0, 2.0)},
	};

	// separated in exactly one dimension -> no intersection
	for (int dim = 0; dim < DIMgeom; ++dim) {
		DomainBox other = baseBox;
		other.min[dim] = 3.0;
		other.max[dim] = 4.0;
		cases.push_back({"separated in dimension " + std::to_string(dim), other, std::nullopt});
	}

	for (const IntersectionCase& testCase : cases) {
		const std::string& msg = testCase.name;
		const std::optional<DomainBox> result = baseBox.intersection(testCase.other);
		const std::optional<DomainBox> reverse = testCase.other.intersection(baseBox);

		// intersection()
		ASSERT_EQUAL_MSG(msg, testCase.expectedIntersection.has_value(), result.has_value());
		ASSERT_EQUAL_MSG(msg, result.has_value(), reverse.has_value());
		if (testCase.expectedIntersection.has_value()) {
			assertBoxEquals(*testCase.expectedIntersection, *result, 0.0);
			assertBoxEquals(*result, *reverse, 0.0);
		}

		// intersectionVolume() has to match intersection()->volume()
		const double expectedVolume = result.has_value() ? result->volume() : 0.0;
		ASSERT_DOUBLES_EQUAL_MSG(msg, expectedVolume, baseBox.intersectionVolume(testCase.other), 1e-12);
		ASSERT_DOUBLES_EQUAL_MSG(msg, expectedVolume, testCase.other.intersectionVolume(baseBox), 1e-12);
	}
}

void DomainBoxTest::testIsNearEqual() {
	const DomainBox a(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);

	ASSERT_TRUE(a.isNearEqual(a));
	ASSERT_TRUE(a.isNearEqual(DomainBox(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)));

	const DomainBox zero(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
	ASSERT_TRUE(zero.isNearEqual(DomainBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)));

	ASSERT_TRUE(!a.isNearEqual(DomainBox(7.0, 8.0, 9.0, 10.0, 11.0, 12.0)));

	std::array<double, 6> values = a.toArray();
	for (std::size_t i = 0; i < values.size(); ++i) {
		// deviation inside the relative tolerance (1e-9) -> equal
		std::array<double, 6> inside = values;
		inside[i] *= 1.0 + 1e-12;
		const DomainBox b = DomainBox::fromArray(inside);
		ASSERT_TRUE(a.isNearEqual(b));
		ASSERT_TRUE(b.isNearEqual(a));

		// deviation outside the relative tolerance (1e-9) -> not equal
		std::array<double, 6> outside = values;
		outside[i] *= 1.0 + 1e-6;
		const DomainBox c = DomainBox::fromArray(outside);
		ASSERT_TRUE(!a.isNearEqual(c));
		ASSERT_TRUE(!c.isNearEqual(a));
	}
}

void DomainBoxTest::testLatchToGridSize() {
	const DomainPoint domainLength{10.0, 10.0, 10.0};

	// rounding down and up to a grid of size 1
	{
		DomainBox box(0.4, 0.6, 2.0, 4.4, 4.6, 7.0);
		box.latchToGridSize({1.0, 1.0, 1.0}, domainLength);
		assertBoxEquals(DomainBox(0.0, 1.0, 2.0, 4.0, 5.0, 7.0), box);
	}

	// different grid size per dimension
	{
		DomainBox box(1.3, 1.1, 0.9, 3.7, 4.4, 5.6);
		box.latchToGridSize({0.5, 2.0, 0.25}, domainLength);
		assertBoxEquals(DomainBox(1.5, 2.0, 1.0, 3.5, 4.0, 5.5), box);
	}

	// values that are already on the grid are not changed
	{
		const DomainBox original(2.0, 4.0, 6.0, 4.0, 6.0, 8.0);
		DomainBox box = original;
		box.latchToGridSize({2.0, 2.0, 2.0}, domainLength);
		assertBoxEquals(original, box);
	}

	// max equal to the domain length is not changed, even if it is not on the grid
	{
		const DomainPoint oddDomainLength{10.3, 10.3, 10.3};
		DomainBox box(2.2, 2.2, 2.2, 10.3, 10.3, 10.3);
		box.latchToGridSize({1.0, 1.0, 1.0}, oddDomainLength);
		assertBoxEquals(DomainBox(2.0, 2.0, 2.0, 10.3, 10.3, 10.3), box, 0.0);
	}
}
