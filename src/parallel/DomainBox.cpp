/**
 * @file DomainBox.h
 * @author Georg von Bismarck
 * @date 27.10.2026
 */

#include "DomainBox.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include "utils/mardyn_assert.h"

DomainBox::DomainBox(const DomainPoint& min_, const DomainPoint& max_) noexcept
    : min(min_), max(max_) {}

DomainBox::DomainBox(double xmin, double ymin, double zmin,
                     double xmax, double ymax, double zmax) noexcept
    : min{xmin, ymin, zmin},
      max{xmax, ymax, zmax} {}


DomainBox DomainBox::fromArray(const std::array<double, 6>& data) noexcept {
    return DomainBox{
        data[0], data[1], data[2],
        data[3], data[4], data[5]
    };
}

std::vector<DomainBox> DomainBox::fromFlatArray(const std::vector<double>& data) {
    if (data.size() % 6 != 0){
        std::ostringstream error_message;
        error_message << "DomainBox: Conversion from a Double Vector to a DomainBox Vector is not possible." << std::endl;
        MARDYN_EXIT(error_message.str());
    }
    std::vector<DomainBox> boxes;
    boxes.reserve(data.size() / 6);

    for (std::size_t i = 0; i < data.size(); i += 6) {
        boxes.emplace_back(
            data[i + 0], data[i + 1], data[i + 2],
            data[i + 3], data[i + 4], data[i + 5]
        );
    }

    return boxes;
}

DomainBox DomainBox::createAsRegularGrid(const DomainPoint& domainLength, const DomainGridPoint& gridSize, const DomainGridPoint& gridCoords) {
	DomainBox box = DomainBox(); 

	// initialize DomainBox as regular grid!
	for (int dim = 0; dim < 3; ++dim) {
		box.min[dim] = gridCoords[dim] * domainLength[dim] / gridSize[dim];
		box.max[dim] = (gridCoords[dim] + 1) * domainLength[dim] / gridSize[dim];
		if (gridCoords[dim] == gridSize[dim] - 1) {
			// ensure that the upper domain boundaries match.
			// lower domain boundaries always match, because they are 0.
			box.max[dim] = domainLength[dim];
		}
	}

    return box;
}

double DomainBox::width() const noexcept {
    return max[0] - min[0];
}

double DomainBox::height() const noexcept {
    return max[1] - min[1];
}

double DomainBox::depth() const noexcept {
    return max[2] - min[2];
}

double DomainBox::volume() const noexcept {
    if (isValid()) {
        return width() * height() * depth();
    }
    return 0;
}

bool DomainBox::isValid() const noexcept {
    return min[0] <= max[0] &&
           min[1] <= max[1] &&
           min[2] <= max[2];
}

std::optional<DomainBox> DomainBox::intersection(const DomainBox& other) const noexcept {
    const DomainPoint newMin{
        std::max(min[0], other.min[0]),
        std::max(min[1], other.min[1]),
        std::max(min[2], other.min[2])
    };

    const DomainPoint newMax{
        std::min(max[0], other.max[0]),
        std::min(max[1], other.max[1]),
        std::min(max[2], other.max[2])
    };

    if (newMin[0] > newMax[0] ||
        newMin[1] > newMax[1] ||
        newMin[2] > newMax[2]) {
        return std::nullopt;
    }

    return DomainBox{newMin, newMax};
}

double DomainBox::intersectionVolume(const DomainBox& other) const noexcept {
    const double w = std::min(max[0], other.max[0]) - std::max(min[0], other.min[0]);
    const double h = std::min(max[1], other.max[1]) - std::max(min[1], other.min[1]);
    const double d = std::min(max[2], other.max[2]) - std::max(min[2], other.min[2]);

    if (w <= 0 || h <= 0 || d <= 0) return 0.0;
    return w * h * d;
}

std::array<double, 6> DomainBox::toArray() const noexcept {
    return { min[0], min[1], min[2], max[0], max[1], max[2] };
}

void DomainBox::latchToGridSize(const DomainPoint& latchedGridSize, const DomainPoint& domainLength) noexcept {
    for (std::size_t ind = 0; ind < 3; ++ind) {
        const double currentGridSize = latchedGridSize[ind];

        min[ind] = std::round(min[ind] / currentGridSize) * currentGridSize;

        if (max[ind] != domainLength[ind]) {
            max[ind] = std::round(max[ind] / currentGridSize) * currentGridSize;
        }
    }
}

std::ostream& operator<<(std::ostream& os, const DomainBox& box)
{
    os << "DomainBox{"
       << "min=[" << box.min[0] << ", " << box.min[1] << ", " << box.min[2] << "], "
       << "max=[" << box.max[0] << ", " << box.max[1] << ", " << box.max[2] << "]"
       << "}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<DomainBox>& boxes) {
    os << "[";
    bool first = true;

    for (const auto& box : boxes) {
        if (!first) {
            os << ", ";
        } else {
            first = false;
        }
        os << box;   
    }

    os << "]";
    return os;
}