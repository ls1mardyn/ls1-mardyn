/**
 * @file DomainBox.h
 * @author Georg von Bismarck
 * @date 27.10.2026
 */

#pragma once

#include <array>
#include <optional>
#include <ostream>
#include <vector>

using DomainGridPoint = std::array<int, 3>;
using DomainPoint = std::array<double, 3>;

// A 3D bounding box defined by two corners
class DomainBox {
public:
    //Lower corner of the DomainBox.
    DomainPoint min{};
    //Upper corner of the DomainBox.
    DomainPoint max{};

    DomainBox() = default;
    
    /**
     * Constructs from two corner points.
     * @param min_ Lower corner
     * @param max_ Upper corner
     */
    DomainBox(const DomainPoint& min_, const DomainPoint& max_) noexcept;
    
    /**
     * Constructs from raw coordinate values.
     * @param xmin Minimum X coordinate
     * @param ymin Minimum Y coordinate
     * @param zmin Minimum Z coordinate
     * @param xmax Maximum X coordinate
     * @param ymax Maximum Y coordinate
     * @param zmax Maximum Z coordinate
     */
    DomainBox(double xmin, double ymin, double zmin,
              double xmax, double ymax, double zmax) noexcept;
    

    /**
     * Constructs from raw array.
     * @param data must be {min.x, min.y, min.z, max.x, max.y, max.z}
     */
    static DomainBox fromArray(const std::array<double, 6>& data) noexcept;
    
    /**
     * Constructs as regular grid!
     * @param domainLength
     * @param gridSize
     * @param gridCoords 
     */
    static DomainBox createAsRegularGrid(const DomainPoint& domainLength, const DomainGridPoint& gridSize, const DomainGridPoint& gridCoords);

    /**
     * Converts a flat list of coordinates into multiple boxes.
     * The input must contain a multiple of 6 values.
     * @param data Flat coordinate list. Every 6 values have to be {min.x, min.y, min.z, max.x, max.y, max.z}
     * @return Vector of DomainBox objects.
     * @throws MARDYN_EXIT if data.size() is not a multiple of 6
     */
    static std::vector<DomainBox> fromFlatArray(const std::vector<double>& data);
    
    /**
     * @return the DomainBox width in x direction.
     */
    double width() const noexcept;
    
    /**
     * @return the DomainBox height in y direction.
     */
    double height() const noexcept;
    
    /**
     * @return the DomainBox depth in z direction.
     */
    double depth() const noexcept;

    /**
     * @return the DomainBox volume, or 0 if the DomainBox is not valid.
     */
    double volume() const noexcept;

    /**
     * Checks whether the DomainBox is valid.
     * @return true if min <= max in all dimensions.
     */
    bool isValid() const noexcept;

    /**
     * Computes the intersection with another DomainBox.
     * @param other The other DomainBox
     * @return The intersection DomainBox, or std::nullopt if there is no overlap.
     */
    std::optional<DomainBox> intersection(const DomainBox& other) const noexcept;
    
    /**
     * Computes the intersection Volume with another DomainBox.
     * @param other The other DomainBox
     * @return The intersection Volume, or 0 if there is no overlap.
     */
    double intersectionVolume(const DomainBox& other) const noexcept;

    /**
     * Converts the DomainBox to a flat array.
     * @return {min.x, min.y, min.z, max.x, max.y, max.z}.
     */
    std::array<double, 6> toArray() const noexcept;

    /**
     * Latches domain boundaries to a grid, which is defined by latchedGridSize.
	 * If boxMax matches the top boundary, it is not changed.
     * @param latchedGridSize the grid to be latched to
     * @param domainLength The total size of the domain
     */
    void latchToGridSize(const DomainPoint& latchedGridSize, const DomainPoint& domainLength) noexcept;

    
};

std::ostream& operator<<(std::ostream& os, const DomainBox& box);
std::ostream& operator<<(std::ostream& os, const std::vector<DomainBox>& boxes);