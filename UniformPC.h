#ifndef UNIFORM_PC_H
#define UNIFORM_PC_H

#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class UniformPC : public PointCloud {
public:

    UniformPC(PointCloud const&p, int voxel_scale);

    std::vector<size_t> planePoints(const Eigen::Hyperplane<double, 3> &thisPlane);
    //std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane);

protected:

    int voxel_size;
    size_t x_voxels, y_voxels, z_voxels;
    std::vector<size_t> limits;
    std::vector<std::vector<std::vector<std::vector<size_t>>>> cells;

    std::vector<size_t> hashCell(const Eigen::Vector3d &p);

    void UniformPC::cleary(std::vector<size_t>& points, double next_yz, std::vector<size_t> cell, const Eigen::Vector3d &dir, const Eigen::Vector3d &norm,
        std::vector<std::vector<std::vector<bool>>>& visited, Eigen::Hyperplane<double, 3> thisPlane, Eigen::Vector3d fpoint);

    void UniformPC::checkcell(size_t x, size_t y, size_t z, std::vector<size_t>& points, std::vector<std::vector<std::vector<bool>>>& visited,
        const Eigen::Hyperplane<double, 3> &thisPlane);

    void UniformPC::addPoints(std::vector<size_t> indexes, std::vector<size_t> &thisPoints,
        const Eigen::Hyperplane<double, 3> &thisPlane);

};
#endif