#ifndef UNIFORM_PC_H
#define UNIFORM_PC_H

#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class UniformPC : public PointCloud {
public:

    UniformPC(PointCloud const&p, float voxel_scale);

    std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane, std::vector<size_t> remainingPoints, unsigned int trial, int plane);

protected:

    double voxel_size;
    int x_voxels, y_voxels, z_voxels;
    std::vector<std::vector<std::vector<std::vector<size_t>>>> cells;

    std::vector<size_t> hashCell(Eigen::Vector3d p);

    std::vector<size_t> UniformPC::cleary(std::vector<size_t> &points, std::vector<size_t> remainingPoints, Eigen::ParametrizedLine<double, 3> ray,
        std::vector<std::vector<std::vector<bool>>>& visited, Eigen::Hyperplane<double, 3> thisPlane, int plane);

    void UniformPC::addPoints(std::vector<size_t> indexes, std::vector<size_t> &thisPoints, std::vector<size_t> remainingPoints,
        Eigen::Hyperplane<double, 3> thisPlane, int plane);

};
#endif