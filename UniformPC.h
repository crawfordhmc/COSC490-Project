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

    std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane);

protected:

    double voxel_size;
    int x_voxels, y_voxels, z_voxels;
    std::vector<std::vector<std::vector<std::vector<size_t>>>> cells;

    std::vector<size_t> hashCell(Eigen::Vector3d p);

    void UniformPC::cleary(std::vector<size_t>& points, std::vector<size_t> cell, Eigen::Vector3d dir,
        std::vector<std::vector<std::vector<bool>>>& visited, Eigen::Hyperplane<double, 3> thisPlane);

    void UniformPC::addPoints(std::vector<size_t> indexes, std::vector<size_t> &thisPoints,
        Eigen::Hyperplane<double, 3> thisPlane);

};
#endif