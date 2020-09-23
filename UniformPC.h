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
    //void removePoints(std::vector<size_t>& planePoints, int plane);

protected:

    int voxel_size;
    size_t x_voxels, y_voxels, z_voxels;
    std::vector<size_t> limits;
    std::vector<std::vector<std::vector<std::vector<size_t>>>> cells;

    std::vector<size_t> hashCell(const Eigen::Vector3d &p);

};
#endif