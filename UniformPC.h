#ifndef UNIFORM_PC_H
#define UNIFORM_PC_H

#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class UniformPC : public PointCloud {
protected:
    double voxel_size;
    int x_voxels, y_voxels, z_voxels;

public:

    UniformPC(const std::string& filepath);
    virtual Point getPoint(int index);
    virtual void setPointPlane(int index, int planeID);
    virtual void setPointColour(int index, Eigen::Vector3i colour);
    virtual std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane, unsigned int trial, float threshold, int plane);
    std::vector<size_t> UniformPC::checkPoints(std::vector<int> indexes, Eigen::Hyperplane<double, 3> thisPlane, unsigned int trial, float threshold, int plane);
    virtual void writeToPly(const std::string& filename);
};
#endif