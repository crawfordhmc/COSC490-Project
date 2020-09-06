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
    std::vector<std::vector<std::vector<std::vector<size_t>>>> cells;
    std::vector<size_t> hashCell(Eigen::Vector3d p);
    void addPoints(std::vector<size_t> indexes, std::vector<size_t> thisPoints, Eigen::Hyperplane<double, 3> plane);
    std::vector<size_t> cleary(std::vector<size_t> points, Eigen::ParametrizedLine<double, 3> ray, std::vector<std::vector<std::vector<bool>>> visited, Eigen::Hyperplane<double, 3> plane);

public:

    UniformPC(const std::string& filepath, float scale_parameter);
    virtual Point getPoint(int index);
    virtual void setPointPlane(int index, int planeID);
    virtual void setPointColour(int index, Eigen::Vector3i colour);
    virtual std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane, unsigned int trial, int plane);
    virtual void writeToPly(const std::string& filename);
};
#endif