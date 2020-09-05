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
    std::vector<std::vector<std::vector<std::vector<size_t>*>>> cells;
    Eigen::Matrix<size_t, 3, 1> hashCell(Eigen::Vector3d p);
    void addPoints(std::vector<size_t>* indexes, std::vector<size_t> thisPoints, Eigen::ParametrizedLine<double, 3> ray);
    std::vector<size_t> cleary(std::vector<size_t> points, Eigen::ParametrizedLine<double, 3> ray);

public:

    UniformPC(const std::string& filepath);
    virtual Point getPoint(int index);
    virtual void setPointPlane(int index, int planeID);
    virtual void setPointColour(int index, Eigen::Vector3i colour);
    virtual std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane);
    virtual void writeToPly(const std::string& filename);
};
#endif