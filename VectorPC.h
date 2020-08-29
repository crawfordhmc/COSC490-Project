#ifndef VECTOR_PC_H
#define VECTOR_PC_H

#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class VectorPC : public PointCloud {
public:

    VectorPC(const std::string& filepath);
    virtual Point getPoint(int index);
    virtual void setPointPlane(int index, int planeID);
    virtual void setPointColour(int index, Eigen::Vector3i colour);
    virtual std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane, unsigned int trial, float threshold, int plane);
    virtual void writeToPly(const std::string& filename);
};
#endif