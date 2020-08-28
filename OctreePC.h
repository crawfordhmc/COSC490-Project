#ifndef OCTREE_PC_H
#define OCTREE_PC_H

#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class OctreePC : public PointCloud {
protected:
    std::vector<Point> pc;

public:

    OctreePC(const std::string& filepath);
    virtual Point getPoint(int index);
    virtual void setPointPlane(int index, int planeID);
    virtual void setPointColour(int index, Eigen::Vector3i colour);
    virtual std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane, unsigned int trial, float threshold, int plane);
    virtual float threshold(float scale_parameter);
    virtual void writeToPly(const std::string& filename);
};
#endif