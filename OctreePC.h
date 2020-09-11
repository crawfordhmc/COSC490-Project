#ifndef OCTREE_PC_H
#define OCTREE_PC_H

#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class OctreePC : public PointCloud {
public:

    OctreePC(const std::string& filepath, float scale_parameter);
    std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane, std::vector<size_t> remainingPoints, unsigned int trial, int plane);
};
#endif