#ifndef OCTREE_PC_H
#define OCTREE_PC_H

#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class OctreePC : public PointCloud {
public:

    OctreePC(PointCloud &p);
    std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane, std::vector<size_t> remainingPoints, unsigned int trial, int plane);
};
#endif