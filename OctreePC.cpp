#include "OctreePC.h"

#include <iostream>
#include <fstream>
#include <sstream>


//key idea for octree:
// construct uniform space subdivision with voxels being the threshold size, using a dense tree instead of a 3d array
// if voxel + neighbours total points < point threshold, merge into a larger voxel and append point vectors to parent node
// going up the tree
// can be done in parallel at each depth
// for ease of programming, accessing voxel dimensions should be a method that returns them calulated from tree position
OctreePC::OctreePC(const std::string& filepath, float scale_parameter) : PointCloud(filepath, scale_parameter) {}


// Returns a vector of points within the threshold to the given hyperplane
// (also prints the number of threads being used for the calculations)
std::vector<size_t> OctreePC::planePoints(Eigen::Hyperplane<double, 3> thisPlane, std::vector<size_t> remainingPoints, unsigned int trial, int plane) {
    std::vector<size_t> thisPoints;
    int threads = 0;
    //OpenMP requires signed integrals for its loop variables... interesting
    signed long long i = 0;
#pragma omp parallel for shared(thisPoints) private (i)
    for (i = 0; i < remainingPoints.size(); ++i) {
        if (thisPlane.absDistance(pc[remainingPoints[i]].location) < threshold)
#pragma omp critical
            thisPoints.push_back(remainingPoints[i]);
        if (omp_get_thread_num() == 0 && trial == 0 && plane == 0)
            threads = omp_get_max_threads();
    }
    if (trial == 0 && plane == 0)
        std::cout << threads << " threads are being used" << std::endl;
    return thisPoints;
}
