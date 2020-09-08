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


// these three may not need to be virtual functions depending on access methods
PointCloud::Point OctreePC::getPoint(int index) { return pc[index]; }
void OctreePC::setPointPlane(int index, int planeID) { pc[index].planeIx = planeID; }
void OctreePC::setPointColour(int index, Eigen::Vector3i colour) { pc[index].colour = colour; }


// Returns a vector of points within the threshold to the given hyperplane
// (also prints the number of threads being used for the calculations)
std::vector<size_t> OctreePC::planePoints(Eigen::Hyperplane<double, 3> thisPlane, std::vector<size_t> removedPoints, unsigned int trial, int plane) {
    std::vector<size_t> thisPoints;
    int threads = 0;
    //OpenMP requires signed integrals for its loop variables... interesting
    signed long long i = 0;
#pragma omp parallel for shared(thisPoints) private (i)
    for (i = 0; i < pc.size(); ++i) {
        if (plane == 0 || !std::binary_search(removedPoints.begin(), removedPoints.end(), i) && thisPlane.absDistance(pc[i].location) < threshold)
#pragma omp critical
            thisPoints.push_back(i);
        if (omp_get_thread_num() == 0 && trial == 0 && plane == 0)
            threads = omp_get_max_threads();
    }
    if (trial == 0 && plane == 0)
        std::cout << threads << " threads are being used" << std::endl;
    return thisPoints;
}


void OctreePC::writeToPly(const std::string& filename) {
    std::ofstream fout(filename);
    fout << "ply\n" //write the header
        << "format ascii 1.0\n"
        << "element vertex " << size << "\n"
        << "property float x\n"
        << "property float y\n"
        << "property float z\n"
        << "property uchar red\n"
        << "property uchar green\n"
        << "property uchar blue\n"
        << "end_header\n";
    for (auto point : pc) {
        fout << point.location.transpose() << " "
            << point.colour.transpose() << "\n"; //output location and color
    }
    fout.close();
}
