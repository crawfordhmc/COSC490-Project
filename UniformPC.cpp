#include "UniformPC.h"

#include <iostream>
#include <fstream>
#include <sstream>

UniformPC::UniformPC(const std::string& filepath) : PointCloud(filepath) {
    // find bounding coordinates
    // assign the given volumes of voxels to the model dimensions
    // each cell is defined by the positive octrant corner point + the voxel size
    // for each point, hash their index in the vector into a cell

    // cell are progressed through by a plane with cleary's algorithm
    Point defaultPoint;
    defaultPoint.location = Eigen::Vector3d::Zero();
    defaultPoint.colour = Eigen::Vector3i::Zero();
    defaultPoint.planeIx = -1;
    pc.resize(verts.size(), defaultPoint);
    signed long long p;
#pragma omp parallel for private(p)
    for (p = 0; p < pc.size(); ++p) {
        pc[p].location = verts[p].cast<double>();
        pc[p].colour = cols[p].cast<int>();
    }

    size = pc.size();
}

// these three may not need to be virtual functions depending on access methods
PointCloud::Point UniformPC::getPoint(int index) { return pc[index]; }
void UniformPC::setPointPlane(int index, int planeID) { pc[index].planeIx = planeID; }
void UniformPC::setPointColour(int index, Eigen::Vector3i colour) { pc[index].colour = colour; }


// Returns a vector of points within the threshold to the given hyperplane
// (also prints the number of threads being used for the calculations)
std::vector<size_t> UniformPC::planePoints(Eigen::Hyperplane<double, 3> thisPlane, unsigned int trial, float threshold, int plane) {
    std::vector<size_t> thisPoints;
    int threads = 0;
    //OpenMP requires signed integrals for its loop variables... interesting
    signed long long i = 0;
#pragma omp parallel for shared(thisPoints) private (i)
    for (i = 0; i < pc.size(); ++i) {
        if (thisPlane.absDistance(pc[i].location) < threshold)
#pragma omp critical
            thisPoints.push_back(i);
        if (omp_get_thread_num() == 0 && trial == 0 && plane == 0)
            threads = omp_get_max_threads();
    }
    if (trial == 0 && plane == 0)
        std::cout << threads << " threads are being used" << std::endl;
    return thisPoints;
}


float UniformPC::threshold(float scale_parameter) {
    // Determine the threshold as a % of model size
    // (coordinate center is all over the place, so biggest/smallest signed point difference gives bounding box)
    double xs = pc[0].location[0];
    double xl = pc[0].location[0];
    double ys = pc[0].location[1];
    double yl = pc[0].location[1];
    double zs = pc[0].location[2];
    double zl = pc[0].location[2];
    // chunk parallelize this if its slow?
    for (size_t i = 1; i < pc.size(); i++) {
        xs = std::min(xs, pc[i].location[0]);
        xl = std::max(xl, pc[i].location[0]);
        ys = std::min(ys, pc[i].location[1]);
        yl = std::max(yl, pc[i].location[1]);
        zs = std::min(zs, pc[i].location[2]);
        zl = std::max(zl, pc[i].location[2]);
    }
    // get x/y/z difference and compute average scale factor for the model
    double scale = (xl - xs + yl - ys + zl - zs) / 3;
    // apply a small % to the value to get a sensible threshold
    return scale_parameter * scale;
}


void UniformPC::writeToPly(const std::string& filename) {
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
