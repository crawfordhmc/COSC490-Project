#include "UniformPC.h"

#include <iostream>
#include <fstream>
#include <sstream>


UniformPC::UniformPC(const std::string& filepath) : PointCloud(filepath) {

    // assign the given volumes of voxels to the model dimensions
    voxel_size = PointCloud::threshold;
    x_voxels = ceil((XL - XS) / voxel_size);
    y_voxels = ceil((YL - YS) / voxel_size);
    z_voxels = ceil((ZL - ZS) / voxel_size);
    // a pointer to a vector is used so the point vectors can be stored elsewhere
    // check how best to order these
    cells = std::vector<std::vector<std::vector<std::vector<size_t>*>>>(x_voxels, std::vector<std::vector<std::vector<size_t>*>>(y_voxels, std::vector<std::vector<size_t>*>(z_voxels)));
    //std::vector<std::vector<std::vector<std::vector<size_t>*>>> cells(x_voxels, std::vector<std::vector<std::vector<size_t>*>>(y_voxels, std::vector<std::vector<size_t>*>(z_voxels)));

    // for each point, hash their index in the vector into a cell
    for (size_t i = 0; i < pc.size(); ++i) {
        Eigen::Matrix<size_t, 3, 1> cell = hashCell(pc[i].location);
        cells[cell[0]][cell[1]][cell[2]]->push_back(i);
    }

}


Eigen::Matrix<size_t, 3, 1> UniformPC::hashCell(Eigen::Vector3d p) {
    size_t x = floor((p[0] - XS) / x_voxels);
    size_t y = floor((p[1] - YS) / y_voxels);
    size_t z = floor((p[2] - ZS) / z_voxels);
    //CHECK
    if (p[0] < XS + x * voxel_size || XS + (x + 1) * voxel_size <= p[0])
        std::cout << "uh oh spaghettios" << std::endl;
    if (p[1] < YS + y * voxel_size || YS + (y + 1) * voxel_size <= p[1])
        std::cout << "uh oh spaghettios" << std::endl;
    if (p[2] < ZS + z * voxel_size || ZS + (z + 1) * voxel_size <= p[2])
        std::cout << "uh oh spaghettios" << std::endl;
    return {x, y, z};
}

// these three may not need to be virtual functions depending on access methods
PointCloud::Point UniformPC::getPoint(int index) { return pc[index]; }
void UniformPC::setPointPlane(int index, int planeID) { pc[index].planeIx = planeID; }
void UniformPC::setPointColour(int index, Eigen::Vector3i colour) { pc[index].colour = colour; }


// Returns a vector of points within the threshold to the given hyperplane
// (also prints the number of threads being used for the calculations)
std::vector<size_t> UniformPC::planePoints(Eigen::Hyperplane<double, 3> thisPlane) {
    std::vector<size_t> indexes;
    // if an origin intersecting with another axis is found, this could be a nice origin to do the points from
    Eigen::ParametrizedLine<double, 3>* start_line = PointCloud::intersectPlanes(thisPlane, Eigen::Hyperplane<double, 3>(thisPlane));
    Eigen::Vector3d tangent; // unit vector tangent to start line pointing towards bounding box
    Eigen::Vector3d p = start_line->origin();
    do {
        cleary(indexes, Eigen::ParametrizedLine<double, 3>(p, tangent));
        p += threshold * start_line->direction();
    } while (); //p stays within the other bounds
    return indexes;

}


void UniformPC::addPoints(std::vector<size_t>* indexes, std::vector<size_t> thisPoints, Eigen::ParametrizedLine<double, 3> ray) {
    //OpenMP requires signed integrals for its loop variables... interesting
    signed long long i = 0;
#pragma omp parallel for shared(thisPoints) private (i)
    for (i = 0; i < indexes->size(); ++i) {
        if (abs(ray.distance(pc[indexes->at(i)].location)) < threshold) //absolute?
#pragma omp critical
            thisPoints.push_back(indexes->at(i));
    }
}


//Returns the points within the threshold of a ray in a 3D bounding box
std::vector<size_t> UniformPC::cleary(std::vector<size_t> points, Eigen::ParametrizedLine<double, 3> ray) {
    Eigen::Vector3d p = ray.origin();
    //find current cell
    Eigen::Matrix<size_t, 3, 1> cell = hashCell(p);
    //check these are correct?
    double theta_x = voxel_size / ray.direction()[0];
    double theta_y = voxel_size / ray.direction()[1];
    double theta_z = voxel_size / ray.direction()[2];
    double dx = theta_x;
    double dy = theta_y;
    double dz = theta_z;

    do {
        if (!cells[cell[0]][cell[1]][cell[2]]->empty())
            //push thresholded points from cell onto points vector
            addPoints(cells[cell[0]][cell[1]][cell[2]], points, ray); // does adding the points need to be with a whole plane or just the ray?
        //work out next cell
        if (abs(dx) < abs(dy) && abs(dx) < abs(dz)) {
            dx += theta_x;
            (theta_x > 0) ? cell[0] += 1 : cell[0] -= 1;
        }
        else if (abs(dy) < abs(dx) && abs(dy) < abs(dz)) {
            dy += theta_y;
            (theta_y > 0) ? cell[1] += 1 : cell[1] -= 1;
        }
        else {
            dz += theta_z;
            (theta_z > 0) ? cell[2] += 1 : cell[2] -= 1;
        }
    } while (cell[0] >= 0 && cell[0] <= x_voxels || cell[1] >= 0 && cell[1] <= y_voxels || cell[2] >= 0 && cell[2] <= z_voxels);
    return points;
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
