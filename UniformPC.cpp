#include "UniformPC.h"

#include <iostream>
#include <fstream>
#include <sstream>


UniformPC::UniformPC(PointCloud const&p, float voxel_scale) : PointCloud(p) {

    // assign the given volumes of voxels to the model dimensions
    voxel_size = voxel_scale*threshold;
    x_voxels = ceil((XL - XS) / voxel_size);
    y_voxels = ceil((YL - YS) / voxel_size);
    z_voxels = ceil((ZL - ZS) / voxel_size);
    // a pointer to a vector is used so the point list vectors can be stored elsewhere
    // check how best to order these for access time
    //cells = std::vector<std::vector<std::vector<std::vector<size_t>*>>>(x_voxels, std::vector<std::vector<std::vector<size_t>*>>(y_voxels, std::vector<std::vector<size_t>*>(z_voxels)));
    cells = std::vector<std::vector<std::vector<std::vector<size_t>>>>(x_voxels, std::vector<std::vector<std::vector<size_t>>>(y_voxels, std::vector<std::vector<size_t>>(z_voxels)));

    std::cout << "Doing uniform space subdivision of " << x_voxels << " by " << y_voxels << " by " << z_voxels << " voxels..." << std::endl;

    std::vector<size_t> cell;
    // for each point, hash their index in the vector into a cell
    for (size_t i = 0; i < pc.size(); ++i) {
        cell = hashCell(pc[i].location);
        cells[cell[0]][cell[1]][cell[2]].push_back(i);
    }
}


std::vector<size_t> UniformPC::hashCell(Eigen::Vector3d p) {
    size_t x = floor((p[0] - XS) / voxel_size);
    size_t y = floor((p[1] - YS) / voxel_size);
    size_t z = floor((p[2] - ZS) / voxel_size);
    //CHECK
    if (x > x_voxels)
        std::cout << "uh oh spaghettios" << std::endl;
    if (y > y_voxels)
        std::cout << "uh oh spaghettios" << std::endl;
    if (z > z_voxels)
        std::cout << "uh oh spaghettios" << std::endl;
    return {x, y, z};
}


// Returns a vector of points within the threshold to the given hyperplane
std::vector<size_t> UniformPC::planePoints(Eigen::Hyperplane<double, 3> thisPlane, std::vector<size_t> remainingPoints, unsigned int trial, int plane) {
    // indexes of points on the plane to be returned
    std::vector<size_t> indexes;
    // 3D truth array of visited voxels
    std::vector<std::vector<std::vector<bool>>> visited;
    visited = std::vector<std::vector<std::vector<bool>>>(x_voxels, std::vector<std::vector<bool>>(y_voxels, std::vector<bool>(z_voxels, false)));

    Eigen::ParametrizedLine<double, 3>* start_line = PointCloud::intersectPlanes(thisPlane, Eigen::Hyperplane<double, 3>(thisPlane));

    // direction of rays to be cast into the bounding box
    Eigen::Vector3d norm = { 0, 0, 0 };
    if (start_line != NULL && abs(start_line->direction()[0]) < 0.001) {
        if (abs(start_line->origin()[0] - XL) < 0.001) // positive side of bounding box
            norm[0] = -1;
        else norm[0] = 1; // negative side of bounding box
    }
    else if (start_line != NULL && abs(start_line->direction()[1]) < 0.001) {
        if (abs(start_line->origin()[0] - XL) < 0.001) // positive side of bounding box
            norm[1] = -1;
        else norm[1] = 1; // negative side of bounding box
    }
    else if (start_line != NULL && abs(start_line->direction()[2]) < 0.001) {
        if (abs(start_line->origin()[0] - XL) < 0.001) // positive side of bounding box
            norm[2] = -1;
        else norm[2] = 1; // negative side of bounding box
    }

    // point along the line to start with
    Eigen::Vector3d p = start_line->origin();
    Eigen::Vector3d step = voxel_size * start_line->direction(); //should this be threshold or voxel size?
    // cast first ray one step in?
    p += step;
    // make sure correct direction to move p along the line into the bounding box
    if (p[0] < XS || p[0] > XL || p[1] < YS || p[1] > YL || p[2] < ZS || p[2] > ZL) { // assuming p's location on the starting axis isn't out by prescision?
        step *= -1;
        p += 2 * step;
    }

    do {
        indexes = cleary(indexes, remainingPoints, Eigen::ParametrizedLine<double, 3>(p, norm), visited, thisPlane, plane);
        p += step;
    } while (p[0] < XS || p[0] > XL || p[1] < YS || p[1] > YL || p[2] < ZS || p[2] > ZL); //p stays within the other bounds

    return indexes;

}


//Returns the points within the threshold of a ray in a 3D bounding box
std::vector<size_t> UniformPC::cleary(std::vector<size_t> &points, std::vector<size_t> remainingPoints, Eigen::ParametrizedLine<double, 3> ray, 
    std::vector<std::vector<std::vector<bool>>> &visited, Eigen::Hyperplane<double, 3> thisPlane, int plane) {
    
    Eigen::Vector3d p = ray.origin();
    //find current cell
    std::vector<size_t> cell = hashCell(p);
    //check these are correct?
    double theta_x = voxel_size / ray.direction()[0];
    double theta_y = voxel_size / ray.direction()[1];
    double theta_z = voxel_size / ray.direction()[2];
    double dx = theta_x;
    double dy = theta_y;
    double dz = theta_z;

    do {
        if (!visited[cell[0]][cell[1]][cell[2]] && !cells[cell[0]][cell[1]][cell[2]].empty()) {
            //push thresholded points from cell onto points vector
            addPoints(cells[cell[0]][cell[1]][cell[2]], points, remainingPoints, thisPlane, plane);
            visited[cell[0]][cell[1]][cell[2]] = true;
        }
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


void UniformPC::addPoints(std::vector<size_t> indexes, std::vector<size_t> &thisPoints, std::vector<size_t> remainingPoints,
    Eigen::Hyperplane<double, 3> thisPlane, int plane) {
    
    //OpenMP requires signed integrals for its loop variables... interesting
    signed long long i = 0;
#pragma omp parallel for shared(thisPoints) private (i)
    for (i = 0; i < remainingPoints.size(); ++i) {
        if (thisPlane.absDistance(pc[remainingPoints[i]].location) < threshold)
#pragma omp critical
            thisPoints.push_back(remainingPoints[i]);
    }
}
