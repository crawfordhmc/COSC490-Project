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
    //theta value for floating-point comparisons NEEDED?
    double t = 1 + 0.01 * threshold;
    std::vector<Eigen::ParametrizedLine<double, 3>> edges = {
        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZS }, { 1, 0, 0 }),
        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZS }, { 0, 1, 0 }),
        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZS }, { 0, 0, 1 }),

        Eigen::ParametrizedLine<double, 3>({ XL, YS, ZS }, { 0, 1, 0 }),
        Eigen::ParametrizedLine<double, 3>({ XL, YS, ZS }, { 0, 0, 1 }),

        Eigen::ParametrizedLine<double, 3>({ XS, YL, ZS }, { 1, 0, 0 }),
        Eigen::ParametrizedLine<double, 3>({ XS, YL, ZS }, { 0, 0, 1 }),

        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZL }, { 1, 0, 0 }),
        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZL }, { 0, 1, 0 }),

        Eigen::ParametrizedLine<double, 3>({ XL, YL, ZL}, { -1, 0, 0 }),
        Eigen::ParametrizedLine<double, 3>({ XL, YL, ZL}, { 0, -1, 0 }),
        Eigen::ParametrizedLine<double, 3>({ XL, YL, ZL}, { 0, 0, -1 })
    };
    // indexes of points on the plane to be returned
    std::vector<size_t> indexes;
    // 3D truth array of visited voxels
    std::vector<std::vector<std::vector<bool>>> visited;
    visited = std::vector<std::vector<std::vector<bool>>>(x_voxels, std::vector<std::vector<bool>>(y_voxels, std::vector<bool>(z_voxels, false)));

    //get first edge intersection
    int edge = 0;
    Eigen::Vector3d p1 = edges[edge].intersectionPoint(thisPlane);
    while (p1[0] < XS * t || p1[0] > XL * t || p1[1] < YS * t || p1[1] > YL * t || p1[2] < ZS * t || p1[2] > ZL * t) p1 = edges[++edge].intersectionPoint(thisPlane);
    //get second edge intersection
    int i = edge + 1;
    Eigen::Vector3d p2 = edges[i].intersectionPoint(thisPlane);
    while (p2[0] < XS * t || p2[0] > XL * t || p2[1] < YS * t || p2[1] > YL * t || p2[2] < ZS * t || p2[2] > ZL * t) {  //end loop once a intersection within the bounding box is found
        i += 1;
        //make sure the lines are adjacent
        if () continue; //FIX THIS LINE
        p2 = edges[i].intersectionPoint(thisPlane);
    }

    //line to cast rays from
    Eigen::ParametrizedLine<double, 3> start_line = Eigen::ParametrizedLine<double, 3>::Through(p1, p2);
    // direction of rays to be cast into the bounding box, following plane
    Eigen::Vector3d norm = { 1, 1, 1 };
    // if plane intersects on the far corner
    if (edge > 8) norm = norm * -1;
    norm = norm - edges[edge].direction();
    if (edges[edge].direction() == edges[i].direction()) //if lines are parallel
        norm = norm - (edges[i].origin() - edges[edge].origin()).normalized();
    else norm = norm - edges[i].direction();
    // if plane intersects on one of the far sides
    if (edge > 2) norm = norm * -1;

    Eigen::Vector3d step = start_line.direction() * voxel_size;
    Eigen::Vector3d dir = thisPlane.projection(norm).normalized(); //look into fixing a direction for 2d clearys
    do {
        cleary(indexes, remainingPoints, Eigen::ParametrizedLine<double, 3>(p1, dir), visited, thisPlane, plane);
        p1 += step;
    } while (start_line.projection(p1).norm() < start_line.projection(p2).norm()); //p stays within the other bounds

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

    while (cell[0] >= 0 && cell[0] < x_voxels && cell[1] >= 0 && cell[1] < y_voxels && cell[2] >= 0 && cell[2] < z_voxels) {
        if (!visited[cell[0]][cell[1]][cell[2]] && !cells[cell[0]][cell[1]][cell[2]].empty()) {
            //push thresholded points from cell onto points vector
            addPoints(cells[cell[0]][cell[1]][cell[2]], points, remainingPoints, thisPlane, plane);
        }
        visited[cell[0]][cell[1]][cell[2]] = true;
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
    }
    return points;
}


void UniformPC::addPoints(std::vector<size_t> indexes, std::vector<size_t> &thisPoints, std::vector<size_t> remainingPoints,
    Eigen::Hyperplane<double, 3> thisPlane, int plane) {
    
    //OpenMP requires signed integrals for its loop variables... interesting
    signed long long i = 0;
#pragma omp parallel for shared(thisPoints) private (i)
    for (i = 0; i < indexes.size(); ++i) {
        if (thisPlane.absDistance(pc[indexes[i]].location) < threshold && std::binary_search(remainingPoints.begin(), remainingPoints.end(), indexes[i]))
#pragma omp critical
            thisPoints.push_back(indexes[i]);
    }
}
