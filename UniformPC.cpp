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
    size_t x = (p[0] - XS < voxel_size) ? 0 : floor((p[0] - XS) / voxel_size);
    size_t y = (p[1] - YS < voxel_size) ? 0 : floor((p[1] - YS) / voxel_size);
    size_t z = (p[2] - ZS < voxel_size) ? 0 : floor((p[2] - ZS) / voxel_size);
    //CHECK
    if (x > x_voxels)
        x -= 1;
        //std::cout << "uh oh spaghettios" << std::endl;
    if (y > y_voxels)
        y -= 1;
        //std::cout << "uh oh spaghettios" << std::endl;
    if (z > z_voxels)
        z -= 1;
        //std::cout << "uh oh spaghettios" << std::endl;
    return {x, y, z};
}


// Returns a vector of points within the threshold to the given hyperplane
std::vector<size_t> UniformPC::planePoints(Eigen::Hyperplane<double, 3> thisPlane, std::vector<size_t> remainingPoints, unsigned int trial, int plane) {
    //theta value for floating-point comparisons NEEDED?
    double t = 1;// +0.001 * threshold;
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

        Eigen::ParametrizedLine<double, 3>({ XS, YL, ZL}, { 1, 0, 0 }),
        Eigen::ParametrizedLine<double, 3>({ XL, YS, ZL}, { 0, 1, 0 }),
        Eigen::ParametrizedLine<double, 3>({ XL, YL, ZS}, { 0, 0, 1 })
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
        //make sure the lines are on the same plane
        if (edges[edge].origin()[0] != edges[i].origin()[0] && edges[edge].origin()[1] != edges[i].origin()[1] && edges[edge].origin()[2] != edges[i].origin()[2]) continue; //FIX THIS LINE
        p2 = edges[i].intersectionPoint(thisPlane);
    }

    //line to cast rays from
    Eigen::ParametrizedLine<double, 3> start_line = Eigen::ParametrizedLine<double, 3>::Through(p1, p2);
    // direction of rays to be cast into the bounding box, following plane
    Eigen::Vector3d dir = { 1, 1, 1 };
    // if plane intersects on the far corner
    if (edge > 8) dir = dir * -1;
    dir = dir - edges[edge].direction();
    if (edges[edge].direction() == edges[i].direction()) //if lines are parallel
        dir = dir - (edges[i].origin() - edges[edge].origin()).normalized();
    else dir = dir - edges[i].direction();
    // if plane intersects on one of the far sides
    if (edge > 2) dir = dir * -1;
    // set ray direction to be moving as opposed to fixed for 2D clearys as the slower-changing one of the start line
    int moving_axis;
    if (abs(dir[0]) == 1) moving_axis = (abs(start_line.direction()[1]) < abs(start_line.direction()[2])) ? 1 : 2;
    else if (abs(dir[1]) == 1) moving_axis = (abs(start_line.direction()[0]) < abs(start_line.direction()[2])) ? 0 : 2;
    else moving_axis = (abs(start_line.direction()[0]) < abs(start_line.direction()[1])) ? 0 : 1;

    Eigen::Vector3d step = start_line.direction() * voxel_size;
    dir[moving_axis] -= thisPlane.normal()[moving_axis];
    //for (int i = 0; i < 3; i++) if (abs(dir[i]) != 1) dir[i] -= thisPlane.normal()[i];
    do { //p stays within the other bounds
        cleary(indexes, remainingPoints, Eigen::ParametrizedLine<double, 3>(p1, dir.normalized()), visited, thisPlane, plane);
        p1 += step;
    } while (start_line.projection(p1).norm() < start_line.projection(p2).norm());
    return indexes;

}


//Returns the points within the threshold of a 2D ray in a 3D bounding box
std::vector<size_t> UniformPC::cleary(std::vector<size_t> &points, std::vector<size_t> remainingPoints, Eigen::ParametrizedLine<double, 3> ray, 
    std::vector<std::vector<std::vector<bool>>> &visited, Eigen::Hyperplane<double, 3> thisPlane, int plane) {
    
    std::vector<int> limits = { x_voxels, y_voxels, z_voxels };
    Eigen::Vector3d p = ray.origin();
    int x1 = ray.direction()[0] == 0 ? 1 : 0;
    int x2 = ray.direction()[2] == 0 ? 1 : 2;
    //int x3 = (x1 == 0) ? (x2 == 1) ? 2 : 1 : 0;
    //double t = 1;// +0.001 * threshold;
    //if (p[0] < XS * t || p[0] > XL * t || p[1] < YS * t || p[1] > YL * t || p[2] < ZS * t || p[2] > ZL * t)
    //    t = 1;

    //find current cell
    std::vector<size_t> cell = hashCell(p);
    //visit x2 adjacent cells just for the first step to catch scragglers
    std::vector<size_t> cell2 = cell;
    cell2[x2] += 1;
    if (cell2[x2] < limits[x2]) {
        if (!cells[cell2[0]][cell2[1]][cell2[2]].empty())
            addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, remainingPoints, thisPlane, plane);
        visited[cell2[0]][cell2[1]][cell2[2]] = true;
    }
    if (cell2[x2] > 1) {
        cell2[x2] -= 2;
        if (!cells[cell2[0]][cell2[1]][cell2[2]].empty())
            addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, remainingPoints, thisPlane, plane);
        visited[cell2[0]][cell2[1]][cell2[2]] = true;
    }

    //size_t up = cell[x3] + 1;
    //size_t down = cell[x3] - 1;
    //check these are correct?
    double theta_x1 = voxel_size / ray.direction()[x1];
    double theta_x2 = voxel_size / ray.direction()[x2];
    double dx1 = theta_x1;
    double dx2 = theta_x2;

    // the >= 0 check isn't needed because a size_t will simply overflow to greater than the limit when negative anyway!
    while (cell[x1] < limits[x1] && cell[x2] < limits[x2]) {
        if (!visited[cell[0]][cell[1]][cell[2]] && !cells[cell[0]][cell[1]][cell[2]].empty()) {
            //push thresholded points from cell onto points vector
            addPoints(cells[cell[0]][cell[1]][cell[2]], points, remainingPoints, thisPlane, plane);
        }
        visited[cell[0]][cell[1]][cell[2]] = true;

        //visit x1 adjacent cells to catch scragglers
        cell2 = cell;
        cell2[x1] += 1;
        if (cell2[x1] < limits[x1]) {
            if (!visited[cell2[0]][cell2[1]][cell2[2]] && !cells[cell2[0]][cell2[1]][cell2[2]].empty())
                addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, remainingPoints, thisPlane, plane);
            visited[cell2[0]][cell2[1]][cell2[2]] = true;
        }
        if (cell2[x1] > 1) {
            cell2[x1] -= 2;
            if (!visited[cell2[0]][cell2[1]][cell2[2]] && !cells[cell2[0]][cell2[1]][cell2[2]].empty())
                addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, remainingPoints, thisPlane, plane);
            visited[cell2[0]][cell2[1]][cell2[2]] = true;
        }

        //work out next cell
        if (abs(dx1) < abs(dx2)) {
            dx1 += theta_x1;
            (theta_x1 > 0) ? cell[x1] += 1 : cell[x1] -= 1;
        }
        else {
            dx2 += theta_x2;
            (theta_x2 > 0) ? cell[x2] += 1 : cell[x2] -= 1;
        }
    }

    //wind back last out of bounds step
    if (cell[x1] > limits[x1]) cell[x1] = 0;
    if (cell[x1] == limits[x1]) cell[x1] = limits[x1] - 1;
    if (cell[x2] > limits[x2]) cell[x2] = 0;
    if (cell[x2] == limits[x2]) cell[x2] = limits[x2] - 1;

    //visit x2 adjacent cells just for the last step to catch scragglers
    cell2 = cell;
    cell2[x2] += 1;
    if (cell2[x2] < limits[x2]) {
        if (!visited[cell2[0]][cell2[1]][cell2[2]] && !cells[cell2[0]][cell2[1]][cell2[2]].empty())
            addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, remainingPoints, thisPlane, plane);
        visited[cell2[0]][cell2[1]][cell2[2]] = true;
    }
    if (cell2[x2] > 1) {
        cell2[x2] -= 2;
        if (!visited[cell2[0]][cell2[1]][cell2[2]] && !cells[cell2[0]][cell2[1]][cell2[2]].empty())
            addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, remainingPoints, thisPlane, plane);
        visited[cell2[0]][cell2[1]][cell2[2]] = true;
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
