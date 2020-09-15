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

    //extend bounding box slightly to fit voxels perfectly
    double x_diff = x_voxels * voxel_size - (XL - XS);
    XS -= x_diff / 2;
    XL += x_diff / 2;
    double y_diff = y_voxels * voxel_size - (YL - YS);
    YS -= y_diff / 2;
    YL += y_diff / 2;
    double z_diff = z_voxels * voxel_size - (ZL - ZS);
    YS -= z_diff / 2;
    YL += z_diff / 2;

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
    //account for exact divisions
    if (x == x_voxels)
        x -= 1;
    if (y == y_voxels)
        y -= 1;
    if (z == z_voxels)
        z == 1;
    return {x, y, z};
}


std::vector<size_t> UniformPC::planePoints(Eigen::Hyperplane<double, 3> thisPlane) {
    //get first edge intersection
    int edge = 0;
    Eigen::Vector3d p1 = edges[edge].intersectionPoint(thisPlane);
    // while the found point is out of bounds find another
    while (p1[0] < XS || p1[0] > XL || p1[1] < YS || p1[1] > YL || p1[2] < ZS || p1[2] > ZL) p1 = edges[++edge].intersectionPoint(thisPlane);
    //get second edge intersection
    int edge2 = edge + 1;
    Eigen::Vector3d p2 = edges[edge2].intersectionPoint(thisPlane);
    Eigen::Vector3d raydir = { 0, 0, 0 };
    // while the found point is out of bounds or the line between them will be on a wall with an x normal, find another point
    while (abs(raydir[0]) == 1 || p2[0] < XS || p2[0] > XL || p2[1] < YS || p2[1] > YL || p2[2] < ZS || p2[2] > ZL) {
        p2 = edges[++edge2].intersectionPoint(thisPlane);
        if (edges[edge].direction() == edges[edge2].direction()) // edges are parallel
            raydir = edges[edge].direction().cross((edges[edge2].origin() - edges[edge].origin()).normalized());
        else raydir = edges[edge].direction().cross(edges[edge2].direction());
        if (std::min(edges[edge].origin()[0], edges[edge2].origin()[0]) == XL 
            || std::min(edges[edge].origin()[1], edges[edge2].origin()[1]) == YL 
            || std::min(edges[edge].origin()[2], edges[edge2].origin()[2]) == ZL) raydir = -1*raydir;
    }

    //line to cast rays from
    Eigen::ParametrizedLine<double, 3> start_line = Eigen::ParametrizedLine<double, 3>::Through(p1, p2);
    // rays are cast one voxel size apart to ensure with the overlap in cleary's, all nearby voxels will be traversed
    Eigen::Vector3d step = start_line.direction() * voxel_size;
    // these rays should be fixed in the x direction, and vary in the y and z directions to fit the plane
    Eigen::Vector3d fixed = {1, 0, 0};
    fixed = thisPlane.normal().cross(fixed);
    // make sure cross product is oriented into the bounding box
    if ((raydir[1] == 1 && fixed[1] < 0) || raydir[1] == -1 && fixed[1] > 0) fixed[1] = -fixed[1];
    else if ((raydir[2] == 1 && fixed[2] < 0) || raydir[2] == -1 && fixed[2] > 0) fixed[2] = -fixed[2];
    raydir = fixed.normalized();
    //test
    double testy = thisPlane.absDistance(p2 + raydir);

    // indexes of points on the plane to be returned
    std::vector<size_t> indexes;
    // 3D truth array of visited voxels
    std::vector<std::vector<std::vector<bool>>> visited;
    visited = std::vector<std::vector<std::vector<bool>>>(x_voxels, std::vector<std::vector<bool>>(y_voxels, std::vector<bool>(z_voxels, false)));
    //move point inwards a little to avoid prescision error in the out-of-bounds checking
    //p1 = p1 + threshold * start_line.direction();
    //p1 = p1 + threshold * raydir;

    while (p1[0] >= XS && p1[0] <= XL && p1[1] >= YS && p1[1] <= YL && p1[2] >= ZS && p1[2] <= ZL) {
        //testing
        //std::vector<size_t> c = hashCell(p1);
        cleary(indexes, p1, raydir, visited, thisPlane);
        p1 += step;
    } 

    // testing
    signed long long i = 0;
#pragma omp parallel for
    for (i = 0; i < pc.size(); ++i) {
        if (thisPlane.absDistance(pc[i].location) < threshold) {
            std::vector<size_t> loc = hashCell(pc[i].location);
            if (!visited[loc[0]][loc[1]][loc[2]])
#pragma omp critical
                std::cout << loc[0] << loc[1] << loc[2] << " distance is " << thisPlane.absDistance(pc[i].location) << std::endl;
        }
    }

    return indexes;

}


//Adds points in place within the threshold of a 2D ray in a 3D bounding box
void UniformPC::cleary(std::vector<size_t> &points, Eigen::Vector3d p, Eigen::Vector3d dir, 
    std::vector<std::vector<std::vector<bool>>> &visited, Eigen::Hyperplane<double, 3> thisPlane) {

    std::vector<size_t> cell = hashCell(p);
    double theta_y = voxel_size / dir[1];
    double theta_z = voxel_size / dir[2];
    // set distance already travelled to how far along p is in the starting cell
    double dy = p[1] - (cell[1]*y_voxels + YS);
    double dz = p[2] - (cell[2]*z_voxels + ZS); // p distance minus smaller wall of cell
    bool x_up = cell[0] < x_voxels - 1;
    bool x_down = cell[0] > 0;
    size_t y = 0;
    size_t z = 0;

    //check the first cell
    if (!cells[cell[0]][cell[1]][cell[2]].empty()) //push thresholded points from cell onto points vector
            addPoints(cells[cell[0]][cell[1]][cell[2]], points, thisPlane);
    visited[cell[0]][cell[1]][cell[2]] = true;
    //visit cells to the left and right
    if (x_down) { //x negative
        if (!visited[cell[0] - 1][cell[1]][cell[2]]) {
            if (!cells[cell[0] - 1][cell[1]][cell[2]].empty())
                addPoints(cells[cell[0] - 1][cell[1]][cell[2]], points, thisPlane);
            visited[cell[0] - 1][cell[1]][cell[2]] = true;
        }
    } 
    if (x_up) { // x positive
        if (!visited[cell[0] + 1][cell[1]][cell[2]]) {
            if (!cells[cell[0] + 1][cell[1]][cell[2]].empty())
                addPoints(cells[cell[0] + 1][cell[1]][cell[2]], points, thisPlane);
            visited[cell[0] + 1][cell[1]][cell[2]] = true;
        }
    }

    // the >= 0 check isn't needed because a size_t will simply overflow to greater than the limit when negative anyway!
    while (cell[1] < y_voxels && cell[2] < z_voxels) {
        //note that due to the adjacent cell checking getting ahead of itself, actually checking the current cell is not needed.

        //z negative
        if (cell[2] > 0) {
            z = cell[2] - 1;
            if (!visited[cell[0]][cell[1]][z]) {
                if (!cells[cell[0]][cell[1]][z].empty())
                    addPoints(cells[cell[0]][cell[1]][z], points, thisPlane);
                visited[cell[0]][cell[1]][z] = true;
                //-x-z
                if (x_down) {
                    if (!visited[cell[0] - 1][cell[1]][z]) {
                        if (!cells[cell[0] - 1][cell[1]][z].empty())
                            addPoints(cells[cell[0] - 1][cell[1]][z], points, thisPlane);
                        visited[cell[0] - 1][cell[1]][z] = true;
                    }
                } //+x-z
                if (x_up) {
                    if (!visited[cell[0] + 1][cell[1]][z]) {
                        if (!cells[cell[0] + 1][cell[1]][z].empty())
                            addPoints(cells[cell[0] + 1][cell[1]][z], points, thisPlane);
                        visited[cell[0] + 1][cell[1]][z] = true;
                    }
                }
            }
        } //z positive
        if (cell[2] < z_voxels - 1) {
            z = cell[2] + 1;
            if (!visited[cell[0]][cell[1]][z]) {
                if (!cells[cell[0]][cell[1]][z].empty())
                    addPoints(cells[cell[0]][cell[1]][z], points, thisPlane);
                visited[cell[0]][cell[1]][z] = true;
                //-x+z
                if (x_down) {
                    if (!visited[cell[0] - 1][cell[1]][z]) {
                        if (!cells[cell[0] - 1][cell[1]][z].empty())
                            addPoints(cells[cell[0] - 1][cell[1]][z], points, thisPlane);
                        visited[cell[0] - 1][cell[1]][z] = true;
                    }
                } //+x+z
                if (x_up) {
                    if (!visited[cell[0] + 1][cell[1]][z]) {
                        if (!cells[cell[0] + 1][cell[1]][z].empty())
                            addPoints(cells[cell[0] + 1][cell[1]][z], points, thisPlane);
                        visited[cell[0] + 1][cell[1]][z] = true;
                    }
                }
            }
        }
        //y negative
        if (cell[1] > 0) {
            y = cell[1] - 1;
            if (!visited[cell[0]][y][cell[2]]) {
                if (!cells[cell[0]][y][cell[2]].empty())
                    addPoints(cells[cell[0]][y][cell[2]], points, thisPlane);
                visited[cell[0]][y][cell[2]] = true;
                //-x-y
                if (x_down) {
                    if (!visited[cell[0] - 1][y][cell[2]]) {
                        if (!cells[cell[0] - 1][y][cell[2]].empty())
                            addPoints(cells[cell[0] - 1][y][cell[2]], points, thisPlane);
                        visited[cell[0] - 1][y][cell[2]] = true;
                    }
                } //+x-y
                if (x_up) {
                    if (!visited[cell[0] + 1][y][cell[2]]) {
                        if (!cells[cell[0] + 1][y][cell[2]].empty())
                            addPoints(cells[cell[0] + 1][y][cell[2]], points, thisPlane);
                        visited[cell[0] + 1][y][cell[2]] = true;
                    }
                }
            }
        } //y positive
        if (cell[1] < y_voxels - 1) {
            y = cell[1] + 1;
            if (!visited[cell[0]][y][cell[2]]) {
                if (!cells[cell[0]][y][cell[2]].empty())
                    addPoints(cells[cell[0]][y][cell[2]], points, thisPlane);
                visited[cell[0]][y][cell[2]] = true;
                //-x+y
                if (x_down) {
                    if (!visited[cell[0] - 1][y][cell[2]]) {
                        if (!cells[cell[0] - 1][y][cell[2]].empty())
                            addPoints(cells[cell[0] - 1][y][cell[2]], points, thisPlane);
                        visited[cell[0] - 1][y][cell[2]] = true;
                    }
                } //+x+y
                if (x_up) {
                    if (!visited[cell[0] + 1][y][cell[2]]) {
                        if (!cells[cell[0] + 1][y][cell[2]].empty())
                            addPoints(cells[cell[0] + 1][y][cell[2]], points, thisPlane);
                        visited[cell[0] + 1][y][cell[2]] = true;
                    }
                }
            }
        }

        //work out next cell
        if (abs(dy) < abs(dz)) {
            dy += theta_y;  // going to the y adjacent cell
            (theta_y > 0) ? cell[1] += 1 : cell[1] -= 1;
        }
        else {
            dz += theta_z;  // going to the z adjacent cell
            (theta_z > 0) ? cell[2] += 1 : cell[2] -= 1;
        }
    }

    return;
}


void UniformPC::addPoints(std::vector<size_t> indexes, std::vector<size_t> &thisPoints, Eigen::Hyperplane<double, 3> thisPlane) {
    
    //OpenMP requires signed integrals for its loop variables... interesting
    signed long long i = 0;
#pragma omp parallel for shared(thisPoints) private (i)
    for (i = 0; i < indexes.size(); ++i) {
        if (thisPlane.absDistance(pc[indexes[i]].location) < threshold && std::binary_search(remainingPoints.begin(), remainingPoints.end(), indexes[i]))
#pragma omp critical
            thisPoints.push_back(indexes[i]);
    }
    comparisons += indexes.size();
}
