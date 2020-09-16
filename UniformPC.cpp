#include "UniformPC.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>


UniformPC::UniformPC(PointCloud const&p, int voxel_scale) : PointCloud(p) {

    XS = floor(XS);//-= x_diff / 2;
    XL = ceil(XL);//+= x_diff / 2;
    YS = floor(YS);//-= y_diff / 2;
    YL = ceil(YL);//+= y_diff / 2;
    ZS = floor(ZS);// -= z_diff / 2;
    ZL = ceil(ZL);// += z_diff / 2;

    // assign the given volumes of voxels to the model dimensions
    voxel_size = voxel_scale;
    bool big = true;
    while ((int)(XL - XS) % voxel_size != 0) {
        if (big) XL += 1, big = false;
        else XS -= 1, big = true;
    }
    while ((int)(YL - YS) % voxel_size != 0) {
        if (big) YL += 1, big = false;
        else YS -= 1, big = true;
    }
    while ((int)(ZL - ZS) % voxel_size != 0) {
        if (big) ZL += 1, big = false;
        else ZS -= 1, big = true;
    }

    x_voxels = (int)(XL - XS) / voxel_size;
    y_voxels = (int)(YL - YS) / voxel_size;
    z_voxels = (int)(ZL - ZS) / voxel_size;


    edges = {
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
    size_t x = (p[0] - XS) / voxel_size;
    size_t y = (p[1] - YS) / voxel_size;
    size_t z = (p[2] - ZS) / voxel_size;
    if (x > x_voxels - 1) x = x_voxels - 1;
    if (y > y_voxels - 1) y = y_voxels - 1;
    if (z > z_voxels - 1) z = z_voxels - 1;

    if (p[0] - XS > (x + 1) * voxel_size)
        std::cout << p[0] << std::endl;
    if (p[0] - XS < x * voxel_size)
        std::cout << p[0] << std::endl;

    if (p[1] - YS > (y + 1) * voxel_size)
        std::cout << p[1] << std::endl;
    if (p[1] - YS < y * voxel_size)
        std::cout << p[1] << std::endl;

    if (p[2] - ZS > (z + 1) * voxel_size)
        std::cout << p[2] << std::endl;
    if (p[2] - ZS < z * voxel_size)
        std::cout << p[2] << std::endl;

    return {x, y, z};
}


std::vector<size_t> UniformPC::planePoints(Eigen::Hyperplane<double, 3> thisPlane) {
    //get first edge intersection
    int edge = 0;
    Eigen::Vector3d p1 = edges[edge].intersectionPoint(thisPlane);
    // while the found point is out of bounds find another
    while (p1[0] < XS || p1[0] > XL || p1[1] < YS || p1[1] > YL || p1[2] < ZS || p1[2] > ZL) 
        p1 = edges[++edge].intersectionPoint(thisPlane);
    //get second edge intersection
    int edge2 = edge + 1;
    Eigen::Vector3d p2 = edges[edge2].intersectionPoint(thisPlane);
    //line to cast rays from
    Eigen::ParametrizedLine<double, 3> start_line = Eigen::ParametrizedLine<double, 3>::Through(p1, p2);
    Eigen::Vector3d norm = edges[edge].direction().cross(start_line.direction());
    while (norm[0] != 0 || p2[0] < XS || p2[0] > XL || p2[1] < YS || p2[1] > YL || p2[2] < ZS || p2[2] > ZL) { //point is out of bounds
        p2 = edges[++edge2].intersectionPoint(thisPlane);
        start_line = Eigen::ParametrizedLine<double, 3>::Through(p1, p2);
        norm = edges[edge].direction().cross(start_line.direction());
    }
    norm.normalize();
    //make sure norm is oriented into the bounding box
    if (abs(norm[1]) == 1 && (norm[1] + p1[1] < YS || norm[1] + p1[1] > YL)) norm[1] = -1 * norm[1];
    else if (norm[2] + p1[2] < ZS || norm[2] + p1[2] > ZL) norm[2] = -1 * norm[2];

    // these rays should be fixed in the x direction, and vary in the y and z directions to fit the plane
    Eigen::Vector3d raydir = {1, 0, 0};
    raydir = thisPlane.normal().cross(raydir);
    // make sure cross product is oriented into the bounding box
    if ((norm[1] == 1 && raydir[1] < 0) || norm[1] == -1 && raydir[1] > 0) raydir = -raydir;
    else if ((norm[2] == 1 && raydir[2] < 0) || norm[2] == -1 && raydir[2] > 0) raydir = -raydir;

    // indexes of points on the plane to be returned
    std::vector<size_t> indexes;
    // 3D truth array of visited voxels
    std::vector<std::vector<std::vector<bool>>> visited;
    visited = std::vector<std::vector<std::vector<bool>>>(x_voxels, std::vector<std::vector<bool>>(y_voxels, std::vector<bool>(z_voxels, false)));

    // calculate exits
    std::vector<size_t> cell = hashCell(p1);
    double theta_x = voxel_size / abs(start_line.direction()[0]);
    double dx = (start_line.direction()[0] > 0 ? (cell[0] + 1) * voxel_size - (p1[0] - XS) : p1[0] - XS - cell[0] * voxel_size) / abs(start_line.direction()[0]);
    double theta_yz = voxel_size / (norm[1] == 0 ? abs(start_line.direction()[1]) : abs(start_line.direction()[2]));
    double dyz = norm[1] == 0 ? (start_line.direction()[1] > 0 ? (cell[1]+1)*voxel_size - (p1[1] - YS) : p1[1] - YS - cell[1]*voxel_size) / abs(start_line.direction()[1])
        : (start_line.direction()[2] > 0 ? (cell[2] + 1) * voxel_size - (p1[2] - ZS) : p1[2] - ZS - cell[2] * voxel_size) / abs(start_line.direction()[2]);

    while (p1[0] >= XS && p1[0] <= XL && p1[1] >= YS && p1[1] <= YL && p1[2] >= ZS && p1[2] <= ZL) { //keep the equals!
        cleary(indexes, p1, raydir, norm, visited, thisPlane, true);
        if (dx < dyz) {
            p1 = start_line.pointAt(dx);
            dx += theta_x;
        }
        else {
            p1 = start_line.pointAt(dyz);
            dyz += theta_yz;
        }
    }

    signed long long i = 0;
    size_t pointys = 0;
//#pragma omp parallel for
        for (i = 0; i < remainingPoints.size(); ++i) {
            if (thisPlane.absDistance(pc[remainingPoints[i]].location) < threshold) {
//#pragma omp critical
                //indexes.push_back(remainingPoints[i]);
                Eigen::Vector3d proj = thisPlane.projection(pc[remainingPoints[i]].location);
                std::vector<size_t> loc = hashCell(pc[remainingPoints[i]].location);
                std::vector<size_t> poc = hashCell(proj);
                if (!visited[loc[0]][loc[1]][loc[2]])
                    std::cout << loc[0] << loc[1] << loc[2] << " location is " << pc[remainingPoints[i]].location[0] << "/" << pc[remainingPoints[i]].location[1] << "/" << pc[remainingPoints[i]].location[2] << ", projection in " << poc[0] << poc[1] << poc[2] << " location is " << proj[0] << "/" << proj[1] << "/" << proj[2] << std::endl;
            }
        }
    //comparisons += remainingPoints.size();
    //if (pointys != indexes.size())
        //std::cout <<"a" << std::endl;
    return indexes;

}


//Adds points in place within the threshold of a 2D ray in a 3D bounding box
void UniformPC::cleary(std::vector<size_t>& points, Eigen::Vector3d p, Eigen::Vector3d dir, Eigen::Vector3d norm,
    std::vector<std::vector<std::vector<bool>>> &visited, Eigen::Hyperplane<double, 3> thisPlane, bool unpadded) {

    std::vector<size_t> cell = hashCell(p);
    double theta_y = voxel_size / abs(dir[1]);
    double theta_z = voxel_size / abs(dir[2]);
    double dy = theta_y;
    double dz = theta_z;
    bool up = dir[1] > 0;
    bool forward = dir[2] > 0;

    // either y or z will start at the maximum/minimum, the other's starting distance will be changed by p's position
    if (abs(norm[1]) == 1) {
        if (up) dy = ((cell[1] + 1) * voxel_size - (p[1] - YS)) / abs(dir[1]);
        else dy = (p[1] - YS - (cell[1]) * voxel_size) / abs(dir[1]);
    }
    else {
        if (forward) dz = ((cell[2] + 1) * voxel_size - (p[2] - ZS)) / abs(dir[2]);
        else dz = (p[2] - ZS - (cell[2]) * voxel_size) / abs(dir[2]);
    }

    bool right = cell[0] < x_voxels - 1;
    bool left = cell[0] > 0;
    size_t y = 0;
    size_t z = 0;

    //check the first cell (might already be visited due to previous rays)
    if (!visited[cell[0]][cell[1]][cell[2]]) {
        if (!cells[cell[0]][cell[1]][cell[2]].empty()) //push thresholded points from cell onto points vector
            addPoints(cells[cell[0]][cell[1]][cell[2]], points, thisPlane);
        visited[cell[0]][cell[1]][cell[2]] = true;
        padX(cell[0], cell[1], cell[2], points, visited, thisPlane, left, right);
    }
    //std::cout << cell[0] << cell[1] << cell[2] << std::endl;
    
    // the >= 0 check isn't needed because a size_t will simply overflow to greater than the limit when negative anyway!
    while (cell[1] < y_voxels && cell[2] < z_voxels) {
        //note that due to the adjacent cell checking getting ahead of itself, actually checking the current cell is not needed.
        //std::cout << cell[0] << cell[1] << cell[2] << std::endl;
        //z negative
        if (cell[2] > 0) {
            z = cell[2] - 1;
            if (!visited[cell[0]][cell[1]][z]) {
                if (!cells[cell[0]][cell[1]][z].empty())
                    addPoints(cells[cell[0]][cell[1]][z], points, thisPlane);
                visited[cell[0]][cell[1]][z] = true;
            }
            padX(cell[0], cell[1], z, points, visited, thisPlane, left, right);
            ////visit the diagonal cells
            if (cell[1] > 0) {
                if (!visited[cell[0]][cell[1] - 1][z]) {
                    if (!cells[cell[0]][cell[1] - 1][z].empty())
                        addPoints(cells[cell[0]][cell[1] - 1][z], points, thisPlane);
                    visited[cell[0]][cell[1] - 1][z] = true;
                    padX(cell[0], cell[1] - 1, z, points, visited, thisPlane, left, right);
                }
            }
            if (cell[1] < y_voxels - 1) {
                if (!visited[cell[0]][cell[1] + 1][z]) {
                    if (!cells[cell[0]][cell[1] + 1][z].empty())
                        addPoints(cells[cell[0]][cell[1] + 1][z], points, thisPlane);
                    visited[cell[0]][cell[1] + 1][z] = true;
                    padX(cell[0], cell[1] + 1, z, points, visited, thisPlane, left, right);
                }
            }
        } //z positive
        if (cell[2] < z_voxels - 1) {
            z = cell[2] + 1;
            if (!visited[cell[0]][cell[1]][z]) {
                if (!cells[cell[0]][cell[1]][z].empty())
                    addPoints(cells[cell[0]][cell[1]][z], points, thisPlane);
                visited[cell[0]][cell[1]][z] = true;
            }
            padX(cell[0], cell[1], z, points, visited, thisPlane, left, right);
            ////visit the diagonal cells
            if (cell[1] > 0) {
                if (visited[cell[0]][cell[1] - 1][z]) {
                    if (!cells[cell[0]][cell[1] - 1][z].empty())
                        addPoints(cells[cell[0]][cell[1] - 1][z], points, thisPlane);
                    visited[cell[0]][cell[1] - 1][z] = true;
                    padX(cell[0], cell[1] - 1, z, points, visited, thisPlane, left, right);
                }
            }
            if (cell[1] < y_voxels - 1) {
                if (!visited[cell[0]][cell[1] + 1][z]) {
                    if (!cells[cell[0]][cell[1] + 1][z].empty())
                        addPoints(cells[cell[0]][cell[1] + 1][z], points, thisPlane);
                    visited[cell[0]][cell[1] + 1][z] = true;
                    padX(cell[0], cell[1] + 1, z, points, visited, thisPlane, left, right);
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
            }
            padX(cell[0], y, cell[2], points, visited, thisPlane, left, right);
        } //y positive
        if (cell[1] < y_voxels - 1) {
            y = cell[1] + 1;
            if (!visited[cell[0]][y][cell[2]]) {
                if (!cells[cell[0]][y][cell[2]].empty())
                    addPoints(cells[cell[0]][y][cell[2]], points, thisPlane);
                visited[cell[0]][y][cell[2]] = true;
            }
            padX(cell[0], y, cell[2], points, visited, thisPlane, left, right);
        }

        //work out next cell
        if (dy < dz) {
            dy += theta_y;  // going to the y adjacent cell
            (up) ? cell[1] += 1 : cell[1] -= 1;
        }
        else {
            dz += theta_z;  // going to the z adjacent cell
            (forward) ? cell[2] += 1 : cell[2] -= 1;
        }
    }

    return;
}


void UniformPC::padX(size_t x, size_t y, size_t z, std::vector<size_t>& points, std::vector<std::vector<std::vector<bool>>>& visited, 
    Eigen::Hyperplane<double, 3> thisPlane, bool left, bool right) {
    if (left) {
        if (!visited[x - 1][y][z]) {
            if (!cells[x - 1][y][z].empty())
                addPoints(cells[x - 1][y][z], points, thisPlane);
            visited[x - 1][y][z] = true;
        }
    }
    if (right) {
        if (!visited[x + 1][y][z]) {
            if (!cells[x + 1][y][z].empty())
                addPoints(cells[x + 1][y][z], points, thisPlane);
            visited[x + 1][y][z] = true;
        }
    }
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
    //std::cout << comparisons << std::endl;
}
