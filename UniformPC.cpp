#include "UniformPC.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>


UniformPC::UniformPC(PointCloud const&p, int voxel_scale) : PointCloud(p) {

    XS = floor(XS);
    XL = ceil(XL);
    YS = floor(YS);
    YL = ceil(YL);
    ZS = floor(ZS);
    ZL = ceil(ZL);

    // assign the given volumes of voxels to the model dimensions
    if (voxel_scale < threshold)
        std::cout << "Warning, voxel size given is smaller than threshold, some points will be missed" << std::endl;
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
    limits = { x_voxels, y_voxels, z_voxels };

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
    signed long long i = 0;
//#pragma omp parallel for
    for (i = 0; i < pc.size(); ++i) {
        cell = hashCell(pc[i].location);
//#pragma omp critical
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

    //if (p[0] - XS > (x + 1) * voxel_size)
    //    std::cout << p[0] << std::endl;
    //if (p[0] - XS < x * voxel_size)
    //    std::cout << p[0] << std::endl;

    //if (p[1] - YS > (y + 1) * voxel_size)
    //    std::cout << p[1] << std::endl;
    //if (p[1] - YS < y * voxel_size)
    //    std::cout << p[1] << std::endl;

    //if (p[2] - ZS > (z + 1) * voxel_size)
    //    std::cout << p[2] << std::endl;
    //if (p[2] - ZS < z * voxel_size)
    //    std::cout << p[2] << std::endl;

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
    Eigen::Vector3d raydir = { 1, 0, 0 };
    if (raydir.isApprox(thisPlane.normal()))
        raydir = norm;
    else
        raydir = thisPlane.normal().cross(raydir);
    raydir.normalize();
    // make sure cross product is oriented into the bounding box (need to reverse the whole vector)
    if ((norm[1] == 1 && raydir[1] < 0) || norm[1] == -1 && raydir[1] > 0) raydir = -raydir;
    else if ((norm[2] == 1 && raydir[2] < 0) || norm[2] == -1 && raydir[2] > 0) raydir = -raydir;
    //double raytest = thisPlane.absDistance(p1 + raydir);

    // indexes of points on the plane to be returned
    std::vector<size_t> indexes;
    // 3D truth array of visited voxels
    std::vector<std::vector<std::vector<bool>>> visited;
    //idea - set all empty cells to visited?
    visited = std::vector<std::vector<std::vector<bool>>>(x_voxels, std::vector<std::vector<bool>>(y_voxels, std::vector<bool>(z_voxels, false)));

    // calculate exits?
    std::vector<size_t> cell = hashCell(p1);
    int yz = norm[1] == 0 ? 1 : 2;
    bool right = start_line.direction()[0] > 0;
    bool up = start_line.direction()[yz] > 0;

    //absolute distance travelled in the y/z direction
    double pyz = p1[yz] - (yz == 1 ? YS : ZS);
    double theta_x = voxel_size / abs(start_line.direction()[0]);
    // the x axis distance to the next cell (will be voxel_size if the lines does not have an x offset)
    double next_x = right ? (cell[0] + 1) * voxel_size - (p1[0] - XS) : p1[0] - XS - cell[0] * voxel_size;
    double dx = next_x / abs(start_line.direction()[0]);
    double theta_yz = voxel_size / abs(start_line.direction()[yz]);
    // the y/z axis distance to the next cell (will be voxel_size if the lines does not have an y/z offset)
    double next_yz = (up ? (cell[yz] + 1) * voxel_size - pyz : pyz - cell[yz] * voxel_size);
    double dyz = next_yz / abs(start_line.direction()[yz]);
    //now change the distance to the next cell to reflect the ray direction
    bool rayup = raydir[yz] > 0;
    if ((up && !rayup) || (!up && rayup))
        next_yz = voxel_size - next_yz;

    //std::cout << thisPlane.normal()[0] << thisPlane.normal()[1] << thisPlane.normal()[2] << ", d = " << thisPlane.coeffs()[3] << std::endl;
    while (cell[0] < x_voxels && cell[yz] < limits[yz]) {
        //pass the indexes to append to in place, the y or z intercept, the starting cell, the ray direction, the plane norm, 
        //the visited array to edit in place and the plane to measure distance from
        cleary(indexes, next_yz, cell, raydir, norm, visited, thisPlane);//, p1);
        if (dx < dyz) {
            (right) ? cell[0] += 1 : cell[0] -= 1;
            dx += theta_x;
            pyz = dx * abs(start_line.direction()[yz]);
            next_yz = (rayup ? (cell[yz] + 1) * voxel_size - pyz : pyz - cell[yz] * voxel_size);
            //p1 = start_line.pointAt(dx);
        }
        else {
            (up) ? cell[yz] += 1 : cell[yz] -= 1;
            dyz += theta_yz;
            pyz += voxel_size;
            next_yz = voxel_size;
            //p1 = start_line.pointAt(dyz);
        }
    }

//    signed long long i = 0;
//    size_t pointys = 0;
////#pragma omp parallel for
//        for (i = 0; i < remainingPoints.size(); ++i) {
//            if (thisPlane.absDistance(pc[remainingPoints[i]].location) < threshold) {
////#pragma omp critical
//                //indexes.push_back(remainingPoints[i]);
//                Eigen::Vector3d proj = thisPlane.projection(pc[remainingPoints[i]].location);
//                std::vector<size_t> loc = hashCell(pc[remainingPoints[i]].location);
//                std::vector<size_t> poc = hashCell(proj);
//                if (!visited[loc[0]][loc[1]][loc[2]]) {
//                    std::cout << loc[0] << loc[1] << loc[2] << " distance is " << thisPlane.absDistance(pc[remainingPoints[i]].location) << std::endl;
//                    std::cout << "projection " << poc[0] << poc[1] << poc[2] << ", covered: " << visited[poc[0]][poc[1]][poc[2]] << ", location is " << proj[0] << "/" << proj[1] << "/" << proj[2] << std::endl;
//                }
//            }
//        }
    //comparisons += remainingPoints.size();
    //if (pointys != indexes.size())
        //std::cout <<"a" << std::endl;
    return indexes;

}


//Adds points in place within the threshold of a 2D ray in a 3D bounding box
void UniformPC::cleary(std::vector<size_t>& points, double next_yz, std::vector<size_t> cell, Eigen::Vector3d dir, Eigen::Vector3d norm,
    std::vector<std::vector<std::vector<bool>>> &visited, Eigen::Hyperplane<double, 3> thisPlane) {//, Eigen::Vector3d fuckingpoint) {

    double theta_y, theta_z, dy, dz;
    theta_y = voxel_size / abs(dir[1]);
    theta_z = voxel_size / abs(dir[2]);
    if (norm[1] != 1) { // if y changed in the starting line
        dy = next_yz / abs(dir[1]);
        dz = theta_z;
    }
    else { // if z changed in the starting line
        dy = theta_y;
        dz = next_yz  / abs(dir[2]);
    }
    bool up = dir[1] > 0;
    bool forward = dir[2] > 0;
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
    //std::cout << "start of ray:" << std::endl;
    
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
            //std::cout << "distance from plane: " << thisPlane.absDistance(fuckingpoint + dir * dy) << std::endl;
        }
        else {
            dz += theta_z;  // going to the z adjacent cell
            (forward) ? cell[2] += 1 : cell[2] -= 1;
            //std::cout << "distance from plane: " << thisPlane.absDistance(fuckingpoint + dir * dz) << std::endl;
        }
        //std::cout << cell[0] << cell[1] << cell[2] << std::endl;
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
