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

    x_voxels = (size_t)(XL - XS) / voxel_size;
    y_voxels = (size_t)(YL - YS) / voxel_size;
    z_voxels = (size_t)(ZL - ZS) / voxel_size;
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

    // for each point, hash their index in the vector into a cell
    signed long long i = 0;
#pragma omp parallel for
    for (i = 0; i < pc.size(); ++i) {
        std::vector<size_t> cell = hashCell(pc[i].location);
#pragma omp critical
        cells[cell[0]][cell[1]][cell[2]].push_back(i);
    }

    //sort the indexes of the populated cells
    for (size_t a = 0; a < x_voxels; ++a) {
        for (size_t b = 0; b < y_voxels; ++b) {
            for (size_t c = 0; c < z_voxels; ++c) {
                if (cells[a][b][c].size() > 1) std::sort(cells[a][b][c].begin(), cells[a][b][c].end());
            }
        }
    }
}


std::vector<size_t> UniformPC::hashCell(const Eigen::Vector3d& p) {
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


std::vector<size_t> UniformPC::planePoints(const Eigen::Hyperplane<double, 3> &thisPlane) {
    std::vector<size_t> indexes;
    int d1, d2, d3;

    if (abs(thisPlane.coeffs()[0]) > abs(thisPlane.coeffs()[1]) && abs(thisPlane.coeffs()[0]) > abs(thisPlane.coeffs()[2])) {
        d3 = 0;
        d2 = 2;
        d1 = 1;
    }
    else if (abs(thisPlane.coeffs()[1]) > abs(thisPlane.coeffs()[0]) && abs(thisPlane.coeffs()[1]) > abs(thisPlane.coeffs()[2])) {
        d3 = 1;
        d2 = 2;
        d1 = 0;
    }
    else {
        d3 = 2;
        d2 = 1;
        d1 = 0;
    }

    Eigen::Vector3d minima = { XS, YS, ZS };
    signed long long i = 0;
    //std::cout << "new planey waney" << std::endl;
#pragma omp parallel for
    for (i = 0; i < limits[d1]; i++) {

        Eigen::Vector3d point = minima;
        size_t thread_comparisons = 0;
        double d1_min = (-thisPlane.coeffs()[d1] * point[d1] - thisPlane.coeffs()[3]) / thisPlane.coeffs()[d3];
        double d1_max = (-thisPlane.coeffs()[d1] * (point[d1] + voxel_size) - thisPlane.coeffs()[3]) / thisPlane.coeffs()[d3];
        std::vector<size_t> cell = { 0, 0, 0 };
        cell[d1] = i;

        //d3 location of the plane in 4 corners of the voxel
        double t1 = d1_min - (thisPlane.coeffs()[d2] * point[d2] / thisPlane.coeffs()[d3]); //bottom left corner
        double t2 = d1_max - (thisPlane.coeffs()[d2] * point[d2] / thisPlane.coeffs()[d3]); //+ d1
        double t3 = t1 - voxel_size * thisPlane.coeffs()[d2] / thisPlane.coeffs()[d3]; // + d2
        double t4 = t2 - voxel_size * thisPlane.coeffs()[d2] / thisPlane.coeffs()[d3]; // + d1 and d2

        for (size_t j = 0; j < limits[d2]; j++) {

            cell[d2] = j;
            point[d1] = point[d1] + i * voxel_size;

            double lower_lim = std::min(t1, t2);
            lower_lim = std::min(lower_lim, t3);
            lower_lim = std::min(lower_lim, t4);
            double upper_lim = std::max(t1, t2);
            upper_lim = std::max(upper_lim, t3);
            upper_lim = std::max(upper_lim, t4);

            lower_lim -= threshold;
            upper_lim += threshold;
            signed long long lower = (lower_lim - minima[d3]) / voxel_size - 1;
            signed long long upper = (upper_lim - minima[d3]) / voxel_size + 1;

            cell[d3] = std::max((long long) 0, lower);

            while (cell[d3] < limits[d3] && (signed long long) cell[d3] <= upper) {
                //std::cout << cell[0] << cell[1] << cell[2] << std::endl;
                size_t olde = indexes.size();
                //initialize iterator to progress during the for loop
                std::vector<size_t>::iterator remain = remainingPoints.begin();
                for (size_t index = 0; index < cells[cell[0]][cell[1]][cell[2]].size(); index++) {
                    // find point in remainingPoints vector
                    remain = std::lower_bound(remain, remainingPoints.end(), cells[cell[0]][cell[1]][cell[2]][index]);
                    if (thisPlane.absDistance(pc[ cells[cell[0]][cell[1]][cell[2]][index] ].location) < threshold
                        && remain != remainingPoints.end())
#pragma omp critical
                        indexes.push_back(cells[cell[0]][cell[1]][cell[2]][index]);
                }
                //if (indexes.size() > olde) std::cout << cell[0] << cell[1] << cell[2] << std::endl;
                thread_comparisons += cells[cell[0]][cell[1]][cell[2]].size();
                cell[d3] += 1;

            }
            point[d2] += voxel_size;
            t1 = t3;
            t2 = t4;
            t3 = t1 - voxel_size * thisPlane.coeffs()[d2] / thisPlane.coeffs()[d3];
            t4 = t2 - voxel_size * thisPlane.coeffs()[d2] / thisPlane.coeffs()[d3];

        }
#pragma omp critical
        comparisons += thread_comparisons;

    }

    //std::cout << "testy westy" << std::endl;
    ////testy
    //std::vector<size_t> thisPoints;
    //for (size_t a = 0; a < x_voxels; a++) {
    //    for (size_t b = 0; b < y_voxels; b++) {
    //        for (size_t c = 0; c < z_voxels; c++) {
    //            size_t old = thisPoints.size();
    //            for (size_t ugh = 0; ugh < cells[a][b][c].size(); ugh++) {
    //                if (thisPlane.absDistance(pc[ cells[a][b][c][ugh] ].location) < threshold && std::binary_search(remainingPoints.begin(), remainingPoints.end(), cells[a][b][c][ugh]))
    //                    thisPoints.push_back(cells[a][b][c][ugh]);
    //            }
    //            if (thisPoints.size() > old) std::cout << a << b << c << ", " << thisPoints.size() - old << std::endl;
    //        }
    //    }
    //}

    //if (thisPoints.size() != indexes.size())
    //    size_t wah = std::max(thisPoints.size(), indexes.size());
    return indexes;
}


//std::vector<size_t> UniformPC::planePoints(Eigen::Hyperplane<double, 3> thisPlane) {
//
//    //get first edge intersection
//    int edge = 0;
//    Eigen::Vector3d p1 = edges[edge].intersectionPoint(thisPlane);
//    // while the found point is out of bounds find another
//    while (p1[0] < XS || p1[0] > XL || p1[1] < YS || p1[1] > YL || p1[2] < ZS || p1[2] > ZL)
//        p1 = edges[++edge].intersectionPoint(thisPlane);
//    //get second edge intersection
//    int edge2 = edge + 1;
//    Eigen::Vector3d p2 = edges[edge2].intersectionPoint(thisPlane);
//    //line to cast rays from
//    Eigen::ParametrizedLine<double, 3> start_line = Eigen::ParametrizedLine<double, 3>::Through(p1, p2);
//    Eigen::Vector3d norm = edges[edge].direction().cross(start_line.direction());
//    while (norm[0] != 0 || p2[0] < XS || p2[0] > XL || p2[1] < YS || p2[1] > YL || p2[2] < ZS || p2[2] > ZL) { //point is out of bounds
//        p2 = edges[++edge2].intersectionPoint(thisPlane);
//        start_line = Eigen::ParametrizedLine<double, 3>::Through(p1, p2);
//        norm = edges[edge].direction().cross(start_line.direction());
//    }
//    norm.normalize();
//    //make sure norm is oriented into the bounding box
//    if (abs(norm[1]) == 1 && (norm[1] + p1[1] < YS || norm[1] + p1[1] > YL)) norm[1] = -1 * norm[1];
//    else if (norm[2] + p1[2] < ZS || norm[2] + p1[2] > ZL) norm[2] = -1 * norm[2];
//
//    // these rays should be fixed in the x direction, and vary in the y and z directions to fit the plane
//    Eigen::Vector3d raydir = { 1, 0, 0 };
//    if (raydir.isApprox(thisPlane.normal()))
//        raydir = norm;
//    else
//        raydir = thisPlane.normal().cross(raydir);
//    raydir.normalize();
//    // make sure cross product is oriented into the bounding box (need to reverse the whole vector)
//    if ((norm[1] == 1 && raydir[1] < 0) || norm[1] == -1 && raydir[1] > 0) raydir = -raydir;
//    else if ((norm[2] == 1 && raydir[2] < 0) || norm[2] == -1 && raydir[2] > 0) raydir = -raydir;
//    double raytest = thisPlane.absDistance(p1 + raydir);
//
//    // indexes of points on the plane to be returned
//    std::vector<size_t> indexes;
//    // 3D truth array of visited voxels
//    std::vector<std::vector<std::vector<bool>>> visited;
//    //idea - set all empty cells to visited?
//    visited = std::vector<std::vector<std::vector<bool>>>(x_voxels, std::vector<std::vector<bool>>(y_voxels, std::vector<bool>(z_voxels, false)));
//
//    // calculate exits?
//    std::vector<size_t> cell = hashCell(p1);
//    std::vector<size_t> cell2 = cell;
//    int yz = norm[1] == 0 ? 1 : 2;
//    int padding_dir = abs(start_line.direction()[0]) < abs(start_line.direction()[yz]) ? 0 : yz;
//    bool right = start_line.direction()[0] > 0;
//    bool up = start_line.direction()[yz] > 0;
//
//    //absolute distance travelled in the y/z direction
//    double pyz = p1[yz] - (yz == 1 ? YS : ZS);
//    double theta_x = voxel_size / abs(start_line.direction()[0]);
//    // the x axis distance to the next cell (will be voxel_size if the lines does not have an x offset)
//    double next_x = right ? (cell[0] + 1) * voxel_size - (p1[0] - XS) : p1[0] - XS - cell[0] * voxel_size;
//    double dx = next_x / abs(start_line.direction()[0]);
//    double theta_yz = voxel_size / abs(start_line.direction()[yz]);
//    // the y/z axis distance to the next cell (will be voxel_size if the lines does not have an y/z offset)
//    double next_yz = (up ? (cell[yz] + 1) * voxel_size - pyz : pyz - cell[yz] * voxel_size);
//    double dyz = next_yz / abs(start_line.direction()[yz]);
//    //now change the distance to the next cell to reflect the ray direction
//    bool rayup = raydir[yz] > 0;
//    if ((up && !rayup) || (!up && rayup))
//        next_yz = voxel_size - next_yz;
//
//    //pad other direction of first cell
//    if (padding_dir == 1) {
//        if (right && cell[0] > 0) cell2[0] -= 1;
//        else if (cell[0] < x_voxels - 1) cell2[0] += 1;
//    }
//    else {
//        if (up && cell[yz] > 0) cell2[yz] -= 1;
//        else if (cell[yz] < limits[yz] - 1) cell2[yz] += 1;
//    }
//    cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//
//    std::cout << thisPlane.normal()[0] << thisPlane.normal()[1] << thisPlane.normal()[2] << ", d = " << thisPlane.coeffs()[3] << std::endl;
//    while (cell[0] < x_voxels && cell[yz] < limits[yz]) {
//        //pass the indexes to append to in place, the y or z intercept, the starting cell, the ray direction, the plane norm, 
//        //the visited array to edit in place and the plane to measure distance from
//        cleary(indexes, next_yz, cell, raydir, norm, visited, thisPlane, p1);
//        if (dx < dyz) {
//            //pad adjacent voxels
//            if (padding_dir != 0) {
//                //pad both yz cells
//                if (cell[yz] > 0) {
//                    cell2[yz] -= 1;
//                    cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//                }
//                if (cell[yz] < limits[yz] - 1) {
//                    cell2[yz] += 1;
//                    cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//                }
//            }
//            else {
//                //only pad x cell behind
//                if (right && cell[0] > 0) {
//                    cell2[0] -= 1;
//                    cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//                }
//                else if (cell[0] < x_voxels - 1){
//                    cell2[0] += 1;
//                    cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//                }
//            }
//            (right) ? cell[0] += 1 : cell[0] -= 1;
//            dx += theta_x;
//            pyz = dx * abs(start_line.direction()[yz]);
//            next_yz = (rayup ? (cell[yz] + 1) * voxel_size - pyz : pyz - cell[yz] * voxel_size);
//            p1 = start_line.pointAt(dx);
//        }
//        else {
//            //pad adjacent voxels
//            if (padding_dir != yz) {
//                //pad both x cells
//                if (cell[0] > 0) {
//                    cell2[0] -= 1;
//                    cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//                }
//                if (cell[0] < limits[yz] - 1) {
//                    cell2[0] += 1;
//                    cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//                }
//            }
//            else {
//                //only pad yz cell behind
//                if (up && cell[yz] > 0) {
//                    cell2[yz] -= 1;
//                    cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//                }
//                else if (cell[yz] < limits[yz] - 1) {
//                    cell2[yz] += 1;
//                    cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//                }
//            }
//            (up) ? cell[yz] += 1 : cell[yz] -= 1;
//            dyz += theta_yz;
//            pyz += voxel_size;
//            next_yz = voxel_size;
//            p1 = start_line.pointAt(dyz);
//        }
//        cell2 = cell;
//    }
//    //pad normally from outside the box if padding direction limit has been reached
//    if (padding_dir == 0 && cell[yz] < limits[yz]) {
//        (right) ? cell[0] -= 1 : cell[0] += 1;
//        cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//    }
//    else if (cell[0] < x_voxels) {
//        (up) ? cell[yz] -= 1 : cell[yz] += 1;
//        cleary(indexes, next_yz, cell2, raydir, norm, visited, thisPlane, p1);
//    }
//
//    signed long long i = 0;
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
//                    //std::cout << "projection " << poc[0] << poc[1] << poc[2] << ", covered: " << visited[poc[0]][poc[1]][poc[2]] << ", location is " << proj[0] << "/" << proj[1] << "/" << proj[2] << std::endl;
//                }
//            }
//        }
//    //comparisons += remainingPoints.size();
//    return indexes;
//
//}


//Adds points in place within the threshold of a 2D ray in a 3D bounding box
void UniformPC::cleary(std::vector<size_t>& points, double next_yz, std::vector<size_t> cell, const Eigen::Vector3d& dir, const Eigen::Vector3d& norm,
    std::vector<std::vector<std::vector<bool>>>& visited, Eigen::Hyperplane<double, 3> thisPlane, Eigen::Vector3d fpoint) {

    int padding_dir = abs(dir[1]) < abs(dir[2]) ? 1 : 2;
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
    size_t y = 0;
    size_t z = 0;

    std::cout << "start of ray:" << std::endl;
    
    // the >= 0 check isn't needed because a size_t will simply overflow to greater than the limit when negative anyway!
    while (cell[1] < y_voxels && cell[2] < z_voxels) {
        //check the current cell
        checkcell(cell[0], cell[1], cell[2], points, visited, thisPlane);
        //work out next cell
        if (dy < dz) {
            dy += theta_y;  // going to the y adjacent cell
            //pad adjacent voxels
            if (padding_dir == 2) {
                //pad both z cells
                if (cell[2] > 0) checkcell(cell[0], cell[1], cell[2] - 1, points, visited, thisPlane);
                if (cell[2] < z_voxels - 1) checkcell(cell[0], cell[1], cell[2] + 1, points, visited, thisPlane);
            }
            else {
                //only pad y cell behind
                if (up && cell[1] > 0) checkcell(cell[0], cell[1] - 1, cell[2], points, visited, thisPlane);
                else if (cell[1] < y_voxels - 1) checkcell(cell[0], cell[1] + 1, cell[2], points, visited, thisPlane);
            }
            (up) ? cell[1] += 1 : cell[1] -= 1;
            std::cout << "distance from plane: " << thisPlane.absDistance(fpoint + dir * dy) << std::endl;
        }
        else {
            dz += theta_z;  // going to the z adjacent cell
            //pad adjacent voxels
            if (padding_dir == 1) {
                //pad both y cells
                if (cell[1] > 0) checkcell(cell[0], cell[1] - 1, cell[2], points, visited, thisPlane);
                if (cell[1] < y_voxels - 1) checkcell(cell[0], cell[1] + 1, cell[2], points, visited, thisPlane);
            }
            else {
                //only pad z cell behind
                if (forward && cell[2] > 0) checkcell(cell[0], cell[1], cell[2] - 1, points, visited, thisPlane);
                else if (cell[2] < z_voxels - 1) checkcell(cell[0], cell[1], cell[2] + 1, points, visited, thisPlane);
            }
            (forward) ? cell[2] += 1 : cell[2] -= 1;
            std::cout << "distance from plane: " << thisPlane.absDistance(fpoint + dir * dz) << std::endl;
        }
        std::cout << cell[0] << cell[1] << cell[2] << std::endl;
    }

    return;
}


void UniformPC::checkcell(size_t x, size_t y, size_t z, std::vector<size_t>& points, std::vector<std::vector<std::vector<bool>>>& visited,
    const Eigen::Hyperplane<double, 3>& thisPlane) {
    if (!visited[x][y][z]) {
        if (!cells[x][y][z].empty()) //push thresholded points from cell onto points vector
            addPoints(cells[x][y][z], points, thisPlane);
        visited[x][y][z] = true;
    }
}


void UniformPC::addPoints(std::vector<size_t> indexes, std::vector<size_t>& thisPoints,
    const Eigen::Hyperplane<double, 3>& thisPlane) {
    
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

//void UniformPC::removePoints(std::vector<size_t>& planePoints, int plane) {
//    signed long long i = 0;
//#pragma omp parallel for
//    for (i = 0; i < x_voxels; ++i) {
//        for (size_t j = 0; j < y_voxels; ++j) {
//            for (size_t k = 0; k < z_voxels; ++k) {
//                
//            }
//        }
//        std::vector<size_t> cell = hashCell(pc[planePoints[i]].location);
//        std::vector<size_t>::iterator spot = std::lower_bound(cells[cell[0]][cell[1]][cell[2]].begin(), cells[cell[0]][cell[1]][cell[2]].end())
//    }
//}
