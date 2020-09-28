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

    // check how best to order these for access time
    cells = std::vector<std::vector<std::vector<std::vector<size_t>>>>(x_voxels, std::vector<std::vector<std::vector<size_t>>>(y_voxels, std::vector<std::vector<size_t>>(z_voxels)));

    // for each point, hash their index in the vector into a cell
    signed long long i = 0;
#pragma omp parallel for num_threads(num_threads)
    for (i = 0; i < pc.size(); ++i) {
        std::vector<size_t> cell = hashCell(pc[i].location);
#pragma omp critical
        cells[cell[0]][cell[1]][cell[2]].push_back(i);
    }

    //sort the indexes of the populated cells
    signed long long a = 0;
#pragma omp parallel for num_threads(num_threads)
    for (a = 0; a < x_voxels; ++a) {
        for (size_t b = 0; b < y_voxels; ++b) {
            for (size_t c = 0; c < z_voxels; ++c) {
                if (cells[a][b][c].size() > 1) std::sort(cells[a][b][c].begin(), cells[a][b][c].end());
            }
        }
    }
    remainingCells = cells;
    std::cout << "Done uniform space subdivision of " << x_voxels << " by " << y_voxels << " by " << z_voxels << " voxels..." << std::endl;
}


std::vector<size_t> UniformPC::hashCell(const Eigen::Vector3d& p) {
    size_t x = (p[0] - XS) / voxel_size;
    size_t y = (p[1] - YS) / voxel_size;
    size_t z = (p[2] - ZS) / voxel_size;
    if (x > x_voxels - 1) x = x_voxels - 1;
    if (y > y_voxels - 1) y = y_voxels - 1;
    if (z > z_voxels - 1) z = z_voxels - 1;

    return {x, y, z};
}


//std::vector<size_t> UniformPC::planePoints(const Eigen::Hyperplane<double, 3> &thisPlane) {
std::vector<bool> UniformPC::planePoints(const Eigen::Hyperplane<double, 3>& thisPlane) {
    //std::vector<size_t> indexes;
    std::vector<bool> indexes = std::vector<bool>(size, false);
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
#pragma omp parallel for num_threads(num_threads)
    for (i = 0; i < limits[d1]; i++) {

        Eigen::Vector3d point = minima;
        //initial bottom left voxel corner of the cell
        point[d1] += i * voxel_size;
        size_t thread_comparisons = 0;
        double d1_min = (-thisPlane.coeffs()[d1] * point[d1] - thisPlane.coeffs()[3]) / thisPlane.coeffs()[d3];
        double d1_max = (-thisPlane.coeffs()[d1] * (point[d1] + voxel_size) - thisPlane.coeffs()[3]) / thisPlane.coeffs()[d3];
        std::vector<size_t> cell = { 0, 0, 0 };
        cell[d1] = i;

        //d3 location of the plane in 4 corners of the voxel
        double t1 = d1_min - (thisPlane.coeffs()[d2] * point[d2] / thisPlane.coeffs()[d3]); //bottom left corner
        double t2 = d1_max - (thisPlane.coeffs()[d2] * point[d2] / thisPlane.coeffs()[d3]); //+ d1
        double shift = voxel_size * thisPlane.coeffs()[d2] / thisPlane.coeffs()[d3];
        double t3 = t1 - shift; // + d2
        double t4 = t2 - shift; // + d1 and d2

        for (size_t j = 0; j < limits[d2]; j++) {

            cell[d2] = j;

            double lower_lim = std::min(t1, t2);
            lower_lim = std::min(lower_lim, t3);
            lower_lim = std::min(lower_lim, t4);
            double upper_lim = std::max(t1, t2);
            upper_lim = std::max(upper_lim, t3);
            upper_lim = std::max(upper_lim, t4);

            lower_lim -= threshold;
            upper_lim += threshold;
            signed long long lower = (lower_lim - minima[d3]) / voxel_size;
            signed long long upper = (upper_lim - minima[d3]) / voxel_size;

            cell[d3] = std::max((long long) 0, lower);

            while (cell[d3] < limits[d3] && (signed long long) cell[d3] <= upper) {
                for (size_t index = 0; index < remainingCells[cell[0]][cell[1]][cell[2]].size(); index++) {
                    if (thisPlane.absDistance(pc[remainingCells[cell[0]][cell[1]][cell[2]][index]].location) < threshold)
                        //#pragma omp critical
                                                //indexes.push_back(remainingCells[cell[0]][cell[1]][cell[2]][index]);
                        indexes[remainingCells[cell[0]][cell[1]][cell[2]][index]] = true;
                }
                thread_comparisons += remainingCells[cell[0]][cell[1]][cell[2]].size();
                cell[d3] += 1;

            }
            t1 = t3;
            t2 = t4;
            t3 -= shift;
            t4 -= shift;

        }
#pragma omp critical
        comparisons += thread_comparisons;

    }
    return indexes;
}


// Sets points plane ID and removes then from the lists of remaining points
void UniformPC::removePoints(const std::vector<size_t>& planePoints, int plane) {
    //remove from the remainingPoints vector
    std::vector<size_t> diff;
    for (size_t i = 0; i < size; ++i) {
        if (planePoints[i]) {
            setPointPlane(i, plane);
        }
        else if (std::binary_search(remainingPoints.begin(), remainingPoints.end(), i)) {
            diff.push_back(i);
        }
    }
    remainingPoints = diff;

    //remove from the remaining cells
    signed long long a = 0;
#pragma omp parallel for num_threads(num_threads)
    for (a = 0; a < x_voxels; ++a) {
        for (size_t b = 0; b < y_voxels; ++b) {
            for (size_t c = 0; c < z_voxels; ++c) {
                // done using erase because i don't want to change the 4D array
                size_t d = 0;
                while (d < remainingCells[a][b][c].size()) {
                    if (planePoints[remainingCells[a][b][c][d]]) { //if not found advance
                        ++d;
                    }
                    else { //erase from remaining cells and keep index the same to move forward
                        remainingCells[a][b][c].erase(remainingCells[a][b][c].begin() + d);
                    }
                }
            }
        }
    }
}


// Helper method for testing when rerunning RANSAC without remaking the point cloud is desired
void UniformPC::resetRemaining() {
    remainingPoints.resize(size);
    comparisons = 0;
    signed long long i = 0;
#pragma omp parallel for num_threads(num_threads)
    for (i = 0; i < size; ++i) {
        remainingPoints[i] = i;
    }
    remainingCells = cells;
}

