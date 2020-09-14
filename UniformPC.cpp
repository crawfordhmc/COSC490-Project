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
    if (x > x_voxels)
        x -= 1;
    if (y > y_voxels)
        y -= 1;
    if (z > z_voxels)
        z -= 1;
    //checking the assignment is correct
    if (p[0] > (x + 1) * voxel_size + XS) std::cout << p[0] << x << voxel_size << XS << std::endl;
    if (p[1] > (y + 1) * voxel_size + YS) std::cout << p[1] << y << voxel_size << YS << std::endl;
    if (p[2] > (z + 1) * voxel_size + ZS) std::cout << p[2] << z << voxel_size << ZS << std::endl;
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
    if (raydir[1] == 1) raydir[2] = -thisPlane.coeffs()[1] / thisPlane.coeffs()[2];
    else raydir[1] = -thisPlane.coeffs()[2] / thisPlane.coeffs()[1];

    // indexes of points on the plane to be returned
    std::vector<size_t> indexes;
    // 3D truth array of visited voxels
    std::vector<std::vector<std::vector<bool>>> visited;
    visited = std::vector<std::vector<std::vector<bool>>>(x_voxels, std::vector<std::vector<bool>>(y_voxels, std::vector<bool>(z_voxels, false)));
    
    do {
        size_t x_vox = hashCell(p1)[0];
        cleary(indexes, Eigen::ParametrizedLine<double, 3>(p1, raydir), visited, thisPlane);
        p1 += step;
    } while (start_line.projection(p1).norm() < start_line.projection(p2).norm());
    
    return indexes;

}


//// Returns a vector of points within the threshold to the given hyperplane
//std::vector<size_t> UniformPC::planePoints(Eigen::Hyperplane<double, 3> thisPlane, unsigned int trial, int plane) {
//    //theta value for floating-point comparisons NEEDED?
//    double t = 1 + 0.01 * voxel_size;
//    std::vector<Eigen::ParametrizedLine<double, 3>> edges = {
//        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZS }, { 1, 0, 0 }),
//        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZS }, { 0, 1, 0 }),
//        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZS }, { 0, 0, 1 }),
//
//        Eigen::ParametrizedLine<double, 3>({ XL, YS, ZS }, { 0, 1, 0 }),
//        Eigen::ParametrizedLine<double, 3>({ XL, YS, ZS }, { 0, 0, 1 }),
//
//        Eigen::ParametrizedLine<double, 3>({ XS, YL, ZS }, { 1, 0, 0 }),
//        Eigen::ParametrizedLine<double, 3>({ XS, YL, ZS }, { 0, 0, 1 }),
//
//        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZL }, { 1, 0, 0 }),
//        Eigen::ParametrizedLine<double, 3>({ XS, YS, ZL }, { 0, 1, 0 }),
//
//        Eigen::ParametrizedLine<double, 3>({ XS, YL, ZL}, { 1, 0, 0 }),
//        Eigen::ParametrizedLine<double, 3>({ XL, YS, ZL}, { 0, 1, 0 }),
//        Eigen::ParametrizedLine<double, 3>({ XL, YL, ZS}, { 0, 0, 1 })
//    };
//    // indexes of points on the plane to be returned
//    std::vector<size_t> indexes;
//    // 3D truth array of visited voxels
//    std::vector<std::vector<std::vector<bool>>> visited;
//    visited = std::vector<std::vector<std::vector<bool>>>(x_voxels, std::vector<std::vector<bool>>(y_voxels, std::vector<bool>(z_voxels, false)));
//    // Normalized plane coefficients
//    Eigen::Vector3d coeffs = { thisPlane.coeffs()[0], thisPlane.coeffs()[1], thisPlane.coeffs()[2] };
//    coeffs.normalize();
//
//    //get first edge intersection
//    int edge = 0;
//    Eigen::Vector3d p1 = edges[edge].intersectionPoint(thisPlane);
//    while (p1[0] < XS * t || p1[0] > XL * t || p1[1] < YS * t || p1[1] > YL * t || p1[2] < ZS * t || p1[2] > ZL * t) p1 = edges[++edge].intersectionPoint(thisPlane);
//    //get second edge intersection
//    int edge2 = edge + 1;
//    Eigen::Vector3d p2 = edges[edge2].intersectionPoint(thisPlane);
//    while (p2[0] < XS * t || p2[0] > XL * t || p2[1] < YS * t || p2[1] > YL * t || p2[2] < ZS * t || p2[2] > ZL * t) {  //end loop once a intersection within the bounding box is found
//        edge2 += 1;
//        //make sure the lines are on the same plane
//        if (edges[edge].origin()[0] != edges[edge2].origin()[0] && edges[edge].origin()[1] != edges[edge2].origin()[1] && edges[edge].origin()[2] != edges[edge2].origin()[2]) continue;
//        p2 = edges[edge2].intersectionPoint(thisPlane);
//    }
//
//    //line to cast rays from
//    Eigen::ParametrizedLine<double, 3> start_line = Eigen::ParametrizedLine<double, 3>::Through(p1, p2);
//    // direction of rays to be cast into the bounding box, following plane
//    Eigen::Vector3d dir = { 1, 1, 1 };
//    // if plane intersects on the far corner
//    if (edge > 8) dir = dir * -1;
//    dir = dir - edges[edge].direction();
//    if (edges[edge].direction() == edges[edge2].direction()) //if lines are parallel
//        dir = dir - (edges[edge2].origin() - edges[edge].origin()).normalized();
//    else dir = dir - edges[edge2].direction();
//    // if plane intersects on one of the far sides
//    if (edge > 2) dir = dir * -1;
//
//    Eigen::Vector3d step = start_line.direction() * voxel_size;
//
//    // set ray direction to be moving as opposed to fixed for 2D clearys as the slower-changing one of the start line
//    if (dir.dot(thisPlane.normal()) < 0.01 * voxel_size) { //if ray direction is not already on plane (dot product with plane normal close to zero)
//        int moving_axis;
//        if (abs(dir[0]) == 1) moving_axis = (abs(start_line.direction()[1]) < abs(start_line.direction()[2])) ? 1 : 2;
//        else if (abs(dir[1]) == 1) moving_axis = (abs(start_line.direction()[0]) < abs(start_line.direction()[2])) ? 0 : 2;
//        else moving_axis = (abs(start_line.direction()[0]) < abs(start_line.direction()[1])) ? 0 : 1;
//        dir[moving_axis] += coeffs[moving_axis]; // add normalized coefficient in the unfixed direction to fit ray to plane
//    }
//    //if plane is completely flat, just do all the cells in that row and call it a day
//    else {
//        std::vector<size_t> cell = hashCell(p1);
//        if (dir[0] == 1 && start_line.direction()[0] == 1) { //variance in y and z planes
//            bool up = cell[0] < x_voxels - 1;
//            bool down = cell[0] > 0;
//            for (size_t i = 0; i < y_voxels; i++) {
//                for (size_t j = 0; j < z_voxels; j++) {
//                    if (up) addPoints(cells[cell[0]+1][i][j], indexes, thisPlane, plane);
//                    addPoints(cells[cell[0]][i][j], indexes, thisPlane, plane);
//                    if (down) addPoints(cells[cell[0] - 1][i][j], indexes, thisPlane, plane);
//                }
//            }
//        }
//        else if (dir[1] == 1) { //variance in x and z planes
//            bool up = cell[1] < y_voxels - 1;
//            bool down = cell[1] > 0;
//            for (size_t i = 0; i < x_voxels; i++) {
//                for (size_t j = 0; j < z_voxels; j++) {
//                    if (up) addPoints(cells[i][cell[1] + 1][j], indexes, thisPlane, plane);
//                    addPoints(cells[i][cell[1]][j], indexes, thisPlane, plane);
//                    if (down) addPoints(cells[i][cell[1] - 1][j], indexes, thisPlane, plane);
//                }
//            }
//        }
//        else {// variance in x and y planes
//            bool up = cell[2] < z_voxels - 1;
//            bool down = cell[2] > 0;
//            for (size_t i = 0; i < x_voxels; i++) {
//                for (size_t j = 0; j < y_voxels; j++) {
//                    if (up) addPoints(cells[i][j][cell[2] + 1], indexes, thisPlane, plane);
//                    addPoints(cells[i][j][cell[2]], indexes, thisPlane, plane);
//                    if (down) addPoints(cells[i][j][cell[2] - 1], indexes, thisPlane, plane);
//                }
//            }
//        }
//        return indexes;
//    }
//
//    do { //p stays within the other bounds
//        cleary(indexes, Eigen::ParametrizedLine<double, 3>(p1, dir.normalized()), visited, thisPlane, plane);
//        p1 += step;
//    } while (start_line.projection(p1).norm() < start_line.projection(p2).norm());
//    return indexes;
//
//}


//Adds points in place within the threshold of a 2D ray in a 3D bounding box
void UniformPC::cleary(std::vector<size_t> &points, Eigen::ParametrizedLine<double, 3> ray, 
    std::vector<std::vector<std::vector<bool>>> &visited, Eigen::Hyperplane<double, 3> thisPlane) {

    //std::vector<int> limits = { x_voxels, y_voxels, z_voxels };
    Eigen::Vector3d p = ray.origin();
    //int x1 = ray.direction()[0] == 0 ? 1 : 0;
    //int x2 = ray.direction()[2] == 0 ? 1 : 2;

    //find current cell
    std::vector<size_t> cell = hashCell(p);

    ////if plane is completely flat, just do all the cells in that row and call it a day
    //if (x1 == x2) { //how to prevent overlaps?
    //    x2 = (x2 + 1) % 3;
    //    for (size_t i = 0; i < limits[x2]; i++) {
    //        cell[x2] = i;
    //        addPoints(cells[cell[0]][cell[1]][cell[2]], points, thisPlane, plane);
    //        cell[x1] += 1;
    //        if (cell[x1] < limits[x1]) addPoints(cells[cell[0]][cell[1]][cell[2]], points, thisPlane, plane);
    //        cell[x1] -= 2;
    //        if (cell[x1] < limits[x1]) addPoints(cells[cell[0]][cell[1]][cell[2]], points, thisPlane, plane);
    //        cell[x1] += 1;
    //    }
    //    return points;
    //}

    ////visit x2 adjacent cells just for the first step to catch scragglers
    //std::vector<size_t> cell2 = cell;
    //cell2[x2] += 1;
    //if (cell2[x2] < limits[x2]) {
    //    if (!cells[cell2[0]][cell2[1]][cell2[2]].empty())
    //        addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, thisPlane, plane);
    //    visited[cell2[0]][cell2[1]][cell2[2]] = true;
    //}
    //if (cell2[x2] > 1) {
    //    cell2[x2] -= 2;
    //    if (!cells[cell2[0]][cell2[1]][cell2[2]].empty())
    //        addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, thisPlane, plane);
    //    visited[cell2[0]][cell2[1]][cell2[2]] = true;
    //}

    //size_t up = cell[x3] + 1;
    //size_t down = cell[x3] - 1;
    //check these are correct?
    double theta_y = voxel_size / ray.direction()[1];
    double theta_z = voxel_size / ray.direction()[2];
    double dy = theta_y;
    double dz = theta_z;
    size_t y = 0;
    size_t z = 0;

    // the >= 0 check isn't needed because a size_t will simply overflow to greater than the limit when negative anyway!
    while (cell[1] < y_voxels && cell[2] < z_voxels) {
        if (!visited[cell[0]][cell[1]][cell[2]]) {
            if (!cells[cell[0]][cell[1]][cell[2]].empty()) //push thresholded points from cell onto points vector
                addPoints(cells[cell[0]][cell[1]][cell[2]], points, thisPlane);
            visited[cell[0]][cell[1]][cell[2]] = true;
        }

        //visit y adjacent cells to catch nearby points
        if (cell[1] < y_voxels - 1) {
            y = cell[1] + 1;
            if (!visited[cell[0]][y][cell[2]]) {
                if (!cells[cell[0]][y][cell[2]].empty())
                    addPoints(cells[cell[0]][y][cell[2]], points, thisPlane);
                visited[cell[0]][y][cell[2]] = true;
            }
        }
        if (cell[1] > 0) {
            y = cell[1] - 1;
            if (!visited[cell[0]][y][cell[2]]) {
                if (!cells[cell[0]][y][cell[2]].empty())
                    addPoints(cells[cell[0]][y][cell[2]], points, thisPlane);
                visited[cell[0]][y][cell[2]] = true;
            }
        }
        //visit z adjacent cells to catch nearby points
        if (cell[2] < z_voxels - 1) {
            z = cell[2] + 1;
            if (!visited[cell[0]][cell[1]][z]) {
                if (!cells[cell[0]][cell[1]][z].empty())
                    addPoints(cells[cell[0]][cell[1]][z], points, thisPlane);
                visited[cell[0]][cell[1]][z] = true;
            }
        }
        if (cell[2] > 0) {
            z = cell[2] - 1;
            if (!visited[cell[0]][cell[1]][z]) {
                if (!cells[cell[0]][cell[1]][z].empty())
                    addPoints(cells[cell[0]][cell[1]][z], points, thisPlane);
                visited[cell[0]][cell[1]][z] = true;
            }
        }

        //work out next cell
        if (abs(dy) < abs(dz)) {
            dy += theta_y;
            (theta_y > 0) ? cell[1] += 1 : cell[1] -= 1;
        }
        else {
            dz += theta_z;
            (theta_z > 0) ? cell[2] += 1 : cell[2] -= 1;
        }
    }

    ////wind back last out of bounds step
    //if (cell[x1] > limits[x1]) cell[x1] = 0;
    //if (cell[x1] == limits[x1]) cell[x1] = limits[x1] - 1;
    //if (cell[x2] > limits[x2]) cell[x2] = 0;
    //if (cell[x2] == limits[x2]) cell[x2] = limits[x2] - 1;

    ////visit x2 adjacent cells just for the last step to catch scragglers
    //cell2 = cell;
    //cell2[x2] += 1;
    //if (cell2[x2] < limits[x2]) {
    //    if (!visited[cell2[0]][cell2[1]][cell2[2]] && !cells[cell2[0]][cell2[1]][cell2[2]].empty())
    //        addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, thisPlane, plane);
    //    visited[cell2[0]][cell2[1]][cell2[2]] = true;
    //}
    //if (cell2[x2] > 1) {
    //    cell2[x2] -= 2;
    //    if (!visited[cell2[0]][cell2[1]][cell2[2]] && !cells[cell2[0]][cell2[1]][cell2[2]].empty())
    //        addPoints(cells[cell2[0]][cell2[1]][cell2[2]], points, thisPlane, plane);
    //    visited[cell2[0]][cell2[1]][cell2[2]] = true;
    //}

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
}
