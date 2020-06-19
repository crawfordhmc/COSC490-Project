#pragma once

#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <string>
#include <vector>
#include <random>

struct Point {
  Eigen::Vector3d location;
  Eigen::Vector3i colour;
  int planeIx;
};

typedef std::vector<Point> PointCloud;

PointCloud readFromPly(const std::string& filename);

void writeToPly(const PointCloud& pc, const std::string& filename);

#endif
