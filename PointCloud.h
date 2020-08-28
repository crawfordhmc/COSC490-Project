#pragma once

#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <omp.h>

class PointCloud {

public:
	struct Point {
		Eigen::Vector3d location;
		Eigen::Vector3i colour;
		int planeIx;
	};

	size_t size = 0;

	PointCloud(const std::string& filepath);
	Eigen::ParametrizedLine<double, 3> PointCloud::intersectPlanes(Eigen::Hyperplane<double, 3> p1, Eigen::Hyperplane<double, 3> p2,
		double xs, double xl, double ys, double yl, double zs, double zl);

	virtual Point getPoint(int index) = 0;
	virtual void setPointPlane(int index, int planeID) = 0;
	virtual void setPointColour(int index, Eigen::Vector3i colour) = 0;
	virtual std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane, unsigned int trial, float threshold, int plane) = 0;
	virtual float threshold(float scale_parameter) = 0;
	virtual void writeToPly(const std::string& filename) = 0;

protected:
	std::vector<Eigen::Vector3f> verts;
	std::vector<Eigen::Matrix<uint8_t, 3, 1>> cols;
};

#endif
