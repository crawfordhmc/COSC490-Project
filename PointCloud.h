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

	//number of points
	size_t size = 0;

	PointCloud(const std::string& filepath);
	Eigen::ParametrizedLine<double, 3> PointCloud::intersectPlanes(Eigen::Hyperplane<double, 3> p1, Eigen::Hyperplane<double, 3> p2,
		double xs, double xl, double ys, double yl, double zs, double zl);
	Eigen::ParametrizedLine<double, 3> PointCloud::intersectPlanes(Eigen::Hyperplane<double, 3> p1, Eigen::Hyperplane<double, 3> p2);

	virtual Point getPoint(int index) = 0;
	virtual void setPointPlane(int index, int planeID) = 0;
	virtual void setPointColour(int index, Eigen::Vector3i colour) = 0;
	virtual std::vector<size_t> planePoints(Eigen::Hyperplane<double, 3> thisPlane, unsigned int trial, float threshold, int plane) = 0;
	float threshold(float scale_parameter);
	virtual void writeToPly(const std::string& filename) = 0;

protected:
	std::vector<Eigen::Vector3f> verts;
	std::vector<Eigen::Matrix<uint8_t, 3, 1>> cols;

	std::vector<Point> pc;
	// interesting point here - a unique_pointer can make a fixed but dynamic array when we know we won't be expanding after creation
	// could improve memory usage
	// auto ints = std::make_unique<int[]>(10);

	// bounding box corners
	double XS;
	double XL;
	double YS;
	double YL;
	double ZS;
	double ZL;
};

#endif
