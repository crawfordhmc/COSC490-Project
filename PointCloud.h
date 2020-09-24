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
	double threshold;
	std::vector<size_t> remainingPoints;
	size_t comparisons = 0;
	unsigned int num_threads = 1;

	PointCloud(const std::string& filepath, float scale_parameter, unsigned int threads);

	virtual std::vector<size_t> planePoints(const Eigen::Hyperplane<double, 3> &thisPlane);

	Point getPoint(size_t index);
	void setPointPlane(size_t index, int planeID);
	void setPointColour(size_t index, const Eigen::Vector3i &colour);
	void writeToPly(const std::string& filename);
	virtual void removePoints(std::vector<size_t>& planePoints, int plane);
	virtual void resetRemaining();

protected:
	std::vector<Eigen::Vector3f> verts;
	std::vector<Eigen::Matrix<uint8_t, 3, 1>> cols;

	std::vector<Point> pc;

	// bounding box corners
	double XS;
	double XL;
	double YS;
	double YL;
	double ZS;
	double ZL;
};

#endif
