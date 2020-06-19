#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>

PointCloud readFromPly(const std::string& filename) {
  std::ifstream fin(filename);
  std::string token = "";
 
  // read up to the vertex block, counting number of lines to skip
  size_t linesToSkip = 0;
  fin >> token;
  if (token != "ply") {
    std::cerr << "ERROR: " << filename << " is not a PLY file" << std::endl;
    exit(-1);
  }
  fin >> token >> token;
  if (token != "ascii") {
    std::cerr << "ERROR: " << filename << " is not an ASCII PLY file" << std::endl;
    exit(-1);
  }
  while (fin >> token && token != "end_header") {
    if (token == "element") {
      fin >> token;
      if (token == "vertex") {
	break;
      } else {
	size_t nElements;
	fin >> nElements;
	linesToSkip += nElements;
      }
    }
  }
  int xIx = -1;
  int yIx = -1;
  int zIx = -1;
  int rIx = -1;
  int gIx = -1;
  int bIx = -1;

  PointCloud pc;
  Point defaultPoint;
  defaultPoint.location = Eigen::Vector3d::Zero();
  defaultPoint.colour = Eigen::Vector3i::Zero();
  defaultPoint.planeIx = -1;
  if (token == "vertex") {
    size_t nVertices;
    fin >> nVertices;
    pc.resize(nVertices, defaultPoint);
    fin >> token;
    int pIx = 0;
    while (token == "property") {
      fin >> token >> token;
      if (token == "x") {
	xIx = pIx;
      }
      if (token == "y") {
	yIx = pIx;
      }
      if (token == "z") {
	zIx = pIx;
      }
      if (token == "red") {
	rIx = pIx;
      }
      if (token == "green") {
	gIx = pIx;
      }
      if (token == "blue") {
	bIx = pIx;
      }
      ++pIx;
      fin >> token;
    } 
  } else {
    std::cerr << "ERROR: no vertex information in file" << std::endl;
    exit(-1);
  }

  if (xIx < 0 || yIx < 0 || zIx < 0) {
    std::cerr << "ERROR: vertex must have x, y, and z properties" << std::endl;
  }

  while (token != "end_header") {
    fin >> token;
  }

  std::string line;
  std::getline(fin, line); // skip past new line

  // Skip over items before vertices
  for (size_t i = 0; i < linesToSkip; ++i) {
    std::getline(fin, line);
  }

  for (size_t p = 0; p < pc.size(); ++p) {
    std::getline(fin, line);
    std::stringstream ss(line);
    double d;
    int ix = 0;
    while (ss >> d) {
      if (ix == xIx) {
	pc[p].location(0) = d;
      } else if (ix == yIx) {
	pc[p].location(1) = d;	
      } else if (ix == zIx) {
	pc[p].location(2) = d;
      } else if (ix == rIx) {
	pc[p].colour(0) = d;
      } else if (ix == gIx) {
	pc[p].colour(1) = d;
      } else if (ix == bIx) {
	pc[p].colour(2) = d;
      }

      ++ix;
    }
  }

  fin.close();

  return pc;
}


void writeToPly(const PointCloud& pc, const std::string& filename) {
  std::ofstream fout(filename);
  fout << "ply\n"
       << "format ascii 1.0\n"
       << "element vertex " << pc.size() << "\n"
       << "property float x\n"
       << "property float y\n"
       << "property float z\n"
       << "property uchar red\n"
       << "property uchar green\n"
       << "property uchar blue\n"
       << "end_header\n";
  for (auto point: pc) {
    fout << point.location.transpose() << " "
	 << point.colour.transpose() << "\n";
  }
  fout.close();
}
