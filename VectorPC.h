#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class VectorPC : public PointCloud {
private:
    std::vector<Point> pc;

public:

    VectorPC(const std::string& filepath);
    Point getPoint(int index);
    void setPointPlane(int index, int planeID);
    void setPointColour(int index, Eigen::Vector3i colour);
    float threshold(float scale_parameter);
    void writeToPly(const std::string& filename);
};