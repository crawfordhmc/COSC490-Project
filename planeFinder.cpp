#include "PointCloud.h"

void initialiseColours(std::vector<Eigen::Vector3i> *colours) {
  colours->resize(6);
  
  (*colours)[0](0) = 255;    (*colours)[0](1) = 0;      (*colours)[0](2) = 0;
  (*colours)[1](0) = 0;      (*colours)[1](1) = 255;    (*colours)[1](2) = 0;
  (*colours)[2](0) = 0;      (*colours)[2](1) = 0;      (*colours)[2](2) = 255;
  (*colours)[3](0) = 255;    (*colours)[3](1) = 255;    (*colours)[3](2) = 0;
  (*colours)[4](0) = 0;      (*colours)[4](1) = 255;     (*colours)[4](2) = 255;
  (*colours)[5](0) = 255;    (*colours)[5](1) = 0;      (*colours)[5](2) = 255;
}


PointCloud ransac(PointCloud pointCloud, std::mt19937 gen, float successProb, float inlierRatio, int numPlanes, float threshold) {

    int numTrials = log(1 - successProb) / log(1 - pow(inlierRatio, 3));
    std::vector<size_t> removedPoints;
    //planes cannot take points from each other in parallel - implement auto plane detection then parallelize
    for (size_t plane = 0; plane < numPlanes; ++plane) { //NOTE: should plane be size_t or int??

        // Create random distribution for the point cloud
        std::uniform_int_distribution<size_t> distr(0, pointCloud.size() - 1);
        Eigen::Hyperplane<double, 3> bestPlane;
        std::vector<size_t> bestPoints;
        //#pragma omp parallel for each trial, queue for points, calculate, queue to compare and update to bestPlane?
        for (size_t trial = 0; trial < numTrials; ++trial) {
            // If not enough points remaining not on a plane, continue to next trial
            if (pointCloud.size() - removedPoints.size() < 3) {
                std::cout << "RANSAC trial " << (trial + 1) << " failed" << std::endl;
                continue;
            }
            // For each trial, generate a plane from 3 random point cloud indexes
            std::vector<size_t> foundPoints;
            while (foundPoints.size() < 3) {
                size_t index = distr(gen);
                // Save random point if not already part of a plane
                if (removedPoints.end() == find(removedPoints.begin(), removedPoints.end(), index))
                    foundPoints.push_back(index);
            }
            Eigen::Hyperplane<double, 3> thisPlane = Eigen::Hyperplane<double, 3>::Through(
                pointCloud[foundPoints[0]].location,
                pointCloud[foundPoints[1]].location,
                pointCloud[foundPoints[2]].location);
            // Add points closer than threshold to this plane
            std::vector<size_t> thisPoints;
            //OpenMP requires signed integrals for its loop variables... interesting
            signed long long i = 0;
#pragma omp parallel for shared(thisPoints) private (i)
            for (i = 0; i < pointCloud.size(); ++i) {
                if (thisPlane.absDistance(pointCloud[i].location) < threshold)
#pragma omp critical
                    thisPoints.push_back(i);
            }
            // Update plane with the most points
            if (thisPoints.size() > bestPoints.size()) {
                bestPlane = thisPlane;
                bestPoints = thisPoints;
            }
        }
        // Save best plane from all trials and its point indexes to be removed
        for (size_t j = 0; j < bestPoints.size(); ++j) {
            pointCloud[bestPoints[j]].planeIx = plane;
            removedPoints.push_back(bestPoints[j]);
        }

    }
    return pointCloud;

}


int main(int argc, char* argv[]) {

    // Command line arguments
    std::string inputFile = "";
    std::string outputFile = "";
    int numPlanes = 0;
    double threshold = 0;

    // Parse the command line
    if (argc != 5) {
        std::cout << "Usage: planeFinder <input file> <output file> <number of planes> <distance threshold, -1 for automatic>" << std::endl;
        exit(-2);
    }
    else {
        inputFile = argv[1];
        outputFile = argv[2];
        numPlanes = atoi(argv[3]);
        threshold = atof(argv[4]);
    }

    // Set up random seed
    // (random size generation taken from https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution)
    std::random_device rd;
    std::mt19937 gen(rd());

    // Read a point cloud from a file
    std::cout << "Reading points from " << inputFile << ", planes = " << numPlanes << ", threshold = " << threshold << std::endl;
    PointCloud pointCloud = readFromPly(inputFile);

    // Checking if number of points is too big for signed long long type (this aint gonna happen lmao)
    if (pointCloud.size() > LLONG_MAX) {
        std::cout << "Model is too big - reduce points to " << LLONG_MAX << " or less" << std::endl;
        return 1;
    }

    // Determine the threshold as a % of model size if auto detect enabled
    if (threshold == -1) {
        // (coordinate center is all over the place, so biggest/smallest signed point difference gives bounding box)
        double xs = pointCloud[0].location[0];
        double xl = pointCloud[0].location[0];
        double ys = pointCloud[0].location[1];
        double yl = pointCloud[0].location[1];
        double zs = pointCloud[0].location[2];
        double zl = pointCloud[0].location[2];
        // chunk parallelize this if its slow?
        for (size_t i = 1; i < pointCloud.size(); i++) {
            xs = std::min(xs, pointCloud[i].location[0]);
            xl = std::max(xl, pointCloud[i].location[0]);
            ys = std::min(ys, pointCloud[i].location[1]);
            yl = std::max(yl, pointCloud[i].location[1]);
            zs = std::min(zs, pointCloud[i].location[2]);
            zl = std::max(zl, pointCloud[i].location[2]);
        }
        // get x/y/z difference and compute average scale factor for the model??
        double scale = (xl - xs + yl - ys + zl - zs) / 3;
        // apply a small % to the value to get a sensible threshold
        threshold = 0.08*scale;
        std::cout << "Auto-generated threshold is " << threshold << std::endl;
    }

    // Set up some colours to assign to the planes that are found
    std::vector<Eigen::Vector3i> colours;
    initialiseColours(&colours);

    // prob could be input or hardcoded. HOW TO GET OUTLIER RATIO
    pointCloud = ransac(pointCloud, gen, 0.99, 0.8, numPlanes, threshold);

    // Recolour points according to their plane then save the results
    std::cout << "Writing points to " << outputFile << std::endl;
    for (auto& point : pointCloud) {
        if (point.planeIx >= 0) {
            point.colour = colours[point.planeIx % colours.size()];
        }
    }
    writeToPly(pointCloud, outputFile);
}
