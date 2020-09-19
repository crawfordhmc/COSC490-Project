#include "UniformPC.h"
#include "OctreePC.h"

void initialiseColours(std::vector<Eigen::Vector3i>* colours) {
    colours->resize(11);

    (*colours)[0](0) = 255;    (*colours)[0](1) = 0;      (*colours)[0](2) = 0; //red
    (*colours)[1](0) = 0;      (*colours)[1](1) = 255;    (*colours)[1](2) = 0; //green
    (*colours)[2](0) = 0;      (*colours)[2](1) = 0;      (*colours)[2](2) = 255; //blue
    (*colours)[3](0) = 255;    (*colours)[3](1) = 255;    (*colours)[3](2) = 0; //yellow
    (*colours)[4](0) = 0;      (*colours)[4](1) = 255;     (*colours)[4](2) = 255; //cyan
    (*colours)[5](0) = 255;    (*colours)[5](1) = 0;      (*colours)[5](2) = 255; //magenta
    (*colours)[6](0) = 0;      (*colours)[6](1) = 0;      (*colours)[6](2) = 0; //black
    (*colours)[7](0) = 255;    (*colours)[7](1) = 255;    (*colours)[7](2) = 255; //white
    (*colours)[8](0) = 127;    (*colours)[8](1) = 0;      (*colours)[8](2) = 0; //dark red
    (*colours)[9](0) = 0;      (*colours)[9](1) = 127;    (*colours)[9](2) = 0; //dark green
    (*colours)[10](0) = 0;     (*colours)[10](1) = 0;     (*colours)[10](2) = 127; //dark blue
}



std::vector<Eigen::Hyperplane<double, 3>> ransac(PointCloud& pointCloud, std::mt19937 gen, double successProb, double noise, double threshold, unsigned int maxTrials) {

    unsigned int plane = 0;
    std::vector<Eigen::Hyperplane<double, 3>> planes;

    //TEST
    Eigen::Vector3d norm = { 0, 1, 0 };
    Eigen::Hyperplane<double, 3> test = Eigen::Hyperplane<double, 3>(norm.normalized(), {0, 0, 0});
    std::vector<size_t> testy = pointCloud.planePoints(test);

    do {
        // Initial number of trials, very high from lowball initial inlier ratio
        double inlierRatio = 0.01;
        unsigned int numTrials = log(1 - successProb) / log(1 - pow(inlierRatio, 3));
        // Create random distribution for the point cloud, with other planes removed
        std::uniform_int_distribution<size_t> distr(0, pointCloud.remainingPoints.size() - 1);
        Eigen::Hyperplane<double, 3> bestPlane;
        std::vector<size_t> bestPoints;

        unsigned int trial = 0;
        while (trial < numTrials && trial < maxTrials) {

            // For each trial, generate a plane from 3 random point cloud indexes
            std::vector<size_t> foundPoints;
            while (foundPoints.size() < 3) {
                size_t index = distr(gen);
                foundPoints.push_back(pointCloud.remainingPoints[index]);
            }
            Eigen::Hyperplane<double, 3> thisPlane = Eigen::Hyperplane<double, 3>::Through(
                pointCloud.getPoint(foundPoints[0]).location,
                pointCloud.getPoint(foundPoints[1]).location,
                pointCloud.getPoint(foundPoints[2]).location);
            // Add points closer than threshold to this plane
            std::vector<size_t> thisPoints = pointCloud.planePoints(thisPlane);

            // Update plane with the most points
            if (thisPoints.size() > bestPoints.size()) {
                bestPlane = thisPlane;
                bestPoints = thisPoints;
                inlierRatio = (float)bestPoints.size() / (pointCloud.remainingPoints.size());
                numTrials = log(1 - successProb) / log(1 - pow(inlierRatio, 3));

            }
            trial++;
        }
        // Other operations with bestPlane can be done here
        std::cout << trial << " RANSAC trials run for plane " << plane + 1 << ", equation: " <<
            bestPlane.coeffs()[0] << "x + " << bestPlane.coeffs()[1] << "y + " << bestPlane.coeffs()[2] << "z + " << bestPlane.coeffs()[3] << " = 0" << std::endl;
        // Remove point indexes of the best plane from all trials
        for (size_t j = 0; j < bestPoints.size(); ++j) {
            pointCloud.setPointPlane(bestPoints[j], plane);
            //if bestpoints is ordered, this iterator can be saved between for loops
            std::vector<size_t>::iterator position = std::lower_bound(pointCloud.remainingPoints.begin(), pointCloud.remainingPoints.end(), bestPoints[j]);
            if (position != pointCloud.remainingPoints.end())
                pointCloud.remainingPoints.erase(position);
            // take iterator back by one to consider the removal?
        }
        plane++;
        planes.push_back(bestPlane);

    } while ((float)pointCloud.remainingPoints.size()/pointCloud.size > noise);

    return planes;

}


void recolor(PointCloud pointCloud, std::string outputFile, std::vector<Eigen::Vector3i> colours) {
    // Recolour points according to their plane then save the results
    std::cout << "Writing points to " << outputFile << std::endl;
    signed long long val = 0;
#pragma omp parallel for
    for (val = 0; val < pointCloud.size; ++val) {
        if (pointCloud.getPoint(val).planeIx >= 0) {
            pointCloud.setPointColour(val, colours[pointCloud.getPoint(val).planeIx % colours.size()]);
            // idea: color planes with limited colors based on avoiding intersecting plane's colors
        }
    }
    pointCloud.writeToPly(outputFile);
}


int main(int argc, char* argv[]) {

    // Command line arguments
    std::string inputFile = "";
    std::string outputFile = "";
    float success = 0.99;
    float noise = 0.01;
    float threshold = -1;
    int maxTrials = 1000;
    float scale_parameter = 0.01;
    std::string structure = "";

    // Parse the command line
    if (argc < 7 || argc > 8) {
        std::cout << "Usage: planeFinder <input file> <output file> <probability of success> <fraction of points to leave out> <max RANSAC trials> <scale factor>" << std::endl;
        std::cout << "optional: \"uniform\" or \"vector\" space subdivision" << std::endl;
        exit(-2);
    }
    else {
        inputFile = argv[1];
        outputFile = argv[2];
        success = atof(argv[3]);
        noise = atof(argv[4]);
        maxTrials = atoi(argv[5]);
        scale_parameter = atof(argv[6]);
        if (argc == 8)
            structure = argv[7];
    }

    // Set up random seed
    // (random size generation taken from https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution)
    std::random_device rd;
    std::mt19937 gen(rd());
    // Set up some colours to assign to the planes that are found
    std::vector<Eigen::Vector3i> colours;
    initialiseColours(&colours);

    std::vector<Eigen::Hyperplane<double, 3>> planes;

    PointCloud pointCloud = PointCloud(inputFile, scale_parameter);
    // Checking if number of points is too big for signed long long type
    if (pointCloud.size > LLONG_MAX) {
        std::cout << "Model is too big - reduce points to " << LLONG_MAX << " or less" << std::endl;
        return 1;
    }
    threshold = pointCloud.threshold;
    std::cout << "Auto-generated threshold is " << threshold << std::endl;
    if (planes.size() > colours.size()) std::cout << "Warning: more planes than colours" << std::endl;

    if (structure == "uniform") {
        UniformPC u = UniformPC(pointCloud, 10);
        std::vector<Eigen::Hyperplane<double, 3>> planes = ransac(u, gen, success, noise, threshold, maxTrials);
        std::cout << "Total point distance calculations made: " << u.comparisons << std::endl;
        recolor(u, outputFile, colours);
    }
    else {
        std::vector<Eigen::Hyperplane<double, 3>> planes = ransac(pointCloud, gen, success, noise, threshold, maxTrials);
        std::cout << "Total point distance calculations made: " << pointCloud.comparisons << std::endl;
        recolor(pointCloud, outputFile, colours);
    }
}

//    if (true) { //fix for structure
//        for (int p1 = 0; p1 < planes.size(); p1++) {
//            for (int p2 = 0; p2 < planes.size(); p2++) {
//                if (p1 == p2) continue; 
//                Eigen::ParametrizedLine<double, 3> *intersection = pointCloud->intersectPlanes(planes[p1], planes[p2]);
//                if (intersection == NULL) continue;
//                signed long long i = 0;
//#pragma omp parallel for
//                for (i = 0; i < pointCloud.size; ++i) {
//                    if (intersection->distance(pointCloud->getPoint(i).location) > threshold) continue;
//                    if (pointCloud->getPoint(i).planeIx == p1 && planes[p2].absDistance(pointCloud->getPoint(i).location) < planes[p1].absDistance(pointCloud.getPoint(i).location)) {
//                        pointCloud->setPointPlane(i, p2);
//                        std::cout << "reassigned some pointy bois" << std::endl;
//                    }
//                    if (pointCloud->getPoint(i).planeIx == p2 && planes[p1].absDistance(pointCloud.getPoint(i).location) < planes[p2].absDistance(pointCloud.getPoint(i).location)) {
//                        pointCloud.setPointPlane(i, p1);
//                        std::cout << "reassigned some pointy bois" << std::endl;
//                    }
//                }
//            }
//        }
//    }