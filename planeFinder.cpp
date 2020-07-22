#include "PointCloud.h"

void initialiseColours(std::vector<Eigen::Vector3i>* colours) {
    colours->resize(11);

    (*colours)[0](0) = 255;    (*colours)[0](1) = 0;      (*colours)[0](2) = 0;
    (*colours)[1](0) = 0;      (*colours)[1](1) = 255;    (*colours)[1](2) = 0;
    (*colours)[2](0) = 0;      (*colours)[2](1) = 0;      (*colours)[2](2) = 255;
    (*colours)[3](0) = 255;    (*colours)[3](1) = 255;    (*colours)[3](2) = 0;
    (*colours)[4](0) = 0;      (*colours)[4](1) = 255;     (*colours)[4](2) = 255;
    (*colours)[5](0) = 255;    (*colours)[5](1) = 0;      (*colours)[5](2) = 255;
    (*colours)[6](0) = 0;      (*colours)[6](1) = 0;      (*colours)[6](2) = 0;
    (*colours)[7](0) = 255;    (*colours)[7](1) = 255;    (*colours)[7](2) = 255;
    (*colours)[8](0) = 127;    (*colours)[8](1) = 0;      (*colours)[8](2) = 0;
    (*colours)[9](0) = 0;      (*colours)[9](1) = 127;    (*colours)[9](2) = 0;
    (*colours)[10](0) = 0;     (*colours)[10](1) = 0;     (*colours)[10](2) = 127;
}


//idea: get user estimate of minimum number of planes
// When running loads of trials, keep that many of the next largest planes
// When a plane is found, run through the next largest planes and reassign points by distance, then re-rank and save
// Then run more trials if % of the scene is not explained

PointCloud ransac(PointCloud pointCloud, std::mt19937 gen, double successProb, double explained, double threshold, unsigned int maxTrials) {

    unsigned int plane = 0;
    std::vector<size_t> removedPoints;
    int threads;
    omp_set_num_threads(4);

    do {
        // Initial number of trials, very high from lowball initial inlier ratio
        double inlierRatio = 0.1;
        unsigned int numTrials = log(1 - successProb) / log(1 - pow(inlierRatio, 3));
        // Create random distribution for the point cloud, with other planes removed
        std::uniform_int_distribution<size_t> distr(0, pointCloud.size() - 1);
        Eigen::Hyperplane<double, 3> bestPlane;
        std::vector<size_t> bestPoints;

        int trial = 0;
        while (trial < numTrials && trial < maxTrials) {

            // If not enough points remaining not on a plane, continue to next trial
            if (pointCloud.size() - removedPoints.size() < 3) {
                std::cout << "RANSAC trial " << (trial + 1) << " failed" << std::endl;
                trial++;
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
                // limit a single thread on the first trial of the first plane to printing the number of threads being used
                if (trial == 0 && omp_get_thread_num() == 0 && plane == 0)
                    threads = omp_get_num_threads();
            }
            // Update plane with the most points
            if (thisPoints.size() > bestPoints.size()) {
                bestPlane = thisPlane;
                bestPoints = thisPoints;
                inlierRatio = (float)bestPoints.size() / (pointCloud.size() - removedPoints.size());
                numTrials = log(1 - successProb) / log(1 - pow(inlierRatio, 3));

            }
            trial++;
        }
        // Other operations with bestPlane can be done here
        std::cout << trial << " RANSAC trials run for plane " << plane + 1 << ", equation: " <<
            bestPlane.coeffs()[0] << "x + " << bestPlane.coeffs()[1] << "y + " << bestPlane.coeffs()[2] << "z + " << bestPlane.coeffs()[3] << " = 0" << std::endl;
        // Save point indexes of the best plane from all trials to be removed
        for (size_t j = 0; j < bestPoints.size(); ++j) {
            pointCloud[bestPoints[j]].planeIx = plane;
            removedPoints.push_back(bestPoints[j]); // opportunity to do reclaiming here?
            // could only compare planes that intersect close to the bounding box?
        }
        plane++;

    } while ((float)removedPoints.size() / pointCloud.size() < explained);
    std::cout << omp_get_max_threads() << " threads were used" << std::endl;

    return pointCloud;

}


int main(int argc, char* argv[]) {

    // Command line arguments
    std::string inputFile = "";
    std::string outputFile = "";
    float success = 0.99;
    float explained = 0.99;
    float threshold = -1;
    int maxTrials = 1000;
    float scale_parameter = 0.01;

    // Parse the command line
    if (argc != 7) {
        std::cout << "Usage: planeFinder <input file> <output file> <probability of success> <ratio of scene to be explained by planes> <max RANSAC trials> <scale factor>" << std::endl;
        exit(-2);
    }
    else {
        inputFile = argv[1];
        outputFile = argv[2];
        success = atof(argv[3]);
        explained = atof(argv[4]);
        maxTrials = atoi(argv[5]);
        scale_parameter = atof(argv[6]);
    }

    // Set up random seed
    // (random size generation taken from https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution)
    std::random_device rd;
    std::mt19937 gen(rd());

    PointCloud pointCloud = readFromPly(inputFile);

    // Checking if number of points is too big for signed long long type (this aint gonna happen lmao)
    if (pointCloud.size() > LLONG_MAX) {
        std::cout << "Model is too big - reduce points to " << LLONG_MAX << " or less" << std::endl;
        return 1;
    }

    // Determine the threshold as a % of model size
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
    // get x/y/z difference and compute average scale factor for the model
    double scale = (xl - xs + yl - ys + zl - zs) / 3;
    // apply a small % to the value to get a sensible threshold
    threshold = scale_parameter * scale;
    std::cout << "Auto-generated threshold is " << threshold << std::endl;

    // Set up some colours to assign to the planes that are found
    std::vector<Eigen::Vector3i> colours;
    initialiseColours(&colours);

    pointCloud = ransac(pointCloud, gen, success, explained, threshold, maxTrials);

    // Recolour points according to their plane then save the results
    std::cout << "Writing points to " << outputFile << std::endl;
    for (auto& point : pointCloud) {
        if (point.planeIx >= 0) {
            point.colour = colours[point.planeIx % colours.size()];
            // idea: color planes with limited colors based on avoiding intersecting plane's colors
        }
    }
    writeToPly(pointCloud, outputFile);
}