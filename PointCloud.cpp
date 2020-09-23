#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <tinyply.h>
using namespace tinyply;


PointCloud::PointCloud(const std::string& filepath, float scale_parameter) {
	std::cout << "........................................................................\n";
	std::cout << "Now Reading: " << filepath << std::endl;
	std::unique_ptr<std::istream> file_stream;
	std::vector<uint8_t> byte_buffer;

	try
		{
		file_stream.reset(new std::ifstream(filepath, std::ios::binary));

		if (!file_stream || file_stream->fail()) throw std::runtime_error("file_stream failed to open " + filepath);

		file_stream->seekg(0, std::ios::end);
		const float size_mb = file_stream->tellg() * float(1e-6);
		file_stream->seekg(0, std::ios::beg);

		PlyFile file;
		file.parse_header(*file_stream);

		// Because most people have their own mesh types, tinyply treats parsed data as structured/typed byte buffers. 
		// See examples below on how to marry your own application-specific data structures with this one. 
		std::shared_ptr<PlyData> vertices, normals, colors, texcoords, faces, tripstrip;

		// The header information can be used to programmatically extract properties on elements
		// known to exist in the header prior to reading the data. For brevity of this sample, properties 
		// like vertex position are hard-coded: 
		try { vertices = file.request_properties_from_element("vertex", { "x", "y", "z" }); }
		catch (const std::exception& e) { //std::cerr << "tinyply exception: " << e.what() << std::endl; 
		} try { normals = file.request_properties_from_element("vertex", { "nx", "ny", "nz" }); }
		catch (const std::exception& e) { //std::cerr << "tinyply exception: " << e.what() << std::endl; 
		} try { colors = file.request_properties_from_element("vertex", { "red", "green", "blue" }); }
		catch (const std::exception& e) { //std::cerr << "tinyply exception: " << e.what() << std::endl; 
		} try { colors = file.request_properties_from_element("vertex", { "r", "g", "b", "a" }); }
		catch (const std::exception& e) { //std::cerr << "tinyply exception: " << e.what() << std::endl; 
		} try { texcoords = file.request_properties_from_element("vertex", { "u", "v" }); }
		catch (const std::exception& e) { //std::cerr << "tinyply exception: " << e.what() << std::endl; 
		}
		// Providing a list size hint (the last argument) is a 2x performance improvement. If you have 
		// arbitrary ply files, it is best to leave this 0. 
		try { faces = file.request_properties_from_element("face", { "vertex_indices" }, 3); }
		catch (const std::exception& e) { //std::cerr << "tinyply exception: " << e.what() << std::endl; 
		}
		// Tristrips must always be read with a 0 list size hint (unless you know exactly how many elements
		// are specifically in the file, which is unlikely); 
		try { tripstrip = file.request_properties_from_element("tristrips", { "vertex_indices" }, 0); }
		catch (const std::exception& e) { //std::cerr << "tinyply exception: " << e.what() << std::endl; 
		}

		file.read(*file_stream);
		std::cout << "PLY file read by tinyply: converting to PointCloud" << std::endl;

		// copy locations to vector
		const size_t numVerticesBytes = vertices->buffer.size_bytes();
		verts = std::vector<Eigen::Vector3f>(vertices->count); // assuming input values are floats
		std::memcpy(verts.data(), vertices->buffer.get(), numVerticesBytes);
		// copy colors to vector
		const size_t numColorsBytes = colors->buffer.size_bytes();
		cols = std::vector<Eigen::Matrix<uint8_t, 3, 1>>(colors->count); // assuming input values are 8 bit unsigned ints
		std::memcpy(cols.data(), colors->buffer.get(), numColorsBytes);

		Point defaultPoint;

		defaultPoint.location = Eigen::Vector3d::Zero();
		defaultPoint.colour = Eigen::Vector3i::Zero();
		defaultPoint.planeIx = -1;
		pc.resize(verts.size(), defaultPoint);
		remainingPoints.resize(pc.size(), 0);

		signed long long i = 0;
		int threads = 1;
#pragma omp parallel for
		for (i = 0; i < pc.size(); i++) {
			pc[i].location = verts[i].cast<double>();
			pc[i].colour = cols[i].cast<int>();
			remainingPoints[i] = i;
			if (omp_get_thread_num() == 0) threads = omp_get_num_threads();
		}
		std::cout << threads << " threads being used" << std::endl;

		XS = pc[0].location[0];
		XL = pc[0].location[0];
		YS = pc[0].location[1];
		YL = pc[0].location[1];
		ZS = pc[0].location[2];
		ZL = pc[0].location[2];

		for (size_t p = 1; p < pc.size(); ++p) {
			XS = std::min(XS, pc[p].location[0]);
			XL = std::max(XL, pc[p].location[0]);
			YS = std::min(YS, pc[p].location[1]);
			YL = std::max(YL, pc[p].location[1]);
			ZS = std::min(ZS, pc[p].location[2]);
			ZL = std::max(ZL, pc[p].location[2]);
		}

		size = pc.size();
		// get x/y/z difference and compute average scale factor for the model
		double scale = (XL - XS + YL - YS + ZL - ZS) / 3;
		// apply a small % to the value to get a sensible threshold
		threshold = scale_parameter * scale;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
	}
 }

PointCloud::Point PointCloud::getPoint(size_t index) { return pc[index]; }
void PointCloud::setPointPlane(size_t index, int planeID) { pc[index].planeIx = planeID; }
void PointCloud::setPointColour(size_t index, const Eigen::Vector3i& colour) { pc[index].colour = colour; }


// Returns a vector of points within the threshold to the given hyperplane
// (also prints the number of threads being used for the calculations)
std::vector<size_t> PointCloud::planePoints(const Eigen::Hyperplane<double, 3>& thisPlane) {
	std::vector<size_t> thisPoints;
	//OpenMP requires signed integrals for its loop variables... interesting
	signed long long i = 0;
#pragma omp parallel for shared(thisPoints) private (i)
	for (i = 0; i < remainingPoints.size(); ++i) {
		if (thisPlane.absDistance(pc[remainingPoints[i]].location) < threshold)
#pragma omp critical
			thisPoints.push_back(remainingPoints[i]);
	}
	comparisons += remainingPoints.size();
	return thisPoints;
}


//Returns the line intersection of the given planes within the given bounding box, if any
//args: plane 1, plane 2, lower x boundary, upper x boundary, lower y boundary etc....
Eigen::ParametrizedLine<double, 3>* PointCloud::intersectPlanes(const Eigen::Hyperplane<double, 3>& p1, const Eigen::Hyperplane<double, 3>& p2,
	double xs, double xl, double ys, double yl, double zs, double zl) {
	// direction of the intersection line
	Eigen::Vector3d vec = p1.normal().cross(p2.normal());
	vec.normalize();
	double a1 = p1.coeffs()[0];
	double b1 = p1.coeffs()[1];
	double c1 = p1.coeffs()[2];
	double d1 = p1.coeffs()[3];
	double a2 = p2.coeffs()[0];
	double b2 = p2.coeffs()[1];
	double c2 = p2.coeffs()[2];
	double d2 = p2.coeffs()[3];
	double x, y, z;

	// assuming none of the coefficients are 0, unlikely with double prescision coefficients?

	// assuming line intersects with x = xs
	z = ((b2 / b1) * (a1 * xs + d1) - a2 * xs - d2) / (c2 - c1 * b2 / b1);
	if (z > zs && z < zl) {
		y = (-c1 * z - a1 * xs - d1) / b1;
		if (y > ys && y < yl) {
			return &Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(xs, y, z), vec);
		}
	}
	// assuming line intersects with x = xl
	z = ((b2 / b1) * (a1 * xl + d1) - a2 * xl - d2) / (c2 - c1 * b2 / b1);
	if (z > zs && z < zl) {
		y = (-c1 * z - a1 * xl - d1) / b1;
		if (y > ys && y < yl)
			return &Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(xl, y, z), vec);
	}

	// assuming line intersects with y = ys
	x = ((c2 / c1) * (b1 * ys + d1) - b2 * ys - d2) / (a2 - a1 * c2 / c1);
	if (x > xs && x < xl) {
		z = (-a1 * x - b1 * ys - d1) / c1;
		if (z > zs && z < zl)
			return &Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(x, ys, z), vec);
	}
	// assuming line intersects with y = yl
	x = ((c2 / c1) * (b1 * yl + d1) - b2 * yl - d2) / (a2 - a1 * c2 / c1);
	if (x > xs && x < xl) {
		z = (-a1 * x - b1 * yl - d1) / c1;
		if (z > zs && z < zl)
			return &Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(x, ys, z), vec);
	}

	// assuming line intersects with z = zs
	y = ((a2 / a1) * (c1 * zs + d1) - c2 * zs - d2) / (b2 - b1 * a2 / a1);
	if (y > ys && y < yl) {
		x = (-b1 * y - c1 * zs - d1) / a1;
		if (x > xs && x < xl)
			return &Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(x, y, zs), vec);
	}
	// assuming line intersects with z = zl
	y = ((a2 / a1) * (c1 * zl + d1) - c2 * zl - d2) / (b2 - b1 * a2 / a1);
	if (y > ys && y < yl) {
		x = (-b1 * y - c1 * zl - d1) / a1;
		if (x > xs && x < xl)
			return &Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(x, y, zs), vec);
	}
	// the two planes do not intersect anywhere within the bounding box
	return NULL;
}

//Returns the line intersection of the given planes within the model's bounding box, if any
//args: plane 1, plane 2, lower x boundary, upper x boundary, lower y boundary etc....
Eigen::ParametrizedLine<double, 3>* PointCloud::intersectPlanes(const Eigen::Hyperplane<double, 3>& p1, const Eigen::Hyperplane<double, 3>& p2) {
	return PointCloud::intersectPlanes(p1, p2, XS - threshold, XL + threshold, YS - threshold, YL + threshold, ZS - threshold, ZL + threshold);
	//return PointCloud::intersectPlanes(p1, p2, XS, XL, YS, YL, ZS, ZL);

}


void PointCloud::writeToPly(const std::string& filename) {
	std::ofstream fout(filename);
	fout << "ply\n" //write the header
		<< "format ascii 1.0\n"
		<< "element vertex " << size << "\n"
		<< "property float x\n"
		<< "property float y\n"
		<< "property float z\n"
		<< "property uchar red\n"
		<< "property uchar green\n"
		<< "property uchar blue\n"
		<< "end_header\n";
	for (auto point : pc) {
		fout << point.location.transpose() << " "
			<< point.colour.transpose() << "\n"; //output location and color
	}
	fout.close();
}

void PointCloud::removePoints(std::vector<size_t> &planePoints, int plane) {
	std::vector<size_t> diff;
	std::sort(planePoints.begin(), planePoints.end());
	std::set_difference(remainingPoints.begin(), remainingPoints.end(), 
		planePoints.begin(), planePoints.end(), std::inserter(diff, diff.begin()));
	remainingPoints = diff;

	signed long long j = 0;
#pragma omp parallel for
	for (j = 0; j < planePoints.size(); ++j) {
		setPointPlane(planePoints[j], plane);
	}
}