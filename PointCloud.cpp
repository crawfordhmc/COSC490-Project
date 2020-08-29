#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <tinyply.h>
using namespace tinyply;


PointCloud::PointCloud(const std::string& filepath) {
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

		}
	catch (const std::exception& e)
	{
		std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
	}
 }


//Returns the line intersection of the given planes within the given bounding box, if any
//args: plane 1, plane 2, lower x boundary, upper x boundary, lower y boundary etc....
Eigen::ParametrizedLine<double, 3> PointCloud::intersectPlanes(Eigen::Hyperplane<double, 3> p1, Eigen::Hyperplane<double, 3> p2,
	double xs, double xl, double ys, double yl, double zs, double zl) {
	Eigen::Vector3d vec = p1.normal().cross(p2.normal());
	double a1 = p1.coeffs()[0];
	double b1 = p1.coeffs()[1];
	double c1 = p1.coeffs()[2];
	double d1 = p1.coeffs()[3];
	double a2 = p1.coeffs()[0];
	double b2 = p1.coeffs()[1];
	double c2 = p1.coeffs()[2];
	double d2 = p1.coeffs()[3];

	//assuming none of the coefficients = 0, unlikely with hyperplane?

	// assuming line intersects with x = xs
	double z = ((b2 / b1)*(a1*xs + d1) - a2*xs - d2) / (c2 - c1*b2 / b1);
	if (z > zs && z < zl) {
		double y = (-c1*z - a1*xs - d1)/b1;
		if (y > ys && y < yl)
			return Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(xs, y, z), vec);
	}
	// assuming line intersects with x = xl
	double z = ((b2/b1)*(a1*xl + d1) - a2*xl - d2)/(c2 - c1*b2 / b1);
	if (z > zs && z < zl) {
		double y = (-c1*z - a1*xl - d1)/b1;
		if (y > ys && y < yl)
			return Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(xl, y, z), vec);
	}

	// assuming line intersects with y = ys
	double x = ((c2/c1)*(b1*ys + d1) - b2*ys - d2)/(a2 - a1*c2/c1);
	if (x > xs && x < xl) {
		double z = (-a1*x - b1*ys - d1)/c1;
		if (z > zs && z < zl)
			return Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(x, ys, z), vec);
	}
	// assuming line intersects with y = yl
	double x = ((c2/c1)*(b1*yl + d1) - b2*yl - d2)/(a2 - a1*c2/c1);
	if (x > xs && x < xl) {
		double z = (-a1*x - b1*yl - d1)/c1;
		if (z > zs && z < zl)
			return Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(x, ys, z), vec);
	}

	// assuming line intersects with z = zs
	double y = ((a2/a1)*(c1*zs + d1) - c2*zs - d2)/(b2 - b1*a2 / a1);
	if (y > ys && y < yl) {
		double x = (-b1*y - c1*zs - d1)/a1;
		if (x > xs && x < xl)
			return Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(x, y, zs), vec);
	}
	// assuming line intersects with z = zl
	double y = ((a2/a1)*(c1*zl + d1) - c2*zl - d2)/(b2 - b1*a2 / a1);
	if (y > ys && y < yl) {
		double x = (-b1*y - c1*zl - d1)/a1;
		if (x > xs && x < xl)
			return Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(x, y, zs), vec);
	}

	// if the plane intersection line does intersect with the lower x, y or z plane of the bounding box,
	// the two planes do not intersect anywhere within the bounding box
	return Eigen::ParametrizedLine<double, 3>(Eigen::Vector3d(NAN, NAN, NAN), Eigen::Vector3d(NAN, NAN, NAN));
}
