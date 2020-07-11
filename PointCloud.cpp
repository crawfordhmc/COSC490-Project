#include "PointCloud.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <tinyply.h>
using namespace tinyply;


PointCloud tinyReadFromPly(const std::string& filepath)
{
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

        /*std::cout << "\t[ply_header] Type: " << (file.is_binary_file() ? "binary" : "ascii") << std::endl;
        for (const auto& c : file.get_comments()) std::cout << "\t[ply_header] Comment: " << c << std::endl;
        for (const auto& c : file.get_info()) std::cout << "\t[ply_header] Info: " << c << std::endl;

        for (const auto& e : file.get_elements())
        {
            std::cout << "\t[ply_header] element: " << e.name << " (" << e.size << ")" << std::endl;
            for (const auto& p : e.properties)
            {
                std::cout << "\t[ply_header] \tproperty: " << p.name << " (type=" << tinyply::PropertyTable[p.propertyType].str << ")";
                if (p.isList) std::cout << " (list_type=" << tinyply::PropertyTable[p.listType].str << ")";
                std::cout << std::endl;
            }
        }*/

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

        /*if (vertices)   std::cout << "\tRead " << vertices->count << " total vertices " << std::endl;
        if (normals)    std::cout << "\tRead " << normals->count << " total vertex normals " << std::endl;
        if (colors)     std::cout << "\tRead " << colors->count << " total vertex colors " << std::endl;
        if (texcoords)  std::cout << "\tRead " << texcoords->count << " total vertex texcoords " << std::endl;
        if (faces)      std::cout << "\tRead " << faces->count << " total faces (triangles) " << std::endl;
        if (tripstrip)  std::cout << "\tRead " << (tripstrip->buffer.size_bytes() / tinyply::PropertyTable[tripstrip->t].stride) << " total indicies (tristrip) " << std::endl;
        */

        // copy locations to vector
        const size_t numVerticesBytes = vertices->buffer.size_bytes();
        std::vector<Eigen::Vector3f> verts(vertices->count); // assuming input values are floats
        std::memcpy(verts.data(), vertices->buffer.get(), numVerticesBytes);
        // copy colors to vector
        const size_t numColorsBytes = colors->buffer.size_bytes();
        std::vector<Eigen::Matrix<uint8_t, 3, 1>> cols(colors->count); // assuming input values are 8 bit unsigned ints
        std::memcpy(cols.data(), colors->buffer.get(), numColorsBytes);

        // Convert to PointCloud
        PointCloud pc;
        Point defaultPoint;
        defaultPoint.location = Eigen::Vector3d::Zero();
        defaultPoint.colour = Eigen::Vector3i::Zero();
        defaultPoint.planeIx = -1;
        pc.resize(vertices->count, defaultPoint);
        signed long long p;
#pragma omp parallel for private(p)
        //omp_set_num_threads(4); for profiling
        for (p = 0; p < pc.size(); ++p) {
            pc[p].location = verts[p].cast<double>();
            pc[p].colour = cols[p].cast<int>();
        }
        return pc;

    }
    catch (const std::exception& e)
    {
        std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
    }
}


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
            }
            else {
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
    if (token == "vertex") { // reading the header
        size_t nVertices;
        fin >> nVertices; // read number of points
        pc.resize(nVertices, defaultPoint);
        fin >> token;
        int pIx = 0;
        while (token == "property") { //read order of point array
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
    }
    else {
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
        int ix = 0; //how many properties read of this point
        while (ss >> d) {
            if (ix == xIx) {
                pc[p].location(0) = d;
            }
            else if (ix == yIx) {
                pc[p].location(1) = d;
            }
            else if (ix == zIx) {
                pc[p].location(2) = d;
            }
            else if (ix == rIx) {
                pc[p].colour(0) = d;
            }
            else if (ix == gIx) {
                pc[p].colour(1) = d;
            }
            else if (ix == bIx) {
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
    fout << "ply\n" //write the header
        << "format ascii 1.0\n"
        << "element vertex " << pc.size() << "\n"
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
