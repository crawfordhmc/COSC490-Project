# COSC490-Project
Finding geometric planes in dense point cloud models can give a useful and concise geometric representation of the object being modelled. 
In this project, accurate and computationally efficient methods of detecting these planes are explored. 
A parallel implementation of Random Sample and Consensus (RANSAC) has been developed and optimized for an 8.5x speedup, scaling up to 12 threads across large and small point clouds. 
Uniform space subdivision is used for an additional 6x serial speedup over plane detection on an undivided point cloud. 
The plane detector presented is capable of reading from file andclassifying 52 million points in just under a minute, with performance benefits on smaller models and thread-limited hardware as well.


This project uses the Eigen library, and the tinyply library to load in larger non-ASCII PLY files.
MeshLab is recommended for viewing point cloud inputs and outputs.

The master branch was directly overwritten with the final code due to merge issues with IDE files - development history and the Oct-tree skeleton class can be viewed in the "lab" branch.

Usage: planeFinder \input file\ \output file\ \probability of success\ \ratio of scene to be explained by planes\ \max RANSAC trials\ \scale factor\
optional: \"uniform" or "vector" space subdivision\
