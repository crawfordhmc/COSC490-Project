# COSC490-Project
Finding geometric planes in dense point cloud models can give a more useful and compact representation of the object being modelled. 
In this project accurate and computationally efficient methods of detecting these planes are explored. 
A parallel implementation of RANSAC with 3.5x speedup over 4 threads has been created, with reasonable classification of noise, but too few planes detected to completely represent the structure of dense models. 
I aim to further improve the accuracy of this process, and to increase performance with other techniques beyond parallelization.
