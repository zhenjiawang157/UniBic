
=================
Usage
==================
This software provides a biclustering module for microarray data. For a set of genes and a set of conditions, the program outputs a block-like structure which shows uniform trending-preserving pattern within the block, the block would contain only subsets of all given genes under subsets of all given conditions. 

Certain parts of the code uses open-source data structure library codes, including:
- fib <http://resnet.uoregon.edu/~gurney_j/jmpc/fib.html>, copyright information in fib.c
- Mark A. Weiss's data structure codes <http://www.cs.fiu.edu/~weiss/>

==================
Installation
==================
Simply put "unibic1.0.tar.gz" in any directory, 

$ tar zxvf unibic1.0.tar.gz

enter the folder "unibic1.0" and type "make" then the compiled codes are within the same directory as the source.

==================
Inputs and outputs
==================
The major program in the provided package is `unibic`, it can parse two 
formats of files, discrete data and continuous data, and examples for each
are provided. See help and look at all available options.

	$ ./unibic -h

Take a look at `toy_example` (discrete data) first. And try to run clustering 

	$ ./unibic -i toy_example -d

-d is important here since it tells the program that this is discrete data.

Then look at a continuous data "example". Try to run

	$ ./unibic -i example -f .25

This restricts no two blocks overlap more than 0.25 of the size of each one. And the other parameters are default value.

For each input file, our program generates three output files, namely,'.blocks' file, '.chars'file and '.rules' file.

In '.blocks' file, you can see all the biclusters the program found, especially, we use a blank line to separate the positively and the negatively (if any) correlated genes in each bicluster.

As to '.chars' file, it provides the qualitative matrix of the microarray data to usrs with some details of how to discrete the data in '.rules' file. You can find further details about how to represent a microarray dataset with a qualitative matrix in our paper.

==================
Changelog
==================
Version 1.0 
- latest version

==================
Contact
==================
Any questions, problems, bugs are welcome and should be dumped to
Zhenjia Wang<zhenjia.sdu@gmail.com>

Creation: Dec. 22, 2014

#--------------------------------------------------------------------------#

