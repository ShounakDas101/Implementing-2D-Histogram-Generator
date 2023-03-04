# Implementing-2D-Histogram-Generator
Title	:	 MadAnalysis5 - Implementing 2D Histogram Generator

Author	: 	Shounak Das

Date	:	04-03-2023

# Problem Statement and Implementation in the code
1. write a C++ code using the attached data....INCORPORATED IN THE CODE

2. Each data line is a four-vector where the first three columns are the momentum vector's x,y, and z directions, and the last is energy....INCORPORATED IN THE CODE

3. All data is given in the units of MeV....NOT INCORPORATED IN THE CODE 

4. The code should contain a reader, a data class and an output system....INCORPORATED IN THE CODE

5. The reader class should be able to read the data file and prepare the appropriate data class....INCORPORATED IN THE CODE
    
6. The data class should have two layers....INCORPORATED IN THE CODE

7. The first layer is the vector which takes x,y, and z input.... INCORPORATED IN THE CODE

8. This vector class should be able to compute basic vector computation, such as extracting the transverse axis or calculating the angle between two axes depending on the requested output....INCORPORATED IN THE CODE:

The transverse axis of a 3D vector is a line perpendicular to the vector and lying in the plane formed by the vector and a reference plane. It is 	also sometimes referred to as the "normal" or "orthogonal" axis. In other words, if we have a 3D vector with components (x, y, z), then the 		transverse axis would be the line perpendicular to this vector and lying in the plane formed by the vector and a chosen reference plane.

The transverse axis of a 3D vector is the axis perpendicular to both the vector and a given reference axis.For example, suppose we have a 3D vector 	v = (1, 2, 3) and a reference axis r = (0, 0, 1). To find the transverse axis of v with respect to r, we can use the cross product:t = v × rwhere × 	denotes the cross product. Using the cross product formula, we get:t = (2 * 1 - 3 * 0, 3 * 0 - 1 * 0, 1 * 0 - 2 * 0) = (2, 0, 0)Therefore, the 		transverse axis of v with respect to r is the x-axis, represented by the vector (1, 0, 0).

9. The second layer is the FourVector class which should inherit the vector class and add additional capabilities. For instance, FourVector class can access energy to compute the given vector's mass where the vector represents a particle....INCORPORATED IN THE CODE
    
10. The output system should generate 1D histogram data for a given distribution....INCORPORATED IN THE CODE

11. The histogram's x-axis corresponds to the bins of the requested distribution, and the y-axis is the probability of each bin.... INCORPORATED IN THE CODE, BUT COMMENTED, IF required please uncomment Line NO 242 to 252. and recomplie and run. To check the result from data file, it is commented.

12. The output should be a two-column text file that one can use to draw a histogram later on, where the first column is the x-axis and the second column is the y-axis of the histogram....INCORPORATED IN THE CODE
    
13. The main program should be able to take the following arguments: datafile, distribution name, number of bins, min value for the x-axis, and max value for the x-axis. INCORPORATED IN THE CODE

14. The code should be able to prepare the following distributions: pT, pX, pY, pZ, energy, and mass.INCORPORATED IN THE CODE

15. If one asks, any other distribution program should not crash but give an error saying that the requested distribution is unavailable. INCORPORATED IN THE CODE



Program Code :

    Histo.cpp

Compilation:

	$ g++ Histo.cpp -o main.o

Execution:

	$ ./main.o input.dat distribution_index num_bins X_min X_max
  
	distribution_index (integer value)
	-----------------------
	 0   for     px
	 1   for     py
	 2   for     pz
	 3   for     E
	 4   for     pT
	 5   for     M         \
	-----------------------
	num_bins   ( integer value)
	X_min       (double value)
	X_max       (double value)
	
	example 
  
  	$ ./main.o input.dat 1 100 -2000.87 10000.9234

Output:
	output file name: output.dat

# Program Code Structure
Vector3 and FourVector are the data classes.

Vector3 is a class that represents a 3-dimensional vector. 

    It has the following public member variables:
    
    double x: represents the x-component of the vector.
    
    double y: represents the y-component of the vector.
    
    double z: represents the z-component of the vector.

Vector operations are listed below:

	Vector3(): a default constructor that initializes all components to 0.

	Vector3(double _x, double _y, double _z): a constructor that takes in the x, y, and z components as arguments and initializes the corresponding member variables.

	Vector3(const Vector3& other): a copy constructor that initializes the new vector with the same components as another vector.

	Vector3& operator=(const Vector3& other): an assignment operator that assigns the components of another vector to the current vector.

	Vector3 operator+(const Vector3& other) const: a binary operator that returns the sum of the current vector and another vector.

	Vector3 operator-(const Vector3& other) const: a binary operator that returns the difference between the current vector and another vector.

	Vector3 operator*(double scalar) const: a binary operator that returns the current vector multiplied by a scalar value.
	    double dot(const Vector3& other) const: a function that returns the dot product between the current vector and another vector.

	Vector3 cross(const Vector3& other) const: a function that returns the cross product between the current vector and another vector.
	    double length() const: a function that returns the length of the current vector.

	Vector3 normalize() const: a function that returns the normalized vector of the current vector.

	Vector3 transverse() const: a function that returns the transverse vector of the current vector.
	    double angle(const Vector3& other) const: a function that returns the angle between the current vector and another vector.

FourVector is a class named FourVector which is publicly inherited from another class Vector3. 

	The FourVector class has a constructor that takes four double inputs representing momentum components (px, py, pz) 
	and energy (E), and initializes the Vector3 part of the object with the momentum components and the E_ data member with the energy component.

Reader class is for the reading input.dat file. It has a public constructor that takes a string filename as input and initializes the filename_ data member with it.

The class also has two public member functions: read3() and read4().

read3() takes a pointer to a System object as input and reads three momentum components from the filename_ file corresponding to each particle in the system, creates a ThreeVector object for each set of momentum components, and adds it to the System object. 
The function returns the number of particles read.

read4() takes a pointer to a System object as input and reads four momentum components from the filename_ file corresponding to each particle in the system, creates a FourVector object for each set of momentum components,and Energy. The function returns the vector of FourVector objects.

The Reader class has a private data member filename_ which holds the name of the file being read.

After the Reader class definition, there is a global function System::read() which takes a character pointer filename as input. The function creates a Reader object using the filename input, reads three momentum components using the read3() member function of the Reader object and adds them to the System object using the this pointer. It then reads momentum components using the read4() member function of the Reader object, assigns the vector of FourVector objects returned to the fourvectors_ data member of the System object using the this pointer.

The System class contains methods for manipulating Vector3 and FourVector objects, and generating histograms.

Member variables:

    std::vector<Vector3> vectors_: a vector of Vector3 objects.
    std::vector<FourVector>

Methods:
    
    System(const std::vector<Vector3>& vectors, const std::vector<FourVector>& fourvectors): a constructor that initializes the vectors_ and fourvectors_ member variables with Vector3 and FourVector objects, respectively.
    const std::vector<Vector3>& vectors() const: a getter method for vectors_.
    const std::vector<FourVector>& fourvectors() const: a getter method for fourvectors_.
    void push3(Vector3& data): appends a Vector3 object to the vectors_ vector.
    void push4(FourVector& data): appends a FourVector object to the fourvectors_ vector.
    void read(char* filename): reads in data from a file.
    void histogram(const std::vector<double>& data, double num_bins, double x_min, double x_max): generates a histogram of the input data.
    std::vector<double> distribution(int index): returns a vector of values from either the x, y, or z component of a Vector3 object, the energy of a FourVector object, the magnitude of the spatial components of a FourVector object, or the spatial component of a FourVector object normalized by its energy.
    int save_histo(double num_bins, double x_min, double x_max): saves histogram data to a file.

# Further Improvements in the code
1. Some of the member variable are in Public, better to put them in Private
2. if any other Distribution should be added into the distribution function, therefore histogram function will be same far all distribution.
