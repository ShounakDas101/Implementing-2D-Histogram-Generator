/*
Title	:	 MadAnalysis5 - Implementing 2D Histogram Generator
Author	: 	Shounak Das 
Date	:	04-03-2023

Program Code : Histo.cpp

Compilation:
	g++ Histo.cpp -o main.o
Execution:
	./main.o input.dat distribution_index num_bins X_min X_max
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
	
	example: ./main.o input.dat 1 100 -2000.87 10000.9234

*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

/*
    This is a class called Vector3 that represents a 3-dimensional vector. 
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

*/
class Vector3 {
public:
    double x, y, z;

    Vector3() : x(0.0), y(0.0), z(0.0) {}//a default constructor that initializes all components to 0.

    Vector3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}//a constructor that takes in the x, y, and z components.

    Vector3(const Vector3& other) : x(other.x), y(other.y), z(other.z) {}// a copy constructor

    Vector3& operator=(const Vector3& other) { //an assignment operator
        if (this != &other) {
            x = other.x;
            y = other.y;
            z = other.z;
        }
        return *this;
    }

    Vector3 operator+(const Vector3& other) const { // sum of the current vector and another vector.
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    Vector3 operator-(const Vector3& other) const { // difference between the current vector and another vector.
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    Vector3 operator*(double scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar); //vector multiplied by a scalar value.
    }

    double dot(const Vector3& other) const { // Dot Product
        return x * other.x + y * other.y + z * other.z;
    }

    Vector3 cross(const Vector3& other) const { // Cross Product
        return Vector3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }

    double length() const { // Length from origin
        return std::sqrt(x * x + y * y + z * z);
    }

    Vector3 normalize() const { //unit vector
        double len = length();
        if (len > 0.0) {
            return Vector3(x / len, y / len, z / len);
        } else {
            return Vector3();
        }
    }
    /*
    	The transverse axis of a 3D vector is a line perpendicular to the vector and lying in the plane formed by the vector and a reference plane. 
    	It is also sometimes referred to as the "normal" or "orthogonal" axis. 
    	In other words, if we have a 3D vector with components (x, y, z), then the transverse axis would be the line perpendicular to this vector 
    	and lying in the plane formed by the vector and a chosen reference plane.

	The transverse axis of a 3D vector is the axis perpendicular to both the vector and a given reference axis.For example, 
	suppose we have a 3D vector v = (1, 2, 3) and a reference axis r = (0, 0, 1). 
	To find the transverse axis of v with respect to r, we can use the cross product:t = v × rwhere × denotes the cross product. 
	Using the cross product formula, we get:t = (2 * 1 - 3 * 0, 3 * 0 - 1 * 0, 1 * 0 - 2 * 0) = (2, 0, 0).
	Therefore, the transverse axis of v with respect to r is the x-axis, represented by the vector (1, 0, 0).
    */
    Vector3 transverse(const Vector3& other) const { //finding transverse axis w.r.t a reference vector
    	Vector3 temp;
    	temp=this->cross(other);
    	return temp.normalize();
    }

    double angle(const Vector3& other) const { //angle between the current vector and another vector.
        double len1 = length();
        double len2 = other.length();
        if (len1 > 0.0 && len2 > 0.0) {
            double dot_prod = dot(other) / (len1 * len2);
            return std::acos(std::max(-1.0, std::min(dot_prod, 1.0)));
        } else {
            return 0.0;
        }
    }
};
/*
This is a class named FourVector which is publicly inherited from another class Vector3. 
The FourVector class has a constructor that takes four double inputs representing momentum components (px, py, pz) 
and energy (E), and initializes the Vector3 part of the object with the momentum components and the E_ data member with the energy component.

The class also has a public member function named 

mass() that calculates and returns the mass from energy and resultant mementum.

*/

class FourVector : public Vector3 {
public:
    // Constructor that takes energy, momentum x, y, and z inputs
    FourVector(double px, double py, double pz, double E)
        : Vector3(px, py, pz), E_(E) {}

    // Additional capabilities of the FourVector class
    double mass() const {
        return std::sqrt(length() * length()/(2*E_));
    }
    
    double E_;
private:
    
};

/*
The System class contains methods for manipulating Vector3 and FourVector objects, and generating histograms.
Member variables:

    std::vector<Vector3> vectors_: a vector of Vector3 objects.
    `std::vector<FourVector>

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

*/

class System {
public:
    // Constructor that takes vectors of Vector and FourVector objects as inputs
    /*
    System(const std::vector<Vector3>& vectors, const std::vector<FourVector>& fourvectors)
        : vectors_(vectors), fourvectors_(fourvectors) {}
	*/

    // Accessors for the vectors and four-vectors
    const std::vector<Vector3>& vectors() const { return vectors_; }
    const std::vector<FourVector>& fourvectors() const { return fourvectors_; }


    void push3(Vector3& data)
    {
        vectors_.push_back(data);
    }

    void push4(FourVector& data)
    {
        fourvectors_.push_back(data);
    }
    void read(char* filename);
    
	
void histogram(const std::vector<double>& data, double num_bins, double x_min, double x_max) {
    // Compute number of bins
    int bin_size = static_cast<int>((x_max - x_min) / num_bins);
    double val;
    
    // check bin_size
    cout<<endl;
    cout << "X_max :" <<x_max<<endl;
    cout << "X_min :" <<x_min<<endl;
    cout << "bin_size ..:" <<bin_size<<endl;
    
    for (std::vector<double>::const_iterator it = data.begin(); it != data.end(); ++it) {
    	val = *it;
    	// check the actual data is read
    	//cout<<val<<endl;
    }


    	
    // Initialize histogram to all zeros
    hist.resize(num_bins, 0.0);

    // Compute total number of data points
    int n = data.size();
	cout<<endl<<"data size :"<<n<<endl;
    
    // Create the bins
    for (int i = 0; i < num_bins; i++) {
    	this->bins.push_back(x_min +i*bin_size);
    	//cout<< this->bins[i] <<endl;
        }
		
    // Loop over data and fill in histogram
    for (int i = 0; i < n; i++) {
        double val = data[i];
        if (val < x_min || val >= x_max) {
            continue;
        }
        int bin = static_cast<int>((val - x_min) / bin_size);
        this->hist[bin]++;
        //cout<<bin<<endl;
    }

// To get probability Distribution, Uncomment following code 
/*
    // Normalize histogram by total counts and bin size
    double total_counts = 0.0;
    for (int i = 0; i < num_bins; i++) {
        total_counts += this->hist[i];
    }
    double norm_factor = bin_size * n;
    for (int i = 0; i < num_bins; i++) {
        this->hist[i] = this->hist[i] / total_counts;
    }
*/
// To get probability Distribution, Uncomment above block of codes   
 
// To print the histogram data on the screen, Uncomment the following codes	
    /*
    for (int i = 0; i < num_bins; i++) {
    	//Test code
        cout<<this->bins[i]<<" "<<this->hist[i] <<endl;
    }
    */
}


std::vector<double> distribution(int index)
{
	vector<double> data;
	double val;
	    if(index==0){    
	   
		for (std::vector<Vector3>::iterator it = vectors_.begin(); it != vectors_.end(); ++it) {
        	val = it->x;
        	data.push_back(val);
		}
	} 
    
	else if(index==1){    
	   
		for (std::vector<Vector3>::iterator it = vectors_.begin(); it != vectors_.end(); ++it) {
        	val = it->y;
        	//cout<<val<<endl;
		data.push_back(val);
		}
	}
	else if(index==2){    
	   
		for (std::vector<Vector3>::iterator it = vectors_.begin(); it != vectors_.end(); ++it) {
        	val = it->z;
		data.push_back(val);
		}
	}
	else if(index==3){    
	   
		for (std::vector<FourVector>::iterator it = fourvectors_.begin(); it != fourvectors_.end(); ++it) {
        	val = it->E_;
        	
		data.push_back(val);
		}
	}
	else if(index==4){    
	   
		for (std::vector<FourVector>::iterator it = fourvectors_.begin(); it != fourvectors_.end(); ++it) {
        	val = sqrt((it->x)*(it->x) + (it->y)*(it->y) +(it->z)*(it->z));
        	data.push_back(val);
		}
	}
	else if(index==5){    
	   
		for (std::vector<FourVector>::iterator it = fourvectors_.begin(); it != fourvectors_.end(); ++it) {
        	val = ((it->x)*(it->x) + (it->y)*(it->y) +(it->z)*(it->z))/(2*(it->E_)); 
        	data.push_back(val);
		}
	}
return data;
} 
	
	


int save_histo( double num_bins, double x_min, double x_max){
 // Generate some sample 2D vector data
 //   std::vector<std::vector<double>> data = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
 // Open output file for writing
 	//int num_bins = static_cast<int>((x_max - x_min) / bin_size);
    std::ofstream outfile("output.dat");
    if (!outfile.is_open()) {
        std::cerr << "Error: failed to open output file" << std::endl;
        return 1;
    }

    // Write data to output file with 2 columns
    for (int i = 0; i < num_bins; i++) {
        outfile << setprecision(12) << this->bins[i] << " " << this->hist[i] << std::endl;
    }

    // Close output file
    outfile.close();
    return 0;
}// End of save_histo() method

private:
    std::vector<Vector3> vectors_;
    std::vector<FourVector> fourvectors_;
    std::vector<double> hist;
    std::vector<int> bins;
}; //End of System Class


/*
This is a class named Reader. The class has a public constructor that takes a string filename as input and initializes the filename_ data member with it.

The class also has two public member functions: read3() and read4().

read3() takes a pointer to a System object as input and reads three momentum components from the filename_ file corresponding to each particle in the system, 
creates a ThreeVector object for each set of momentum components, and adds it to the System object. 
The function returns the number of particles read.

read4() takes a pointer to a System object as input and reads four momentum components from the filename_ file corresponding to each particle in the system, 
creates a FourVector object for each set of momentum components, 
and adds it to a vector of FourVector objects. The function returns the vector of FourVector objects.

The Reader class has a private data member filename_ which holds the name of the file being read.

After the Reader class definition, there is a global function System::read() which takes a character pointer filename as input. The function creates a Reader object using the filename input, reads three momentum components using the read3() member function of the Reader object and adds them to the System object using the this pointer. It then reads four momentum components using the read4() member function of the Reader object, assigns the vector of FourVector objects returned to the fourvectors_ data member of the System object using the this pointer.


*/
class Reader {
public:
    Reader(std::string filename) {
        filename_ = filename;
    }

    int read3(System* system);
    vector<FourVector> read4(System* system);	
private:
    std::string filename_;
};

// Define System::read function after Reader class definition
void System::read(char* filename)
{
    Reader my_read(filename);
    my_read.read3(this);
    this->fourvectors_=my_read.read4(this);
    
}



// Define Reader::read3 function after System class definition
int Reader::read3(System* system)
{
    std::ifstream infile(filename_);
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string x_str, y_str, z_str, w_str;

        if (!(iss >> x_str >> y_str >> z_str >> w_str)) {
            std::cerr << "Error reading line\n";
            return 0;
        }

        double x = std::stod(x_str);
        double y = std::stod(y_str);
        double z = std::stod(z_str);
        double w = std::stod(w_str);
	
        Vector3 V(x, y, z);
        system->push3(V);
        
    }

    infile.close();
    return 1;
}

// Define Reader::read3 function after System class definition
vector<FourVector> Reader::read4(System* system)
{
    std::ifstream infile(filename_);
    std::string line;
    vector<FourVector> Vlist;
    double x,y,z,w;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string x_str, y_str, z_str, w_str;

	
        if (!(iss >> x_str >> y_str >> z_str >> w_str)) {
            std::cerr << "Error reading line\n";
            
        }
	else
	{
		x = std::stod(x_str);
        	y = std::stod(y_str);
        	z = std::stod(z_str);
        	w = std::stod(w_str);
		
        	FourVector V(x, y, z, w);
       	 	Vlist.push_back(V);
	}
        
        
    }
    
    infile.close();
    return Vlist;
}


//////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
    cout << "You have entered " << argc
         << " arguments:" << endl;
  /*
  argv[0] = main.o
  argv[1] = filename 
  argv[2] = distribution: pT, pX, pY, pZ, E, and M.
  argv[3] = num_bins
  argv[4] = X_min
  argv[5] = X_max
  */
    for (int i = 0; i < argc; ++i)
        cout << argv[i] << endl;
    /* Run Example : $ main.o input.dat px 100 0.0 100.0 */ 
    char* filename = argv[1];
    int distribution = atoi(argv[2]);
    int num_bins = atoi(argv[3]); 
    float X_min = atof(argv[4]);
    float X_max = atof(argv[5]);
    // calculate number of bins
    int bin_size = (int)((X_max - X_min) / (num_bins - 1));
    /*
    distribution_index (integer value)
	-----------------------
	 0   for     px
	 1   for     py
	 2   for     pz
	 3   for     E
	 4   for     pT
	 5   for     M         \
	-----------------------
    distribution index other than 0-5, program should exit
    */
    if((distribution== 0)||(distribution==1)||(distribution==2)||(distribution==3)||(distribution==4)||(distribution==5))
    {
    	
    System Sys;// Instantiates a System object
    Sys.read(filename);// Read data fron input.dat file
    //Sys.distribution(distribution);
    Sys.histogram(Sys.distribution(distribution), num_bins, X_min, X_max);// Generates histogram data
    Sys.save_histo( num_bins, X_min, X_max);// Save the data in a output.dat file
    return 1;
    }
    else
    {
    cout<<"Additional Support needed, Ckeck Distribution Available"<<endl;
    	return 0;
    }
}// End of main()   


