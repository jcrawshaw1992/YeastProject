#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <stdlib.h>
// #include <iostream>
 
/*
 * A class to read data from a csv file.
 */
class CSVReader
{
	std::string fileName;
	std::string delimeter;
 
public:
	CSVReader(std::string filename, std::string delm = ",") :
			fileName(filename), delimeter(delm)
	{ }
 
	// Function to fetch data from a CSV File
	std::vector<std::vector<double > > getData();

	// std::vector<std::vector<double > > InterpolateData(std::vector<std::vector<double> > Vector);
};
 
/*
* Parses through csv file line by line and returns the data
* in vector of vector of strings.
*/

// std::vector<std::vector<double> > CSVReader::InterpolateData(std::vector<std::vector<double> > Vec)
// {
// 	std::vector<std::vector<double> > Interploated;

// 		for (int i = 0  ; i< Vec.size() ; ++i)
// 		{
// 			PRINT_VECTOR(Vec[i]);
// 			PRINT_VECTOR(Vec[i+1]);
// 			int Next = i+1;
// 			std::vector<double>  Upper = Vec[i];
// 			std::vector<double>  Lower = Vec[Next];
// 			std::vector<double>  MidPoint;
// 			MidPoint[0] =  Upper[0] +Lower[0];
// 			//  = Upper +Lower;
// 			Interploated.push_back(Vec[i]);
// 			// Interploated.push_back(MidPoint);



// 			// std::vector<std::vector<double> >::iterator sec = ++iter;
// 			// // PRINT_VECTOR( Vector(*sec));
// 			// PRINT_VECTOR( Vector(*iter));
// 			// 
// 			// 

// 			// double value = std::strtod(iter->c_str(), NULL);
// 			// DoubleVec.push_back(value);

// 			// // PRINT_2_VARIABLES(* iter, value);
			
// 		}

		


// return Interploated;

// }

std::vector<std::vector<double> > CSVReader::getData()
{
	std::ifstream file(fileName);
 
	std::vector<std::vector< double> > dataList;
 
	std::string line = "";
	// Iterate through each line and split the content using delimeter

	


	while (getline(file, line))
	{
		std::vector<std::string> vec;
		std::vector<double> DoubleVec;
		boost::algorithm::split(vec, line, boost::is_any_of(delimeter));

		//  for (std::vector<unsigned>::iterator iter = Neighbourhood.begin();
        // //          iter != Neighbourhood.end();
        // //          ++iter)
        //     {


		for (std::vector<std::string>::iterator iter = vec.begin();
				 iter != vec.end(); ++iter)
				{
					double value = std::strtod(iter->c_str(), NULL);
					DoubleVec.push_back(value);

					// PRINT_2_VARIABLES(* iter, value);
					
				}
		dataList.push_back(DoubleVec);
	}
	// Close the File
	file.close();

	// std::vector<std::vector<double> > Interploated = InterpolateData(dataList);
	
 
	return dataList;
}


// ------------------
// How to use this code
// ------------------
// int main()
// {
// 	// Creating an object of CSVWriter
// 	CSVReader reader("example.csv");
 
// 	// Get the data from CSV File
// 	std::vector<std::vector<double> > dataList = reader.getData();
 
// 	// Print the content of row by row on screen
// 	for(std::vector<double> vec : dataList)
// 	{
// 		for(double data : vec)
// 		{
// 			std::cout<<data << " , ";
// 		}
// 		std::cout<<std::endl;
// 	}
// 	return 0;
 
// }