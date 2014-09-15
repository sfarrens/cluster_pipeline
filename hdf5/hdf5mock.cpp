/*
	Code:			hdf5mock.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	13/09/2012
	Version: 		0.1
	Description: 	Code with utility functions to read and extract hdf5 mock files generated by Alex/Euclid.
					Makes use of hdf5 lib (http://www.hdfgroup.org/)
	Usage:			Include "hdf5mock.h" to use in other C/C++ codes.
					
	Compile (in VS): cl /MD /EHsc hdf5mock.cpp /I"C:\Program Files\HDF Group\HDF5\1.8.9\include" /I"C:\Program Files\HDF Group\HDF5\1.8.9\include\cpp" hdf_lib\hdf5_cppdll.lib hdf_lib\hdf5dll.lib
*/

#include "hdf5mock.h"
#include "H5Cpp.h"

std::string hdf5mock::getType(const std::string filename, const std::string column){
	H5::H5File file( filename, H5F_ACC_RDONLY );
	H5::DataSet dataset = file.openDataSet( column );
	H5T_class_t type_class = dataset.getTypeClass();
	
	size_t size;
	H5::IntType intype;
	H5::FloatType floattype;
		
	switch(type_class){
		case H5T_INTEGER:
			intype = dataset.getIntType();
			size = intype.getSize()*8;
			if (size <= 32)
				return "int";
			else
				return "long long";
			break;
		case H5T_FLOAT:
			floattype = dataset.getFloatType();
			size = floattype.getSize()*8;
			if (size <= 32)
				return "float";
			else
				return "double";
			break;
	}
	return "unknown";
}

int hdf5mock::getNumLines(const H5::H5File file, const H5::DataSpace dataspace){
	hsize_t dims_out[1];
	int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
	int nlines = (int) dims_out[0];
	return nlines;
}

// read double
std::vector<double> hdf5mock::readColumnDouble(const std::string filename, const std::string column){
	H5::H5File file( filename, H5F_ACC_RDONLY );
	H5::DataSet dataset = file.openDataSet( column );
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	int nlines = hdf5mock::getNumLines(file, dataspace);

	hsize_t dims_out[1];
	dims_out[0] = nlines;
	int rank = 1;
	H5::DataSpace mspace1(rank, dims_out);
	
	double *data_out;
	data_out = new double[nlines];
		
	dataset.read( data_out, H5::PredType::NATIVE_DOUBLE, mspace1, dataspace );
	
	std::vector<double> data_vec (&data_out[0], &data_out[nlines]);
	
	delete [] data_out;
	
	return data_vec;
}

// read int
std::vector<int> hdf5mock::readColumnInt(const std::string filename, const std::string column){
	H5::H5File file( filename, H5F_ACC_RDONLY );
	H5::DataSet dataset = file.openDataSet( column );
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	int nlines = hdf5mock::getNumLines(file, dataspace);

	hsize_t dims_out[1];
	dims_out[0] = nlines;
	int rank = 1;
	H5::DataSpace mspace1(rank, dims_out);
	
	int *data_out;
	data_out = new int[nlines];
		
	dataset.read( data_out, H5::PredType::NATIVE_INT, mspace1, dataspace );
	
	std::vector<int> data_vec (&data_out[0], &data_out[nlines]);
	
	delete [] data_out;
	
	return data_vec;
}

// read float
std::vector<float> hdf5mock::readColumnFloat(const std::string filename, const std::string column){
	H5::H5File file( filename, H5F_ACC_RDONLY );
	H5::DataSet dataset = file.openDataSet( column );
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	int nlines = hdf5mock::getNumLines(file, dataspace);

	hsize_t dims_out[1];
	dims_out[0] = nlines;
	int rank = 1;
	H5::DataSpace mspace1(rank, dims_out);
	
	float *data_out;
	data_out = new float[nlines];
		
	dataset.read( data_out, H5::PredType::NATIVE_FLOAT, mspace1, dataspace );
	
	std::vector<float> data_vec (&data_out[0], &data_out[nlines]);
	
	delete [] data_out;
	
	return data_vec;
}

// read long long or int 64bit
std::vector<long long> hdf5mock::readColumnLongLong(const std::string filename, const std::string column){
	H5::H5File file( filename, H5F_ACC_RDONLY );
	H5::DataSet dataset = file.openDataSet( column );
	
	H5::DataSpace dataspace = dataset.getSpace();
	
	int nlines = hdf5mock::getNumLines(file, dataspace);

	hsize_t dims_out[1];
	dims_out[0] = nlines;
	
	
	
	int rank = 1;
	H5::DataSpace mspace1(rank, dims_out);
	
	long long *data_out;
	data_out = new long long[nlines];
	
	//std::cout << nlines << std::endl;
		
	dataset.read( data_out, H5::PredType::NATIVE_LLONG, mspace1, dataspace );
	
	std::vector<long long> data_vec (&data_out[0], &data_out[nlines]);
	
	
	
	delete [] data_out;
	
	
	
	return data_vec;
}





