/*
 This code pull out information of interest from the output of evaluate_gtf.pl (part of eval). 
 It adds a line with the specified information.
 
 The command line arguments are:
 eval output file, 
 file to add to, 
 line where the first piece of information, 
 column (tab, space delimited) where that information is, 
 rows and columns for the rest of the information to be added. 
 
 For example: 
 parseeval evalout_3 exonsensitivity_3 79 2 79 3 11 3 12 3 > garbage
 */

#define MAXCHILDREN 10
#define MAXDEPTH 10
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

template <class T>
bool from_string(T& t, 
                 const std::string& s, 
                 std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

int main(int argc, char *argv[]) {
	//Check command line
	if (argc < 5 || argc % 2 != 1) {
		cerr << "Please execute like " << argv[0] << " EvalOutputFile DataFile Linenumber1 ColumnNumber1 LineNumber2 ColumnNumber2 ..." << endl;
		exit(1);
	}
	
	//Set variables to hold data
	int numarguments=(argc-3)/2;
	string line;
	string target[numarguments];
	float data[numarguments];
	
	//Open file and find entries
	ifstream inFile;
	for (int j=0; j<numarguments; j++) {
		inFile.open(argv[1]);
		if (!inFile) {
			cerr << "Unable to open file " << argv[1] << endl;
			exit(1);
		}
		for (int i=1; i<atoi(argv[3+2*j]); i++) {
			getline(inFile,line);
			cout << line << endl;
		}
		for (int i=0; i<atoi(argv[4+2*j]); i++) {
			inFile >> target[j];
		}
		inFile.close();
	}

	//Turn target into a float and divide by 100 if it is a percent.
	for (int i=0; i<numarguments; i++) {
		if (target[i][target[i].size()-1]=='%') {
			target[i]=target[i].substr(0,target[i].size()-1);
			if(!from_string<float>(data[i],target[i],std::dec)) {
				cerr << "Failed to convert " << target[i] << endl;
			}
			data[i]/=100.0;
		}
		else {
			if(!from_string<float>(data[i],target[i],std::dec)) {
				cerr << "Failed to convert " << target[i] << endl;
			}
		}
	}	
		
	//Add data to DataFile
	ofstream outFile;
	inFile.open(argv[2]);
	int numlines=0;
	if (inFile) {
		while (getline(inFile,line)) {
			numlines++;
		}
	}
	inFile.close();
	string DataFile[numlines];
	if (numlines>0) {
		inFile.open(argv[2]);
		for (int i=0; i<numlines; i++) {
			getline(inFile,DataFile[i]);
		}
		inFile.close();
	}
	outFile.open(argv[2]);
	if (numlines>0) {
		for (int i=0; i<numlines; i++) {
			outFile << DataFile[i] << endl;
		}
	}
	for (int i=0; i<numarguments; i++) {
		outFile<< data[i] << " ";
	}
	outFile << endl;
}