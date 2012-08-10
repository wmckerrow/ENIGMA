/*
 This code changes all the values in the last column in a dot file to 2, getting rid of any information about how many species are involved in that alignment.
 This is so that we can test whether or not adding more aligned species actually helps.
 The only command line argument is the dot file to be edited. The new dot file is sent to standard out. For example:
 all2way SM_25.dot > SM_25_2way.dot
 */

#define NUMSPECIES 2
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

struct dotline {
	int starts[NUMSPECIES];
	int length;
	int directions[NUMSPECIES];
	int numinalignment;
};

dotline getDotLine (string line) { //read space for tab...
	//cout << "Getting dotline from " << line << endl;
	int tab=-1;
	int nexttab;
	string data[2*NUMSPECIES+2];
	for (int i=0; i<2*NUMSPECIES+2; i++) {
		nexttab = tab+1;
		while (line[nexttab] != ' ' && nexttab < line.size()) {
			nexttab++;
		}
		data[i] = line.substr(tab+1,nexttab-tab-1);
		tab=nexttab;
	}
	dotline row;
	
	for (int i=0; i<NUMSPECIES; i++) {
		if (!from_string<int>(row.starts[i],data[i],std::dec)) {
			cerr << "from_string failed on data[" << i << "]" << endl;
		}
	}
	if (!from_string<int>(row.length,data[NUMSPECIES],std::dec)) {
		cerr << "from_string failed on data[" << NUMSPECIES << "]" << endl;
	}
	for (int i=0; i<NUMSPECIES; i++) {
		if (data[NUMSPECIES+1+i][0] == '+') {
			row.directions[i] = 1;
		}
		else {
			row.directions[i]=-1;
		}
	}
	if (!from_string<int>(row.numinalignment,data[2*NUMSPECIES+1],std::dec)) {
		cerr << "from_string failed on data[" << 2*NUMSPECIES+1 << "]" << endl;
	}
	//cout << "Got dotline" << endl;
	return row;
}

int main (int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Execute like " << argv[0] << " dotfile" << endl;
		exit(1);
	}

	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr << "Unable to open file " << argv[1] << endl;
		exit(1);
	}
	string line;
	while (getline(inFile,line)) {
		dotline row=getDotLine(line);
		for (int i=0; i<NUMSPECIES; i++) {
			cout << row.starts[i] << " ";
		}
		cout << row.length << " ";
		for (int i=0; i<NUMSPECIES; i++) {
			if (row.directions[i] == 1) {
				cout << "+ ";
			}
			else {
				cout << "- ";
			}
		}
		cout << 2;
		cout << endl;
	}
	return 0;
}