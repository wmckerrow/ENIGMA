#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 3) {
		cerr << "Please execute like " << " old.stk newtree" << endl;
		exit(1);
	}
	string newtree = argv[2];
	cout << newtree << endl;
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr << "Unable to open stk file " << argv[1] << endl;
		exit(1);
	}
	string line;
	getline(inFile,line);
	while (getline(inFile,line)) {
		cout << line << endl;
	}
}