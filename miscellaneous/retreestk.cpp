/*
 This code removes the tree from an stk and replaces it with a new one.
 The command line arguments are: the old stk, the new tree. For example:
 retreestk old.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1)ABCD:.1)AB:.1)A;"
 */

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