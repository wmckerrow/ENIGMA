#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 5) {
		cerr << "Please execute like " << argv[0] << " expandedgrammar.eg forwardName reverseName NodesToReverse" << endl;
		exit(1);
	}
	string forwardName = argv[2];
	string reverseName = argv[3];
	string line;
	
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr << "Unable to open expanded grammar " << argv[1] << endl;
		exit(1);
	}
	
	while (getline(inFile,line)) {
		string newline = line;
		for (int i=4; i<argc; i++) {
			string temp = argv[i];
			string findthis = "(label " + temp + ")";
			if (line.find(findthis) != string::npos) {
				size_t forwardNamePos = line.find(forwardName);
				if (forwardNamePos == string::npos) {
					cerr << "Can't find " << forwardName << " in " << line << endl;
				}
				else {
					newline = line.substr(0,forwardNamePos) + reverseName + line.substr(forwardNamePos + forwardName.size(), line.size()-(forwardNamePos + forwardName.size()) );
				}
			}
		}
		cout << newline << endl;
	}
}