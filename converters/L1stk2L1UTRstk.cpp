// To convert from L1 to L1utr we only need to replace y with X.
// The code also replaces X with x since in older iterations X was also used for x.

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Please execute like " << argv[0] << " L1stk" << endl;
	}
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr << "Unable to open file " << argv[1] << endl;
		exit(1);
	}
	string line;
	bool space;
	while (getline(inFile,line)) {
		if (line[0]!='#') {
			space=0;
			for (int i=0; i<line.size(); i++) {
				if (line[i]==' ') {
					space=1;
				}
				if (space) {
					switch (line[i]) {
						case 'X':
							cout << 'x';
							break;
						case 'y':
							cout << 'X';
							break;
						default:
							cout << line[i];
							break;
					}
				}
				else {
					cout << line[i];
				}

			}
			cout << endl;
		}
		else {
			cout << line << endl;
		}

	}
	return 0;
}