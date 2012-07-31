/* Reletters from newL0 alphabet to old one
 For each line for each character if we haven't seen a space yet in this line we are still in the label, so we output the character exactly.
 If we have seen a space output the equivalent oldL0 letter.
 */

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
				if (space == 1) {
					switch (line[i]) {
						case 'p':
							cout << 'f';
							break;
						case 'u':
							cout << 'j';
							break;
						case 'E':
							cout << 'g';
							break;
						case 'I':
							cout << 'k';
							break;
						case 'T':
							cout << 'u';
							break;
						case 'S':
							cout << 'v';
							break;
						case 'D':
							cout << 'b';
							break;
						case 'A':
							cout << 'c';
							break;
						case 'P':
							cout << 'h';
							break;
						case 'U':
							cout << 'l';
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