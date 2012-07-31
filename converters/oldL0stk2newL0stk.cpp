/* Reletters from old L0 alphabet to new one.
 For each line for each character if we haven't seen a space yet in this line we are still in the label, so we output the character exactly.
 If we have seen a space output the equivalent newL0 letter.
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
						case 'f':
							cout << 'p';
							break;
						case 'j':
							cout << 'u';
							break;
						case 'g':
							cout << 'E';
							break;
						case 'k':
							cout << 'I';
							break;
						case 'u':
							cout << 'T';
							break;
						case 'v':
							cout << 'S';
							break;
						case 'b':
							cout << 'D';
							break;
						case 'c':
							cout << 'A';
							break;
						case 'h':
							cout << 'P';
							break;
						case 'l':
							cout << 'U';
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