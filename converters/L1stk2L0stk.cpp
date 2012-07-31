/* Removes reading from a L1 stk converting it to L0
 For each line for each character if we haven't seen a space yet in this line we are still in the label, so we output the character exactly.
 If we have seen a space output the equivalent L0 letter.
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
				if (space) {
					switch (line[i]) {
						case 'f':
							cout << 'e';
							break;
						case 'g':
							cout << 'e';
							break;
						case 'j':
							cout << 'i';
							break;
						case 'k':
							cout << 'i';
							break;
						case 'b':
							cout << 'a';
							break;
						case 'c':
							cout << 'a';
							break;
						case 'o':
							cout << 'd';
							break;
						case 'n':
							cout << 'd';
							break;
						case 'q':
							cout << 'p';
							break;
						case 'r':
							cout << 'p';
							break;
						case 'v':
							cout << 'u';
							break;
						case 'w':
							cout << 'u';
							break;
						case 'F':
							cout << 'E';
							break;
						case 'G':
							cout << 'E';
							break;
						case 'J':
							cout << 'I';
							break;
						case 'K':
							cout << 'I';
							break;
						case 'B':
							cout << 'A';
							break;
						case 'C':
							cout << 'A';
							break;
						case 'O':
							cout << 'D';
							break;
						case 'N':
							cout << 'D';
							break;
						case 'Q':
							cout << 'P';
							break;
						case 'R':
							cout << 'P';
							break;
						case 'V':
							cout << 'U';
							break;
						case 'W':
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