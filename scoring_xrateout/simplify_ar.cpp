/* xrate puts ancestral reconstruction results in comments.
 This code creates an stk file containing only the ancestral reconstruction info.
 Use this code before any of the accuracy codes.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Please execute like " << argv[0] << " ar.stk" << endl;
	}
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr << "Unable to open ancestral reconstruction " << argv[1] << endl;
		exit(1);
	}
	string line;
	while (getline(inFile,line)) {
		if (line[0] != '#') {
			if (!getline(inFile,line)) {
				break;
			}
			int endlabel=5;
			while (line[endlabel] != ' ') {
				endlabel++;
			}
			string label=line.substr(5,endlabel-5);
			int lastspace=line.size()-1;
			while (line[lastspace] != ' ') {
				lastspace--;
			}
			cout << label;
			for (int i=0; i<lastspace-label.size(); i++) {
				cout << " ";
			}
			cout << line.substr(lastspace+1,line.size()-lastspace-1) << endl;
		}
	}
}