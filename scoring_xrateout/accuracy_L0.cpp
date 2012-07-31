#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <iomanip>
using namespace std;

template <class T>
bool from_string(T& t, 
                 const std::string& s, 
                 std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}
struct stkline {
	string label;
	string sequence;
};

stkline readstkline (string line) {
	int position=0;
	while (line[position] != ' ') {
		position++;
		if (position >= line.size()) {
			cerr << line << " is not a proper stk line." << endl;
			exit(1);
		}
	}
	stkline result;
	result.label = line.substr(0,position);
	while (line[position] == ' ') {
		position++;
	}
	result.sequence=line.substr(position,line.size()-position);
	return result;
}

bool is_section(char letter) {
	string sectionletters="egikxX";
	if (sectionletters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

bool is_exonsection(char letter) {
	string exonsectionletters="eg";
	if (exonsectionletters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

bool is_nonexonsection(char letter) {
	string nonexonsectionletters="ikxX";
	if (nonexonsectionletters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

int main (int argc, char *argv[]) {
	if (argc < 4) {
		cerr << "Please execute like " << argv[0] << " nameOfTarget xrateout genout" << endl;
		exit(1);
	}
	string line;
	
	string nameOfTarget=argv[1];
	
	ifstream inFile;
	inFile.open(argv[2]);
	if (!inFile) {
		cerr << "Unable to open xrate output" << argv[2] << endl;
		exit(1);
	}
	stkline xrateTarget;
	while (getline(inFile,line)) {
		xrateTarget = readstkline(line);
		if (xrateTarget.label == nameOfTarget) {
			break;
		}
	}
	inFile.close();
	
	inFile.open(argv[3]);
	if (!inFile) {
		cerr << "Unable to open genout file " << argv[3] << endl;
	}
	stkline genoutTarget;
	while (getline(inFile,line)) {
		genoutTarget = readstkline(line);
		if (genoutTarget.label == nameOfTarget) {
			break;
		}
	}
	inFile.close();
	
	int completeagreements=0;
	int exonagreements=0;
	int sections=0;
	for (int i=0; i<xrateTarget.sequence.size(); i++) {
		if (is_section(xrateTarget.sequence[i])) {
			sections++;
			if (is_exonsection(xrateTarget.sequence[i]) && is_exonsection(genoutTarget.sequence[i])) {
				exonagreements++;
			}
			if (is_nonexonsection(xrateTarget.sequence[i]) && is_nonexonsection(genoutTarget.sequence[i])) {
				exonagreements++;
			}
			if (xrateTarget.sequence[i] == genoutTarget.sequence[i]) {
				completeagreements++;
			}
		}
	}
	cout << (float)exonagreements/sections << " " << (float)completeagreements/sections << endl;
}