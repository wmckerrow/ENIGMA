/*
 This code adds reading from to the oldL0 grammar containing only forward letters.
 It takes as input a .stk file in the eix alphabet.
 It reads through the stk input creating a new stk as it goes.
 1. It reads the stk into memory.
 2. It copies the line labels to the new stk.
 3. It goes through each sequence column by column.
	a. If every character in the column is x. Choose a random strand.
	b. If the forward strand is chosen for this column add the column to the new stk.
	c. If the reverse strand is chosed for thsi column add the reverse letter to the new stk.
 4. Once the entire stk input has been read the new stk is printed to stdout.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

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

int main (int argc, char *argv[]) {
	srand(time(NULL));
	if (argc < 2) {
		cerr << "Please execute like " << argv[0] << " unstrandedStkFile" << endl;
		exit(1);
	}
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr  << "Unable to open unstranded stkfile " << argv[1] << endl;
		exit(1);
	}
	string line;
	int numstklines=0;
	while (getline(inFile,line)) {
		if (line[0] != '#') {
			numstklines++;
		}
	}
	inFile.close();
	inFile.open(argv[1]);
	stkline unstrandedStk[numstklines];
	for (int i=0; i<numstklines; i++) {
		getline(inFile,line);
		while (line[0] == '#') {
			getline(inFile,line);
		}
		unstrandedStk[i] = readstkline(line);
	}
	inFile.close();
	stkline strandedStk[numstklines];
	for (int i=0; i<numstklines; i++) {
		strandedStk[i].label = unstrandedStk[i].label;
		strandedStk[i].sequence = "";
	}

	int strand=rand()%2;
	
	for (int j=0; j<unstrandedStk[0].sequence.size(); j++) {
		bool allx=1;
		for (int i=0; i<numstklines; i++) {
			if (unstrandedStk[i].sequence[j] != 'x') {
				allx=0;
				break;
			}
		}
		if (allx) {
			strand=rand()%2;
		}
		if (strand == 1) {
			for (int i=0; i<numstklines; i++) {
				strandedStk[i].sequence += unstrandedStk[i].sequence[j];
			}
		}
		else {
			for (int i=0; i<numstklines; i++) {
				if (unstrandedStk[i].sequence[j] == 'e') {
					strandedStk[i].sequence += 'g';
				}
				if (unstrandedStk[i].sequence[j] == '1') {
					strandedStk[i].sequence += '9';
				}
				if (unstrandedStk[i].sequence[j] == '2') {
					strandedStk[i].sequence += '8';
				}
				if (unstrandedStk[i].sequence[j] == 'i') {
					strandedStk[i].sequence += 'k';
				}
				if (unstrandedStk[i].sequence[j] == 'q') {
					strandedStk[i].sequence += 'w';
				}
				if (unstrandedStk[i].sequence[j] == 'r') {
					strandedStk[i].sequence += 'z';
				}
				if (unstrandedStk[i].sequence[j] == 'x') {
					strandedStk[i].sequence += 'x';
				}
				if (unstrandedStk[i].sequence[j] == 's') {
					strandedStk[i].sequence += 'u';
				}
				if (unstrandedStk[i].sequence[j] == 't') {
					strandedStk[i].sequence += 'v';
				}
				if (unstrandedStk[i].sequence[j] == 'a') {
					strandedStk[i].sequence += 'b';
				}
				if (unstrandedStk[i].sequence[j] == 'd') {
					strandedStk[i].sequence += 'c';
				}
				if (unstrandedStk[i].sequence[j] == 'f') {
					strandedStk[i].sequence += 'h';
				}
				if (unstrandedStk[i].sequence[j] == 'j') {
					strandedStk[i].sequence += 'l';
				}
				if (unstrandedStk[i].sequence[j] == 'y') {
					strandedStk[i].sequence += 'y';
				}
				if (unstrandedStk[i].sequence[j] == '*') {
					strandedStk[i].sequence += '*';
				}
			}
		}
	}
	
	int maxlabellength=0;
	for (int i=0; i<numstklines; i++) {
		if (strandedStk[i].label.size() > maxlabellength) {
			maxlabellength = strandedStk[i].label.size();
		}
	}
	int labelcolumns = max(maxlabellength+3,10);
	for (int i=0; i<numstklines; i++) {
		cout << strandedStk[i].label;
		for (int j=0; j<labelcolumns-strandedStk[i].label.size(); j++) {
			cout << " ";
		}
		cout << strandedStk[i].sequence << endl;
	}
}