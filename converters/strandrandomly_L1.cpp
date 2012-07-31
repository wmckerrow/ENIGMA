/*
 This code adds reading from to the L1 grammar containing only forward letters.
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
					strandedStk[i].sequence += 'E';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'f') {
					strandedStk[i].sequence += 'F';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'g') {
					strandedStk[i].sequence += 'G';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'i') {
					strandedStk[i].sequence += 'I';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'j') {
					strandedStk[i].sequence += 'J';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'k') {
					strandedStk[i].sequence += 'K';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'x') {
					strandedStk[i].sequence += 'x';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 's') {
					strandedStk[i].sequence += 'T';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 't') {
					strandedStk[i].sequence += 'S';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'a') {
					strandedStk[i].sequence += 'D';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'b') {
					strandedStk[i].sequence += 'O';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'c') {
					strandedStk[i].sequence += 'N';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'd') {
					strandedStk[i].sequence += 'A';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'o') {
					strandedStk[i].sequence += 'B';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'n') {
					strandedStk[i].sequence += 'C';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'p') {
					strandedStk[i].sequence += 'P';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'q') {
					strandedStk[i].sequence += 'Q';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'r') {
					strandedStk[i].sequence += 'R';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'u') {
					strandedStk[i].sequence += 'U';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'v') {
					strandedStk[i].sequence += 'V';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'w') {
					strandedStk[i].sequence += 'W';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == 'y') {
					strandedStk[i].sequence += 'y';
					continue;
				}
				if (unstrandedStk[i].sequence[j] == '*') {
					strandedStk[i].sequence += '*';
					continue;
				}
				cerr << "No rule for " << unstrandedStk[i].sequence[j] << endl;
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