/*
 This code adds transions to the L0 grammar containing only the letters e,i and x.
 It takes as input a .stk file in the eix alphabet.
 It reads through the stk input creating a new stk as it goes.
 1. It reads the stk into memory.
 2. It copies the line labels to the new stk.
 3. It goes through each sequence character by character.
	a. For the first character of a sequence it simply add the character to the sequence for the new stk.
	b. For each other character it checks the preceding character and uses that to decide what type of transition occured.
	c. It add the transition letter followed by the current character to the new stk.
 4. Once the entire stk input has been read the new stk is printed to stdout.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

//Store each stk line as two strings: one for the label and one for the sequence.
struct stkline {
	string label;
	string sequence;
};

//Turn single string containing an entire stk line into two strings: one for the label and one for the sequence.
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
	//Check input
	if (argc < 2) {
		cerr << "Please execute like " << argv[0] << " eixStkFile" << endl;
		exit(1);
	}
	
	//Make sure we can open the stk file and check how  many lines it has
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr  << "Unable to open eix alphabet stkfile " << argv[1] << endl;
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
	
	//Read stk file into memory
	inFile.open(argv[1]);
	stkline eixStk[numstklines];
	for (int i=0; i<numstklines; i++) {
		getline(inFile,line);
		while (line[0] == '#') {
			getline(inFile,line);
		}
		eixStk[i] = readstkline(line);
	}
	inFile.close();
	
	//Copy line labels and first letter of sequences into a new stk
	stkline stadnStk[numstklines];
	for (int i=0; i<numstklines; i++) {
		stadnStk[i].label = eixStk[i].label;
		stadnStk[i].sequence = eixStk[i].sequence[0];
	}
	
	//Check add right transition letter to the new stk
	for (int j=1; j<eixStk[0].sequence.size(); j++) {
		for (int i=0; i<numstklines; i++) {
			int added=0;
			if (eixStk[i].sequence[j-1] == 'e' && eixStk[i].sequence[j] == 'e') {
				stadnStk[i].sequence += 'f';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'e' && eixStk[i].sequence[j] == 'i') {
				stadnStk[i].sequence += 'd';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'e' && (eixStk[i].sequence[j] == 'x' || eixStk[i].sequence[j] == 'X')) {
				stadnStk[i].sequence += 't';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'i' && eixStk[i].sequence[j] == 'e') {
				stadnStk[i].sequence += 'a';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'i' && eixStk[i].sequence[j] == 'i') {
				stadnStk[i].sequence += 'j';
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'x' || eixStk[i].sequence[j-1] == 'X') && eixStk[i].sequence[j] == 'e') {
				stadnStk[i].sequence += 's';
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'x' || eixStk[i].sequence[j-1] == 'X') && (eixStk[i].sequence[j] == 'x' || eixStk[i].sequence[j] == 'X')) {
				stadnStk[i].sequence += 'y';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '*' || eixStk[i].sequence[j] == '*') {
				stadnStk[i].sequence += '*';
			added++;
			}
			if (added != 1) {
				stadnStk[i].sequence += '*';
			}
		}
		// add the letter that comes after the transition
		for (int i=0; i<numstklines; i++) {
			stadnStk[i].sequence += eixStk[i].sequence[j];
		}
	}
	
	//print the new stk
	int maxlabellength=0;
	for (int i=0; i<numstklines; i++) {
		if (stadnStk[i].label.size() > maxlabellength) {
			maxlabellength = stadnStk[i].label.size();
		}
	}
	int labelcolumns = max(maxlabellength+3,10);
	for (int i=0; i<numstklines; i++) {
		cout << stadnStk[i].label;
		for (int j=0; j<labelcolumns-stadnStk[i].label.size(); j++) {
			cout << " ";
		}
		cout << stadnStk[i].sequence << endl;
	}
}
