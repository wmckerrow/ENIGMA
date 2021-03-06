/*
 This code adds transions to the L1 grammar containing only the section letters.
 It takes as input a .stk file in the eix alphabet.
 It reads through the stk input creating a new stk as it goes.
 1. It reads the stk into memory.
 2. It copies the line labels to the new stk.
 3. It goes through each sequence character by character.
	a. For the first character of a sequence it simply add the character to the sequence for the new stk.
	b. For each other character it checks the preceding character and uses that to decide what type of transition occured.
	c. It adds the transition character followed by the current character to the new stk.
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
		cerr << "Please execute like " << argv[0] << " eixStkFile" << endl;
		exit(1);
	}
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr  << "Unable to open eix alphabet stkfile " << argv[1] << endl;
		exit(1);
	}
	string line;
	int numstklines=0;
	while (getline(inFile,line)) {
		if (line[0] != '#' && line.size()!=0) {
			numstklines++;
		}
	}
	inFile.close();
	inFile.open(argv[1]);
	stkline eixStk[numstklines];
	for (int i=0; i<numstklines; i++) {
		getline(inFile,line);
		while (line[0] == '#' || line.size()==0) {
			getline(inFile,line);
			cout << line << endl;
		}
		eixStk[i] = readstkline(line);
	}
	inFile.close();
	stkline stadnStk[numstklines];
	for (int i=0; i<numstklines; i++) {
		stadnStk[i].label = eixStk[i].label;
		stadnStk[i].sequence = eixStk[i].sequence[0];
	}
		
	for (int j=1; j<eixStk[0].sequence.size(); j++) {
		for (int i=0; i<numstklines; i++) {
			int added=0;
			if (eixStk[i].sequence[j-1] == 'x' && eixStk[i].sequence[j] == 'x') {
				stadnStk[i].sequence += 'X';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'e' && eixStk[i].sequence[j] == 'e') {
				stadnStk[i].sequence += 'p';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'f' && eixStk[i].sequence[j] == 'f') {
				stadnStk[i].sequence += 'q';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'g' && eixStk[i].sequence[j] == 'g') {
				stadnStk[i].sequence += 'r';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'i' && eixStk[i].sequence[j] == 'i') {
				stadnStk[i].sequence += 'u';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'j' && eixStk[i].sequence[j] == 'j') {
				stadnStk[i].sequence += 'v';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'k' && eixStk[i].sequence[j] == 'k') {
				stadnStk[i].sequence += 'w';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'E' && eixStk[i].sequence[j] == 'E') {
				stadnStk[i].sequence += 'P';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'F' && eixStk[i].sequence[j] == 'F') {
				stadnStk[i].sequence += 'Q';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'G' && eixStk[i].sequence[j] == 'G') {
				stadnStk[i].sequence += 'R';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'I' && eixStk[i].sequence[j] == 'I') {
				stadnStk[i].sequence += 'U';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'J' && eixStk[i].sequence[j] == 'J') {
				stadnStk[i].sequence += 'V';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'K' && eixStk[i].sequence[j] == 'K') {
				stadnStk[i].sequence += 'W';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '5' && eixStk[i].sequence[j] == '5') {
				stadnStk[i].sequence += '6';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'h' && eixStk[i].sequence[j] == 'h') {
				stadnStk[i].sequence += 'l';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '3' && eixStk[i].sequence[j] == '3') {
				stadnStk[i].sequence += '4';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'm' && eixStk[i].sequence[j] == 'm') {
				stadnStk[i].sequence += 'y';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '&' && eixStk[i].sequence[j] == '&') {
				stadnStk[i].sequence += '\'';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'H' && eixStk[i].sequence[j] == 'H') {
				stadnStk[i].sequence += 'L';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '$' && eixStk[i].sequence[j] == '$') {
				stadnStk[i].sequence += '%';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'M' && eixStk[i].sequence[j] == 'M') {
				stadnStk[i].sequence += 'Y';
				added++;
			}
			
			if (eixStk[i].sequence[j-1] == 'x' && eixStk[i].sequence[j] == 'e') {
				stadnStk[i].sequence += 's';;
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'e' || eixStk[i].sequence[j-1] == 'f' || eixStk[i].sequence[j-1] == 'g') && eixStk[i].sequence[j] == 'x') {
				stadnStk[i].sequence += 't';;
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'i' && eixStk[i].sequence[j] == 'e') {
				stadnStk[i].sequence += 'a';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'j' && eixStk[i].sequence[j] == 'f') {
				stadnStk[i].sequence += 'b';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'k' && eixStk[i].sequence[j] == 'g') {
				stadnStk[i].sequence += 'c';
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'e' || eixStk[i].sequence[j-1] == 'f' || eixStk[i].sequence[j-1] == 'g') && eixStk[i].sequence[j] == 'i') {
				stadnStk[i].sequence += 'd';
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'e' || eixStk[i].sequence[j-1] == 'f' || eixStk[i].sequence[j-1] == 'g') && eixStk[i].sequence[j] == 'j') {
				stadnStk[i].sequence += 'o';
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'e' || eixStk[i].sequence[j-1] == 'f' || eixStk[i].sequence[j-1] == 'g') && eixStk[i].sequence[j] == 'k') {
				stadnStk[i].sequence += 'n';
				added++;
			}
			
			if (eixStk[i].sequence[j-1] == 'x' && eixStk[i].sequence[j] == '5') {
				stadnStk[i].sequence += '{';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '5' && eixStk[i].sequence[j] == 'h') {
				stadnStk[i].sequence += '~';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'h' && eixStk[i].sequence[j] == '5') {
				stadnStk[i].sequence += '|';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '5' && eixStk[i].sequence[j] == 'e') {
				stadnStk[i].sequence += 's';
				added++;
			}
			
			if ((eixStk[i].sequence[j-1] == 'e' || eixStk[i].sequence[j-1] == 'f' || eixStk[i].sequence[j-1] == 'g') && eixStk[i].sequence[j] == '3') {
				stadnStk[i].sequence += 't';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '3' && eixStk[i].sequence[j] == 'm') {
				stadnStk[i].sequence += 'z';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'm' && eixStk[i].sequence[j] == '3') {
				stadnStk[i].sequence += '?';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '3' && eixStk[i].sequence[j] == 'x') {
				stadnStk[i].sequence += '}';
				added++;
			}
			
			if (eixStk[i].sequence[j-1] == 'x' && eixStk[i].sequence[j] == '$') {
				stadnStk[i].sequence += ']';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '$' && eixStk[i].sequence[j] == 'M') {
				stadnStk[i].sequence += '0';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'M' && eixStk[i].sequence[j] == '$') {
				stadnStk[i].sequence += 'Z';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '$' && eixStk[i].sequence[j] == 'e') {
				stadnStk[i].sequence += 'T';
				added++;
			}
			
			if ((eixStk[i].sequence[j-1] == 'e' || eixStk[i].sequence[j-1] == 'f' || eixStk[i].sequence[j-1] == 'g') && eixStk[i].sequence[j] == '&') {
				stadnStk[i].sequence += 'S';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '&' && eixStk[i].sequence[j] == 'H') {
				stadnStk[i].sequence += '\\';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'H' && eixStk[i].sequence[j] == '&') {
				stadnStk[i].sequence += '^';
				added++;
			}
			if (eixStk[i].sequence[j-1] == '&' && eixStk[i].sequence[j] == 'x') {
				stadnStk[i].sequence += '[';
				added++;
			}
			
			if (eixStk[i].sequence[j-1] == 'x' && eixStk[i].sequence[j] == 'E') {
				stadnStk[i].sequence += 'T';
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'E' || eixStk[i].sequence[j-1] == 'F' || eixStk[i].sequence[j-1] == 'G') && eixStk[i].sequence[j] == 'x') {
				stadnStk[i].sequence += 'S';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'I' && eixStk[i].sequence[j] == 'E') {
				stadnStk[i].sequence += 'D';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'J' && eixStk[i].sequence[j] == 'F') {
				stadnStk[i].sequence += 'O';
				added++;
			}
			if (eixStk[i].sequence[j-1] == 'K' && eixStk[i].sequence[j] == 'G') {
				stadnStk[i].sequence += 'N';
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'E' || eixStk[i].sequence[j-1] == 'F' || eixStk[i].sequence[j-1] == 'G') && eixStk[i].sequence[j] == 'I') {
				stadnStk[i].sequence += 'A';
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'E' || eixStk[i].sequence[j-1] == 'F' || eixStk[i].sequence[j-1] == 'G') && eixStk[i].sequence[j] == 'J') {
				stadnStk[i].sequence += 'B';
				added++;
			}
			if ((eixStk[i].sequence[j-1] == 'E' || eixStk[i].sequence[j-1] == 'F' || eixStk[i].sequence[j-1] == 'G') && eixStk[i].sequence[j] == 'K') {
				stadnStk[i].sequence += 'C';
				added++;
			}
			
			if (eixStk[i].sequence[j-1] == '*' || eixStk[i].sequence[j] == '*') {
				stadnStk[i].sequence += '*';
				added++;
			}
			if (added != 1) {
				stadnStk[i].sequence += '*';
				cerr << "No rule for " << eixStk[i].sequence[j-1] << " " << eixStk[i].sequence[j] << endl;
			}
		}
		for (int i=0; i<numstklines; i++) {
			stadnStk[i].sequence += eixStk[i].sequence[j];
		}
	}
	
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