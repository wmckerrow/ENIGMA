/* 
 This code generates a valid set of positions for the segments in an L1 stk file.
 1. Read stk file into memory
 2. Start with position=1
 3. For each column
	a. Print position
	b. If the next column contains a transition to intron, add 0, 1, 2 to position so that we are in the right reading frame
	c. If the next column contains a transition to exon, remember which reading frame we are in
	d. If this column is not a transition column add 30 to postion.
 
 Note: There is no guarantee that there exists a valid set of positions for an L1 stk. If there is no valid position this code will happily output one that is not valid.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <math.h>
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

bool is_section(char letter) {
	string sectionletters="egfEFGijkIJKxX";
	if (sectionletters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

bool is_trans2int(char letter) {
	string trans2intronletters="donABC";
	if (trans2intronletters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

int main (int argc, char *argv[]) {
	//make sure we have enough input
	if (argc < 2) {
		cerr << "Please execute like " << argv[0] << "AligmentFileName" << endl;
		exit(1);
	}
	
	//Read stk file
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr  << "Unable to open stkfile " << argv[1] << endl;
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
	stkline L1stk[numstklines];
	for (int i=0; i<numstklines; i++) {
		getline(inFile,line);
		while (line[0] == '#') {
			getline(inFile,line);
		}
		L1stk[i] = readstkline(line);
	}
	inFile.close();
	stkline stadnStk[numstklines];
	for (int i=0; i<numstklines; i++) {
		stadnStk[i].label = L1stk[i].label;
		stadnStk[i].sequence = L1stk[i].sequence[0];
	}
	
	//Generate positions
	int currentpos = 1;
	int readingframe[numstklines];
	for (int i=0; i<numstklines; i++) {
		readingframe[i] = 0;
	}
	for (int j=0; j<L1stk[0].sequence.size(); j++) {
		//output position
		cout << currentpos << " ";
		
		//update reading frame
		if (j != L1stk[0].sequence.size()-1) {
			for (int i=0; i<numstklines; i++) {
				bool changepos=0;
				//if next letter is the start of an intron, adjust the position to match reading frame
				//if the next letter is d or A, we want currentpos%3=readingframe%3
				if (L1stk[i].sequence[j+1] == 'd' || L1stk[i].sequence[j+1] == 'A') {
					if (currentpos%3 == (readingframe[i]+1)%3) {
						currentpos+=2;
						if (changepos) {
							cerr << "a)double change pos: "  << i << " " << j << endl;
							cerr << "currentpos=" << currentpos << " readingframe=" << readingframe[i] << endl;
						}
						changepos=1;
					}
					if (currentpos%3 == (readingframe[i]+2)%3) {
						currentpos++;
						if (changepos) {
							cerr << "b)double change pos: "  << i << " " << j << endl;
							cerr << "currentpos=" << currentpos << " readingframe=" << readingframe[i] << endl;
						}
						changepos=1;
					}
				}
				//if the next letter is o or B we want (currentpos-1)%3 == readingframe[i]%3
				if (L1stk[i].sequence[j+1] == 'o' || L1stk[i].sequence[j+1] == 'B') {
					//cerr << "At " << j << " next is o or B; currentpos=" << currentpos << "; readingframe[i]=" << readingframe[i] << "; (currentpos-readingframe[i])%3=" << (currentpos-readingframe[i])%3 << endl;
					if ((currentpos-1)%3 == (readingframe[i]+1)%3) {
						currentpos+=2;
						if (changepos) {
							cerr << "c)double change pos: " << i << " " << j << endl;
							cerr << "currentpos=" << currentpos << " readingframe=" << readingframe[i] << endl;
						}
						changepos=1;
					}
					if ((currentpos-1)%3 == (readingframe[i]+2)%3) {
						currentpos++;
						if (changepos) {
							cerr << "d)double change pos: " << i << " " << j << endl;
							cerr << "currentpos=" << currentpos << " readingframe=" << readingframe[i] << endl;
						}
						changepos=1;
					}
				}
				//if the next letter is n or C we want (currentpos-2)%3 == readingframe[i]%3
				if (L1stk[i].sequence[j+1] == 'n' || L1stk[i].sequence[j+1] == 'C') {
					if ((currentpos-2)%3 == (readingframe[i]+1)%3) {
						currentpos+=2;
						if (changepos) {
							cerr << "e)double change pos: " << i << " " << j << endl;
							cerr << "currentpos=" << currentpos << " readingframe=" << readingframe[i] << endl;
						}
						changepos=1;
					}
					if ((currentpos-2)%3 == (readingframe[i]+2)%3) {
						currentpos++;
						if (changepos) {
							cerr << "f)double change pos: " << i << " " << j << endl;
							cerr << "currentpos=" << currentpos << " readingframe=" << readingframe[i] << endl;
						}
						changepos=1;
					}
				}
				
				//if start of exon, set reading frame
				if (L1stk[i].sequence[j+1] == 's' || L1stk[i].sequence[j+1] == 'T' || L1stk[i].sequence[j+1] == 'a' || L1stk[i].sequence[j+1] == 'D') {
					readingframe[i] = currentpos%3;
				}
				if (L1stk[i].sequence[j+1] == 'b' || L1stk[i].sequence[j+1] == 'O') {
					readingframe[i] = (currentpos+2) %3;
				}
				if (L1stk[i].sequence[j+1] == 'c' || L1stk[i].sequence[j+1] == 'N') {
					readingframe[i] = (currentpos+1) %3;
				}
			}
		}
		
		//add 30 if not a transition column
		if (is_section(L1stk[0].sequence[j])) {
			currentpos += 30;
		}
	}
}
