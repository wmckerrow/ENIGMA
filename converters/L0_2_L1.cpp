/*
 This takes an L0 stk together with a set of position and generates an L1 stk. It will not give a good answer if there is an * in the input. It will also not give a good answer if the stk does not start x.
 It takes as input a .stk file in the eix alphabet.
 It reads through the stk input creating a new stk as it goes.
 1. It reads the stk into memory.
 2. It reads the pos file.
 3. It copies the line labels to the new stk.
 4. It goes through each sequence character by character.
	a. If the letter is s or T then reading frame is reset to 0
	b. If the letter is a or D, the start of the exon is remembered.
	c. If the letter is d or A, the length of the exon is calculated mod 3 and added to reading frame.
	d. Given the current letter and the reading frame, the corresponding L1 letter is added to the new stk.
 5. Once the entire stk input has been read the new stk is printed to stdout.
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
	//Make sure there are enough arguments
	if (argc < 3) {
		cerr << "Please execute like " << argv[0] << " eixStkFile posFile" << endl;
		exit(1);
	}
	
	//Read stk file
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
	
	//open segment postions file and read into memory
	int segposition[eixStk[0].sequence.size()];
	ifstream segfilein;
	segfilein.open(argv[2]);
	if (!segfilein) {
		cout << "Unable to open segmentation positions file " << argv[2] << endl;
		exit(1);
	}
	for (int i=0; i<eixStk[0].sequence.size(); i++) {
		segfilein >> segposition[i];
	}
	segfilein.close();
	
	//Initialize new stk
	stkline L1stk[numstklines];
	for (int i=0; i<numstklines; i++) {
		L1stk[i].label=eixStk[i].label;
	}
	
	//Update stk;
	for (int i=0; i<numstklines; i++) {
		int readingframe=0;
		int lastexonstart=0;
		for (int j=0; j<eixStk[0].sequence.size(); j++) {
			char eixLetter = eixStk[i].sequence[j];
			if (eixLetter=='s' || eixLetter=='T') {
				readingframe = 0;
				lastexonstart = segposition[j];
			}
			if (eixLetter=='a' || eixLetter=='D') {
				lastexonstart = segposition[j];
			}
			if (eixLetter=='d' || eixLetter=='A') {
				readingframe = (readingframe + segposition[j] - lastexonstart) % 3;
			}
			
			if (readingframe < 0 || readingframe > 2) {
				cerr << "readingframe=" << readingframe << endl;
			}
			
			switch (eixLetter) {
				case 'e':
					if (readingframe == 0) {
						L1stk[i].sequence+='e';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='f';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='g';
					}
					break;
				case 'E':
					if (readingframe == 0) {
						L1stk[i].sequence+='E';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='F';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='G';
					}
					break;
				case 'i':
					if (readingframe == 0) {
						L1stk[i].sequence+='i';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='j';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='k';
					}
					break;
				case 'I':
					if (readingframe == 0) {
						L1stk[i].sequence+='I';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='J';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='K';
					}
					break;
				case 'a':
					if (readingframe == 0) {
						L1stk[i].sequence+='a';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='b';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='c';
					}
					break;
				case 'A':
					if (readingframe == 0) {
						L1stk[i].sequence+='A';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='B';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='C';
					}
					break;
				case 'd':
					if (readingframe == 0) {
						L1stk[i].sequence+='d';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='o';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='n';
					}
					break;
				case 'D':
					if (readingframe == 0) {
						L1stk[i].sequence+='D';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='O';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='N';
					}
					break;
				case 'p':
					if (readingframe == 0) {
						L1stk[i].sequence+='p';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='q';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='r';
					}
					break;
				case 'P':
					if (readingframe == 0) {
						L1stk[i].sequence+='P';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='Q';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='R';
					}
					break;
				case 'u':
					if (readingframe == 0) {
						L1stk[i].sequence+='u';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='v';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='w';
					}
					break;
				case 'U':
					if (readingframe == 0) {
						L1stk[i].sequence+='U';
					}
					if (readingframe == 1) {
						L1stk[i].sequence+='V';
					}
					if (readingframe == 2) {
						L1stk[i].sequence+='W';
					}
					break;
				default:
					L1stk[i].sequence += eixLetter;
					break;
			}
		}
	}
	
	/*
	for (int i=0; i<numstklines; i++) {
		cerr << eixStk[i].sequence.size() << " " << L1stk[i].sequence.size() << endl;
	}
	 */
	
	int maxlabellength=0;
	for (int i=0; i<numstklines; i++) {
		if (L1stk[i].label.size() > maxlabellength) {
			maxlabellength = L1stk[i].label.size();
		}
	}
	int labelcolumns = max(maxlabellength+3,10);
	for (int i=0; i<numstklines; i++) {
		cout << L1stk[i].label;
		for (int j=0; j<labelcolumns-L1stk[i].label.size(); j++) {
			cout << " ";
		}
		cout << L1stk[i].sequence << endl;
	}
}