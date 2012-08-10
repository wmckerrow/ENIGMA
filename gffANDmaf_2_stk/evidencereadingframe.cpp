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
	if (argc < 3) {
		cerr << "Please execute like " << argv[0] << " inputStk posfile lineName evidenceType" << endl;
		cerr << "evidenceTypes are 1=est, 2=full" << endl;
		exit(1);
	}
	
	int numtochange = (argc-3)/2;
	
	string evidenceNames[numtochange];
	int evidenceTypes[numtochange];
	for (int i=0; i<numtochange; i++) {
		evidenceNames[i] = argv[3 + 2*i];
		evidenceTypes[i] = atoi(argv[4 + 2*i]);
	}
	
	//Read stk file
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr  << "Unable to open input alphabet stkfile " << argv[1] << endl;
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
	stkline inputStk[numstklines];
	for (int i=0; i<numstklines; i++) {
		getline(inFile,line);
		while (line[0] == '#' || line.size()==0) {
			cout << line << endl;
			getline(inFile,line);
		}
		inputStk[i] = readstkline(line);
	}
	inFile.close();
	
	//open segment postions file and read into memory
	int segposition[inputStk[0].sequence.size()];
	ifstream segfilein;
	segfilein.open(argv[2]);
	if (!segfilein) {
		cout << "Unable to open segmentation positions file " << argv[2] << endl;
		exit(1);
	}
	for (int i=0; i<inputStk[0].sequence.size(); i++) {
		segfilein >> segposition[i];
	}
	segfilein.close();
	
	//Initialize new stk
	stkline outputStk[numstklines];
	for (int i=0; i<numstklines; i++) {
		outputStk[i].label=inputStk[i].label;
	}
	
	//Update stk;
	for (int i=0; i<numstklines; i++) {
		int type=0;
		for (int j=0; j<numtochange; j++) {
			if (inputStk[i].label == evidenceNames[j]) {
				type = evidenceTypes[j];
			}
		}
		if (type != 1 && type != 2) {
			outputStk[i].sequence = inputStk[i].sequence;
		}
		
		if (type == 1) {
			for (int j=0; j<inputStk[0].sequence.size(); j++) {
				char inputLetter = inputStk[i].sequence[j];
				
				if (inputLetter == 's' || inputLetter == '{') {
					outputStk[i].sequence+='+';
					continue;
				}
				if (inputLetter == '5' || inputLetter == '3' || inputLetter == 'e' || inputLetter == 'f' || inputLetter == 'g') {
					outputStk[i].sequence+='-';
					continue;
				}
				if (inputLetter == '6' || inputLetter == '4' || inputLetter == 'p' || inputLetter == 'q' || inputLetter == 'r') {
					outputStk[i].sequence+=',';
					continue;
				}
				if (inputLetter == '~' || inputLetter == 'z' || inputLetter == 'd' || inputLetter == 'o' || inputLetter == 'n') {
					outputStk[i].sequence+='.';
					continue;
				}
				if (inputLetter == 'h' || inputLetter == 'm' || inputLetter == 'i' || inputLetter == 'j' || inputLetter == 'k') {
					outputStk[i].sequence+='/';
					continue;
				}
				if (inputLetter == 'l' || inputLetter == 'y' || inputLetter == 'u' || inputLetter == 'v' || inputLetter == 'w') {
					outputStk[i].sequence+=':';
					continue;
				}
				if (inputLetter == '|' || inputLetter == '?' || inputLetter == 'a' || inputLetter == 'b' || inputLetter == 'c') {
					outputStk[i].sequence+='<';
					continue;
				}
				if (inputLetter == 't' || inputLetter == '}') {
					outputStk[i].sequence+='=';
					continue;
				}
				
				if (inputLetter == 'T' || inputLetter == ']') {
					outputStk[i].sequence+='>';
					continue;
				}
				if (inputLetter == '$' || inputLetter == '@' || inputLetter == 'E' || inputLetter == 'F' || inputLetter == 'G') {
					outputStk[i].sequence+='@';
					continue;
				}
				if (inputLetter == '%' || inputLetter == '\'' || inputLetter == 'P' || inputLetter == 'Q' || inputLetter == 'R') {
					outputStk[i].sequence+='_';
					continue;
				}
				if (inputLetter == '0' || inputLetter == '\\' || inputLetter == 'A' || inputLetter == 'B' || inputLetter == 'C') {
					outputStk[i].sequence+='1';
					continue;
				}
				if (inputLetter == 'M' || inputLetter == 'H' || inputLetter == 'I' || inputLetter == 'J' || inputLetter == 'K') {
					outputStk[i].sequence+='2';
					continue;
				}
				if (inputLetter == 'Y' || inputLetter == 'L' || inputLetter == 'U' || inputLetter == 'V' || inputLetter == 'W') {
					outputStk[i].sequence+='7';
					continue;
				}
				if (inputLetter == 'Z' || inputLetter == '^' || inputLetter == 'D' || inputLetter == 'O' || inputLetter == 'N') {
					outputStk[i].sequence+='8';
					continue;
				}
				if (inputLetter == 'S' || inputLetter == '[') {
					outputStk[i].sequence+='9';
					continue;
				}
				
				outputStk[i].sequence+=inputLetter;
			}
		}
		
		if (type == 2) {
			int readingframe=0;
			int lastexonstart=0;
			for (int j=0; j<inputStk[0].sequence.size(); j++) {
				char inputLetter = inputStk[i].sequence[j];
				if (inputLetter=='s' || inputLetter=='T') {
					readingframe = 0;
					lastexonstart = segposition[j];
				}
				if (inputLetter=='a' || inputLetter=='D') {
					lastexonstart = segposition[j];
				}
				if (inputLetter=='d' || inputLetter=='A') {
					readingframe = (readingframe + segposition[j] - lastexonstart) % 3;
				}
				
				if (readingframe < 0 || readingframe > 2) {
					cerr << "readingframe=" << readingframe << endl;
				}
				
				switch (inputLetter) {
					case 'e':
						if (readingframe == 0) {
							outputStk[i].sequence+='e';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='f';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='g';
						}
						break;
					case 'E':
						if (readingframe == 0) {
							outputStk[i].sequence+='E';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='F';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='G';
						}
						break;
					case 'i':
						if (readingframe == 0) {
							outputStk[i].sequence+='i';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='j';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='k';
						}
						break;
					case 'I':
						if (readingframe == 0) {
							outputStk[i].sequence+='I';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='J';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='K';
						}
						break;
					case 'a':
						if (readingframe == 0) {
							outputStk[i].sequence+='a';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='b';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='c';
						}
						break;
					case 'A':
						if (readingframe == 0) {
							outputStk[i].sequence+='A';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='B';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='C';
						}
						break;
					case 'd':
						if (readingframe == 0) {
							outputStk[i].sequence+='d';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='o';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='n';
						}
						break;
					case 'D':
						if (readingframe == 0) {
							outputStk[i].sequence+='D';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='O';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='N';
						}
						break;
					case 'p':
						if (readingframe == 0) {
							outputStk[i].sequence+='p';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='q';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='r';
						}
						break;
					case 'P':
						if (readingframe == 0) {
							outputStk[i].sequence+='P';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='Q';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='R';
						}
						break;
					case 'u':
						if (readingframe == 0) {
							outputStk[i].sequence+='u';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='v';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='w';
						}
						break;
					case 'U':
						if (readingframe == 0) {
							outputStk[i].sequence+='U';
						}
						if (readingframe == 1) {
							outputStk[i].sequence+='V';
						}
						if (readingframe == 2) {
							outputStk[i].sequence+='W';
						}
						break;
					default:
						outputStk[i].sequence += inputLetter;
						break;
				}
			}
		}
	}
	
	/*
	 for (int i=0; i<numstklines; i++) {
	 cerr << inputStk[i].sequence.size() << " " << outputStk[i].sequence.size() << endl;
	 }
	 */
	
	int maxlabellength=0;
	for (int i=0; i<numstklines; i++) {
		if (outputStk[i].label.size() > maxlabellength) {
			maxlabellength = outputStk[i].label.size();
		}
	}
	int labelcolumns = max(maxlabellength+3,10);
	for (int i=0; i<numstklines; i++) {
		cout << outputStk[i].label;
		for (int j=0; j<labelcolumns-outputStk[i].label.size(); j++) {
			cout << " ";
		}
		cout << outputStk[i].sequence << endl;
	}
}