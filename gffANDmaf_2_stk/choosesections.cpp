/*
 This program reads through one or more outputs of prepare4stk and gives only those sections that have an exon evidence alignment coverage score (EAC) above a certain threshold.
 The command line arguments are: outputs from prepare4stk, a threshold for EAC. For example:
 ./choosesections alignmentSections1 alignmentSections2 .1 > alignmentSectionsMinIs0.1
*/

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

template <class T>
bool from_string(T& t, 
                 const std::string& s, 
                 std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

class sectionScore {
public:
	int start;
	int end;
	float score;
};

int main (int argc, char *argv[]) {
	if (argc < 4) {
		cerr << "Please execute like " << argv[0] << "numLists sectionScoresLists threshold" << endl;
		exit(1);
	}
	
	string line;
	int numLists=atoi(argv[1]);
	float threshold = atof(argv[argc-1]);
	ifstream inFile;
	int numrawsections=0;
	
	//Read all lists into a single long (raw) list
	for (int i=0; i<numLists; i++) {
		inFile.open(argv[2+i]);
		if (!inFile) {
			cerr << "Unable to open " << argv[2+i] << endl;
			exit(1);
		}
		while (getline(inFile,line)) {
			for (int j=line.size()-2; j>=0; j--) {
				if (line[j]==' ') {
					float score;
					if (!from_string<float>(score,line.substr(j+1,line.size()-(j+1)),std::dec)) {
						cerr << "from_string failed on" << line.substr(j+1,line.size()-(j+1)) << endl;
					}
					if (score >= threshold) {
						numrawsections++;
					}
					break;
				}
			}
		}
		inFile.close();
	}
		
	sectionScore rawlist[numrawsections];
	int whichrawsection=0;
	
	for (int i=0; i<numLists; i++) {
		inFile.open(argv[2+i]);
		while (!inFile.eof()) {
			int start;
			int end;
			float score;
			inFile >> start;
			inFile >> end;
			inFile >> score;
			if (score >= threshold) {
				rawlist[whichrawsection].start=start;
				rawlist[whichrawsection].end=end;
				rawlist[whichrawsection].score=score;
				whichrawsection++;
			}
		}
		inFile.close();
	}
	
	/*
	cout << numrawsections << endl;
	
	for (int i=0; i<numrawsections; i++) {
		cout << rawlist[i].start << " " << rawlist[i].end << " " << endl;
	}
	 //*/
	
	//Bubble sort
	bool switched=1;
	while (switched) {
		switched=0;
		for (int i=0; i<numrawsections-1; i++) {
			if (rawlist[i].start > rawlist[i+1].start) {
				switched=1;
				sectionScore tempSectionScore=rawlist[i];
				rawlist[i]=rawlist[i+1];
				rawlist[i+1]=tempSectionScore;
			}
		}
	}
	
	//Output chosen sections
	int i=0;
	while (i < numrawsections) {
		int j=i+1;
		int end=rawlist[i].end;
		while (rawlist[i].end >= rawlist[j].start && j<numrawsections) {
			if (rawlist[j].end > end) {
				end=rawlist[j].end;
			}
			j++;
		}
		cout << rawlist[i].start << " " << end << " " << endl;
		i=j;
	}
}