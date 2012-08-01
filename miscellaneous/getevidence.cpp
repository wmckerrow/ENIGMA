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
	if (argc < 3) {
		cerr << "Please execute like " << argv[0] << " genout newtree numevidence evidenceLabels geneInclusionProbs numAlignedSpecies alignedSpeciesLabels geneInclusionProbs" << endl;
		exit(1);
	}
	string newtree = argv[2];
	cout << newtree << endl;
	
	int numevidence = atoi(argv[3]);
	string evidenceLabels[numevidence];
	float evidenceGeneInclusion[numevidence];
	for (int i=0; i<numevidence; i++) {
		evidenceLabels[i] = argv[4+i];
		evidenceGeneInclusion[i] = atof(argv[4+numevidence+i]);
	}
	
	int numAlignedSpecies = atoi(argv[4+2*numevidence]);
	string alignedSpeciesLabels[numAlignedSpecies];
	float alignedGeneInclusion[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		alignedSpeciesLabels[i] = argv[5+2*numevidence+i];
		alignedGeneInclusion[i] = atof(argv[5+2*numevidence+numAlignedSpecies+i]);
	}
	
	/*
	for (int i=0; i<numevidence; i++) {
		cout << evidenceLabels[i] << " " << evidenceGeneInclusion[i] << endl;
	}
	for (int i=0; i<numAlignedSpecies; i++) {
		cout << alignedSpeciesLabels[i] << " " << alignedGeneInclusion[i] << endl;
	}
	 //*/
	
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr  << "Unable to open genout " << argv[1] << endl;
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
	stkline genout[numstklines];
	for (int i=0; i<numstklines; i++) {
		getline(inFile,line);
		while (line[0] == '#') {
			getline(inFile,line);
		}
		genout[i] = readstkline(line);
	}
	inFile.close();
	stkline evidenceStk[numevidence+numAlignedSpecies];
	for (int i=0; i<numevidence; i++) {
		evidenceStk[i].label = evidenceLabels[i];
		evidenceStk[i].sequence = "";
	}
	for (int i=0; i<numAlignedSpecies; i++) {
		evidenceStk[numevidence+i].label = alignedSpeciesLabels[i];
		evidenceStk[numevidence+i].sequence = "";
	}
	
	bool evidencehasthisgene[numevidence];
	for (int i=0; i<numevidence; i++) {
		float randreal = (float)rand()/(float)RAND_MAX;
		if (randreal <= evidenceGeneInclusion[i]) {
			evidencehasthisgene[i] = 1;
		}
		else {
			evidencehasthisgene[i] = 0;
		}
	}
	
	bool alignedhasthisgene[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		float randreal = (float)rand()/(float)RAND_MAX;
		if (randreal <= alignedGeneInclusion[i]) {
			alignedhasthisgene[i] = 1;
		}
		else {
			alignedhasthisgene[i] = 0;
		}
	}
	
	int strand=rand()%2;
	
	for (int j=0; j<genout[0].sequence.size(); j++) {
		bool allx=1;
		for (int i=0; i<numstklines; i++) {
			if (genout[i].sequence[j] != 'x') {
				allx=0;
				break;
			}
		}
		if (allx) {
			for (int i=0; i<numevidence; i++) {
				float randreal = (float)rand()/(float)RAND_MAX;
				if (randreal <= evidenceGeneInclusion[i]) {
					evidencehasthisgene[i] = 1;
				}
				else {
					evidencehasthisgene[i] = 0;
				}
			}
			for (int i=0; i<numAlignedSpecies; i++) {
				float randreal = (float)rand()/(float)RAND_MAX;
				if (randreal <= alignedGeneInclusion[i]) {
					alignedhasthisgene[i] = 1;
				}
				else {
					alignedhasthisgene[i] = 0;
				}
			}
		}
		for (int i=0; i<numevidence; i++) {
			for (int k=0; k<numstklines; k++) {
				if (evidenceLabels[i] == genout[k].label) {
					if (evidencehasthisgene[i]) {
						evidenceStk[i].sequence+=genout[k].sequence[j];
					}
					else {
						evidenceStk[i].sequence+='x';
					}
					break;
				}
			}
		}
		for (int i=0; i<numAlignedSpecies; i++) {
			for (int k=0; k<numstklines; k++) {
				if (alignedSpeciesLabels[i] == genout[k].label) {
					if (alignedhasthisgene[i]) {
						evidenceStk[i+numevidence].sequence+=genout[k].sequence[j];
					}
					else {
						evidenceStk[i+numevidence].sequence+='*';
					}
					break;
				}
			}
		}
	}
	
	int maxlabellength=0;
	for (int i=0; i<numevidence+numAlignedSpecies; i++) {
		if (evidenceStk[i].label.size() > maxlabellength) {
			maxlabellength = evidenceStk[i].label.size();
		}
	}
	int labelcolumns = max(maxlabellength+3,10);
	for (int i=0; i<numevidence+numAlignedSpecies; i++) {
		cout << evidenceStk[i].label;
		for (int j=0; j<labelcolumns-evidenceStk[i].label.size(); j++) {
			cout << " ";
		}
		cout << evidenceStk[i].sequence << endl;
	}
}