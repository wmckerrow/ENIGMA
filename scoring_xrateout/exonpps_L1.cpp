/*
 At each column position, gives the posterior probability that there is an exon at the postion, followed by whether or not there actually is an exon there.
 Can be used to generate ROC curves.
 Use with L1 grammars: L1_trans_strands_4transpose.eg
 First it reads through the posterior probability output given by xrate's -pp option to get the posterior probability for each letter.
 For each section column it outputs the total posterior probability of exon.
 Then it checks the correct stk and output 1 if it is exon and 0 if it is not exon.
 */

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

class gffrow {
public:
	string sequence;
	string source;
	string featureType;
	int start;
	int end;
	float score;
	char strand;
	int frame;
	string group;
	
	void print();
};

void gffrow::print() {
	cout << sequence << '\t' << source << '\t' << featureType << '\t' << start << '\t' << end << '\t' << score << '\t';
	if (strand == -1) {
		cout << ".";
	}
	else {
		cout << strand;
	}
	cout << '\t' << frame << '\t' << group << endl;
}

gffrow getGFFrow(string line) {
	int numtabs=0;
	for (int i=0; i<line.length(); i++) {
		if (line[i]=='\t') {
			numtabs++;
		}
	}
	if (numtabs < 8) {
		cerr << line << endl << "is not a proper gff line" << endl;
		exit(1);
	}
	
	int tab=-1;
	int nexttab;
	string data[9];
	for (int i=0; i<9; i++) {
		nexttab = tab+1;
		while (line[nexttab] != '\t' && nexttab < line.size()) {
			nexttab++;
		}
		data[i] = line.substr(tab+1,nexttab-tab-1);
		tab=nexttab;
	}
	gffrow row;
	row.sequence=data[0];
	row.source=data[1];
	row.featureType=data[2];
	if (!from_string<int>(row.start,data[3],std::dec)) {
		cerr << "from_string failed on data[3]" << endl;
	}
	if (!from_string<int>(row.end,data[4],std::dec)) {
		cerr << "from_string failed on data[4]" << endl;
	}
	if (!from_string<float>(row.score,data[5],std::dec)) {
		row.score=0;
	}
	row.strand=data[6][0];
	if (!from_string<int>(row.frame,data[7],std::dec)) {
		row.frame=-1;
	}
	row.group=data[8];
	//cout << "got gffrow" << endl;
	return row;
}

char findletter(int position, gffrow *gff, int numgffrows) {
	if (position < 1) {
		return '*';
	}
	bool geneCDSbefore=0;
	bool geneCDSafter=0;
	bool inthisgene=0;
	for (int i=0; i<numgffrows; i++) {
		if (gff[i].featureType == "gene") {
			inthisgene=0;
		}
		if (gff[i].start <= position && position <= gff[i].end) {
			if (gff[i].featureType == "CDS") {
				return 'e';
			}
			if (gff[i].featureType == "gene") {
				inthisgene = 1;
			}
		}
		if (inthisgene && position > gff[i].end && gff[i].featureType == "CDS") {
			geneCDSbefore=1;
		}
		if (inthisgene && position < gff[i].start && gff[i].featureType == "CDS") {
			geneCDSafter=1;
		}
		if (geneCDSbefore && geneCDSafter) {
			return 'i';
		}
	}
	return 'x';
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


int main (int argc, char *argv[]) {
	if (argc < 4) {
		cerr << "Please execute like " << argv[0] << " nameOfTarget ppsFromXrate genout" << endl;
		exit(1);
	}
	
	string nameOfTarget=argv[1];
	
	ifstream inFile;
	inFile.open(argv[2]);
	if (!inFile) {
		cerr << "Unable to open xrate output" << argv[2] << endl;
		exit(1);
	}
	
	string line;
	int numpplines=0;
	while (getline(inFile,line)) {
		if (line.substr(0,8)=="#=GF GFF") {
			numpplines++;
		}
	}
	//cerr << numpplines << endl;
	gffrow *ppdata;
	ppdata = new gffrow[numpplines];
	
	inFile.close();
	inFile.open(argv[2]);
	
	string one;
	string two;
	int whichline=0;
	
	while (!inFile.eof()) {
		inFile >> one;
		inFile >> two;
		getline(inFile,line);
		if (one == "#=GF" && two == "GFF") {
			if (whichline >= numpplines) {
				cerr << "Miscount" << endl;
				break;
			}
			line = line.substr(19,line.size()-19);
			ppdata[whichline] = getGFFrow(line);
			whichline++;
		}
	}
	
	inFile.close();
		
	/*
	for (int i=0; i<numpplines; i++) {
		ppdata[i].print();
	}
	 */
	
	int numcols;
	for (int i=0; i<numpplines; i++) {
		if (ppdata[i].start==ppdata[i].end && ppdata[i].start != 1) {
			numcols=ppdata[i].end;
			break;
		}
	}

	double **forwardprobs;
	forwardprobs = new double*[numcols];
	for (int i=0; i<numcols; i++) {
		forwardprobs[i] = new double[43];
	}
	
	double **backwardprobs;
	backwardprobs = new double*[numcols];
	for (int i=0; i<numcols; i++) {
		backwardprobs[i] = new double[43];
	}
	
	for (int i=0; i<numcols; i++) {
		for (int j=0; j<43; j++) {
			forwardprobs[i][j]=0;
			backwardprobs[i][j]=0;
		}
	}
	
	for (int i=0; i<numpplines; i++) {
		//ppdata[i].print();
		if (ppdata[i].end > numcols) {
			cerr << ppdata[i].end << " > " << numcols << endl;
		}
		if (ppdata[i].start < 1) {
			cerr << ppdata[i].start << " < 1" << endl;
		}
		if (ppdata[i].featureType=="E0_forward_section" || ppdata[i].featureType=="E_forward_section") {
				if (ppdata[i].start==1) {
					forwardprobs[ppdata[i].end-1][0]=exp(ppdata[i].score);
					continue;
				}
				else {
					backwardprobs[ppdata[i].start-1][0]=exp(ppdata[i].score);
					continue;
				}
		}
		if (ppdata[i].featureType=="E1_forward_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][1]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][1]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="E2_forward_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][2]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][2]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I0_forward_section" || ppdata[i].featureType=="I_forward_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][3]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][3]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I1_forward_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][4]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][4]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I2_forward_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][5]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][5]=exp(ppdata[i].score);
				continue;
			}
		}
		
		if (ppdata[i].featureType=="E0_reverse_section" || ppdata[i].featureType=="E_reverse_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][6]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][6]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="E1_reverse_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][7]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][7]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="E2_reverse_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][8]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][8]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I0_reverse_section" || ppdata[i].featureType=="I_reverse_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][9]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][9]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I1_reverse_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][10]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][10]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I2_reverse_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][11]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][11]=exp(ppdata[i].score);
				continue;
			}
		}
		
		if (ppdata[i].featureType=="X_section") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][12]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][12]=exp(ppdata[i].score);
				continue;
			}
		}
		
		if (ppdata[i].featureType=="S_forward") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][13]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][13]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="T_forward") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][14]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][14]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="A0_forward" || ppdata[i].featureType=="A_forward") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][15]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][15]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="A1_forward") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][16]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][16]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="A2_forward") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][17]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][17]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="D0_forward" || ppdata[i].featureType=="D_forward") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][18]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][18]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="D1_forward") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][19]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][19]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="D2_forward") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][20]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][20]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="E0_forward_trans" || ppdata[i].featureType=="E_forward_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][21]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][21]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="E1_forward_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][22]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][22]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="E2_forward_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][23]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][23]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I0_forward_trans" || ppdata[i].featureType=="I_forward_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][24]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][24]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I1_forward_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][25]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][25]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I2_forward_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][26]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][26]=exp(ppdata[i].score);
				continue;
			}
		}
		
		if (ppdata[i].featureType=="X_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][27]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][27]=exp(ppdata[i].score);
				continue;
			}
		}

		if (ppdata[i].featureType=="T_reverse") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][28]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][28]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="S_reverse") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][29]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][29]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="D0_reverse" || ppdata[i].featureType=="D_reverse") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][30]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][30]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="D1_reverse") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][31]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][31]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="D2_reverse") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][32]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][32]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="A0_reverse" || ppdata[i].featureType=="A_reverse") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][33]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][33]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="A1_reverse") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][34]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][34]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="A2_reverse") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][35]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][35]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="E0_reverse_trans" || ppdata[i].featureType=="E_reverse_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][36]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][36]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="E1_reverse_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][37]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][37]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="E2_reverse_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][38]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][38]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I0_reverse_trans" || ppdata[i].featureType=="I_reverse_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][39]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][39]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I1_reverse_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][40]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][40]=exp(ppdata[i].score);
				continue;
			}
		}
		if (ppdata[i].featureType=="I2_reverse_trans") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][41]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][41]=exp(ppdata[i].score);
				continue;
			}
		}
		
		if (ppdata[i].featureType=="START") {
			if (ppdata[i].start==1) {
				forwardprobs[ppdata[i].end-1][42]=exp(ppdata[i].score);
				continue;
			}
			else {
				backwardprobs[ppdata[i].start-1][42]=exp(ppdata[i].score);
				continue;
			}
		}
		cout << "No case for " << ppdata[i].featureType << endl;
	}	
		
	double **pps;
	pps = new double*[numcols];
	for (int i=0; i<numcols; i++) {
		pps[i] = new double[42];
	}
	
	for (int i=0; i<numcols; i++) {
		double probsum=0;
		for (int j=0; j<42; j++) {
			probsum+=forwardprobs[i][j]*backwardprobs[i][j];
		}
		for (int j=0; j<42; j++) {
			pps[i][j] = forwardprobs[i][j]*backwardprobs[i][j]/probsum;
			//cout << fixed;
			//cout << setprecision(5) << backwardprobs[i][j] << " ";
			//cout << setprecision(5) << pps[i][j] << " ";
		}
		//cout << endl;
	}

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
	
	int numlinestoprint = 0;
	for (int i=1; i<numcols; i++) {
		double largeprob = 0;
		for (int j=0; j<=12; j++) {
			largeprob += pps[i][j];
		}
		if (largeprob > .5) {
			numlinestoprint++;
		}
	}
	
	cout << numlinestoprint << " " << numlinestoprint << endl;
	
	for (int i=1; i<numcols; i++) {
		double largeprob = 0;
		for (int j=0; j<=12; j++) {
			largeprob += pps[i][j];
		}
		if (largeprob > .5) {
			cout << fixed << setprecision(16) << pps[i][0]+pps[i][1]+pps[i][2]+pps[i][6]+pps[i][7]+pps[i][8] << " ";
			if (genoutTarget.sequence[i] == 'e' || genoutTarget.sequence[i] == 'f' || genoutTarget.sequence[i] == 'g' || genoutTarget.sequence[i] == 'E' || genoutTarget.sequence[i] == 'F' || genoutTarget.sequence[i] == 'G') {
				cout << 1;
			}
			else {
				cout << 0;
			}
			cout << endl;
		}
	}

}