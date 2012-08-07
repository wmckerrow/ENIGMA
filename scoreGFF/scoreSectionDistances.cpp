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

int sgn(int input) {
	if (input < 0) {
		return -1;
	}
	if (input > 0) {
		return 1;
	}
	return 0;
}

string convertInt(int number)
{
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

int isminus(int input) {
	if (input < 0) {
		return 1;
	}
	return 0;
}

int isplus(int input) {
	if (input > 0) {
		return 1;
	}
	return 0;
}

struct gffrow {
	string sequence;
	string source;
	string featureType;
	int start;
	int end;
	float score;
	char strand;
	int frame;
	string group;
};

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

bool incds(int position, gffrow *gff, int numgffrows) {
	for (int i=0; i<numgffrows; i++) {
		if (gff[i].start <= position && position <= gff[i].end) {
			if (gff[i].featureType == "CDS") {
				return 1;
			}
		}
	}
	return 0;
}

bool inunaligned(int position, gffrow *gff, int numgffrows) {
	for (int i=0; i<numgffrows; i++) {
		if (gff[i].start <= position && position <= gff[i].end) {
			if (gff[i].featureType == "unaligned") {
				return 1;
			}
		}
	}
	return 0;
}

int *getintsfromstring(string line) {
	int pos=0;
	int *result;
	result = new int[2];
	while (line[pos] != ' ') {
		pos++;
	}
	if (!from_string<int>(result[0],line.substr(0,pos),std::dec)) {
		cerr << line.substr(0,pos) << " is not an integer" << endl;
		exit(1);
	}
	if (!from_string<int>(result[1],line.substr(pos+1,line.size()-pos),std::dec)) {
		cerr << line.substr(0,pos) << " is not an integer" << endl;
		exit(1);
	}
	return result;
}

int main (int argc, char *argv[]) {
	if (argc < 3) {
		cerr << "Please execute like " << argv[0] << " sectionList numGFFs gffFiles table? fullOutput?" << endl;
		exit(1);
	}
	
	string line;
	int numGFFs = atoi(argv[2]);
	bool table;
	if (argv[argc-2][0] == 'y' || argv[argc-2][0] == 'Y') {
		table = 1;
	}
	else {
		table = 0;
	}
	bool fullOutput;
	if (argv[argc-1][0] == 'y' || argv[argc-1][0] == 'Y') {
		fullOutput = 1;
	}
	else {
		fullOutput = 0;
	}
	
	ifstream gffFile;
	gffrow **GFFs;
	GFFs = new gffrow*[numGFFs];
	int GFFsize[numGFFs];
	for (int i=0; i<numGFFs; i++) {
		gffFile.open(argv[3+i]);
		//cout << argv[3+i] << endl;
		if (!gffFile) {
			cerr << "Unable to open gffFile " << argv[3+i] << endl;
			exit(1);
		}
		GFFsize[i]=0;
		while (getline(gffFile,line)) {
			if (line[0] != '#') {
				GFFsize[i]++;
			}
		}
		gffFile.close();
		GFFs[i]=new gffrow[GFFsize[i]];
		gffFile.open(argv[3+i]);
		for (int j=0; j<GFFsize[i]; j++) {
			getline(gffFile,line);
			while (line[0] == '#') {
				getline(gffFile,line);
			}
			GFFs[i][j]=getGFFrow(line);
		}
		gffFile.close();
	}
	
	for (int i=0; i<numGFFs; i++) {
		int TP[numGFFs];
		int TN[numGFFs];
		int FP[numGFFs];
		int FN[numGFFs];
		int unaligned[numGFFs];
		for (int k=0; k<numGFFs; k++) {
			TP[k]=0;
			TN[k]=0;
			FP[k]=0;
			FN[k]=0;
		}
		ifstream sectionList;
		sectionList.open(argv[1]);
		if (!sectionList) {
			cerr << "Unable to open " << argv[1] << endl;
			exit(1);
		}
		while (getline(sectionList,line)) {
			int *startend;
			startend = getintsfromstring(line);
			int start = startend[0];
			int end = startend[1];
			for (int j=start; j<=end; j++) {
				for (int k=0; k<numGFFs; k++) {
					bool kcds = incds(j,GFFs[k],GFFsize[k]);
					bool icds = incds(j,GFFs[i],GFFsize[i]);
					if (icds && kcds) {
						TP[k]++;
					}
					if (icds && !kcds) {
						FP[k]++;
					}
					if (!icds && kcds) {
						FN[k]++;
					}
					if (!icds && !kcds) {
						TN[k]++;
					}
				}
			}
		}
		
		if (!table) {
			cout << "For " << argv[3+i] << ":" << endl;
		}

		for (int k=0; k<numGFFs; k++) {
			float sensitivity = (float)TP[k]/(TP[k]+FN[k]);
			float backwardsensitivity = (float)TP[k]/(TP[k]+FP[k]);
			float distance=1-(sensitivity + backwardsensitivity) / 2;
			if (table) {
				cout << distance << " ";
			}
			else {
				cout << "Distance to " << argv[3+k] << ": " << distance << endl;
			}

		}
		cout << endl;
	}
}