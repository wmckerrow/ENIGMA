#define MAXALIGNMENT 3
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

struct dotline {
	int starts[NUMSPECIES];
	int length;
};

struct partposition {
	string name;
	int chromosome;
	int start;
	char strand;
};

dotline getDotLine (string line) { //read space for tab...
	int tab=-1;
	int nexttab;
	string data[NUMSPECIES+1];
	for (int i=0; i<NUMSPECIES+1; i++) {
		nexttab = tab+1;
		while (line[nexttab] != ' ' && nexttab < line.size()) {
			nexttab++;
		}
		data[i] = line.substr(tab+1,nexttab-tab-1);
		tab=nexttab;
	}
	dotline row;
	
	for (int i=0; i<NUMSPECIES; i++) {
		if (!from_string<int>(row.starts[i],data[i],std::dec)) {
			cerr << "from_string failed on data[" << i << "]" << endl;
		}
	}
	if (!from_string<int>(row.length,data[NUMSPECIES],std::dec)) {
		cerr << "from_string failed on data[" << NUMSPECIES << "]" << endl;
	}
	return row;
}

partposition getPartPosition (string line) {
	int tab=-1;
	int nexttab;
	string data[4];
	for (int i=0; i<4; i++) {
		nexttab = tab+1;
		while (line[nexttab] != ' ' && nexttab < line.size()) {
			nexttab++;
		}
		data[i] = line.substr(tab+1,nexttab-tab-1);
		tab=nexttab;
	}
	partposition row;
	
	row.name=data[0];
	if (!from_string<int>(row.chromosome,data[1],std::dec)) {
		cerr << "from_string failed on data[" << 1 << "]" << endl;
	}
	if (!from_string<int>(row.start,data[2],std::dec)) {
		cerr << "from_string failed on data[" << 2 << "]" << endl;
	}
	row.strand=data[3][0];
	
	return row;
}

int main(int argc, char *argv[]) {
	if (argc < 6) {
		cerr << "Please execute like " argv[0] " dotFile genomepositions arabidopsisInfo medicagoInfo soybeanInfo" << endl;
		exit(1);
	}
	
	ifstream genomeposfile;
	genomeposfile.open(argv[2]);
	if (!genomeposfile) {
		cerr << "Unable to open file " << argv[2] << endl;
		exit(1);
	}
	partposition arabpos[6];
	partposition medipos[13];
	partposition soypos[4];
	string line;
	for (int i=0; i<6; i++) {
		getline(genomeposfile,line);
		arabpos[i] = getPartPosition[i];
	}
	for (int i=0; i<13; i++) {
		getline(genomeposfile,line);
		medipos[i] = getPartPosition[i];
	}
	for (int i=0; i<4; i++) {
		getline(genomeposfile,line);
		soypos[i] = getPartPosition[i];
	}
	genomeposfile.close();
	
	string SpeciesNames[3];
	SpeciesNames[0] = "Arabidopsis";
	SpeciesNames[1] = "Medicago";
	SpeciesNames[2] = "Soybean";
	ifstream InfoFile;
	string line;
	int numsections=3;
	int maxlines=0;
	for (int i=3; i<6; i++) {
		InfoFile.open(argv[i]);
		if (!InfoFile) {
			cerr << "Unable to open file " << argv[i] << endl;
			exit(1);
		}
		int numlines=0;
		while (getline(InfoFile,line)) {
			numlines++;
		}
		if (numlines > maxlines) {
			maxlines = numlines;
		}
		InfoFile.close();
	}
	string SectionNames[3][maxlines/3];
	int SectionStarts[3][maxlines/3];
	int SectionLengths[3][maxlines/3];
	for (int i=0; i<3; i++) {
		for (int j=0; j<maxlines/3; j++) {
			SectionNames[i][j]="";
			SectionStarts[i][j]=0;
		}
	}
	for (int i=0; i<3; i++) {
		InfoFile.open(argv[i+3]);
		int linenum=0;
		while (getline(InfoFile,line)) {
			SectionNames[i][linenum]=line;
			getline(InfoFile,line);
			if (!from_string<int>(SectionStarts[i][linenum],line.substr(10,line.size()-10),std::dec)) {
				cerr << "from_string failed" << endl;
			}
			getline(InfoFile,line);
			if (!from_string<int>(SectionLengths[i][linenum],line.substr(8,line.size()-5),std::dec)) {
				cerr << "from_string failed" << endl;
			}
			linenum++;
		}
		InfoFile.close();
	}
	
	ifstream dotFile;
	while (getline(dotFile,line)) {
		dotline row=getDotLine(line);
		
	}
	
	return 0;
}