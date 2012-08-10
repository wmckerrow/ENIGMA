/*
 The dotplotter program takes a .maf file and creates a .dot file containing alignment information about two of the genomes described in the .maf file.
 Since in my tests I am concatenating a series of medicago, soybean and arabidopsis sections, dotplotter also peroforms that concatenation.
 Information about this concatenation is given by "info files" (see below.)
 Optionally, with -a it can only give those lines that correspond to all way alignments of every species considered.
 The result is written to standard out.
 
 The arguments for dotplotter are as follows:
 The .maf file, a series of info file pairs, the genome that will be in column 1,
 the genome that will be in column 2, an optional -a.
 The info file pairs are the location of an info file, followed by the species that info file is for.
 
 For expample:
 ./dotplotter someMAF.maf ArabInfoFile Arabidopsis MediInfoFile Medicago SoyInfoFile Soybean Arabidopsis Medicago -a > ArabidopsisMedicagoDotFile.dot
 */

#define MAXALIGNMENT 25
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

struct SpeciesStart {
	string species;
	int start;
	int length;
};

struct AlignmentRow {
	string sequencename;
	int startposition;
	int alignmentlength;
	char direction;
	int sequencelength;
	string alignmentsequence;
};

AlignmentRow extractLine(string line) {
	//cout << "called on "  << line << endl;
	int tab=1;
	int nexttab;
	string data[6];
	for (int i=0; i<6; i++) {
		if (line[tab] != '\t') {
			cerr << "line[" << tab << "] is not a tab" << endl; 
		}
		nexttab = tab+1;
		while (line[nexttab] != '\t' && nexttab < line.size()) {
			nexttab++;
		}
		data[i] = line.substr(tab+1,nexttab-tab-1);
		tab=nexttab;
	}
	AlignmentRow row;
	row.sequencename=data[0];
	if (!from_string<int>(row.startposition,data[1],std::dec)) {
		cerr << "from_string failed on data[1]" << endl;
	}
	if (!from_string<int>(row.alignmentlength,data[2],std::dec)) {
		cerr << "from_string failed on data[2]" << endl;
	}
	row.direction=data[3][0];
	if (!from_string<int>(row.sequencelength,data[4],std::dec)) {
		cerr << "from_string failed on data[4]" << endl;
	}
	row.alignmentsequence=data[5];
	
	return row;
}

int main (int argc, char *argv[]) {
	if (argc < 3) {
		cerr << "Please execute like " << argv[0] << " mafFile pairs xAxisSpecies yAxisSpecies" << endl;
		cerr << "Where pairs look like: InfoFile SpeciesName." << endl;
		exit(1);
	}
	bool limit=0;
	int limitstart;
	int limitend;
	for (int i=0; i<argc; i++) {
		string argument=argv[i];
		if (argument == "-p") {
			limit=1;
			limitstart = atoi(argv[i+1]);
			limitend = atoi(argv[i+2]);
			argc-=3;
		}
	}
	bool alignall = 0;
	for (int i=0; i<argc; i++) {
		string argument=argv[i];
		if (argument == "-a") {
			alignall=1;
			argc-=1;
		}
	}
	string xAxisSpecies=argv[argc-2];
	string yAxisSpecies=argv[argc-1];
	//cout << xAxisSpecies << " " << yAxisSpecies << endl;
	argc-=2;
	string SpeciesNames[(argc-2)/2];
	for (int i=0; i<(argc-2)/2; i++) {
		SpeciesNames[i] = argv[2*i+3];
	}
	ifstream InfoFile;
	string line;
	int numsections=(argc-2)/2;
	int maxlines=0;
	for (int i=0; i<(argc-2)/2; i++) {
		InfoFile.open(argv[2*i+2]);
		if (!InfoFile) {
			cerr << "Unable to open file " << argv[2*i+2] << endl;
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
	string SectionNames[(argc-2)/2][maxlines];
	int SectionStarts[(argc-2)/2][maxlines];
	int SectionLengths[(argc-2)/2][maxlines];
	for (int i=0; i<(argc-2)/2; i++) {
		for (int j=0; j<maxlines; j++) {
			SectionNames[i][j]="";
			SectionStarts[i][j]=0;
			SectionLengths[i][j]=0;
		}
	}
	for (int i=0; i<(argc-2)/2; i++) {
		InfoFile.open(argv[2*i+2]);
		int linenum=0;
		while (getline(InfoFile,line)) {
			SectionNames[i][linenum]=line;
			getline(InfoFile,line);
			if (!from_string<int>(SectionStarts[i][linenum],line.substr(10,line.size()-10),std::dec)) {
				cerr << "from_string failed" << endl;
			}
			getline(InfoFile,line);
			if (!from_string<int>(SectionLengths[i][linenum],line.substr(8,line.size()-13),std::dec)) {
				cerr << "from_string failed" << endl;
			}
			linenum++;
		}
		InfoFile.close();
	}	
	
	/*
	for (int i=0; i<(argc-2)/2; i++) {
		for (int j=0; j<maxlines; j++) {
			if (SectionNames[i][j] != "") {
				cout << SectionNames[i][j] << " starts at " << SpeciesNames[i] << " " << SectionStarts[i][j] << endl;
			}
		}
	}
	 */
	
	string alignment[MAXALIGNMENT];
	ifstream mafFile;
	mafFile.open(argv[1]);
	if (!mafFile) {
		cerr << "Unable to open .maf file " << argv[1] << endl;
		exit(1);
	}
	int alignmentpos=0;
	
	while (getline(mafFile,line)) {
		if (line[0] == 'a') {
			alignmentpos=0;
			for (int i=0; i<MAXALIGNMENT; i++) {
				alignment[i]="";
			}
		}
		if (line[0] == 's') {
			if (alignmentpos >= MAXALIGNMENT) {
				cerr << "MAXALIGMNET is too small." <<endl;
			}
			alignment[alignmentpos]=line;
			alignmentpos++;
		}
		if (line == "" && alignment[0] != "") {
			bool gotSpecies[(argc-2)/2];
			for (int i=0; i<(argc-2)/2; i++) {
				gotSpecies[i]=0;
			}
			for (int i=0; i<MAXALIGNMENT; i++) {
				if (alignment[i] != "") {
					AlignmentRow row=extractLine(alignment[i]);
					for (int k=0; k<(argc-2)/2; k++) {
						for (int l=0; l<maxlines; l++) {
							if (row.sequencename == SectionNames[k][l]) {
								gotSpecies[k]=1;
							}
						}
					}
				}
			}
			bool gotAllSpecies = 1;
			int numinalignment = 0;
			for (int i=0; i<(argc-2)/2; i++) {
				if (!gotSpecies[i]) {
					gotAllSpecies=0;
				}
				else {
					numinalignment++;
				}
			}
			if (gotAllSpecies || !alignall) {
				for (int i=0; i<MAXALIGNMENT; i++) {
					for (int j=i; j<MAXALIGNMENT; j++) {
						if (alignment[i] != "" && alignment[j] != "") {
							AlignmentRow row1=extractLine(alignment[i]);
							AlignmentRow row2=extractLine(alignment[j]);
							SpeciesStart ss1;
							for (int k=0; k<(argc-2)/2; k++) {
								for (int l=0; l<maxlines; l++) {
									if (row1.sequencename == SectionNames[k][l]) {
										ss1.species=SpeciesNames[k];
										ss1.start=SectionStarts[k][l];
										ss1.length=SectionLengths[k][l];
									}
								}
							}
							SpeciesStart ss2;
							for (int k=0; k<(argc-2)/2; k++) {
								for (int l=0; l<maxlines; l++) {
									if (row2.sequencename == SectionNames[k][l]) {
										ss2.species=SpeciesNames[k];
										ss2.start=SectionStarts[k][l];
										ss2.length=SectionLengths[k][l];
									}
								}
							}
							if (ss1.species == xAxisSpecies && ss2.species == yAxisSpecies) {
								if (row1.direction == '+') {
									cout << ss1.start + row1.startposition << " ";
								}
								else {
									cout << ss1.start + ss1.length - row1.startposition -1 << " ";
								}
								if (row2.direction == '+') {
									cout << ss2.start + row2.startposition << " ";
								}
								else {
									cout << ss2.start + ss2.length - row2.startposition -1 << " ";
								}
								cout << row1.alignmentlength << " " << row1.direction << " " <<row2.direction << " " << numinalignment << endl;
							}
							if (ss2.species == xAxisSpecies && ss1.species == yAxisSpecies) {
								if (row2.direction == '+') {
									cout << ss2.start + row2.startposition << " ";
								}
								else {
									cout << ss2.start + ss2.length - row2.startposition -1 << " ";
								}
								if (row1.direction == '+') {
									cout << ss1.start + row1.startposition << " ";
								}
								else {
									cout << ss1.start + ss1.length - row1.startposition -1 << " ";
								}
								cout << row1.alignmentlength << " " << row2.direction << " " <<row1.direction << " " << numinalignment << endl;
							}
						}
					}
				}
			}
		}
	}
	
	bool gotSpecies[(argc-2)/2];
	for (int i=0; i<(argc-2)/2; i++) {
		gotSpecies[i]=0;
	}
	for (int i=0; i<MAXALIGNMENT; i++) {
		if (alignment[i] != "") {
			AlignmentRow row=extractLine(alignment[i]);
			for (int k=0; k<(argc-2)/2; k++) {
				for (int l=0; l<maxlines; l++) {
					if (row.sequencename == SectionNames[k][l]) {
						gotSpecies[k]=1;
					}
				}
			}
		}
	}
	bool gotAllSpecies = 1;
	int numinalignment = 0;
	for (int i=0; i<(argc-2)/2; i++) {
		if (!gotSpecies[i]) {
			gotAllSpecies=0;
		}
		else {
			numinalignment++;
		}

	}
	
	if (gotAllSpecies || !alignall) {
		for (int i=0; i<MAXALIGNMENT; i++) {
			for (int j=i; j<MAXALIGNMENT; j++) {
				if (alignment[i] != "" && alignment[j] != "") {
					AlignmentRow row1=extractLine(alignment[i]);
					AlignmentRow row2=extractLine(alignment[j]);
					SpeciesStart ss1;
					for (int k=0; k<(argc-2)/2; k++) {
						for (int l=0; l<maxlines; l++) {
							if (row1.sequencename == SectionNames[k][l]) {
								ss1.species=SpeciesNames[k];
								ss1.start=SectionStarts[k][l];
								ss1.length=SectionLengths[k][l];
							}
						}
					}
					SpeciesStart ss2;
					for (int k=0; k<(argc-2)/2; k++) {
						for (int l=0; l<maxlines; l++) {
							if (row2.sequencename == SectionNames[k][l]) {
								ss2.species=SpeciesNames[k];
								ss2.start=SectionStarts[k][l];
								ss2.length=SectionLengths[k][l];
							}
						}
					}
					if (ss1.species == xAxisSpecies && ss2.species == yAxisSpecies) {
						if (row1.direction == '+') {
							cout << ss1.start + row1.startposition << " ";
						}
						else {
							cout << ss1.start + ss1.length - row1.startposition -1 << " ";
						}
						if (row2.direction == '+') {
							cout << ss2.start + row2.startposition << " ";
						}
						else {
							cout << ss2.start + ss2.length - row2.startposition -1 << " ";
						}
						cout << row1.alignmentlength << " " << row1.direction << " " << row2.direction << " " << numinalignment << endl;
					}
					if (ss2.species == xAxisSpecies && ss1.species == yAxisSpecies) {
						if (row2.direction == '+') {
							cout << ss2.start + row2.startposition << " ";
						}
						else {
							cout << ss2.start + ss2.length - row2.startposition -1 << " ";
						}
						if (row1.direction == '+') {
							cout << ss1.start + row1.startposition << " ";
						}
						else {
							cout << ss1.start + ss1.length - row1.startposition -1 << " ";
						}
						cout << row1.alignmentlength << " " << row2.direction << " " <<row1.direction << " " << numinalignment << endl;
					}
				}
			}
		}
	}
}