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
	//cout << "called on "  << line << endl;
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
	//cout << "returning" << endl;
	return row;
}

void printgffrow(gffrow row) {
	cout << row.sequence << '\t' << row.source << '\t' << row.featureType << '\t' << row.start << '\t' << row.end << '\t' << row.score << '\t' << row.strand << '\t' << row.frame << '\t' << row.group << endl;
}

int main (int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Please execute like " << argv[0] << " gffFile" << endl;
		exit(1);
	}
	
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr << "Unable to open << " << argv[1] << endl;
		exit(1);
	}
	string line;
	int mrnanum=1;
	while (getline(inFile,line)) {
		gffrow row=getGFFrow(line);
		printgffrow(row);
		if (row.featureType == "gene") {
			row.featureType = "mrna";
			row.group = "ID=mrna" + mrnanum;
			printgffrow(row);
			mrnanum++;
		}
	}
	return 0;
}