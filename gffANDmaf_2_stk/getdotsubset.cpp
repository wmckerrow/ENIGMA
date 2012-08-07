#define NUMSPECIES 2
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

struct dotline {
	int starts[NUMSPECIES];
	int length;
	int directions[NUMSPECIES];
	int numinalignment;
};

dotline getDotLine (string line) { //read space for tab...
	//cout << "Getting dotline from " << line << endl;
	int tab=-1;
	int nexttab;
	string data[2*NUMSPECIES+2];
	for (int i=0; i<2*NUMSPECIES+2; i++) {
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
	for (int i=0; i<NUMSPECIES; i++) {
		if (data[NUMSPECIES+1+i][0] == '+') {
			row.directions[i] = 1;
		}
		else {
			row.directions[i]=-1;
		}
	}
	if (!from_string<int>(row.numinalignment,data[2*NUMSPECIES+1],std::dec)) {
		cerr << "from_string failed on data[" << 2*NUMSPECIES+1 << "]" << endl;
	}
	//cout << "Got dotline" << endl;
	return row;
}

int main (int argc, char *argv[]) {
	if (argc < 7) {
		cerr << "Execute like " << argv[0] << " dotfile 1st_start 1st_end 2nd_start 2nd_end minlength" << endl;
		exit(1);
	}
	int first_start=atoi(argv[2]);
	int first_end=atoi(argv[3]);
	int second_start=atoi(argv[4]);
	int second_end=atoi(argv[5]);
	int minlength=atoi(argv[6]);
	ifstream inFile;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr << "Unable to open file " << argv[1] << endl;
		exit(1);
	}
	string line;
	while (getline(inFile,line)) {
		dotline row=getDotLine(line);
		//cout << arab_start << " " << row.starts[0] << " " << arab_end << endl;
		if (row.starts[0] >= first_start && row.starts[0] <= first_end) {
			//cout << "arabmatch" << endl;
			if (row.starts[1] >= second_start && row.starts[1] <= second_end) {
				//cout << "medmatch" << endl;
				if (row.length >= minlength) {
					for (int i=0; i<NUMSPECIES; i++) {
						cout << row.starts[i] << " ";
					}
					cout << row.length << " ";
					for (int i=0; i<NUMSPECIES; i++) {
						if (row.directions[i] == 1) {
							cout << "+ ";
						}
						else {
							cout << "- ";
						}
					}
					cout << row.numinalignment;
					cout << endl;
				}
			}
		}
	}
	return 0;
}