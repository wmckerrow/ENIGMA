#define MAXALIGN 10
#define NUMSPECIES 2
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

int debugint=0;

struct intchar {
	int integer;
	char character;
};

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

class dotline {
public:
	int starts[NUMSPECIES];
	int length;
	int directions[NUMSPECIES];
	int bottom(int);
	int top(int);
	void shortenbottom(int,int);
	void shortentop(int,int);
};

int dotline::bottom(int column) {
	return starts[column] - isminus(directions[column])*(length-1);
}

int dotline::top(int column) {
	return starts[column] + isplus(directions[column])*(length-1);
}

void dotline::shortenbottom(int shortenby, int column) {
	if (directions[column] == -1) {
		for (int i=0; i<NUMSPECIES; i++) {
			starts[i] -= directions[i]*directions[column];
		}
	}
	length=length-shortenby;
}

void dotline::shortentop(int shortenby, int column) {
	if (directions[column] == 1) {
		for (int i=0; i<NUMSPECIES; i++) {
			starts[i] += directions[i]*directions[column];
		}
	}
	length=length-shortenby;
}

string convertInt(int number)
{
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
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

dotline getDotLine (string line) { //read space for tab...
	int numspaces=0;
	for (int i=0; i<line.length(); i++) {
		if (line[i]==' ') {
			numspaces++;
		}
	}
	if (numspaces < 2*NUMSPECIES) {
		cerr << line << endl << "is not a proper dotline" << endl;
		exit(1);
	}
	
	int tab=-1;
	int nexttab;
	string data[2*NUMSPECIES+1];
	for (int i=0; i<2*NUMSPECIES+1; i++) {
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
	//cout << "Got dotline" << endl;
	return row;
}

/*
struct dotlistelem {
	dotline data;
	dotlistelem *next;
};
 */

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

int oldfindalignment(int position, int alignmentFrom, int alignmentTo, dotline dots[], int numdots, int fudgeFactor) { //If we have >1 alignments we are picking the longest one
	int alignment=-1;
	int score=0;
	for (int i=0; i<numdots; i++) {
		if (position >= dots[i].starts[alignmentFrom] - isminus(dots[i].directions[alignmentFrom])*(dots[i].length-1)) {
			if (position <= dots[i].starts[alignmentFrom] + isplus(dots[i].directions[alignmentFrom])*(dots[i].length-1)) {
				if (dots[i].length > score) {
					int distancefromstart = dots[i].directions[alignmentFrom]*(position-dots[i].starts[alignmentFrom]);
					alignment = dots[i].starts[alignmentTo] + dots[i].directions[alignmentTo]*distancefromstart;
					score=dots[i].length;
					if (fudgeFactor == 2) {
						cout << dots[i].starts[0] << " " << dots[i].starts[1] << " " << dots[i].length << " " << dots[i].directions[0] << " " << dots[i].directions[1] << endl;
					}
				}
			}
		}
	}
	return alignment;
}

int findalignment(int position, int alignmentFrom, int alignmentTo, dotline dot) { //If we have >1 alignments we are picking the longest one
	int alignment=-1;
	if (position >= dot.bottom(alignmentFrom)) {
		if (position <= dot.top(alignmentFrom)) {
			int distancefromstart = dot.directions[alignmentFrom]*(position-dot.starts[alignmentFrom]);
			alignment = dot.starts[alignmentTo] + dot.directions[alignmentTo]*distancefromstart;
		}
	}
	return alignment;
}

int finddirection (int position, int alignmentFrom, int alignmentTo, dotline dots[], int numdots, int fudgeFactor) {
	int direction=1;
	int score=0;
	for (int i=0; i<numdots; i++) {
		if (position >= dots[i].starts[alignmentFrom] - isminus(dots[i].directions[alignmentFrom])*dots[i].length) {
			if (position < dots[i].starts[alignmentFrom] + isplus(dots[i].directions[alignmentFrom])*dots[i].length) {
				if (dots[i].length > score) {
					direction=dots[i].directions[alignmentFrom]*dots[i].directions[alignmentTo];
					score=dots[i].length;
				}
			}
		}
	}
	return direction;
}

char findletter(int position, gffrow *gff, int numgffrows) {
	if (position < 1) {
		return '*';
	}
	bool geneCDSbefore=0;
	bool geneCDSafter=0;
	bool inthisgene=0;
	char strand;
	for (int i=0; i<numgffrows; i++) {
		if (gff[i].featureType == "gene") {
			inthisgene=0;
			geneCDSbefore=0;
			geneCDSafter=0;
		}
		if (gff[i].start <= position && position <= gff[i].end) {
			if (gff[i].featureType == "CDS") {
				if (gff[i].strand == '+') {
					return 'e';
				}
				else {
					return 'E';
				}
			}
			if (gff[i].featureType == "gene") {
				inthisgene = 1;
				strand = gff[i].strand;
			}
		}
		if (inthisgene && position > gff[i].end && gff[i].featureType == "CDS") {
			geneCDSbefore=1;
		}
		if (inthisgene && position < gff[i].start && gff[i].featureType == "CDS") {
			geneCDSafter=1;
		}
		if (geneCDSbefore && geneCDSafter) {
			if (strand == '+') {
				return 'i';
			}
			else {
				return 'I';
			}
		}
	}
	return 'x';
}

char findletterfromintchar (int position, intchar *letters, int numintchars, int direction) {
	if (position < 1) {
		return '*';
	}
	for (int i=0; i<numintchars-1; i++) {
		if (position < letters[i+1].integer) {
			char theletter=letters[i].character;
			while (direction == -1) {
				if (theletter == 'e') {
					theletter = 'E';
					break;
				}
				if (theletter == 'f') {
					theletter = 'F';
					break;
				}
				if (theletter == 'g') {
					theletter = 'G';
					break;
				}
				if (theletter == 'i') {
					theletter = 'I';
					break;
				}
				if (theletter == 'j') {
					theletter = 'J';
					break;
				}
				if (theletter == 'k') {
					theletter = 'K';
					break;
				}
				if (theletter == 'E') {
					theletter = 'e';
					break;
				}
				if (theletter == 'F') {
					theletter = 'f';
					break;
				}
				if (theletter == 'G') {
					theletter = 'g';
					break;
				}
				if (theletter == 'I') {
					theletter = 'i';
					break;
				}
				if (theletter == 'J') {
					theletter = 'j';
					break;
				}
				if (theletter == 'K') {
					theletter = 'k';
					break;
				}
				break;
			}
			if (direction == 1 && position==letters[i].integer) {
				if (theletter == 'e' || theletter == 'f' || theletter == 'g') {
					if (i==0) {
						return 's';
					}
					if (letters[i-1].character=='x') {
						return 's';
					}
					if (letters[i-1].character=='i') {
						return 'a';
					}
					if (letters[i-1].character=='j') {
						return 'b';
					}
					if (letters[i-1].character=='k') {
						return 'c';
					}
				}
				if (theletter == 'E' || theletter == 'F' || theletter == 'G') {
					if (i==0) {
						return 'T';
					}
					if (letters[i-1].character=='x') {
						return 'T';
					}
					if (letters[i-1].character=='I') {
						return 'D';
					}
					if (letters[i-1].character=='J') {
						return 'O';
					}
					if (letters[i-1].character=='K') {
						return 'N';
					}
				}
				if (theletter=='i' || theletter=='j' || theletter=='k') {
					if (i == 0 || letters[i-1].character == 'e') {
						return 'd';
					}
					if (i == 0 || letters[i-1].character == 'f') {
						return 'o';
					}
					if (i == 0 || letters[i-1].character == 'g') {
						return 'n';
					}
				}
				if (theletter=='I' || theletter=='J' || theletter=='K') {
					if (i == 0 || letters[i-1].character == 'E') {
						return 'A';
					}
					if (i == 0 || letters[i-1].character == 'F') {
						return 'B';
					}
					if (i == 0 || letters[i-1].character == 'G') {
						return 'C';
					}
				}
				if (theletter == 'x') {
					if (i == 0 || letters[i-1].character == 'e' || letters[i-1].character == 'f' || letters[i-1].character == 'g') {
						return 't';
					}
				}
				if (theletter == 'x') {
					if (i == 0 || letters[i-1].character == 'E' || letters[i-1].character == 'F' || letters[i-1].character == 'G') {
						return 'S';
					}
				}
			}
			if (i != numintchars-1) {
				if (direction == -1 && position==letters[i+1].integer-1) {
					if (theletter == 'e' || theletter == 'f' || theletter == 'g') {
						if (letters[i+1].character=='x') {
							return 's';
						}
						if (letters[i+1].character=='I') {
							return 'a';
						}
						if (letters[i+1].character=='J') {
							return 'b';
						}
						if (letters[i+1].character=='K') {
							return 'c';
						}
					}
					if (theletter == 'E' || theletter == 'F' || theletter == 'G') {
						if (letters[i+1].character=='x') {
							return 'T';
						}
						if (letters[i+1].character=='i') {
							return 'D';
						}
						if (letters[i+1].character=='j') {
							return 'O';
						}
						if (letters[i+1].character=='k') {
							return 'N';
						}
					}
					if (theletter=='i' || theletter=='j' || theletter=='k') {
						if (letters[i+1].character == 'E') {
							return 'd';
						}
					}
					if (theletter=='I') {
						if (letters[i+1].character == 'e' || letters[i+1].character == 'f' || letters[i+1].character == 'g') {
							return 'A';
						}
					}
					if (theletter=='J') {
						if (letters[i+1].character == 'e' || letters[i+1].character == 'f' || letters[i+1].character == 'g') {
							return 'B';
						}
					}
					if (theletter=='K') {
						if (letters[i+1].character == 'e' || letters[i+1].character == 'f' || letters[i+1].character == 'g') {
							return 'C';
						}
					}
					if (theletter == 'x') {
						if (letters[i+1].character == 'E' || letters[i+1].character == 'F' || letters[i+1].character == 'G') {
							return 't';
						}
					}
					if (theletter == 'x') {
						if (letters[i+1].character == 'e' || letters[i+1].character == 'f' || letters[i+1].character == 'g') {
							return 'S';
						}
					}
				}
			}
			return theletter;
		}
	}
	return letters[numintchars-1].character;
}

char findletterfromintchar_nostadn (int position, intchar *letters, int numintchars, int direction) {
	if (position < 1) {
		return '*';
	}
	for (int i=0; i<numintchars-1; i++) {
		if (position < letters[i+1].integer) {
			char theletter=letters[i].character;
			while (direction == -1) {
				if (theletter == 'e') {
					theletter = 'E';
					break;
				}
				if (theletter == 'f') {
					theletter = 'F';
					break;
				}
				if (theletter == 'g') {
					theletter = 'G';
					break;
				}
				if (theletter == 'i') {
					theletter = 'I';
					break;
				}
				if (theletter == 'j') {
					theletter = 'J';
					break;
				}
				if (theletter == 'k') {
					theletter = 'K';
					break;
				}
				if (theletter == 'E') {
					theletter = 'e';
					break;
				}
				if (theletter == 'F') {
					theletter = 'f';
					break;
				}
				if (theletter == 'G') {
					theletter = 'g';
					break;
				}
				if (theletter == 'I') {
					theletter = 'i';
					break;
				}
				if (theletter == 'J') {
					theletter = 'j';
					break;
				}
				if (theletter == 'K') {
					theletter = 'k';
					break;
				}
				break;
			}
			return theletter;
		}
	}
	return letters[numintchars-1].character;
}

bool is_stad (char character) {
	if (character =='s') {
		return 1;
	}
	if (character == 't') {
		return 1;
	}
	if (character == 'a') {
		return 1;
	}
	if (character == 'b') {
		return 1;
	}
	if (character == 'c') {
		return 1;
	}
	if (character == 'd') {
		return 1;
	}
	if (character == 'o') {
		return 1;
	}
	if (character == 'n') {
		return 1;
	}
	if (character =='T') {
		return 1;
	}
	if (character == 'S') {
		return 1;
	}
	if (character == 'D') {
		return 1;
	}
	if (character == 'O') {
		return 1;
	}
	if (character == 'N') {
		return 1;
	}
	if (character == 'A') {
		return 1;
	}
	if (character == 'B') {
		return 1;
	}
	if (character == 'C') {
		return 1;
	}
	return 0;
}

bool is_cds (char character) {
	if (character =='e') {
		return 1;
	}
	if (character == 'f') {
		return 1;
	}
	if (character == 'g') {
		return 1;
	}
	if (character == 'E') {
		return 1;
	}
	if (character == 'F') {
		return 1;
	}
	if (character == 'G') {
		return 1;
	}
	return 0;
}

/*
bool is_sa (char character) {
	if (character =='s') {
		return 1;
	}
	if (character == 'a') {
		return 1;
	}
	return 0;
}
 */

/*
bool is_td (char character) {
	if (character =='t') {
		return 1;
	}
	if (character == 'd') {
		return 1;
	}
	return 0;
}
 */

int nextlettertrans (int position, intchar *letters, int numintchars) {
	for (int i=0; i<numintchars-1; i++) {
		if (position < letters[i+1].integer) {
			return letters[i+1].integer;
		}
	}
	return 0;
}

int prevlettertrans (int position, intchar *letters, int numintchars) {
	for (int i=numintchars-1; i>0; i--) {
		if (position > letters[i].integer-1) {
			return letters[i].integer-1;
		}
	}
	return 0;
}

void removerow(dotline *removefrom, int removethis, int dotsize) {
	for (int i=removethis; i<dotsize-1; i++) {
		removefrom[i]=removefrom[i+1];
	}
}

int removeoverlaps(dotline *sortthis, int dotsize, int targetpos) {
	bool removed=1;
	while (removed) {
		removed=0;
		int i=1;
		while (i<dotsize) {
			if (sortthis[i-1].top(targetpos) >= sortthis[i].bottom(targetpos)) {
				if (sortthis[i-1].top(targetpos) >= sortthis[i].top(targetpos)) {
					removerow(sortthis,i,dotsize);
					dotsize--;
				}
				else {
					if (sortthis[i-1].length > sortthis[i].length) {
						sortthis[i].shortenbottom(sortthis[i-1].top(targetpos)-sortthis[i].bottom(targetpos)+1,targetpos);
					}
					else {
						sortthis[i-1].shortentop(sortthis[i-1].top(targetpos)-sortthis[i].bottom(targetpos)+1,targetpos);
					}
					i++;
				}
			}
			else {
				i++;
			}
		}
	}
	return dotsize;
}

void sortdots(dotline *sortthis, int dotsize, int targetpos, int alignedpos) {
	bool switched=1;
	while (switched) {
		switched=0;
		for (int i=1; i<dotsize; i++) {
			if (sortthis[i-1].bottom(targetpos) > sortthis[i].bottom(targetpos) || (sortthis[i-1].bottom(targetpos) == sortthis[i].bottom(targetpos) && sortthis[i-1].bottom(alignedpos) > sortthis[i].bottom(alignedpos)) ) {
				dotline temp=sortthis[i-1];
				sortthis[i-1]=sortthis[i];
				sortthis[i]=temp;
				switched=1;
			}
		}
	}
}

void sortintchars (intchar *sortthis, int sizetosort) {
	bool switched=1;
	while (switched) {
		switched=0;
		for (int i=1; i<sizetosort; i++) {
			if (sortthis[i-1].integer > sortthis[i].integer) {
				intchar temp=sortthis[i-1];
				sortthis[i-1]=sortthis[i];
				sortthis[i]=temp;
				switched=1;
			}
		}
	}
}

void removeintchar (intchar *removefrom, int removethis, int intcharsize) {
	for (int i=removethis; i<intcharsize-1; i++) {
		removefrom[i]=removefrom[i+1];
	}
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
		cerr << "Please execute like " << argv[0] << " numAlignedSpecies alignedSpeciesInfo targetName numgffs4target gffs4target positionOutputFile sectionsFile tree" << endl;
		cerr << "alignedSpeciesInfo should look like SpeciesName dotFile targetPosInDotFile alignedPosInDotFile gffFile" << endl;
		//cerr << "gffs4targetInfo should look like gffFileName evidenceType." << endl;
		//cerr << "Types are 1: Full evicence; 2: No reading frame." << endl;
		exit(1);
	}
	
	string line;
	int numAlignedSpecies=atoi(argv[1]);
	
	/*
	if (numAlignedSpecies > MAXALIGN) {
		cerr << "MAXALIGN is " << MAXALIGN << ". Recompile with a value of atleast " << numAlignedSpecies << "." << endl;
		exit(1);
	}
	 */
	
	string alignedNames[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		alignedNames[i]=argv[2+5*i];
	}
	
	int targetPosInDotFile[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		targetPosInDotFile[i] = atoi(argv[4+5*i]);
	}
	
	int posInDotFile[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		posInDotFile[i] = atoi(argv[5+5*i]);
	}
	
	cerr << "Reading and sorting MAF blocks." << endl;
	
	ifstream dotFile;
	dotline **alignments = new dotline*[numAlignedSpecies];
	int numdotlines[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		dotFile.open(argv[3+5*i]);
		if (!dotFile) {
			cerr << "Unable to open file " << argv[3+5*i] << endl;
		}
		numdotlines[i]=0;
		while (getline(dotFile,line)) {
			numdotlines[i]++;
		}
		dotFile.close();
		dotFile.open(argv[3+5*i]);
		alignments[i] = new dotline[numdotlines[i]];
		for (int j=0; j<numdotlines[i]; j++) {
			getline(dotFile,line);
			alignments[i][j] = getDotLine(line);
		}
		dotFile.close();
		sortdots(alignments[i],numdotlines[i],targetPosInDotFile[i],posInDotFile[i]);
		numdotlines[i]=removeoverlaps(alignments[i],numdotlines[i],targetPosInDotFile[i]);
		
		
		ofstream sorteddot;
		string filename = "sortedMAFs" + convertInt(i);
		sorteddot.open(filename.c_str());
		for (int j=0; j<numdotlines[i]; j++) {
			sorteddot << alignments[i][j].starts[0] << " " << alignments[i][j].starts[1] << " " << alignments[i][j].length << " " << alignments[i][j].directions[0] << " " << alignments[i][j].directions[1] << endl;
		}
	}
	
	string targetName=argv[2+5*numAlignedSpecies];
	int numgffs4target=atoi(argv[3+5*numAlignedSpecies]);
	
	cerr << "Reading GFF files for related species." << endl;
	
	ifstream gffFile;
	int thisgfflines[numAlignedSpecies];
	gffrow **alignedGFFs;
	intchar *alignedGFFletters[numAlignedSpecies];
	alignedGFFs = new gffrow*[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		thisgfflines[i]=0;
		gffFile.open(argv[6+5*i]);
		while (getline(gffFile,line)) {
			thisgfflines[i]++;
		}
		gffFile.close();
		alignedGFFs[i] = new gffrow[thisgfflines[i]];
		alignedGFFletters[i] = new intchar[2*thisgfflines[i]+1];
		gffFile.open(argv[6+5*i]);
		for (int j=0; j<thisgfflines[i]; j++) {
			getline(gffFile,line);
			alignedGFFs[i][j]=getGFFrow(line);
			alignedGFFletters[i][2*j+1].integer=alignedGFFs[i][j].start;
			alignedGFFletters[i][2*j+2].integer=alignedGFFs[i][j].end+1;
		}
		gffFile.close();
	}
	
	cerr << "Getting letters." << endl;
	for (int i=0; i<numAlignedSpecies; i++) {
		alignedGFFletters[i][0].integer=1;
		alignedGFFletters[i][0].character='x';
		for (int j=1; j<2*thisgfflines[i]+1; j++) {
			alignedGFFletters[i][j].character = findletter(alignedGFFletters[i][j].integer, alignedGFFs[i], thisgfflines[i]);
		}
	}
	
	/*
	for (int i=0; i<numAlignedSpecies; i++) {
		for (int j=0; j<2*thisgfflines[i]+1; j++) {
			cout << alignedGFFletters[i][j].integer << " " << alignedGFFletters[i][j].character << endl;
		}
		cout << endl;
	}
	 */
	
	cerr << "Sorting letters." << endl;
	for (int i=0; i<numAlignedSpecies; i++) {
		sortintchars (alignedGFFletters[i], 2*thisgfflines[i]+1);
	}
	
	/*
	for (int i=0; i<numAlignedSpecies; i++) {
		for (int j=0; j<2*thisgfflines[i]+1; j++) {
			cout << alignedGFFletters[i][j].integer << " " << alignedGFFletters[i][j].character << endl;
		}
		cout << endl;
	}
	 */
	
	cerr << "Removing duplicates." << endl;
	int numlettertrans[numAlignedSpecies];
	for (int i=0; i < numAlignedSpecies; i++) {
		numlettertrans[i] = 2*thisgfflines[i]+1;
		int thisintchar=1;
		while (thisintchar < numlettertrans[i]) {
			while (alignedGFFletters[i][thisintchar].character == alignedGFFletters[i][thisintchar-1].character && thisintchar < numlettertrans[i]) {
				removeintchar (alignedGFFletters[i], thisintchar, numlettertrans[i]);
				numlettertrans[i]--;
			}
			thisintchar++;
		}
	}
	
	cerr << "Copying for reverse alignments." << endl;
	intchar *reversealignedGFFletters[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		reversealignedGFFletters[i] = new intchar[numlettertrans[i]];
		for (int j=0; j<numlettertrans[i]; j++) {
			reversealignedGFFletters[i][j] = alignedGFFletters[i][j];
		}
	}
	
	/*
	for (int i=0; i<numAlignedSpecies; i++) {
		for (int j=0; j<numlettertrans[i]; j++) {
			cout << alignedGFFletters[i][j].integer << " " << alignedGFFletters[i][j].character << endl;
		}
		cout << endl;
	}
	 //*/
	
	cerr << "Adding reading frame." << endl;
	for (int i=0; i<numAlignedSpecies; i++) {
		for (int j=1; j<numlettertrans[i]; j++) {
			if (alignedGFFletters[i][j].character == 'e' && alignedGFFletters[i][j-1].character == 'j') {
				alignedGFFletters[i][j].character = 'f';
			}
			if (alignedGFFletters[i][j].character == 'e' && alignedGFFletters[i][j-1].character == 'k') {
				alignedGFFletters[i][j].character = 'g';
			}
			
			if (alignedGFFletters[i][j].character == 'i') {
				int frame = alignedGFFletters[i][j].integer - alignedGFFletters[i][j-1].integer;
				if (alignedGFFletters[i][j-1].character == 'g') {
					frame--;
				}
				if (alignedGFFletters[i][j-1].character == 'f') {
					frame-=2;
				}
				frame = frame % 3;
				if (frame == 1) {
					alignedGFFletters[i][j].character = 'j';
				}
				if (frame == 2) {
					alignedGFFletters[i][j].character = 'k';
				}
			}
			
			if (alignedGFFletters[i][j].character == 'E' && alignedGFFletters[i][j-1].character == 'J') {
				alignedGFFletters[i][j].character = 'F';
			}
			if (alignedGFFletters[i][j].character == 'E' && alignedGFFletters[i][j-1].character == 'K') {
				alignedGFFletters[i][j].character = 'G';
			}
			
			if (alignedGFFletters[i][j].character == 'I') {
				int frame = alignedGFFletters[i][j].integer - alignedGFFletters[i][j-1].integer;
				if (alignedGFFletters[i][j-1].character == 'G') {
					frame--;
				}
				if (alignedGFFletters[i][j-1].character == 'F') {
					frame-=2;
				}
				frame = frame % 3;
				if (frame == 1) {
					alignedGFFletters[i][j].character = 'J';
				}
				if (frame == 2) {
					alignedGFFletters[i][j].character = 'K';
				}
			}
		}
	}
	
	cerr << "Adding reading frame for reverse aligned." << endl;
	for (int i=0; i<numAlignedSpecies; i++) {
		for (int j=numlettertrans[i]-2; j>=0; j--) {
			if (reversealignedGFFletters[i][j].character == 'e' && reversealignedGFFletters[i][j+1].character == 'j') {
				reversealignedGFFletters[i][j].character = 'f';
			}
			if (reversealignedGFFletters[i][j].character == 'e' && reversealignedGFFletters[i][j+1].character == 'k') {
				reversealignedGFFletters[i][j].character = 'g';
			}
			
			if (reversealignedGFFletters[i][j].character == 'i') {
				int frame = reversealignedGFFletters[i][j+2].integer - reversealignedGFFletters[i][j+1].integer;
				if (reversealignedGFFletters[i][j+1].character == 'g') {
					frame--;
				}
				if (reversealignedGFFletters[i][j+1].character == 'f') {
					frame-=2;
				}
				frame = frame % 3;
				if (frame == 1) {
					reversealignedGFFletters[i][j].character = 'j';
				}
				if (frame == 2) {
					reversealignedGFFletters[i][j].character = 'k';
				}
			}
			
			if (reversealignedGFFletters[i][j].character == 'E' && reversealignedGFFletters[i][j+1].character == 'J') {
				reversealignedGFFletters[i][j].character = 'F';
			}
			if (reversealignedGFFletters[i][j].character == 'E' && reversealignedGFFletters[i][j+1].character == 'K') {
				reversealignedGFFletters[i][j].character = 'G';
			}
			
			if (reversealignedGFFletters[i][j].character == 'I') {
				int frame = reversealignedGFFletters[i][j+2].integer - reversealignedGFFletters[i][j+1].integer;
				if (reversealignedGFFletters[i][j+1].character == 'G') {
					frame--;
				}
				if (reversealignedGFFletters[i][j+1].character == 'F') {
					frame-=2;
				}
				frame = frame % 3;
				if (frame == 1) {
					reversealignedGFFletters[i][j].character = 'J';
				}
				if (frame == 2) {
					reversealignedGFFletters[i][j].character = 'K';
				}
			}
		}
	}
	
	/*
	for (int i=0; i<numAlignedSpecies; i++) {
		for (int j=0; j<numlettertrans[i]; j++) {
			cout << alignedGFFletters[i][j].integer << " " << alignedGFFletters[i][j].character << endl;
		}
		cout << endl;
	}
	for (int i=0; i<numAlignedSpecies; i++) {
		for (int j=0; j<numlettertrans[i]; j++) {
			cout << reversealignedGFFletters[i][j].integer << " " << reversealignedGFFletters[i][j].character << endl;
		}
		cout << endl;
	}
	//*/
	
	cerr << "Reading GFFs for target." << endl;
	gffrow **targetGFFs;
	intchar *targetGFFletters[numgffs4target];
	targetGFFs = new gffrow*[numgffs4target];
	int targetGFFsize[numgffs4target];
	for (int i=0; i<numgffs4target; i++) {
		gffFile.open(argv[4+5*numAlignedSpecies+i]);
		if (!gffFile) {
			cerr << "Unable to open gffFile " << argv[4+5*numAlignedSpecies+i] << endl;
			exit(1);
		}
		targetGFFsize[i]=0;
		while (getline(gffFile,line)) {
			targetGFFsize[i]++;
		}
		gffFile.close();
		targetGFFs[i]=new gffrow[targetGFFsize[i]];
		targetGFFletters[i] = new intchar[2*targetGFFsize[i]+1];
		gffFile.open(argv[4+5*numAlignedSpecies+i]);
		for (int j=0; j<targetGFFsize[i]; j++) {
			getline(gffFile,line);
			targetGFFs[i][j]=getGFFrow(line);
			targetGFFletters[i][2*j+1].integer=targetGFFs[i][j].start;
			targetGFFletters[i][2*j+2].integer=targetGFFs[i][j].end+1;
		}
		gffFile.close();
	}
	
	
	cerr << "Getting letters." << endl;
	for (int i=0; i<numgffs4target; i++) {
		targetGFFletters[i][0].integer=1;
		targetGFFletters[i][0].character='x';
		for (int j=1; j<2*targetGFFsize[i]+1; j++) {
			targetGFFletters[i][j].character = findletter(targetGFFletters[i][j].integer, targetGFFs[i], targetGFFsize[i]);
		}
	}
	
	/*
	 for (int i=0; i<numAlignedSpecies; i++) {
	 for (int j=0; j<2*thisgfflines[i]+1; j++) {
	 cout << alignedGFFletters[i][j].integer << " " << alignedGFFletters[i][j].character << endl;
	 }
	 cout << endl;
	 }
	 */
	
	cerr << "Sorting letters." << endl;
	for (int i=0; i<numgffs4target; i++) {
		sortintchars (targetGFFletters[i], 2*targetGFFsize[i]+1);
	}
	
	/*
	 for (int i=1; i<numAlignedSpecies; i++) {
	 for (int j=0; j<2*thisgfflines[i]+1; j++) {
	 cout << alignedGFFletters[i][j].integer << " " << alignedGFFletters[i][j].character << endl;
	 }
	 cout << endl;
	 }
	 // */
	
	cerr << "Removing duplicates." << endl;
	int numtargetlettertrans[numgffs4target];
	for (int i=0; i < numgffs4target; i++) {
		numtargetlettertrans[i] = 2*targetGFFsize[i]+1;
		int thisintchar=1;
		while (thisintchar < numtargetlettertrans[i]) {
			while (targetGFFletters[i][thisintchar].character == targetGFFletters[i][thisintchar-1].character && thisintchar < numtargetlettertrans[i]) {
				removeintchar (targetGFFletters[i], thisintchar, numtargetlettertrans[i]);
				numtargetlettertrans[i]--;
			}
			thisintchar++;
		}
	}
	
	/*
	for (int i=0; i<numgffs4target; i++) {
		for (int j=0; j<numtargetlettertrans[i]; j++) {
			cout << targetGFFletters[i][j].integer << " " << targetGFFletters[i][j].character << endl;
		}
		cout << endl;
	}
	 //*/
	
	
	/*
	cerr << "Adding reading frame." << endl;
	for (int i=0; i<numgffs4target; i++) {
		for (int j=0; j<numtargetlettertrans[i]; j++) {
			if (targetGFFletters[i][j].character == 'e' && targetGFFletters[i][j-1].character == 'j') {
				targetGFFletters[i][j].character = 'f';
			}
			if (targetGFFletters[i][j].character == 'e' && targetGFFletters[i][j-1].character == 'k') {
				targetGFFletters[i][j].character = 'g';
			}
			
			if (targetGFFletters[i][j].character == 'i') {
				int frame = targetGFFletters[i][j].integer - targetGFFletters[i][j-1].integer;
				if (targetGFFletters[i][j-1].character == 'f') {
					frame--;
				}
				if (targetGFFletters[i][j-1].character == 'g') {
					frame-=2;
				}
				frame = frame % 3;
				if (frame == 1) {
					targetGFFletters[i][j].character = 'k';
				}
				if (frame == 2) {
					targetGFFletters[i][j].character = 'j';
				}
			}
			
			if (targetGFFletters[i][j].character == 'E' && targetGFFletters[i][j-1].character == 'J') {
				targetGFFletters[i][j].character = 'F';
			}
			if (targetGFFletters[i][j].character == 'E' && targetGFFletters[i][j-1].character == 'K') {
				targetGFFletters[i][j].character = 'G';
			}
			
			
			if (targetGFFletters[i][j].character == 'I') {
				int frame = targetGFFletters[i][j].integer - targetGFFletters[i][j-1].integer;
				if (targetGFFletters[i][j-1].character == 'F') {
					frame--;
				}
				if (targetGFFletters[i][j-1].character == 'G') {
					frame-=2;
				}
				frame = frame % 3;
				if (frame == 1) {
					targetGFFletters[i][j].character = 'K';
				}
				if (frame == 2) {
					targetGFFletters[i][j].character = 'J';
				}
			}
		}
	}
	*/
	 
	/*
	 for (int i=0; i<numgffs4target; i++) {
	 for (int j=0; j<numtargetlettertrans[i]; j++) {
	 cout << targetGFFletters[i][j].integer << " " << targetGFFletters[i][j].character << endl;
	 }
	 cout << endl;
	 }
	 //*/
	
	cerr << "Preparing to generate stk." << endl;
	
	string tree=argv[argc-1];
	
	int initialmafblock=alignments[0][0].starts[targetPosInDotFile[0]];
	int lastmafblock=initialmafblock;
	for (int i=0; i<numAlignedSpecies; i++) {
		for (int j=0; j<numdotlines[i]; j++) {
			int thismafblockstart=alignments[i][j].starts[targetPosInDotFile[i]] - isminus(alignments[i][j].directions[targetPosInDotFile[i]])*(alignments[i][j].length-1);
			if (thismafblockstart < initialmafblock) {
				initialmafblock=thismafblockstart;
			}
			int thismafblockend = alignments[i][j].starts[targetPosInDotFile[i]] + isplus(alignments[i][j].directions[targetPosInDotFile[i]])*(alignments[i][j].length-1);
			if (thismafblockend > lastmafblock) {
				lastmafblock=thismafblockend;
			}
		}
	}
	
	ofstream posFile;
	posFile.open(argv[argc-3]);
	posFile << "1 ";
	
	char lastcolumn[numgffs4target+numAlignedSpecies];
	for (int i=0; i<numgffs4target+numAlignedSpecies; i++) {
		lastcolumn[i]='x';
	}
	string stkLines[numgffs4target+numAlignedSpecies];
	for (int i=0; i<numgffs4target+numAlignedSpecies; i++) {
		stkLines[i]="";
	}
	
	//cout << initialmafblock << " " << lastmafblock << endl;
	
	int currentdots[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		currentdots[i]=0;
	}
	int nextdotbottoms[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		nextdotbottoms[i]=alignments[i][1].bottom(targetPosInDotFile[i]);
	}
	int initialdotbottoms[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		initialdotbottoms[i]=alignments[i][0].bottom(targetPosInDotFile[i]);
	}
	
	string sections[numAlignedSpecies+1];
	for (int i=0; i<numAlignedSpecies+1; i++) {
		sections[i]="";
	}
	
	cerr << "Generating stk." << endl;
	
	ifstream sectionsFile;
	int lastiout=0;
	bool notfirst=0;
	
	ofstream forcex;
	forcex.open("forcex");
	
	sectionsFile.open(argv[argc-2]);
	while (getline(sectionsFile,line)) {
		if (stkLines[0] != "") {
			forcex << stkLines[0].size()-1 << " ";
		}
		
		int *firstlast=getintsfromstring(line);
		initialmafblock=firstlast[0];
		lastmafblock=firstlast[1];
		int i=initialmafblock;
		if (notfirst) {
			posFile << initialmafblock << " ";
		}
		bool last_is_stad=0;
		while (i<=lastmafblock) {
			
			///*
			if (i > lastiout + 10000) {
				cerr << i << " ";
				lastiout=i;
			}
			// */
			
			for (int j=0; j<numAlignedSpecies; j++) {
				while (i >= nextdotbottoms[j] && currentdots[j]<numdotlines[j]-1) {
					currentdots[j]++;
					nextdotbottoms[j]=alignments[j][currentdots[j]+1].bottom(targetPosInDotFile[j]);
				}
			}
			
			char thiscolumn[numgffs4target+numAlignedSpecies];
			char thiscolumn_nostadn[numgffs4target+numAlignedSpecies];
			for (int j=0; j<numgffs4target; j++) {
				thiscolumn[j] =  findletterfromintchar (i, targetGFFletters[j], numtargetlettertrans[j],1);
				thiscolumn_nostadn[j] =  findletterfromintchar_nostadn (i, targetGFFletters[j], numtargetlettertrans[j],1);
			}
			for (int j=0; j<numAlignedSpecies; j++) {
				int positionOnAligned=findalignment(i, targetPosInDotFile[j], posInDotFile[j], alignments[j][currentdots[j]]);
				if (positionOnAligned != -1) {
					if (alignments[j][currentdots[j]].directions[0]*alignments[j][currentdots[j]].directions[1] == 1) {
						thiscolumn[numgffs4target+j]= findletterfromintchar(positionOnAligned, alignedGFFletters[j], numlettertrans[j],1);
						thiscolumn_nostadn[numgffs4target+j]= findletterfromintchar_nostadn(positionOnAligned, alignedGFFletters[j], numlettertrans[j],1);
						
					}
					else {
						thiscolumn[numgffs4target+j]= findletterfromintchar(positionOnAligned, reversealignedGFFletters[j], numlettertrans[j],-1);
						thiscolumn_nostadn[numgffs4target+j] = findletterfromintchar_nostadn(positionOnAligned, reversealignedGFFletters[j], numlettertrans[j],-1);
					}
				}
				else {
					int lastaligned=alignments[j][currentdots[j]].top(targetPosInDotFile[j]);
					if (lastaligned > i) {
						lastaligned=initialmafblock;
					}
					
					/*
					 int oldlastaligned=i;
					 while (oldfindalignment(oldlastaligned, targetPosInDotFile[j], posInDotFile[j], alignments[j], numdotlines[j], fudgeFactor) == -1 && oldlastaligned > initialmafblock) {
					 oldlastaligned--;
					 }
					 
					 if (oldlastaligned != lastaligned) {
					 cerr << j << " " << currentdots[j] << " " << lastaligned << " " << i << " " << oldlastaligned << endl;
					 oldfindalignment(oldlastaligned, targetPosInDotFile[j], posInDotFile[j], alignments[j], numdotlines[j], 2);
					 }
					 */
					
					int nextaligned=alignments[j][currentdots[j]+1].bottom(targetPosInDotFile[j]);
					
					int lastalignedto=findalignment(lastaligned, targetPosInDotFile[j], posInDotFile[j], alignments[j][currentdots[j]]);
					int nextalignedto=findalignment(nextaligned, targetPosInDotFile[j], posInDotFile[j], alignments[j][currentdots[j]+1]);
					
					int lastdirection=alignments[j][currentdots[j]].directions[0]*alignments[j][currentdots[j]].directions[1];
					int nextdirection=alignments[j][currentdots[j]+1].directions[0]*alignments[j][currentdots[j]+1].directions[1];
					
					//lastalignedto =  oldfindalignment(lastaligned, targetPosInDotFile[j], posInDotFile[j], alignments[j], numdotlines[j], 1);
					
					char lastalignedletter;
					if (lastdirection == 1) {
						lastalignedletter=findletterfromintchar(lastalignedto, alignedGFFletters[j], numlettertrans[j],1);
					}
					else {
						lastalignedletter=findletterfromintchar(lastalignedto, reversealignedGFFletters[j], numlettertrans[j],-1);
					}

					
					bool allsameletter=1;
					if (lastalignedto < nextalignedto) {
						for (int k=lastalignedto; k<nextalignedto; k++) {
							char thisalignedletter;
							if (lastdirection == 1) {
								thisalignedletter=findletterfromintchar(k, alignedGFFletters[j], numlettertrans[j],1);
							}
							else {
								thisalignedletter=findletterfromintchar(k, reversealignedGFFletters[j], numlettertrans[j],-1);
							}

							if (thisalignedletter != lastalignedletter) {
								allsameletter=0;
								break;
							}
						}
					}
					if (nextalignedto < lastalignedto) {
						for (int k=nextalignedto; k<lastalignedto; k++) {
							char thisalignedletter;
							if (lastdirection == 1) {
								thisalignedletter=findletterfromintchar(k, alignedGFFletters[j], numlettertrans[j],1);
							}
							else {
								thisalignedletter=findletterfromintchar(k, reversealignedGFFletters[j], numlettertrans[j],-1);
							}

							
							if (thisalignedletter != lastalignedletter) {
								allsameletter=0;
								break;
							}
						}
					}
					int distanceontarget=nextaligned-lastaligned;
					distanceontarget=sgn(distanceontarget)*distanceontarget;
					int distanceonaligned=nextalignedto-lastalignedto;
					distanceonaligned=sgn(distanceonaligned)*distanceonaligned;
					if (allsameletter && !is_cds(lastalignedletter) && !is_stad(lastalignedletter) /*&& distanceontarget < 5*distanceonaligned && distanceontarget > 0.2*distanceonaligned*/) {
						thiscolumn[numgffs4target+j]=lastalignedletter;
					}
					else {
						/*
						 cout << "Putting a start at " << i << " because ";
						 if (!allsameletter) {
						 cout << "not all same letter ";
						 }
						 if (lastalignedletter == 'e') {
						 cout << "lastalignedletter is e";
						 }
						 cout << endl;
						 */
						thiscolumn[numgffs4target+j]='*';
					}
				}
			}
			
			bool newstad=0;
			for (int j=0; j<numgffs4target+numAlignedSpecies; j++) {
				if (is_stad(thiscolumn[j])) {
					if (last_is_stad) {
						thiscolumn[j] = thiscolumn_nostadn[j];
					}
					else {
						newstad=1;
					}
				}
			}
			if (newstad) {
				last_is_stad=1;
			}
			else {
				last_is_stad=0;
			}
			
			bool stadcolumn=0;
			for (int j=0; j<numgffs4target+numAlignedSpecies; j++) {
				if (is_stad(thiscolumn[j])) {
					stadcolumn=1;
					break;
				}
			}
			if (stadcolumn) {
				for (int j=0; j<numgffs4target+numAlignedSpecies; j++) {
					if (thiscolumn[j] == 'e') {
						thiscolumn[j]='p';
					}
					if (thiscolumn[j] == 'f') {
						thiscolumn[j]='q';
					}
					if (thiscolumn[j] == 'g') {
						thiscolumn[j]='r';
					}
					if (thiscolumn[j] == 'i') {
						thiscolumn[j]='u';
					}
					if (thiscolumn[j] == 'j') {
						thiscolumn[j]='v';
					}
					if (thiscolumn[j] == 'k') {
						thiscolumn[j]='w';
					}
					if (thiscolumn[j] == 'E') {
						thiscolumn[j]='P';
					}
					if (thiscolumn[j] == 'F') {
						thiscolumn[j]='Q';
					}
					if (thiscolumn[j] == 'G') {
						thiscolumn[j]='R';
					}
					if (thiscolumn[j] == 'I') {
						thiscolumn[j]='U';
					}
					if (thiscolumn[j] == 'J') {
						thiscolumn[j]='V';
					}
					if (thiscolumn[j] == 'K') {
						thiscolumn[j]='W';
					}
					if (thiscolumn[j] == 'x') {
						thiscolumn[j]='y';
					}
				}
			}

			
			for (int j=0; j<numgffs4target+numAlignedSpecies; j++) {
				if (thiscolumn[j] != lastcolumn[j]) {
					for (int k=0; k<numgffs4target+numAlignedSpecies; k++) {
						stkLines[k] += lastcolumn[k];
					}
					posFile << i << " ";
					
					/*
					int positionOnAligned=findalignment(i, targetPosInDotFile[0], posInDotFile[0], alignments[0][currentdots[0]]);
					if (positionOnAligned == -1) {
						sections[0] += "z";
					}
					else if (positionOnAligned < 106839) {
						sections[0] += "a";
					}
					else if (positionOnAligned < 196282) {
						sections[0] += "b";
					}
					else if (positionOnAligned < 262419) {
						sections[0] += "c";
					}
					else if (positionOnAligned < 331199) {
						sections[0] += "d";
					}
					else if (positionOnAligned < 608804) {
						sections[0] += "e";
					}
					else {
						sections[0] += "f";
					}					 
					
					if (i < 126188) {
						sections[1]+="a";
					}
					else if (i < 267216) {
						sections[1]+="b";
					}
					else if (i < 386443) {
						sections[1]+="c";
					}
					else if (i < 498598) {
						sections[1]+="d";
					}
					else if (i < 612656) {
						sections[1]+="e";
					}
					else if (i < 737182) {
						sections[1]+="f";
					}
					else if (i < 836613) {
						sections[1]+="g";
					}
					else if (i < 956143) {
						sections[1]+="h";
					}
					else if (i < 1048243) {
						sections[1]+="i";
					}
					else if (i < 1155095) {
						sections[1]+="j";
					}
					else if (i < 1223555) {
						sections[1]+="k";
					}
					else if (i < 1343613) {
						sections[1]+="l";
					}
					else {
						sections[1]+="m";
					}
					
					if (numAlignedSpecies > 1) {
						positionOnAligned=findalignment(i, targetPosInDotFile[1], posInDotFile[1], alignments[1][currentdots[1]]);
						if (positionOnAligned == -1) {
							sections[2] += "z";
						}
						else if (positionOnAligned < 127198) {
							sections[2] += "a";
						}
						else if (positionOnAligned < 463111) {
							sections[2] += "b";
						}
						else if (positionOnAligned < 813065) {
							sections[2] += "c";
						}
						else {
							sections[2] += "d";
						}

					}
					 //*/

					break;
				}
			}
			
			for (int j=0; j<numgffs4target+numAlignedSpecies; j++) {
				lastcolumn[j]=thiscolumn[j];
			}
			
			int mindistance = lastmafblock - i;
			for (int j=0; j<numAlignedSpecies; j++) {
				int thisdottop = alignments[j][currentdots[j]].top(targetPosInDotFile[j]);
				if (i < thisdottop+1) {
					if (thisdottop+1 - i < mindistance) {
						mindistance = thisdottop+1 - i;
					}
					int alignedto = findalignment(i, targetPosInDotFile[j], posInDotFile[j], alignments[j][currentdots[j]]);
					int alignmentdirection = alignments[j][currentdots[j]].directions[0]*alignments[j][currentdots[j]].directions[1];
					if (alignedto != -1) {
						if (alignmentdirection == 1) {
							int nextfeature = nextlettertrans (alignedto, alignedGFFletters[j], numlettertrans[j]);
							if (nextfeature != 0 && nextfeature-alignedto < mindistance) {
								mindistance = nextfeature-alignedto;
							}
						}
						else {
							int nextfeature = prevlettertrans (alignedto, alignedGFFletters[j], numlettertrans[j]);
							if (nextfeature != 0 && alignedto - nextfeature < mindistance) {
								mindistance = alignedto - nextfeature;
							}
						}
					}
				}
				if (i < nextdotbottoms[j] && nextdotbottoms[j] - i < mindistance) {
					mindistance = nextdotbottoms[j] - i;
				}
				if (i < initialdotbottoms[j] && initialdotbottoms[j] - i < mindistance) {
					mindistance = initialdotbottoms[j] - i;
				}
			}
			
			for (int j=0; j<numgffs4target; j++) {
				int nextfeature = nextlettertrans (i, targetGFFletters[j], numtargetlettertrans[j]);
				if (nextfeature - i < mindistance && nextfeature != 0) {
					mindistance = nextfeature - i;
				}
			}
			
			mindistance=max(mindistance,1);
			//mindistance=1;
			if (last_is_stad) {
				mindistance=0;
			}
			
			i+=mindistance;
			
		}
		
		for (int j=0; j<numAlignedSpecies+numgffs4target; j++) {
			stkLines[j]+='X';
		}
		for (int j=0; j<numAlignedSpecies+1; j++) {
			sections[j]+=' ';
		}
		notfirst=1;
	}
	
	cerr << endl;
	
	for (int k=0; k<numgffs4target+numAlignedSpecies; k++) {
		stkLines[k] += lastcolumn[k];
	}
	
	int namechars=(targetName + convertInt(numgffs4target+1)).size();
	for (int i=0; i<numAlignedSpecies; i++) {
		if (namechars < alignedNames[i].size()) {
			namechars = alignedNames[i].size();
		}
	}
	namechars += 5;
	namechars = max(10,namechars);
	
	cout << "#=GF NH " << tree << ";" << endl;
	for (int i=0; i<numgffs4target; i++) {
		cout << targetName << i+1;
		for (int j=0; j<namechars - (targetName + convertInt(i+1)).size(); j++) {
			cout << " ";
		}
		cout << stkLines[i] << endl;
	}
	for (int i=0; i<numAlignedSpecies; i++) {
		cout << alignedNames[i];
		for (int j=0; j<namechars - alignedNames[i].size(); j++) {
			cout << " ";
		}
		cout << stkLines[numgffs4target + i] << endl;
	}
	
	/*
	cout << "A_sect";
	for (int i=6; i< namechars+1; i++) {
		cout << " ";
	}
	cout << sections[0] << endl;
	
	cout << "M_sect";
	for (int i=6; i< namechars+1; i++) {
		cout << " ";
	}
	cout << sections[1] << endl;
	
	if (numAlignedSpecies > 1) {
		cout << "S_sect";
		for (int i=6; i< namechars+1; i++) {
			cout << " ";
		}
		cout << sections[2] << endl;
	}
	 //*/
	
	
	return 0;
}