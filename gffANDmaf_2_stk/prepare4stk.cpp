//This code is a little sketchy. If you get new data, you might consider running valigrind first to make sure everything is okay.

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

string convertInt(int number)
{
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
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

struct intlistelem {
	int data;
	intlistelem *next;
};

class intlist {
	intlistelem *root;
	intlistelem *last;
public:
	int size;
	void initialize();
	int value(unsigned int);
	void set(unsigned int,int);
	void add(int);
	void remove(int);
	void insert(unsigned int, int);
};

void intlist::initialize() {
	root=new intlistelem;
	root->next=NULL;
	root->data=0;
	size=0;
	last=root;
}

int intlist::value(unsigned int position) {
	intlistelem *conductor=root;
	for (int i=0; i<position; i++) {
		if (conductor != NULL) {
			conductor=conductor->next;
		}
		else {
			return 0;
		}
	}
	return conductor->data;
}

void intlist::set(unsigned int position, int newdata) {
	intlistelem *conductor=root;
	for (int i=0; i<position; i++) {
		if (conductor->next == NULL) {
			conductor->next=new intlistelem;
			conductor=conductor->next;
			conductor->next=NULL;
			conductor->data=0;
			last = conductor;
		}
		else {
			conductor=conductor->next;
		}
	}
	conductor->data=newdata;
	if (size < position + 1) {
		size=position+1;
	}
}

void intlist::add(int addthis) {
	if (size == 0) {
		root->data=addthis;
		root->next=NULL;
		size++;
		return;
	}
	last->next = new intlistelem;
	last = last->next;
	last->data = addthis;
	last->next = NULL;
	size++;
}

void intlist::remove (int removeThisPosition) {
	intlistelem *conductor=root;
	for (int i=0; i<removeThisPosition-1; i++) {
		conductor->next = conductor->next->next;
	}
	size--;
}

void intlist::insert (unsigned int insertHere, int insertThisData) {
	intlistelem *insertion = new intlistelem;
	insertion->data=insertThisData;
	size++;
	if (insertHere == 0) {
		insertion->next=root;
		root = insertion;
		return;
	}
	intlistelem *before=root;
	for (int i=0; i<insertHere-1; i++) {
		before=before->next;
	}
	insertion->next=before->next;
	before->next=insertion;
	return;
}	


struct intarraylistelem {
	int *data;
	int datalength;
	intarraylistelem *next;
};

class intarraylist {
	intarraylistelem *root;
	intarraylistelem *last;
public:
	int size;
	void initialize();
	int *value(unsigned int);
	int arraysize(unsigned int);
	void set(unsigned int,int*,int);
	void add(int*,int);
	void remove(int);
	void insert(unsigned int, int*, int);
};

void intarraylist::initialize() {
	root=new intarraylistelem;
	root->next=NULL;
	root->data=NULL;
	root->datalength=0;
	size=0;
	last=root;
}

int *intarraylist::value(unsigned int position) {
	intarraylistelem *conductor=root;
	for (int i=0; i<position; i++) {
		if (conductor != NULL) {
			conductor=conductor->next;
		}
		else {
			return 0;
		}
	}
	return conductor->data;
}

int intarraylist::arraysize(unsigned int position) {
	intarraylistelem *conductor=root;
	for (int i=0; i<position; i++) {
		if (conductor != NULL) {
			conductor=conductor->next;
		}
		else {
			return 0;
		}
	}
	return conductor->datalength;
}

void intarraylist::set(unsigned int position, int *newdata, int newdatalength) {
	intarraylistelem *conductor=root;
	for (int i=0; i<position; i++) {
		if (conductor->next == NULL) {
			conductor->next=new intarraylistelem;
			conductor=conductor->next;
			conductor->next=NULL;
			conductor->data=0;
			last = conductor;
		}
		else {
			conductor=conductor->next;
		}
	}
	conductor->data=newdata;
	conductor->datalength = newdatalength;
	if (size < position + 1) {
		size=position+1;
	}
}

void intarraylist::add(int *addthis, int length) {
	if (size == 0) {
		root->data=addthis;
		root->datalength=length;
		root->next=NULL;
		size++;
		return;
	}
	last->next = new intarraylistelem;
	last = last->next;
	last->data = addthis;
	last->datalength = length;
	last->next = NULL;
	size++;
}

void intarraylist::remove (int removeThisPosition) {
	intarraylistelem *conductor=root;
	for (int i=0; i<removeThisPosition-1; i++) {
		conductor->next = conductor->next->next;
	}
	size--;
}

void intarraylist::insert (unsigned int insertHere, int *insertThisData, int length) {
	intarraylistelem *insertion = new intarraylistelem;
	insertion->data=insertThisData;
	insertion->datalength=length;
	size++;
	if (insertHere == 0) {
		insertion->next=root;
		root = insertion;
		return;
	}
	intarraylistelem *before=root;
	for (int i=0; i<insertHere-1; i++) {
		before=before->next;
	}
	insertion->next=before->next;
	before->next=insertion;
	return;
}	

class chain {
public:
	int alignFromStart;
	int alignFromEnd;
	int *scores;
	int numscores;
	int *mafslist;
	int nummafs;
};

void bubblesortchains (chain *sortthis, int sizetosort) {
	bool switched=1;
	while (switched) {
		switched=0;
		for (int i=1; i<sizetosort; i++) {
			if (sortthis[i-1].alignFromEnd > sortthis[i].alignFromEnd) {
				chain temp=sortthis[i-1];
				sortthis[i-1]=sortthis[i];
				sortthis[i]=temp;
				switched=1;
			}
		}
	}
}

class dotline {
public:
	int starts[2];
	int length;
	int directions[2];
	int numinalignment;
	int bottom(int);
	int top(int);
	void shortenbottom(int,int);
	void shortentop(int,int);
	void print();
	string makestring();
};

int dotline::bottom(int column) {
	return starts[column] - isminus(directions[column])*(length-1);
}

int dotline::top(int column) {
	return starts[column] + isplus(directions[column])*(length-1);
}

void dotline::shortenbottom(int shortenby, int column) {
	if (directions[column] == -1) {
		for (int i=0; i<2; i++) {
			starts[i] -= directions[i]*directions[column];
		}
	}
	length=length-shortenby;
}

void dotline::shortentop(int shortenby, int column) {
	if (directions[column] == 1) {
		for (int i=0; i<2; i++) {
			starts[i] += directions[i]*directions[column];
		}
	}
	length=length-shortenby;
}

void dotline::print() {
	for (int i=0; i<2; i++) {
		cout << starts[i] << " ";
	}
	cout << length << " ";
	for (int i=0; i<2; i++) {
		if (directions[i] > 0) {
			cout << "+ ";
		}
		else {
			cout << "- ";
		}
	}
	cout << endl;
}

string dotline::makestring() {
	string result="";
	for (int i=0; i<2; i++) {
		result += (convertInt(starts[i]) + " ");
	}
	result += (convertInt(length) + " ");
	for (int i=0; i<2; i++) {
		if (directions[i] > 0) {
			result += "+ ";
		}
		else {
			result += "- ";
		}
	}
	return result;
}

dotline getDotLine (string line) { //read space for tab...
	int numspaces=0;
	for (int i=0; i<line.length(); i++) {
		if (line[i]==' ') {
			numspaces++;
		}
	}
	if (numspaces < 5) {
		cerr << line << endl << "is not a proper dotline" << endl;
		exit(1);
	}

	int tab=-1;
	int nexttab;
	string data[2*2+2];
	for (int i=0; i<2*2+2; i++) {
		nexttab = tab+1;
		while (line[nexttab] != ' ' && nexttab < line.size()) {
			nexttab++;
		}
		data[i] = line.substr(tab+1,nexttab-tab-1);
		tab=nexttab;
	}
	dotline row;

	for (int i=0; i<2; i++) {
		if (!from_string<int>(row.starts[i],data[i],std::dec)) {
			cerr << "from_string failed on data[" << i << "]" << endl;
		}
	}
	if (!from_string<int>(row.length,data[2],std::dec)) {
		cerr << "from_string failed on data[" << 2 << "]" << endl;
	}
	for (int i=0; i<2; i++) {
		if (data[2+1+i][0] == '+') {
			row.directions[i] = 1;
		}
		else {
			row.directions[i]=-1;
		}
	}

	if (!from_string<int>(row.numinalignment,data[2*2+1],std::dec)) {
		cerr << "from_string failed on data[" << 2*2+1 << "]" << endl;
	}

	//cout << "Got dotline" << endl;
	return row;
}

/*
class dotline3 {
public:
	int starts[3];
	int length;
	int directions[3];
	int bottom(int);
	int top(int);
	void shortenbottom(int,int);
	void shortentop(int,int);
};

int dotline3::bottom(int column) {
	return starts[column] - isminus(directions[column])*(length-1);
}

int dotline3::top(int column) {
	return starts[column] + isplus(directions[column])*(length-1);
}

void dotline3::shortenbottom(int shortenby, int column) {
	if (directions[column] == -1) {
		for (int i=0; i<3; i++) {
			starts[i] -= directions[i]*directions[column];
		}
	}
	length=length-shortenby;
}

void dotline3::shortentop(int shortenby, int column) {
	if (directions[column] == 1) {
		for (int i=0; i<3; i++) {
			starts[i] += directions[i]*directions[column];
		}
	}
	length=length-shortenby;
}

dotline3 getdotline3 (string line) { //read space for tab...
	int numspaces=0;
	for (int i=0; i<line.length(); i++) {
		if (line[i]==' ') {
			numspaces++;
		}
	}
	if (numspaces < 2*3) {
		cerr << line << endl << "is not a proper dotline3" << endl;
		exit(1);
	}
	
	int tab=-1;
	int nexttab;
	string data[2*3+1];
	for (int i=0; i<2*3+1; i++) {
		nexttab = tab+1;
		while (line[nexttab] != ' ' && nexttab < line.size()) {
			nexttab++;
		}
		data[i] = line.substr(tab+1,nexttab-tab-1);
		tab=nexttab;
	}
	dotline3 row;
	
	for (int i=0; i<3; i++) {
		if (!from_string<int>(row.starts[i],data[i],std::dec)) {
			cerr << "from_string failed on data[" << i << "]" << endl;
		}
	}
	if (!from_string<int>(row.length,data[3],std::dec)) {
		cerr << "from_string failed on data[" << 3 << "]" << endl;
	}
	for (int i=0; i<3; i++) {
		if (data[3+1+i][0] == '+') {
			row.directions[i] = 1;
		}
		else {
			row.directions[i]=-1;
		}
	}
	//cout << "Got dotline3" << endl;
	return row;
}
*/

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

char findletterfromintchar (int position, intchar *letters, int numintchars) {
	if (position < 1) {
		return '*';
	}
	for (int i=0; i<numintchars-1; i++) {
		if (position < letters[i+1].integer) {
			return letters[i].character;
		}
	}
	return letters[numintchars-1].character;
}

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
 
int main (int argc, char *argv[]) {
	if (argc < 3) {
		cerr << "Please execute like " << argv[0] << " numAlignedSpecies alignedSpeciesInfo targetName numgffs4target gffs4target start stop" << endl;
		cerr << "alignedSpeciesInfo should look like SpeciesName dotFile targetPosInDotFile alignedPosInDotFile gffFile" << endl;
		exit(1);
	}
	
	string line;
	int numAlignedSpecies=atoi(argv[1]);
	
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
	int maxalignment = 0;
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
			if (alignments[i][j].numinalignment > maxalignment) {
				maxalignment = alignments[i][j].numinalignment;
			}
		}
		dotFile.close();
		sortdots(alignments[i],numdotlines[i],targetPosInDotFile[i],posInDotFile[i]);
		//numdotlines[i]=removeoverlaps(alignments[i],numdotlines[i],targetPosInDotFile[i]);
		
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
	
	/*
	for (int i=0; i<numAlignedSpecies; i++) {
		for (int j=0; j<numlettertrans[i]; j++) {
			cout << alignedGFFletters[i][j].integer << " " << alignedGFFletters[i][j].character << endl;
		}
		cout << endl;
	}
	 */
	
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
	 for (int i=0; i<numAlignedSpecies; i++) {
	 for (int j=0; j<2*thisgfflines[i]+1; j++) {
	 cout << alignedGFFletters[i][j].integer << " " << alignedGFFletters[i][j].character << endl;
	 }
	 cout << endl;
	 }
	 */
	
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
	
	
	cerr << "Preparing" << endl;
	
	int start = atoi(argv[argc-2]);
	int end = atoi(argv[argc-1]);
	
	int currentdots[numAlignedSpecies];
	for (int i=0; i<numAlignedSpecies; i++) {
		currentdots[i]=0;
	}
	
	intlist sectionstarts;
	sectionstarts.initialize();
	intlist sectionends;
	sectionends.initialize();
	int thisstart=start;
	int thisend;
	
	int whichlettertrans[numgffs4target];
	for (int i=0; i<numgffs4target; i++) {
		whichlettertrans[i]=0;
		while (targetGFFletters[i][whichlettertrans[i]+1].integer <= start && whichlettertrans[i]+1 < numtargetlettertrans[i]-1) {
			whichlettertrans[i]++;
		}
	}
	
	bool ingene=0;
	for (int j=0; j<numgffs4target; j++) {
		if (targetGFFletters[j][whichlettertrans[j]].character != 'x') {
			ingene=1;
			break;
		}
	}
	
	int dotsforsection[numAlignedSpecies][2];
	for (int i=0; i<numAlignedSpecies; i++) {
		dotsforsection[i][0]=0;
		dotsforsection[i][1]=-1;
	}
		
	ofstream goodmafs[numAlignedSpecies];
	for (int j=0; j<numAlignedSpecies; j++) {
		string filename="goodmafs" + convertInt(j);
		goodmafs[j].open(filename.c_str());
	}
	
	
	cerr << "Beginning..." << endl;
	bool lastallx=0;
	int i=start;
	while (i < end) {
		//cerr << i << " " << end << endl;
		
		//Find the evidence letters at position i and see if they are all x
		for (int j=0; j<numgffs4target; j++) {
			while (targetGFFletters[j][whichlettertrans[j]+1].integer <= i && whichlettertrans[j]+1 < numtargetlettertrans[j]-1) {
				whichlettertrans[j]++;
			}
		}
		
		bool allx=1;
		for (int j=0; j<numgffs4target; j++) {
			if (targetGFFletters[j][whichlettertrans[j]].character != 'x') {
				allx=0;
				break;
			}
		}
		
		//If they are all x and last time they were not all x, then i marks the end of a section.
		bool finishedsection=0;
		if (allx) {
			/*
			for (int j=0; j<numgffs4target; j++) {
				cout << j << " " << whichlettertrans[j] << " " << targetGFFletters[j][whichlettertrans[j]].integer << " " << targetGFFletters[j][whichlettertrans[j]].character << endl;
			}
			 //*/
			if (ingene) {
				thisend=i;
				sectionstarts.add(thisstart);
				sectionends.add(thisend);
				ingene=0;
				finishedsection=1;
			}
			else {
				//cout << "here" << endl;
				thisstart=i;
			}
		}
		else {
			ingene=1;
			if (lastallx) {
				thisstart=i-1;
			}
		}
		
		if (allx) {
			lastallx=1;
		}
		else {
			lastallx=0;
		}

		// At the end of a section we must pick out the best set of alignments for that section.
		if (finishedsection) {
			//cerr << "Trying: " << thisstart << " " << thisend << endl;
			
			//Find all the alignments that are in this section.
			for (int j=0; j<numAlignedSpecies; j++) {
				while (alignments[j][dotsforsection[j][0]].top(targetPosInDotFile[j]) < thisstart && dotsforsection[j][0] < numdotlines[j]-1) {
					dotsforsection[j][0]++;
				}
				while (alignments[j][dotsforsection[j][1]+1].bottom(targetPosInDotFile[j]) < thisend && dotsforsection[j][1]+1  < numdotlines[j]-1) {
					dotsforsection[j][1]++;
				}
			}
			
			/*
			cout << thisstart << " " << thisend << endl;
			for (int j=dotsforsection[0][0]; j<=dotsforsection[0][1]; j++) {
				cout << alignments[0][j].bottom(targetPosInDotFile[0]) << " " << alignments[0][j].top(targetPosInDotFile[0]) << endl;
			}
			//cout << endl;
			 //*/
			
			//Group alignments that align into the same gene on the aligned (non target) species.
			for (int j=0; j<numAlignedSpecies; j++) {
				
				//cout << "new aligned scpecies" << endl;
				intlist aligntostarts;
				aligntostarts.initialize();
				intlist	aligntoends;
				aligntoends.initialize();
				intlist alignfromstarts;
				alignfromstarts.initialize();
				intlist alignfromends;
				alignfromends.initialize();
				intarraylist scores;
				scores.initialize();
				intarraylist mafslist;
				mafslist.initialize();
								
				for (int k=dotsforsection[j][0]; k<=dotsforsection[j][1]; k++) {
					//cout << dotsforsection[j][0] << " " << dotsforsection[j][1] << endl;
					bool havethis=0;
					for (int l=0; l<aligntostarts.size; l++) {
						if (alignments[j][k].top(posInDotFile[j]) > aligntostarts.value(l) && alignments[j][k].bottom(posInDotFile[j]) < aligntoends.value(l)) {
							havethis=1;
							alignfromstarts.set(l,min(alignfromstarts.value(l),alignments[j][k].bottom(targetPosInDotFile[j])));
							alignfromends.set(l,max(alignfromends.value(l),alignments[j][k].top(targetPosInDotFile[j])));
						}
					}
					if (!havethis) {
						int lastxendbefore=0;
						int firstxafter=0;
						//cout << "thenumberis: " << numlettertrans[j] << " "<< j << endl;
						for (int l=0; l<numlettertrans[j]; l++) {
							if (alignedGFFletters[j][l].character == 'x' && alignedGFFletters[j][l].integer < alignments[j][k].bottom(posInDotFile[j]) && l < numlettertrans[j]-1) {
								lastxendbefore = alignedGFFletters[j][l+1].integer;
							}
							if (alignedGFFletters[j][l].character == 'x' && alignedGFFletters[j][l].integer > alignments[j][k].top(posInDotFile[j])) {
								firstxafter = alignedGFFletters[j][l].integer;
								break;
							}
						}
						if (lastxendbefore != 0 && firstxafter != 0 && lastxendbefore < alignments[j][k].top(posInDotFile[j]) && firstxafter > alignments[j][k].bottom(posInDotFile[j])) {
							aligntostarts.add(lastxendbefore);
							aligntoends.add(firstxafter);
							alignfromstarts.add(alignments[j][k].bottom(targetPosInDotFile[j]));
							alignfromends.add(alignments[j][k].top(targetPosInDotFile[j]));
						}
						else {
							aligntostarts.add(alignments[j][k].bottom(posInDotFile[j]));
							aligntoends.add(alignments[j][k].top(posInDotFile[j]));
							alignfromstarts.add(alignments[j][k].bottom(targetPosInDotFile[j]));
							alignfromends.add(alignments[j][k].top(targetPosInDotFile[j]));
						}

					}
				}
				
				//cout << thisstart << " " << thisend << endl;
				for (int k=dotsforsection[0][0]; k<=dotsforsection[0][1]; k++) {
					//alignments[j][k].print();
				}
				
				//Keep track of the alignments in each group and check to see if there are any overlaps.
				for (int k=0; k<aligntostarts.size; k++) {
					//cout << aligntostarts.value(k) << " " << aligntoends.value(k) << " " << alignfromstarts.value(k) << " " << alignfromends.value(k) << endl;
										
					int numforwardmafs=0;
					int numbackwardmafs=0;
					int *forwardmafs = new int[max(1,dotsforsection[0][1]-dotsforsection[0][0])+1];
					int forwardcrossing[max(1,dotsforsection[0][1]-dotsforsection[0][0])];
					int *backwardmafs = new int[max(1,dotsforsection[0][1]-dotsforsection[0][0])+1];
					int backwardcrossing[max(1,dotsforsection[0][1]-dotsforsection[0][0])];
					 
					for (int l=dotsforsection[0][0]; l<=dotsforsection[0][1]; l++) {
						if (alignments[j][l].top(posInDotFile[j]) > aligntostarts.value(k) && alignments[j][l].bottom(posInDotFile[j]) < aligntoends.value(k)) {
							if (alignments[j][l].top(targetPosInDotFile[j]) > alignfromstarts.value(k) && alignments[j][l].bottom(targetPosInDotFile[j]) < alignfromends.value(k)) {
								if (alignments[j][l].directions[0]*alignments[j][l].directions[1] > 0) {
									//alignments[j][l].print();
									
									forwardmafs[numforwardmafs]=l;
									forwardcrossing[numforwardmafs]=0;
									for (int m=0; m<numforwardmafs; m++) {
										if (alignments[j][forwardmafs[numforwardmafs]].bottom(posInDotFile[j]) < alignments[j][forwardmafs[m]].top(posInDotFile[j]) || alignments[j][forwardmafs[numforwardmafs]].bottom(targetPosInDotFile[j]) < alignments[j][forwardmafs[m]].top(targetPosInDotFile[j])) {
											forwardcrossing[m]+=alignments[j][forwardmafs[numforwardmafs]].length;
											forwardcrossing[numforwardmafs]+=alignments[j][forwardmafs[m]].length;
										}
									}
									numforwardmafs++;
								}
								else {
									//alignments[j][l].print();
									
									backwardmafs[numbackwardmafs]=l;
									backwardcrossing[numbackwardmafs]=0;
									for (int m=0; m<numbackwardmafs; m++) {
										if (alignments[j][backwardmafs[numbackwardmafs]].top(posInDotFile[j]) > alignments[j][backwardmafs[m]].bottom(posInDotFile[j]) || alignments[j][backwardmafs[numbackwardmafs]].bottom(targetPosInDotFile[j]) < alignments[j][backwardmafs[m]].top(targetPosInDotFile[j])) {
											backwardcrossing[m]+=alignments[j][backwardmafs[numbackwardmafs]].length;
											backwardcrossing[numbackwardmafs]+=alignments[j][backwardmafs[m]].length;
										}
									}
									numbackwardmafs++;
								}
							}
						}
					}
					
					bool forwardnonzero=0;
					for (int l=0; l<numforwardmafs; l++) {
						if (forwardcrossing[l] != 0) {
							forwardnonzero=1;
							break;
						}
					}
					
					//If there are overlaps check to see if they can be eliminated by splitting the group. This often happens if the two parts of the section align to the same gene. (i.e. a gene duplication that one of the evidences joins into one gene.)
					if (forwardnonzero) {
						bool cut=0;
						for (int l=0; l<numforwardmafs-1; l++) {
							for (int m=l+1; m<numforwardmafs; m++) {
								if (alignments[j][forwardmafs[l]].length == alignments[j][forwardmafs[m]].length) {
									if (alignments[j][forwardmafs[l]].starts[posInDotFile[j]] == alignments[j][forwardmafs[m]].starts[posInDotFile[j]]) {
										aligntostarts.add(aligntostarts.value(k));
										aligntoends.add(aligntoends.value(k));
										alignfromstarts.add( (alignments[j][forwardmafs[m]].bottom(targetPosInDotFile[j])+alignments[j][forwardmafs[m-1]].top(targetPosInDotFile[j]))/2 );
										alignfromends.add(alignfromends.value(k));
										
										alignfromends.set(k,(alignments[j][forwardmafs[m]].bottom(targetPosInDotFile[j])+alignments[j][forwardmafs[m-1]].top(targetPosInDotFile[j]))/2-1);
										numforwardmafs=m;
										
										cut = 1;
										break;
									}
								}
							}
							if (cut) {
								break;
							}
						}
						
						if (cut) {
							for (int l=0; l<numforwardmafs; l++) {
								forwardcrossing[l]=0;
							}
							for (int l=0; l<numforwardmafs-1; l++) {
								for (int m=l+1; m<numforwardmafs; m++) {
									if (alignments[j][forwardmafs[m]].bottom(posInDotFile[j]) < alignments[j][forwardmafs[l]].top(posInDotFile[j]) || alignments[j][forwardmafs[numforwardmafs]].bottom(targetPosInDotFile[j]) < alignments[j][forwardmafs[m]].top(targetPosInDotFile[j])) {
										forwardcrossing[m]+=alignments[j][forwardmafs[l]].length;
										forwardcrossing[l]+=alignments[j][forwardmafs[m]].length;
									}
								}
								
							}
							
							forwardnonzero=0;
							for (int l=0; l<numforwardmafs; l++) {
								if (forwardcrossing[l] != 0) {
									forwardnonzero=1;
									break;
								}
							}
							
						}
					}
					
					bool backwardnonzero=0;
					for (int l=0; l<numbackwardmafs; l++) {
						if (backwardcrossing[l] != 0) {
							backwardnonzero=1;
							break;
						}
					}
					
					if (backwardnonzero) {
						bool cut=0;
						for (int l=0; l<numbackwardmafs-1; l++) {
							for (int m=l+1; m<numbackwardmafs; m++) {
								if (alignments[j][backwardmafs[l]].length == alignments[j][backwardmafs[m]].length) {
									if (alignments[j][backwardmafs[l]].starts[posInDotFile[j]] == alignments[j][backwardmafs[m]].starts[posInDotFile[j]]) {
										aligntostarts.add(aligntostarts.value(k));
										aligntoends.add(aligntoends.value(k));
										alignfromstarts.add( (alignments[j][backwardmafs[m-1]].bottom(targetPosInDotFile[j])+alignments[j][backwardmafs[m]].top(targetPosInDotFile[j]))/2 );
										alignfromends.add(alignfromends.value(k));
										
										alignfromends.set(k,(alignments[j][backwardmafs[m-1]].bottom(targetPosInDotFile[j])+alignments[j][backwardmafs[m]].top(targetPosInDotFile[j]))/2-1);
										numbackwardmafs=m;
										
										cut = 1;
										break;
									}
								}
							}
							if (cut) {
								break;
							}
						}
						
						if (cut) {
							for (int l=0; l<numbackwardmafs; l++) {
								backwardcrossing[l]=0;
							}
							for (int l=0; l<numbackwardmafs-1; l++) {
								for (int m=l+1; m<numbackwardmafs; m++) {
									if (alignments[j][backwardmafs[m]].bottom(posInDotFile[j]) > alignments[j][backwardmafs[l]].top(posInDotFile[j]) || alignments[j][backwardmafs[numbackwardmafs]].bottom(targetPosInDotFile[j]) < alignments[j][backwardmafs[m]].top(targetPosInDotFile[j])) {
										backwardcrossing[m]+=alignments[j][backwardmafs[l]].length;
										backwardcrossing[l]+=alignments[j][backwardmafs[m]].length;
									}
								}
								
							}
							
							backwardnonzero=0;
							for (int l=0; l<numbackwardmafs; l++) {
								if (backwardcrossing[l] != 0) {
									backwardnonzero=1;
									break;
								}
							}
						}
					}
					
					/*
					//Aggressive Cut
					if (forwardnonzero) {
						for (int l=0; l<numforwardmafs; l++) {
							if (forwardcrossing[l] > 0) {
								for (int m=l+1; m<numforwardmafs; m++) {
									forwardmafs[m-1]=forwardmafs[m];
									forwardcrossing[m-1]=forwardcrossing[m];
								}
								numforwardmafs--;
							}
						}
					}
					if (backwardnonzero) {
						for (int l=0; l<numbackwardmafs; l++) {
							if (backwardcrossing[l] > 0) {
								for (int m=l+1; m<numbackwardmafs; m++) {
									backwardmafs[m-1]=backwardmafs[m];
									backwardcrossing[m-1]=backwardcrossing[m];
								}
								numbackwardmafs--;
							}
						}
					}
					// */
					
					///*
					//Soft Cut: Remove overlap mafs one at a time, in order of how much they overlap and how many way the alignment is until there are no more overlaps.
					while (forwardnonzero) {
						int maxcrossing=0;
						int maxcrossingalignmentsize=100;
						int maxcrossingpos=0;
						for (int l=0; l<numforwardmafs; l++) {
							if (forwardcrossing[l] > maxcrossing && alignments[j][forwardmafs[l]].numinalignment < maxcrossingalignmentsize) {
								maxcrossing = forwardcrossing[l];
								maxcrossingpos=l;
								maxcrossingalignmentsize = alignments[j][forwardmafs[l]].numinalignment;
							}
						}
						for (int l=maxcrossingpos+1; l<numforwardmafs; l++) {
							forwardmafs[l-1]=forwardmafs[l];
						}
						numforwardmafs--;
						for (int l=0; l<numforwardmafs; l++) {
							forwardcrossing[l]=0;
						}
						for (int l=0; l<numforwardmafs-1; l++) {
							for (int m=l+1; m<numforwardmafs; m++) {
								if (alignments[j][forwardmafs[m]].bottom(posInDotFile[j]) < alignments[j][forwardmafs[l]].top(posInDotFile[j]) || alignments[j][forwardmafs[numforwardmafs]].bottom(targetPosInDotFile[j]) < alignments[j][forwardmafs[m]].top(targetPosInDotFile[j])) {
									forwardcrossing[m]+=alignments[j][forwardmafs[l]].length;
									forwardcrossing[l]+=alignments[j][forwardmafs[m]].length;
								}
							}
						}
						forwardnonzero=0;
						for (int l=0; l<numforwardmafs; l++) {
							if (forwardcrossing[l] != 0) {
								forwardnonzero=1;
								break;
							}
						}
					}
					while (backwardnonzero) {
						int maxcrossing=0;
						int maxcrossingalignmentsize=100;
						int maxcrossingpos=0;
						for (int l=0; l<numbackwardmafs; l++) {
							if (backwardcrossing[l] > maxcrossing && alignments[j][backwardmafs[l]].numinalignment < maxcrossingalignmentsize) {
								maxcrossing = backwardcrossing[l];
								maxcrossingpos=l;
							}
						}
						for (int l=maxcrossingpos+1; l<numbackwardmafs; l++) {
							backwardmafs[l-1]=backwardmafs[l];
						}
						numbackwardmafs--;
						for (int l=0; l<numbackwardmafs; l++) {
							backwardcrossing[l]=0;
						}
						for (int l=0; l<numbackwardmafs-1; l++) {
							for (int m=l+1; m<numbackwardmafs; m++) {
								if (alignments[j][backwardmafs[m]].bottom(posInDotFile[j]) < alignments[j][backwardmafs[l]].top(posInDotFile[j]) || alignments[j][backwardmafs[numbackwardmafs]].bottom(targetPosInDotFile[j]) < alignments[j][backwardmafs[m]].top(targetPosInDotFile[j])) {
									backwardcrossing[m]+=alignments[j][backwardmafs[l]].length;
									backwardcrossing[l]+=alignments[j][backwardmafs[m]].length;
								}
							}
							
						}
						backwardnonzero=0;
						for (int l=0; l<numbackwardmafs; l++) {
							if (backwardcrossing[l] != 0) {
								backwardnonzero=1;
								break;
							}
						}
					}
					// */
					
					//Now we want to score each group. Each group receives a score for the amount 2+way, 3+way, etc alignment that it contains.
					int *forwardscore;
					forwardscore = new int[maxalignment-1];
					for (int l=0; l<maxalignment-1; l++) {
						forwardscore[l]=0;
					}
					for (int l=0; l<numforwardmafs; l++) {
						for (int m=0; m<maxalignment-1; m++) {
							if (alignments[j][forwardmafs[l]].numinalignment >= m+2) {
								forwardscore[m]+=alignments[j][forwardmafs[l]].length;
							}
						}
					}
					
					int *backwardscore;
					backwardscore = new int[maxalignment-1];
					for (int l=0; l<maxalignment-1; l++) {
						backwardscore[l]=0;
					}
					for (int l=0; l<numbackwardmafs; l++) {
						for (int m=0; m<maxalignment; m++) {
							if (alignments[j][backwardmafs[l]].numinalignment >= m+2) {
								backwardscore[m]+=alignments[j][backwardmafs[l]].length;
							}
						}
					}
					
					bool chooseforward = 1;
					for (int l=maxalignment-2; l >= 0; l--) {
						if (backwardscore[l] > forwardscore[l]) {
							chooseforward = 0;
							break;
						}
						if (forwardscore[l] > backwardscore[l]) {
							chooseforward = 1;
							break;
						}
					}
					
					if (chooseforward) {
						mafslist.add(forwardmafs,numforwardmafs);
						scores.add(forwardscore,maxalignment-1);
					}
					else {
						mafslist.add(backwardmafs,numbackwardmafs);
						scores.add(backwardscore,maxalignment-1);
					}

					/*
					//if (forwardnonzero || backwardnonzero) {
						//cout << thisstart << " " << thisend << endl;
						cout << aligntostarts.value(k) << " " << aligntoends.value(k) << " " << alignfromstarts.value(k) << " " << alignfromends.value(k) << " " << scores.value(k) << endl;
						for (int l=0; l<numforwardmafs; l++) {
							//alignments[j][forwardmafs[l]].print();
						}
						for (int l=0; l<numbackwardmafs; l++) {
							//alignments[j][backwardmafs[l]].print();
						}
						//cout << endl;
					//}
					// */
				}
				
				//For ease we copy the relevent intformation about the groups into "chains" and sort them by end position.
				chain chains[aligntostarts.size];
				for (int k=0; k<aligntostarts.size; k++) {
					chains[k].alignFromStart=alignfromstarts.value(k);
					chains[k].alignFromEnd=alignfromends.value(k);
					chains[k].scores=scores.value(k);
					chains[k].numscores=scores.arraysize(k);
					chains[k].mafslist=mafslist.value(k);
					chains[k].nummafs=mafslist.arraysize(k);
				}
				bubblesortchains(chains,aligntostarts.size);
				
				/*
				for (int k=0; k<aligntostarts.size; k++) {
					cout << chains[k].alignFromStart << " ";
					cout << chains[k].alignFromEnd << " ";
					cout << chains[k].score << endl;
				}
				 */
				
				//Find the set of chains that give the best score without overlapping.
				if (aligntostarts.size > 0) {
					chain bestcombination[aligntostarts.size][aligntostarts.size];
					bestcombination[0][0]=chains[0];
					int numinbestcombination[aligntostarts.size];
					numinbestcombination[0]=1;
					int bestcombinationscore[aligntostarts.size][maxalignment-1];
					for (int m=0; m<maxalignment-1; m++) {
						bestcombinationscore[0][m]=chains[0].scores[m];
					}
					
					for (int k=1; k<aligntostarts.size; k++) {
						//cout << "for " << k << endl;
						bool tackedon=0;
						for (int l=k-1; l>=0; l--) {
							//cout << "comparing to " << l << endl;
							if (bestcombination[l][numinbestcombination[l]-1].alignFromEnd <= chains[k].alignFromStart) {
								//cout << "can add " << bestcombination[l][numinbestcombination[l]-1].alignFromEnd << endl;
								int sumscore[maxalignment-1];
								for (int m=0; m<maxalignment-1; m++) {
									sumscore[m] = bestcombinationscore[l][m] + chains[k].scores[m];
								}
								
								bool usenewscore = 0;
								for (int m=maxalignment-2; m >= 0; m--) {
									if (sumscore[m] > bestcombinationscore[k-1][m]) {
										usenewscore = 1;
										break;
									}
									if (sumscore[m] < bestcombinationscore[k-1][m]) {
										usenewscore = 0;
										break;
									}
								}
								
								if (usenewscore) {
									//cout << "adding" << endl;
									for (int m=0; m<numinbestcombination[l]; m++) {
										bestcombination[k][m]=bestcombination[l][m];
									}
									bestcombination[k][numinbestcombination[l]]=chains[k];
									numinbestcombination[k]=numinbestcombination[l]+1;
									for (int m=0; m<maxalignment-1; m++) {
										bestcombinationscore[k][m] = sumscore[m];
									}
								}
								else {
									//cout << "not adding" << endl;
									for (int m=0; m<numinbestcombination[k-1]; m++) {
										bestcombination[k][m]=bestcombination[k-1][m];
									}
									numinbestcombination[k]=numinbestcombination[k-1];
									for (int m=0; m<maxalignment-1; m++) {
										bestcombinationscore[k][m]=bestcombinationscore[k-1][m];
									}
								}
								tackedon=1;
								break;
							}
						}
						if (!tackedon) {
							//cout << "new" << endl;
							
							bool usenewscore = 0;
							for (int m=maxalignment-2; m >= 0; m--) {
								if (chains[k].scores[m] > bestcombinationscore[k-1][m]) {
									usenewscore = 1;
									break;
								}
								if (chains[k].scores[m] < bestcombinationscore[k-1][m]) {
									usenewscore = 0;
									break;
								}
							}
							if (usenewscore) {
								//cout << "replacing" << endl;
								bestcombination[k][0]=chains[k];
								numinbestcombination[k]=1;
								for (int m=0; m<maxalignment-1; m++) {
									bestcombinationscore[k][m]=chains[k].scores[m];
								}
							}
							else {
								//cout << "not replacing" << endl;
								for (int m=0; m<numinbestcombination[k-1]; m++) {
									bestcombination[k][m]=bestcombination[k-1][m];
								}
								numinbestcombination[k]=numinbestcombination[k-1];
								for (int m=0; m<maxalignment-1; m++) {
									bestcombinationscore[k][m]=bestcombinationscore[k-1][m];
								}
							}
						}
					}
					
					//Save the chosen alignments to disk.
					int nummafsinblocks=0;
					for (int k=0; k<numinbestcombination[aligntostarts.size-1]; k++) {
						nummafsinblocks+=bestcombination[aligntostarts.size-1][k].nummafs;
					}
					dotline alignmentsforthissection[nummafsinblocks];
					int position=0;
					for (int k=0; k<numinbestcombination[aligntostarts.size-1]; k++) {
						for (int l=0; l<bestcombination[aligntostarts.size-1][k].nummafs; l++) {
							alignmentsforthissection[position]=alignments[j][bestcombination[aligntostarts.size-1][k].mafslist[l]];
							string temp = alignmentsforthissection[position].makestring();
							goodmafs[j] << temp << endl;
							position++;
						}
					}
					
					int firstintchar;
					for (int k=0; k<numlettertrans[j]; k++) {
						if (alignedGFFletters[j][k].integer <= thisstart) {
							firstintchar=k;
						}
						else {
							break;
						}
					}
					int lastintchar;
					for (int k=firstintchar; k<numlettertrans[j]; k++) {
						if (alignedGFFletters[j][k].integer <= thisend) {
							lastintchar=k;
						}
						else {
							break;
						}
					}
					
					//See how much the alignment there is for this section so that we can decide if it is worth keeping.
					int numalignede=0;
					int numcoveredalignede=0;
					for (int k=firstintchar; k<=lastintchar; k++) {
						if (alignedGFFletters[j][k].character == 'e') {
							numalignede+=alignedGFFletters[j][k+1].integer-alignedGFFletters[j][k].integer;
							for (int l=0; l<nummafsinblocks; l++) {
								if (alignmentsforthissection[l].top(targetPosInDotFile[j]) >= alignedGFFletters[j][k].integer && alignmentsforthissection[l].bottom(targetPosInDotFile[j]) <= alignedGFFletters[j][k+1].integer) {
									numcoveredalignede+=min(alignmentsforthissection[l].top(targetPosInDotFile[j])+1,alignedGFFletters[j][k+1].integer) - max(alignedGFFletters[j][k].integer,alignmentsforthissection[l].bottom(targetPosInDotFile[j]));
								}
							}
						}
					}
					
					int numevidencee=0;
					int numcoveredevidencee=0;
					for (int m=0; m<numgffs4target; m++) {
						//cout << j << endl;
						//cout << numlettertrans[j] << endl;
						for (int k=0; k<numlettertrans[j]; k++) {
							//cout << targetGFFletters[m][k].integer << " " << thisstart << endl;
							if (targetGFFletters[m][k].integer <= thisstart) {
								firstintchar=k;
							}
							else {
								break;
							}
						}
						for (int k=firstintchar; k<numlettertrans[j]; k++) {
							if (targetGFFletters[m][k].integer <= thisend) {
								lastintchar=k;
							}
							else {
								break;
							}
						}
						for (int k=firstintchar; k<=lastintchar; k++) {
							if (targetGFFletters[j][k].character == 'e') {
								numevidencee+=targetGFFletters[m][k+1].integer-targetGFFletters[m][k].integer;
								//cout << "adding " << targetGFFletters[m][k+1].integer-targetGFFletters[m][k].integer << " to denominator" << endl;
								for (int l=0; l<nummafsinblocks; l++) {
									if (alignmentsforthissection[l].top(targetPosInDotFile[j]) >= targetGFFletters[m][k].integer && alignmentsforthissection[l].bottom(targetPosInDotFile[j]) <= targetGFFletters[m][k+1].integer) {
										numcoveredevidencee+=min(alignmentsforthissection[l].top(targetPosInDotFile[j])+1,targetGFFletters[m][k+1].integer) - max(targetGFFletters[m][k].integer,alignmentsforthissection[l].bottom(targetPosInDotFile[j]));
										//cout << "adding " << min(alignmentsforthissection[l].top(targetPosInDotFile[j])+1,targetGFFletters[m][k+1].integer) - max(targetGFFletters[m][k].integer,alignmentsforthissection[l].bottom(targetPosInDotFile[j])) << " to numerator" << endl;
									}
								}
							}
						}
					}
					
					//Output the beginning and end of the section followed by what fraction of the evidence that says "exon" is covered by an alignment.
					cout << thisstart << " " << thisend << " " /*<< bestcombinationscore[aligntostarts.size-1] << " " << (float)numcoveredalignede/numalignede << " "*/ << (float)numcoveredevidencee/numevidencee << endl;
					//cout << numcoveredevidencee << " " << numevidencee << endl;
					
				}
				else {
					//cout << thisstart << " " << thisend << " " << 0 << endl;
				}
			}
			
			thisstart = thisend;
			
		}
		
		//Choose the next value of i so that we don't skip any intchars.
		//cerr << "choosing nexti" << endl;
		int nexti=end;
		for (int j=0; j<numgffs4target; j++) {
			if (whichlettertrans[j]+1 <= numtargetlettertrans[j]-1) {
				if (targetGFFletters[j][whichlettertrans[j]+1].integer < nexti) {
					nexti = targetGFFletters[j][whichlettertrans[j]+1].integer;
				}
			}
		}
		nexti = max(nexti,i+1);
		//cerr << "chose nexti" << endl;
		
		i = nexti;
		//i++;
	}
	
	for (int j=0; j<numAlignedSpecies; j++) {
		goodmafs[j].close();
	}
	//cout << endl;

}