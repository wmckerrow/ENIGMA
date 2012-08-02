#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <math.h>
using namespace std;

//finds out if a letter is a cds letter
bool iscds(char letter) {
	string cdsletters="efgEFGpqrPQRabcDONsT";
	if (cdsletters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

//finds out if a letter is an intergenic letter
bool isx(char letter) {
	string cdsletters="xXtS";
	if (cdsletters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

//finds out if a letter is a forward letter
bool isforward(char letter) {
	string forwardletters="{56~hl|sefgpqrdonijkuvwabct34zmy?}";
	if (forwardletters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

//finds out if a letter is a backward letter
bool isbackward(char letter) {
	string backwardletters="]$%0MYZTEFGPQRABCIJKUVWDONS&'\\HL^[";
	if (backwardletters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

bool is5prime(char letter) {
	string prime5letters="{56|&'^";
	if (prime5letters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

bool is3prime(char letter) {
	string prime3letters="34?]$%Z";
	if (prime3letters.find(letter) != string::npos) {
		return 1;
	}
	else {
		return 0;
	}
}

int main (int argc, char *argv[]) {
	//make sure we have enough input
	if (argc < 5) {
		cerr << "Please execute like " << argv[0] << " NonalignmentRows AligmentFileName NameForgffFiles PostionFile" << endl;
		exit(1);
	}
		
	//open stk file
	ifstream inFile;
	inFile.open(argv[2]);
	if (!inFile) {
		cerr << "Unable to open .stk file: " << argv[2] << endl;
		exit(1);
	}

	//read stk file into memory
	int UnusedRows=atoi(argv[1]);
	int numlines=0;
	string line;
	while (getline(inFile,line)) {
		numlines++;
	}
	inFile.close();
	numlines-=UnusedRows;
	string stockholmfile[numlines];
	inFile.open(argv[2]);
	if (UnusedRows>0) {
		for (int i=0; i<UnusedRows; i++) {
			getline(inFile,line);
		}
	}
	for (int i=0; i < numlines; i++) {
		getline(inFile,stockholmfile[i]);
		//cout << stockholmfile[i] << endl;
	}
	inFile.close();
	
	//find number of label columns and length of sequence
	int labelcolumns=0;
	while (stockholmfile[0][labelcolumns] != ' ') {
		labelcolumns++;
	}
	while (stockholmfile[0][labelcolumns] == ' ') {
		labelcolumns++;
	}
	int numsegpositions=stockholmfile[0].length()-labelcolumns+1;
	int segposition[numsegpositions];
	
	//open segment postions file and read into memory
	ifstream segfilein;
	segfilein.open(argv[4]);
	if (!segfilein) {
		cout << "Unable to open segmentation positions file " << argv[4] << endl;
		exit(1);
	}
	for (int i=0; i<numsegpositions; i++) {
		segfilein >> segposition[i];
	}
	segfilein.close();

	//output for debugging
	/*
	for (int i=0; i<stockholmfile[0].length()-labelcolumns+1; i++) {
		cout << segposition[i] << " ";
	}
	
	cout << stockholmfile[0] << endl;
	*/	
	 
	//write gff files
	ofstream outFile;
	string filenname;
	int labelchars;
	string label;
	int startgene;
	int endgene;
	bool foundcds;
	int startcds;
	int endcds;
	bool found5prime;
	int start5prime;
	int end5prime;
	bool found3prime;
	int start3prime;
	int end3prime;
	bool foundunaligned;
	int startunaligned;
	int endunaligned;
	int mrnanum;
	int times=0; //debugging line
	for (int i=0; i<numlines; i++) {
		labelchars=0;
		while (stockholmfile[i][labelchars]!=' ') {
			labelchars++;
		}
		label=stockholmfile[i].substr(0,labelchars);
		filenname = argv[3];
		filenname += "_" + label + ".gff";
		outFile.open(filenname.c_str());
		//outFile << "##gff-version 3" << endl;
		mrnanum=0;
		startgene=labelcolumns+1;
		while (startgene < stockholmfile[i].length()-1) {
			while (!(isx(stockholmfile[i][startgene-1]) && !isx(stockholmfile[i][startgene])) && startgene < stockholmfile[i].length()-1) {
				startgene++;
			}
			endgene=startgene;
			while (!isx(stockholmfile[i][endgene]) && endgene < stockholmfile[i].length()-1) {
				endgene++;
			}
			/*
			if (times < 5) {
				cout << startgene << " " << endgene << endl;
				times++;
			}
			 */
			if (startgene < stockholmfile[i].length()-1) {
				char strand = 'x';
				if (isforward(stockholmfile[i][startgene])) {
					strand = '+';
				}
				if (isbackward(stockholmfile[i][startgene])) {
					strand = '-';
				}
				mrnanum++;
				outFile << label << "	" << "." << '\t' << "gene" << '\t' << segposition[startgene-labelcolumns] << '\t' << segposition[endgene-labelcolumns]-1 << '\t' << "." << '\t' << strand << '\t' << "." << '\t' << "ID=gene" << mrnanum << endl;
				outFile << label << "	" << "." << '\t' << "mRNA" << '\t' << segposition[startgene-labelcolumns] << '\t' << segposition[endgene-labelcolumns]-1 << '\t' << "." << '\t' << strand << '\t' << "." << '\t' << "ID=mRNA" << mrnanum << ";parent=gene" << mrnanum << endl;
				foundcds=0;
				found5prime=0;
				found3prime=0;
				foundunaligned=0;
				for (int j=startgene; j<endgene; j++) {
					if (iscds(stockholmfile[i][j]) && !foundcds) {
						foundcds=1;
						startcds=j;
					}
					if (is5prime(stockholmfile[i][j]) && !found5prime) {
						found5prime=1;
						start5prime=j;
					}
					if (is3prime(stockholmfile[i][j]) && !found3prime) {
						found3prime=1;
						start3prime=j;
					}
					if (stockholmfile[i][j] == '*' && !foundunaligned) {
						foundunaligned=1;
						startunaligned=j;
					}
				}
				int fiveprimenum=0;
				int threeprimenum=0;
				int cdsnum=0;
				int unalignednum=0;
				if (found5prime) {
					while (start5prime < endgene) {
						end5prime=start5prime;
						while (stockholmfile[i][end5prime] == '5' && end5prime <= endgene -1 ) {
							end5prime++;
						}
						end5prime--;
						fiveprimenum++;
						outFile << label << "	" << "." << '\t' << "exon" << '\t' << segposition[start5prime-labelcolumns] << '\t' << segposition[end5prime-labelcolumns+1]-1 << '\t' << "." << '\t' << strand << '\t' << "." << '\t' << "ID=five_prime" << fiveprimenum << ";Parent=mRNA" << mrnanum << endl;
						outFile << label << "	" << "." << '\t' << "five_prime_UTR" << '\t' << segposition[start5prime-labelcolumns] << '\t' << segposition[end5prime-labelcolumns+1]-1 << '\t' << "." << '\t' << strand << '\t' << "." << '\t' << "ID=five_prime" << fiveprimenum << ";Parent=mRNA" << mrnanum << endl;
						start5prime=end5prime+1;
						while (stockholmfile[i][start5prime] != '5' && start5prime < endgene ) {
							start5prime++;
						}
					}
				}
				if (foundcds) {
					while (startcds < endgene) {
						endcds=startcds;
						while (iscds(stockholmfile[i][endcds]) && endcds <= endgene -1 ) {
							endcds++;
						}
						endcds--;
						cdsnum++;
						outFile << label << "	" << "." << '\t' << "exon" << '\t' << segposition[startcds-labelcolumns] << '\t' << segposition[endcds-labelcolumns+1]-1 << '\t' << "." << '\t' << strand << '\t' << "." << '\t' << "ID=exon" << cdsnum << ";Parent=mRNA" << mrnanum << endl;
						outFile << label << "	" << "." << '\t' << "CDS" << '\t' << segposition[startcds-labelcolumns] << '\t' << segposition[endcds-labelcolumns+1]-1 << '\t' << "." << '\t' << strand << '\t' << "." << '\t' << "ID=CDS" << cdsnum << ";Parent=mRNA" << mrnanum << endl;
						startcds=endcds+1;
						while (!iscds(stockholmfile[i][startcds]) && startcds < endgene ) {
							startcds++;
						}
					}
				}
				if (found3prime) {
					while (start3prime < endgene) {
						end3prime=start3prime;
						while (stockholmfile[i][end3prime] == '3' && end3prime <= endgene -1 ) {
							end3prime++;
						}
						end3prime--;
						threeprimenum++;
						outFile << label << "	" << "." << '\t' << "exon" << '\t' << segposition[start3prime-labelcolumns] << '\t' << segposition[end3prime-labelcolumns+1]-1 << '\t' << "." << '\t' << strand << '\t' << "." << '\t' << "ID=three_prime" << threeprimenum << ";Parent=mRNA" << mrnanum << endl;
						outFile << label << "	" << "." << '\t' << "three_prime_UTR" << '\t' << segposition[start3prime-labelcolumns] << '\t' << segposition[end3prime-labelcolumns+1]-1 << '\t' << "." << '\t' << strand << '\t' << "." << '\t' << "ID=three_prime" << threeprimenum << ";Parent=mRNA" << mrnanum << endl;
						start3prime=end3prime+1;
						while (stockholmfile[i][start3prime] != '3' && start3prime < endgene ) {
							start3prime++;
						}
					}
				}
				if (foundunaligned) {
					while (startunaligned < endgene) {
						endunaligned=startunaligned;
						while (stockholmfile[i][endunaligned] == '*' && endunaligned <= endgene -1 ) {
							endunaligned++;
						}
						endunaligned--;
						unalignednum++;
						outFile << label << "	" << "." << '\t' << "unaligned" << '\t' << segposition[startunaligned-labelcolumns] << '\t' << segposition[endunaligned-labelcolumns+1]-1 << '\t' << "." << '\t' << strand << '\t' << "." << '\t' << "ID=unaligned" << unalignednum << ";Parent=mRNA" << mrnanum << endl;
						startunaligned=endunaligned+1;
						while (stockholmfile[i][startunaligned] != '*' && startunaligned < endgene ) {
							startunaligned++;
						}
					}
				}
			}
			startgene=endgene;
		}
		outFile.close();
	}
	
	return 0;
	
}