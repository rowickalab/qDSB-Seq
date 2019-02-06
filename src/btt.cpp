#include "btt.h"
//#define FAKE_ELAND 1
int main (int argc, char**argv){
   int argcnum=1;
   string mapFileName="";
   string argshelp="usage:> btt <bowtie output file>";
   if (argc <= argcnum){
     cout << endl << argshelp;
   //  cout << endl << "Please Input File Name which contains mapped addresses: ";
   //  cin >> mapFileName;
   }
   else
     mapFileName=argv[argcnum++];
   if (checkFile(mapFileName)!=1) terminate(1);
   
   bool verifyWithGenome=0;
   string genomeFileName="";
   string fastaData="";
   unsigned int verificationChrNum=0;
   if (argc <= argcnum){
     cout << endl << "Verification Disabled";
   }
   else{
     verifyWithGenome=1;
     genomeFileName=argv[argcnum++];
     if (checkFile(genomeFileName)!=1) terminate(1);
     boost::iostreams::stream_buffer<boost::iostreams::file_source> fastaFileSource(genomeFileName.c_str());
     istream fastaFile(&fastaFileSource);
     bool headerRead=0;
     string f_header="";
     while (!fastaFile.eof()){
       string f_line;
       getline (fastaFile, f_line);
       if (f_line.size() == 0)
	  break;
       if (f_line[0] == '>'){ //Header Located
	     f_header = f_line;
	     unsigned c=extract_chrnum(f_header,0);
	     //cout << endl << c << f_header; WAITUSER;
	     if (c>0){
	       if (headerRead==0){
		  cout << endl << "Using chromosome " << c << " as a reference for verification";
		  headerRead=1;
		  verificationChrNum=c;
		  //WAITUSER;
	       }
	       else if (headerRead==1)
		 break;
	     }
       }
       else{
	 if (headerRead==1)
	   fastaData+=f_line;
       }
     }
     if (headerRead==0){
       cout << endl << "Given fasta file could not be parsed successfully; Verification Disabled";
       verifyWithGenome=0;
     }
   }
//   cout << endl << fastaData.size(); WAITUSER;
   bool verified=0;
   boost::iostreams::stream_buffer<boost::iostreams::file_source> mapFileMapSource(mapFileName.c_str());
   istream mapFile (&mapFileMapSource);
   string transMapFileName=mapFileName+".btt";
   boost::iostreams::stream_buffer<boost::iostreams::file_sink> transFileMapSink(transMapFileName.c_str());
   ostream mapFileTrans(&transFileMapSink);
   #ifdef FAKE_ELAND
   string elandMapFileName=mapFileName+".eland";
   boost::iostreams::stream_buffer<boost::iostreams::file_sink> elandMapFileSink(elandMapFileName.c_str());
   ostream mapFileEland(&elandMapFileSink);
   #endif
   while (!mapFile.eof()){
      string line;
      getline(mapFile,line);
      //cout << endl << line; 
      unsigned int columns=0;
      size_t tab;
      string lineNum, strand, chrNum, chrPos, sequence, quality, repeats, changeString;
      do{
	tab =  line.find_first_of('\t');
	columns++;
	switch (columns){
	    case 1: lineNum = line.substr(0, tab);
	    break;
	    case 2: strand = line.substr(0, tab);
	    break;
	    case 3: chrNum = line.substr(0, tab);
	    break;
	    case 4: chrPos= line.substr(0, tab);
	    break;
	    case 5: sequence = line.substr(0, tab);
	    break;
	    case 6: quality = line.substr(0, tab);
	    break;
	    case 7: repeats = line.substr(0, tab);
	    break;
	    case 8: changeString = line.substr(0,tab);
	    break;
	    default:
	    break;
	  }
	  //cout << endl << "TAB pos = " << tab;
	  if (tab!=ULONG_MAX) line = line.substr(tab+1, line.size());
      }while(tab!=ULONG_MAX); 
      //WAITUSER;
      /*
	Base	Name	Bases Represented	Complementary Base
	A	Adenine			A	T
	T	Thymidine		T	A
	U	Uridine(RNA only)	U 	A
	G	Guanidine		G	C
	C	Cytidine		C	G
	Y	pYrimidine		C T	R
	R	puRine			A G	Y
	S	Strong(3Hbonds)		G C	S*
	W	Weak(2Hbonds)		A T	W*
	K	Keto			T/U G	M
	M	aMino			A C	K
	B	not A			C G T	V
	D	not C			A G T	H
	H	not G			A C T	D
	V	not T/U			A C G	B
	N	Unknown			A C G T	N
	*/
      if (columns==8 || columns==7){
	#ifdef FAKE_ELAND
	  string elandName, elandSeq, elandType, elandNum0, elandNum1, elandNum2, elandGenomeFile, elandPos;
	    //		1	2	3		4	5	6		7		8
	  string elandDir, elandInterpret, elandSubst1, elandSubst2;
	    //		9	10		11		12
	  elandName="FAKE:"+lineNum;
	  elandSeq=sequence;
	  (repeats=="0"? elandType="U": elandType="R");
	  if (changeString.size() == 0)
	      elandType+="0";
	  else{
	      size_t commas = changeString.find_first_of(',');
	      if (commas==ULONG_MAX) elandType+="1";
	      if (commas!=ULONG_MAX) elandType+="2";
	  }
	  elandNum0 = "0";
	  elandNum1 = "0";
	  elandNum2 = "0";
	  if (elandType=="U1") elandNum1="1";
	  else if (elandType=="U2") elandNum2="1";
	  else if (elandType=="R0") elandNum0=repeats; 
	  else if (elandType=="R1") elandNum1=repeats;
	  else if (elandType=="R2") elandNum2=repeats;
	  elandGenomeFile=chrNum;
	  elandPos=chrPos;
	  switch ((char)strand[0]){
	      case '-':elandDir="R";
	      break;
	      case '+':elandDir="F";
	      break;
	      default: cout << endl << "Strand information not found";
	      elandDir="";
	      break;
	  }
	  #endif 
	  if (strand=="-") {
	    //cout << endl << lineNum << "\t" << strand << "\t" << sequence ; WAITUSER;
	    string sequenceComplement;
	    for (unsigned int n=sequence.size();n>0;--n){
	      char bp;
	      //cout << endl << sequence[n-1];
	      switch(sequence[n-1]){
		  case 'A':bp='T';
		  break;
		  case 'C':bp='G';
		  break;
		  case 'G':bp='C';
		  break;
		  case 'T':bp='A';
		  break;
  // 		case 'U':bp='A';
  // 		break;
		  case 'Y':bp='R';
		  break;
		  case 'R':bp='Y';
		  break;
		  case 'K':bp='M';
		  break;
		  case 'M':bp='K';
		  break;
		  case 'B':bp='V';
		  break;
		  case 'D':bp='H';
		  break;
		  case 'H':bp='D';
		  break;
		  case 'V':bp='B';
		  break;
		  case 'N':bp='N';
		  break;
		  case 'S':bp='S';
		  break;
		  case 'W':bp='W';
		  break;
		  default:bp=sequence[n]; cout << endl << "Not expecting bp " << bp; WAITUSER;
		  break;
	      }
	      //cout << "\t complement = " << bp; WAITUSER;
	      sequenceComplement+=bp;
	    }
	    unsigned int newChrPos=0;
	    stringstream ss_chrPos(chrPos);
	    ss_chrPos >> newChrPos;
	    newChrPos+=sequence.size()-1;
	    mapFileTrans << lineNum << "\t"<< strand << "\t" << chrNum << "\t" << \
	    newChrPos << "\t" << sequenceComplement << "\t" << quality << "\t" << repeats << changeString << endl;
	    #ifdef FAKE_ELAND
	      mapFileEland << elandName << "\t" << sequenceComplement << "\t" <<  elandType << "\t" << elandNum0 << "\t" << \
	      elandNum1 << "\t" << elandNum2 << "\t" << elandGenomeFile << "\t" << newChrPos << "\t" << elandDir << endl;
	    #endif
	    if (verifyWithGenome==1&&verified==0){
	      if (verificationChrNum==extract_chrnum(chrNum,0)){
		unsigned int i_chrPos=0;
		stringstream ss_chrPos(chrPos);
		ss_chrPos >> i_chrPos;
		string sequenceSegmentW = fastaData.substr(i_chrPos,sequence.size());
		cout << endl << i_chrPos << "\t" << sequenceSegmentW << "\t" << sequence << "\t" << sequenceComplement ;
		if (StringToUpper(sequenceSegmentW)==sequence){
		   cerr << endl << "File looks untranslated";
		   verified=1;
		}
	      }
	    }
	  }
	  else if (strand=="+"){
	    mapFileTrans << lineNum << "\t"<< strand << "\t" << chrNum << "\t" << \
	    chrPos << "\t" << sequence << "\t" << quality << "\t" << repeats << changeString << endl;
   	    #ifdef FAKE_ELAND
	      mapFileEland << elandName << "\t" << sequence << "\t" <<  elandType << "\t" << elandNum0 << "\t" << \
	      elandNum1 << "\t" << elandNum2 << "\t" << elandGenomeFile << "\t" << chrPos << "\t" << elandDir << endl;
	    #endif
	  }
      }
      else
	mapFileTrans << line << endl;
   }
  if (transFileMapSink.is_open())transFileMapSink.close();
  #ifdef FAKE_ELAND
    if (elandMapFileSink.is_open())elandMapFileSink.close();
  #endif
  cout << endl;
}
