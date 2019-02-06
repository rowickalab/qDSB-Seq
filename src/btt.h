#include<iostream>
#include<fstream>
#include<climits>
#include<vector>
#include<string>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<cmath>
#include<ctype.h>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#define WAITUSER cout << endl << "Press Return";cin.get();
#define X_CHROMOSOME 1000
#define Y_CHROMOSOME 1001
using namespace std;

void terminate(int retval)
{
  cout << endl << "Terminating! because of Error " << retval << endl;
  exit(retval);
}

int checkFile(string fileName)
{
  ifstream fileOpenTest;
  int isFileOpen=-1;
  fileOpenTest.open(fileName.c_str());
  if (fileOpenTest.is_open())
  {
    isFileOpen=1;
    fileOpenTest.close();
  }
  return isFileOpen;
}

unsigned int extract_chrnum(string s_chrNum, bool XY){
  unsigned int len=s_chrNum.size();
  unsigned int count =0;
  unsigned int chrNum=0;
  //cout << endl << "LEN" << len;WAITUSER;
  for (unsigned int n=len; n>0; --n)
  {
    char c = s_chrNum[n-1];
    if (c >= '0' && c<='9'){
	++count;
	if (count==1)
	{
	  chrNum+= (c-'0');
	}
	else if (count==2)
	{
	  chrNum+= (c-'0')*10;
	}
	else if (count==3)
	{
	  chrNum+= (c-'0')*100;
	}
	else if (count > 3){
	  chrNum=0;
	  count=0;
	}
    }
    //cout << endl << "char:" << c << " len:" << n << " chrnum" << chrNum;WAITUSER;
  }
    
  if (chrNum==0 && XY==0)
  { 
    bool counting=0;
    unsigned int count=0;
    for (unsigned int n=len; n>0; --n)
    {
      char c = s_chrNum[n-1];
      if (c == 'X' || c== 'V' || c=='I'){ 
	if (counting==0) counting = 1;
	++count;
	if (count==1)
	{
	  if (c=='X')
	    chrNum+=10;
	  else if (c=='V')
	    chrNum+=5;
	  else if (c=='I')
	    chrNum+=1;
	}
	else if (count>1 && count <=4)
	{
	 // cout << endl << " " << count << " " << char(s_chrNum[n]); WAITUSER;
	  if (c < char(s_chrNum[n])){	  
	    if (c=='X')
	      chrNum-=10;
	    else if (c=='V')
	      chrNum-=5;
	    else if (c=='I')
	      chrNum-=1;
	  }
	  else if (c >= char(s_chrNum[n])){
	    if (c=='X')
	    chrNum+=10;
	  else if (c=='V')
	    chrNum+=5;
	  else if (c=='I')
	    chrNum+=1;
	  }
	}
	else if (count >4){
	  chrNum=0;
	  counting=0;
	  count=0;
	}
      }
      else{ 
	if (chrNum==0){
	  counting=0;
	  chrNum=0;
	}
	else
	  break;
      }
    }
    //cout << endl << "Roman Numeral " << chrNum; WAITUSER;
  }
  else if (chrNum==0 && XY==1)//XY chromosome
    for (unsigned int n=len; n>0; --n){
      char c=s_chrNum[n-1];
      if (c=='X')
	chrNum=X_CHROMOSOME;
      else if (c=='Y')
	chrNum=Y_CHROMOSOME;
    }
  return chrNum;
}

string convert_chrnum_to_string(unsigned int chrNum){
  string s_chrNum="";
  stringstream ss_chrNum;
  if (chrNum <10){
    string s_chrNumTemp="";
    ss_chrNum << chrNum;
    ss_chrNum >> s_chrNumTemp;
    s_chrNum = "0" + s_chrNumTemp;
  }
  else if (chrNum <1000){
    ss_chrNum << chrNum;
    ss_chrNum >> s_chrNum;
  }
  else {
    char c_chrNum=chrNum-1000+'X';
    s_chrNum=c_chrNum;
  }
    
  return s_chrNum;
}


string StringToUpper(string strToConvert) {//change each element of the string to upper case
  for(unsigned int i=0;i<strToConvert.length();i++){
    strToConvert[i] = toupper(strToConvert[i]);
  }
  return strToConvert;//return the converted string
}

string StringToLower(string strToConvert)
 
      {//change each element of the string to lower case
  
         for(unsigned int i=0;i<strToConvert.length();i++)
  
         {
 
            strToConvert[i] = tolower(strToConvert[i]);
 
         }
 
         return strToConvert;//return the converted string
 
      } 