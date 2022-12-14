/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Parser class                                                            //
//                                                                         //
// Burkhard Militzer                                    Urbana 4-1-99      //
// based on earlier version by John Shumway                                //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "Standard.h"
#include "Parser.h"

Parser::Parser(const string & fn):
  ifs(fn.c_str()),
  is(ifs),
  filename(fn),
  ignoreComments(true),
  ignoreEmptyLines(true),
  nLinesRead(0),
  nNonEmptyLinesRead(0)
{
  if (!is)
    error("Parser: could not open file",filename);
}
 
Parser::Parser(istream& is1): 
  is(is1),
  ignoreComments(true),
  ignoreEmptyLines(true),
  nLinesRead(0),
  nNonEmptyLinesRead(0)
{}

// There can be a problem with the tyep of npos
// int          using KCC-3.4d
// unsigned int using gcc
// Earlier code used KCC defines "_CFE", use this as indicator
/*
#ifndef _CFE
unsigned int find_first_not_of(const string & s, const char c) {
  unsigned int i=0;
  while (i<s.length() && s[i]==c)
    i++;
  if (i==s.length()) i=s.npos;
  return i;
}
#else
int find_first_not_of(const string & s, const char c) {
  int i=0;
  while (i<s.length() && s[i]==c)
    i++;
  if (i==s.length()) i=s.npos;
  return i;
}
#endif
*/
string::size_type find_first_not_of(const string & s, const char c) {
  string::size_type i=0;
  while (i<s.length() && s[i]==c)
    i++;
  if (i==s.length()) i=s.npos;
  return i;
}

bool Parser::ReadLine() {
  words.clear();
  string word; 
  do {
    while(getline(is, lineString)) {
      nLinesRead++;
      if (!ignoreEmptyLines) break;
      //      if ( lineString.length() != 0 && find_first_not_of(lineString,' ')!=lineString.npos ) {
      // This following line should be used but does not work with KCC earlier than 3.4f
      // but must use if for PGI 5.1 !!!
      if ( lineString.length() != 0 && lineString.find_first_not_of(' ')!=lineString.npos ) {
	break;
      }
    }
    if (!is) return is;

    istringstream line_stream(lineString);
    line_stream >> word;
  } while (word[0]=='#' && ignoreComments);

  nNonEmptyLinesRead++;

  istringstream line_stream(lineString);
  line_stream >> word;
  while (word.length()!=0 && line_stream) {
    words.push_back(word);
    //    Write2(words.size(),words[words.size()-1]);
    line_stream >> word;
  }  
  // if a line only contains a "tab" then CUPID will have problems, fix later
  return is; //Return boolean status of the input stream.
}

const string& Parser::GetString(const int i) const {
  if (i<0 || i>=int(words.size())) 
    error("Index out of bounds in Parser",i,words.size(),lineString,filename);
  return words[i];
}

void Parser::CheckNotNumber(const string & s) const {
  const string numbers = "0123456789";
  if (s.find_first_not_of(numbers)==string::npos) 
    error("Parser: Name string contains only numbers",s,filename);
}

// Check if it contains not just numbers
const string& Parser::GetNameString(const int i) const {
  const string & s=GetString(i);
  CheckNotNumber(s);
  return s;
}

const string Parser::GetStringLowerCase(const int i) const {
  string s=GetString(i);
  string sl;
  for(string::const_iterator p=s.begin();p!=s.end();++p) {
    //    s.replace(*p,*p+1,tolower(*p));
    sl += tolower(*p);
  }
  //  Write2(s,sl);
  return sl;
}

const string Parser::GetNameStringLowerCase(const int i) const {
  string s=GetStringLowerCase(i);
  CheckNotNumber(s);
  return s;
}

const string Parser::GetStringUpperCase(const int i) const {
  string s=GetString(i);
  return UpperCase(s);
}

const string Parser::GetNameStringUpperCase(const int i) const {
  string s=GetStringUpperCase(i);
  CheckNotNumber(s);
  return s;
}

int Parser::GetInt(const int i) const {
  if (i<0 || i>=int(words.size())) 
    error("Index out of bounds in Parser",i,words.size(),lineString,filename);

  const string intChars(" +-0123456789");
  //  if (size_type j=words[i].find_first_not_of(intChars) != string::npos)
  if (int j=words[i].find_first_not_of(intChars) != string::npos)
    error("Parser: Not an integer ","\""+words[i]+"\"",i,j+1,"\""+lineString+"\"","\""+filename+"\"");
  
  istringstream wordstream(words[i]);
  int value;
  wordstream >> value;
  return value;
}

double Parser::GetDouble(const int i) const {
  if (i<0 || i>=int(words.size())) 
    error("Index out of bounds in Parser",i,words.size(),lineString,filename);

  const string doubleChars(" +-0123456789.edDE");
  //  if (size_type j=words[i].find_first_not_of(doubleChars) != string::npos)
  if (int j=words[i].find_first_not_of(doubleChars) != string::npos)
    error("Parser: Not a double ",words[i],j+1,filename);
  
  istringstream wordstream(words[i]);
  double value;
  wordstream >> value;
  return value;
}

void Parser::ProcessKeyword(const string & keyword) {
  ReadLineSafely();

  string kw=GetNameStringLowerCase(0);
  if (kw!=LowerCase(keyword))
    error("Parser: Wrong keyword found",keyword,GetLineString(),filename);
}

void Parser::ProcessKeyword(const string & keyword, const int n) {
  ProcessKeyword(keyword);

  if (n!=GetNWords())
    error("Parser: Incorrect number of arguments",n,GetNWords(),GetLineString(),filename);

}

void Parser::ProcessLine(const int n) {
  ReadLineSafely();

  if (n!=GetNWords())
    error("Parser: Incorrect number of arguments",n,GetNWords(),GetLineString(),filename);

}

void Parser::PutLineBack() {
  is.putback('\n');
  for(int i=lineString.length()-1; i>=0; i--) {
    is.putback(lineString[i]);
  }
  nLinesRead--;
  if (words.size()>0) nNonEmptyLinesRead--; 
}
