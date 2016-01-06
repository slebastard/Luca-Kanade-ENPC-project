#ifndef OPTIONS_H
#define OPTIONS_H

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <ctime>
using namespace std;

namespace command
{
	class Option 
	{
	public:
		int *arg = NULL;
		virtual void filter(int& argc, char *argv[]) = 0;
	};

	template <typename T>
	class OptionField: public Option
	{
	public:
		T *arg = NULL;
		char *opt_name[];
		OptionField<T>(char *name[]);
		void filter(int& argc, char *argv[]);
		T* argument(){ return arg };
	};

	class CmdLine
	{
	public:
		vector<Option*> opt;
		CmdLine();
		CmdLine(vector<Option*> opt0);
		void add(Option* opt0);
		void process(int& argc, char *argv[]);
	};

	template <typename T>
	Option* make_option(char c, string s1, T var);
}

string char2string(char* c);
string last_part(char *c, char limit);

#endif