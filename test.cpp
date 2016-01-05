#include <iostream>
using namespace std;

#include <stdio.h>
#include <stdlib.h>

#include "options.h"

time_t rawtime;

int main(int argc, char *argv[])
{
	bool help, version, seconds, time;
	int year = 2016;
	command::CmdLine cmd;
	cmd.add( command::make_option('h', "help", help) );
	cmd.add( command::make_option('v', "version", version) );
	cmd.add( command::make_option('s', "seconds", seconds) );
	cmd.add( command::make_option('t', "time", time) );
	cmd.add( command::make_option('y', "year", year) );
	try
	{
		cmd.process(argc, argv);
	}
	catch (exception &e)
	{
		cerr << e.what();
	}
	if (cmd.opt[0]->arg != NULL)
		printf("Help command activated.");
	if (cmd.opt[1]->arg != NULL)
		printf("This is version 4.3.5");
	if (cmd.opt[2]->arg != NULL)
	{
		std::time(&rawtime);
		printf("Time (secondes since 1970) is %s", ctime(&rawtime));
	}
	if (cmd.opt[3]->arg != NULL)
	{
		struct tm *timeinfo = localtime(&rawtime);
		printf("Current time is %s", asctime(timeinfo));
	}
	if (cmd.opt[4]->arg != NULL)
	{
		struct tm y2k = {0};
		int seconds;
		int year = int(cmd.opt[4]->arg);

		y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
		y2k.tm_year = year - 1900; y2k.tm_mon = 0; y2k.tm_mday = 1;

		seconds = int(difftime(rawtime, mktime(&y2k)));
		cout << "Time in seconds since January, 1st of " << year << " is " << seconds << endl;
	}


	return(EXIT_SUCCESS);
}