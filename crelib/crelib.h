#pragma once
#ifndef CRELIB_H
#define CRELIB_H
#include <stdio.h>
#include <stdlib.h>
#include "defines.h"
namespace cre
{
	void die(const char * Format=NULL, ...);
	void updateTime(const char * What, const char * OtherMessage=NULL, FILE * LogFile=stderr);
	unsigned long long getTimeInMuse();//get time in microseconds. Muse no minasan, daisuki!^.^

	void warn(const char * Format, ...);

	class Logger
	{
		public:
		FILE* LogFile;
		int Verbosity;//0:info, 1:verbose, 3: debug
		void log(int Level, const char * Format, ...);//Level: 0:info, 1:verbose, 3: debug, 4: warn, 5: critical, 6: error
		void vlog(int Level, const char * Format, va_list argptr);//Level: 0:info, 1:verbose, 3: debug, 4: warn, 5: critical, 6: error

		void error(const char * Format, ...);//After this, the program shall die.
		void critical(const char * Format, ...);//Like error, but the program don't need to die.
		void warn(const char * Format, ...);//Waring.
		void debug(const char * Format, ...);//Information for debugging.
		void verbose(const char * Format, ...);//More information
		void info(const char * Format, ...);//Information for users.

		Logger(FILE* LogFile=stderr, int Verbosity=0);
	};
}
#endif
