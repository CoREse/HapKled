/* File: crelib.cpp
 * Author: CRE
 * Last Edited: Fri Oct 28 16:18:03 2016
 */

#include "crelib.h"
#include <stdlib.h>
#include <stdarg.h>
#include <sys/time.h>
#include <time.h>

void cre::die(const char * Format, ...)
{
	va_list argptr;
	va_start(argptr, Format);
	vfprintf(stderr,Format,argptr);
	fprintf(stderr, "\n");
	exit(1);
}

void cre::warn(const char * Format, ...)
{
	va_list argptr;
	va_start(argptr, Format);
	vfprintf(stderr,Format,argptr);
	fprintf(stderr, "\n");
}

unsigned long long LastTime=0;

void cre::updateTime(const char * What, const char * OtherMessage, FILE* LogFile)
{
	unsigned long long ThisTime;
	if (LastTime==0)
	{
		LastTime=getTimeInMuse();
		fprintf(LogFile, "Performed time initial.\n");
	}
	ThisTime=getTimeInMuse();
	fprintf(LogFile,"%s cost %lld ms.",What, (ThisTime-LastTime)/1000);
	LastTime=ThisTime;
	if (OtherMessage!=NULL) fprintf(LogFile," %s",OtherMessage);
	fprintf(LogFile,"\n");
}

unsigned long long cre::getTimeInMuse()
{
	timeval VTime;
	gettimeofday(&VTime, NULL);
	return VTime.tv_sec*1000000ull+VTime.tv_usec;
}

cre::Logger::Logger(FILE* LogFile, int Verbosity)
:LogFile(LogFile), Verbosity(Verbosity)
{}

void cre::Logger::vlog(int Level, const char * Format, va_list argptr)
{
	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];
	timeval VTime;

	time (&rawtime);
	gettimeofday(&VTime, NULL);
	timeinfo = localtime (&(rawtime));

	strftime (buffer,80,"%Y-%m-%d %H:%M:%S",timeinfo);
	fprintf(LogFile, "%s.%d,%d",buffer,VTime.tv_usec/1000,VTime.tv_usec%1000);
	if (Level<2) fprintf(LogFile," [INFO] ");
	else if (Level==2) fprintf(LogFile," [DEBUG] ");
	else if (Level==3) fprintf(LogFile," [WARNING] ");
	else if (Level==4) fprintf(LogFile," [CRITICAL] ");
	else if (Level==5) fprintf(LogFile," [ERROR] ");

	// va_list argptr;
	// va_start(argptr, Format);
	vfprintf(stderr,Format,argptr);
	fprintf(stderr, "\n");
}

void cre::Logger::log(int Level, const char * Format, ...)
{
	va_list argptr;
	va_start(argptr, Format);
	vlog(Level,Format,argptr);
}

void cre::Logger::error(const char * Format, ...)
{
	va_list argptr;
	va_start(argptr, Format);
	vlog(5,Format,argptr);
}
void cre::Logger::critical(const char * Format, ...)
{
	va_list argptr;
	va_start(argptr, Format);
	vlog(4,Format,argptr);
}
void cre::Logger::warn(const char * Format, ...)
{
	va_list argptr;
	va_start(argptr, Format);
	vlog(3,Format,argptr);
}
void cre::Logger::debug(const char * Format, ...)
{
	if (Verbosity<2) return;
	va_list argptr;
	va_start(argptr, Format);
	vlog(2,Format,argptr);
}
void cre::Logger::verbose(const char * Format, ...)
{
	if (Verbosity<1) return;
	va_list argptr;
	va_start(argptr, Format);
	vlog(1,Format,argptr);
}
void cre::Logger::info(const char * Format, ...)
{
	va_list argptr;
	va_start(argptr, Format);
	vlog(0,Format,argptr);
}