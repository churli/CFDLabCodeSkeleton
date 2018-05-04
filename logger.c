//
// Created by tommaso on 05/05/18.
//

#include "logger.h"
#include <stdio.h>
#include <stdarg.h>

/*
 * Logging machinery - how to:
 * 1) Make sure to call openLogFile() once at the beginning of main(), before calling any logmsg().
 * 2) To log a message both to log file and to console, call logmsg() as follows:
 *      logEvent(time, message, [arg1, arg2, ...]);
 *   time must be the current simulation time, message a printf-friendly string,
 *   then optional args in prtf style can be passed.
 * 3) Make sure to call closeLogfile() before exiting the main (this closes the file at OS level).
 */
static char* LOG_FILE_NAME = "sim.log";
static FILE* LOG_FILE;

void openLogFile()
{
    LOG_FILE = fopen(LOG_FILE_NAME, "w");
}

void logEvent(double t, char *fmt, ...)
{
    // Newline at the end of the message is included.
    va_list args;
    va_start(args,fmt);
    printf("[%12.9f] ", t);
    vprintf(fmt, args);
    printf("\n");
    va_end(args);
    va_start(args,fmt);
    fprintf(LOG_FILE, "[%12.9f] ", t);
    vfprintf(LOG_FILE, fmt, args);
    fprintf(LOG_FILE, "\n");
    va_end(args);
}

void logMsg(char *fmt, ...)
{
    // Newline at the end of the message is included.
    va_list args;
    va_start(args,fmt);
    printf("---> ");
    vprintf(fmt, args);
    printf("\n");
    va_end(args);
    va_start(args,fmt);
    fprintf(LOG_FILE, "---> ");
    vfprintf(LOG_FILE, fmt, args);
    fprintf(LOG_FILE, "\n");
    va_end(args);
}

void closeLogFile()
{
    fclose(LOG_FILE);
}
