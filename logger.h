//
// Created by tommaso on 05/05/18.
//

#ifndef SIM_LOGGER_H
#define SIM_LOGGER_H

void openLogFile();
void logEvent(double t, char *fmt, ...);
void logMsg(char *fmt, ...);
void closeLogFile();

#endif //SIM_LOGGER_H
