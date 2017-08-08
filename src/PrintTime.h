/*
 * PrintTime.h
 *
 *  Created on: Mar 7, 2017
 *      Author: cochez
 */

#ifndef PRINTTIME_H_
#define PRINTTIME_H_

//taken from by http://stackoverflow.com/a/16358264

#include <iostream>
#include <ctime>

using namespace std;

inline string currentTime() {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, sizeof(buffer), "%Y-%m-%d %I:%M:%S ", timeinfo);
	std::string str(buffer);

	return str;

}

#endif /* PRINTTIME_H_ */
