/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONSOLE_COLOR_H_
#define CONSOLE_COLOR_H_

#include <iostream>

#if defined(WIN32)
#include <windows.h>
#endif

inline std::ostream& blue(std::ostream &s) {
#if defined(WIN32)
  HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
  SetConsoleTextAttribute(hStdout,
                          FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
#else
  s << "\x1b[0;34m";
#endif
  return s;
}

inline std::ostream& red(std::ostream &s) {
#if defined(WIN32)
  HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
  SetConsoleTextAttribute(hStdout, FOREGROUND_RED | FOREGROUND_INTENSITY);
#else
  s << "\x1b[0;31m";
#endif
  return s;
}

inline std::ostream& green(std::ostream &s) {
#if defined(WIN32)
  HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
  SetConsoleTextAttribute(hStdout, FOREGROUND_GREEN | FOREGROUND_INTENSITY);
#else
  s << "\x1b[0;32m";
#endif
  return s;
}

inline std::ostream& yellow(std::ostream &s) {
#if defined(WIN32)
  HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
  SetConsoleTextAttribute(hStdout,
                          FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_INTENSITY);
#else
  s << "\x1b[0;33m";
#endif
  return s;
}

inline std::ostream& white(std::ostream &s) {
#if defined(WIN32)
  HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
  SetConsoleTextAttribute(hStdout,
                          FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
#else
  s << "\x1b[0;37m";
#endif
  return s;
}

inline std::ostream& black_on_white(std::ostream &s) {
#if defined(WIN32)
  HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
  SetConsoleTextAttribute(hStdout,
                          BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE);
#else
  s << "\x1b[0;33m";
#endif
  return s;
}


#endif  // CONSOLE_COLOR_H_
