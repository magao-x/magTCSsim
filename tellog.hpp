/** \file tellog.hpp
  * \author Jared R. Males
  * \brief Declaration and definition of logging utilities for magTCSsim
  */

//***********************************************************************//
// Copyright 2018 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of magTCSsim.
//
// magTCSsim is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// magTCSsim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with magTCSsim.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef tellog_hpp
#define tellog_hpp

#define lmsg_size 100

///The common logger for magtelescope applications, to use for tracking and debugging the AO software.
/** Holds a static FILE * to the open log file to write messages during execution.
  * To set the log filename, call with msg == 0, and pass the name as setlog.  This opens the file.
  * Also uses a static char * to hold the application name.  Set the same way as the log name.
  * If msg != 0, then stapp and setlog are ignored.
  * Call with all params == 0 to close the log file.
  * 
  * \retval 0 = success. (always).
  * 
  */
inline
int tellog( const char *msg,    ///< [in] the message to log.
              const char *fname,  ///< [in] the filename from which the log entry is generated.
              int lineno,         ///< [in] the line number from which the log entry is generated.
              const char *setapp, ///< [in] the application name to use in the log entry.  Ignored unless msg==0.
              const char *setlog  ///< [in] the log file name to open and use.  Ignored unless msg == 0.
            )
{
   static FILE * lf = 0;
   static char app[20];
   
   char logtime[100]; 

   time_t tl;
   tl = time(0);
   strftime(logtime, 100, "%m/%d/%Y %H:%M:%S (%Z)", localtime(&tl));

   if(msg == 0)
   {
      if(setlog != 0)
      {
         if(lf != 0)
         {
            fclose(lf);
         }
         lf = fopen(setlog, "a");
      }
      if(setlog == 0 && lf != 0)
      {
         fclose(lf);
         lf = 0;
      }
      if(setapp !=0)
      {
         strncpy(app, setapp, 20);
      }
   }
   else
   {
      if(lf != 0)
      {
         fprintf(lf, "%s [%s]: %s (in %s at line %i)\n", app, logtime, msg, fname, lineno);
         fflush(lf);
      }
      else
      {
         fprintf(stderr, "%s log [%s]: %s (in %s at line %i)\n", app, logtime, msg, fname, lineno);
      }
   }
   return 0;
}//__tellog(char *msg, char *fname, int lineno, char *setlog)

#ifndef TELLOG
///Macro for making a log entry, calls \ref tellog.
#define TELLOG(msg) tellog(msg, __FILE__, __LINE__, 0, 0)
#endif


#endif //tellog_hpp
