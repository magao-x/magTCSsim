/** \file magtelUtils.hpp
  * \author Jared R. Males
  * \brief Declaration of some utilities for magTCSsim
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

#ifndef magtelUtils_hpp
#define magtelUtils_hpp

#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/time.h>


inline
int tcs_fputs( char *string,
               FILE *fp
             )
{
   if (fseek (fp, 0L, SEEK_CUR) != 0) 
   {
       ;
   }
   if (fputs (string, fp) < 0) 
   {
       return (-1);
   }
   if (fputc ('\n', fp) < 0) 
   {
       return (-2);
   }
   if (fflush (fp) != 0) 
   {
       return (-3);
   }
   return 0;
}

int tcs_fgets( char *string, 
               int length,
               FILE *fp
             )
{
   if (fseek (fp, 0L, SEEK_CUR) != 0)              
   {
      //Do nothing?
      ;
   }
   
   if ( fgets (string, length, fp) == NULL ) 
   {
       return -1;
   }
   
   if (string[strlen(string)-1] == '\n') string[strlen(string)-1] = 0;
   
   return 0;
}

inline
int deg2dms(double hms[3], double ang)
{
   int tmp;

   tmp = (int) ang;

   hms[0] = tmp;

   tmp = (int) ((ang-hms[0])*60.0);

   hms[1] = tmp;

   hms[2] = ((ang-hms[0]-hms[1]/60.0)*3600.0);

   if(hms[2] >= 59.995)
   {
      hms[1] += 1;
      hms[2] = 0;
   }

   return 0;

}

inline
int mts_format_date(char dtstr[11], int dtarr[3])
{
   int i;
   char tmp[5];
   strncpy(dtstr, "0000-00-00",11);

   //Year first
   if(dtarr[0] > 9999 || dtarr[0] <1000)
   {
      fprintf(stderr, "Attempt to format invalid year, must be 4 digits.\n");
      return -1;
   }
   //First turn it into a string.
   snprintf(tmp, 5,"%i", dtarr[0]);

   //then copy
   for(i=0;i<4; i++) dtstr[i] = tmp[i];

   for(i=1; i<3; i++)
   {
      if(dtarr[i] > 99 || dtarr[i] < 0)
      {
         fprintf(stderr, "Attempt to format day or month value longer than 2 digits.\n");
         return -1;
      }

      //First turn it into a string.
      snprintf(tmp, 3,"%i", dtarr[i]);

      //Now copy the digits
      if(dtarr[i] > 9) 
      {  //2 Digits
         dtstr[3*i+2] = tmp[0];
         dtstr[3*i+3] = tmp[1];
      }
      else dtstr[3*i + 3] = tmp[0]; //only 1 digit
   }

   return 0;
}//int mts_format_date(char dtstr[11], int dtarr[3])

inline
int mts_format_time_uint(char tmstr[10], double tmarr[3])
{
   int i, negval = 0;

   if(tmarr[0] <=0 && tmarr[1] <= 0 && tmarr[2] <= 0) negval = 1;
   else if(!(tmarr[0] >=0 && tmarr[1] >= 0 && tmarr[2] >= 0))
   {
      fprintf(stderr, "Inconsistent signs in time/angle array.  Arith will fail.\n");
      return -1;
   }

   if(tmarr[2] >= 59.5)
   {
      tmarr[1] += 1;
      tmarr[2] = 0;
   }
   for(i=0; i<3; i++)
   {
      if(fabs(tmarr[i]) > 99)
      {
         fprintf(stderr, "Time value longer than 3 digits. %i %f\n", i, tmarr[i]);
         return -1;
      }
   }
   if(negval) sprintf(tmstr, "-%02.0f:%02.0f:%02.0f", fabs(tmarr[0]), fabs(tmarr[1]), fabs(tmarr[2]));
   else sprintf(tmstr, "%02.0f:%02.0f:%02.0f", fabs(tmarr[0]), fabs(tmarr[1]), fabs(tmarr[2]));

   return 0;
}//int mts_format_time_uint(char tmstr[9], int tmarr[3])

inline
int mts_format_time_doub(char tmstr[13], double tmarr[3])
{
   int i, negval = 0;

   if(tmarr[0] <=0 && tmarr[1] <= 0 && tmarr[2] <= 0) negval = 1;
   else if(!(tmarr[0] >=0 && tmarr[1] >= 0 && tmarr[2] >= 0))
   {
      fprintf(stderr, "Inconsistent signs in time/angle array.  Arith will fail.\n");
      return -1;
   }

   if(negval) strncpy(tmstr, "-00:00:00.00",13);
   else strncpy(tmstr, "00:00:00.00", 13);

   for(i=0; i<3; i++)
   {
      if(fabs(tmarr[i]) > 99)
      {
         fprintf(stderr, "Time value longer than 3 digits. %i, %f\n", i, tmarr[i]);
         return -1;
      }
   }

   if(negval) sprintf(tmstr, "-%02.0f:%02.0f:%02.2f", fabs(tmarr[0]), fabs(tmarr[1]), fabs(tmarr[2]));
   else sprintf(tmstr, "%02.0f:%02.0f:%02.2f", fabs(tmarr[0]), fabs(tmarr[1]), fabs(tmarr[2]));

   return 0;
}//int mts_format_time_doub(char tmstr[13], double tmarr[3])


#endif //magtelUtils_hpp
