/** \file magtelescope.hpp
  * \author Jared R. Males
  * \brief Declaration of the magtelescope class.
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

#ifndef magtelescope_hpp
#define magtelescope_hpp

#include "magtelUtils.hpp"


#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#include <stdlib.h>

#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
#include <list>

#include <map>

#include <boost/math/constants/constants.hpp>
#define PI (boost::math::constants::pi<double>())

#include <mx/timeUtils.hpp>
#include <mx/astro/astroDynamics.hpp>


#include "cmdstr.hpp"
#include "tellog.hpp"



#define dms2dd(d, m, s) (d + (m/60.0) + (s/3600.0))

#define NASW 0
#define NASE 1
#define CASS 2
#define AUX1 3
#define AUX2 4
#define AUX3 5



class magtelescope
{
   
protected:
   double m_lon;      ///< The telescope longitude [degrees].
   double m_lat;      ///< The telescope latitude [degrees].
   int m_dateobs[3];  ///< The broken down universal date [YYYY, MM, DD].
   double m_ut {0};   ///< The current universal time [hrs].
   double m_ra {0};   ///< The current right ascension, hours
   double m_dec {0};  ///< The current declination, degrees
   
   double next_ra;  ///< The next value of RA, which will be slewed to.
   double next_dec; ///< The next value of Dec, which will be slewed to.
   
   double off_ra;   ///< The offset in RA
   double off_dec;  ///< The offset in Dec
   
   double az;       ///< azimuth, north referenced, degrees
   double el;       ///< elevation, degrees
   double next_az;  ///< The next value of azimuth, which will be slewed to.
   double next_el;  ///< The next value of elevatoin, which will be slewed to.
   
   double epoch;    ///< coordinate epoch
   double airmass;  ///< the current airmass

   int telfocus;    //= 001725
   int telx;
   int tely;
   double telh;
   double telv;
   
   //int secx;
   //int secy;

   int roi ;//the rotator of interest
   double rotangle; // =  318.4265
   double next_rotangle;
   
   double rotenc; //The rotator encoder angle
   
   

   //Simulation flags
   int simulating;
   int tracking;
   int rottracking;

   int slewing;///<1 if currently slewing
   int was_tracking;///<to reset tracking after slewing

   int guiding;///<1 if guiding - just a flag, doesn't do anything here.
   int was_guiding; ///<to reset to guiding after slewing
   
   int rotmode; ///<1 if following
   int rotwas_following; ///<to reset to following after slewing

   //Simulation time hacks
   double simstart;
   double trackstart;
   double slewstart;
   double rotstart;

   pthread_mutex_t comMutex;
   
   pthread_mutex_t logmutex;
   
   std::string logstr;
   char logval[50];
   
   //Simulation motion
   double az_accel;
   double az_rate;
   double curr_az_rate;
   double az_lastt;
   double el_accel;
   double el_rate;
   double curr_el_rate;
   double el_lastt;

   double rotslewing;
   double rotrate;
   double init_rotangle;
   double rotmovesz;
   int rotmovedir;

   char tgtname[25];
   double cat_ra;
   double cat_dc;
   double cat_ep;
   double cat_ro;
   int cat_rm;
   
   char lmsg[lmsg_size];
   
   
   /// Map of commands to command number, for fast lookup.
   std::map<std::string, cmdidx> commands;
   
public:
   
   /// Default c'tor
   magtelescope();

protected:
   
   /// Load the commands map.
   void loadCommands();

public:
      
   /// Set the typical values of parameters.
   void set_typval();

   /// Get the telescope latitude.
   /** Units: degrees.
     * 
     * \returns the current value of m_lat.
     */
   double lat();
   
   /// Get the telescope longitude
   /** Units: degrees.
     * 
     * \returns the current value of m_lon.
     */
   double lon();
   
   /// Set dateobs to the current date.
   /** Gets current UT date from the O/S.
     * 
     * \returns 0 on success.
     * \returns -1 on error.
     */
   int dateobsNow();
   
   /// Get the current value of dateobs.
   /** Call dateobsNow, then fills in the array with YYYY,MM,DD.
     */
   void dateobs(int dobs[3] /**< [out] The current value of dateobs [YYY,MM,DD].*/);
   
   /// Get the current UT time.
   /** Sets m_ut to the current UT time and returns it.
     *
     * \returns the current value of m_ut
     */
   double ut();
   
   /// Get the current sidereal time.
   /** Calculates LMST and returns it.
     *
     * \returns current value of LMST 
     */
   double st();

   /// Get the current value of Right Ascension.
   /** If nocalc==false and the telescope is not tracking, then RA is re-calculated for the current time, az, and el.
     * 
     * \returns the current value of m_ra.
     */
   double ra(bool nocalc = false /**< [in] [optional] flag to control whether the RA is recalculated.*/);
   
   /// Set the current value of Right Ascension.
   /**
     * \returns 0 on success.
     * \returns -1 on error.
     */ 
   int ra(double ra /**< [in] the new value of Right Ascension*/);
   
   /// Get the current value of Declination.
   /** If nocalc==false and the telescope is not tracking, then Dec is re-calculated for the current time, az, and el.
     * 
     * \returns the current value of m_dec.
     */
   double dec(bool nocalc = false /**< [in] [optional] flag to control whether the Dec is recalculated.*/);
   
   /// Set the current value of Declination.
   /**
     * \returns 0 on success.
     * \returns -1 on error.
     */
   int dec(double d /**< [in] the new value of Declination*/);

      
   void set_next_ra(double r);
   double get_next_ra(){return next_ra;}
   void set_off_ra(double r);
   double get_off_ra();
   
   

   
   void set_next_dec(double d);
   double get_next_dec(){return next_dec;}
   void set_off_dec(double d);   double get_off_dec();
	
  

   void set_az(double a);
   double get_az(bool nocalc = false);

   void set_el(double e);
   double get_el(bool nocalc = false);

   int set_ra_dec(double nra, double ndec);
   int set_az_el(double naz, double nel);
   double get_ha(bool nocalc = false);
   double get_next_ha();

   double get_epoch();    // = 2000.00
   double get_airmass();  // =  1.000
   double get_zd();       // = 0.040398
   double get_pa();
   double get_next_pa();


   double get_wxtemp();
   double get_wxpres();
   double get_wxhumid();
   double get_wxwind();
   double get_wxdir();
   
   
   int get_telfocus(){return telfocus;}    //= 001725
   int get_telz(){return get_telfocus();}
   void set_telfocus(int tf);
   
   int get_telx(){return telx;}
   void set_telx(int tx);
   
   int get_tely(){return tely;}
   void set_tely(int ty);
   
   double get_telh(){return telh;}
   void set_telh(double th);
   
   double get_telv(){return telv;}
   void set_telv(double tv);
   
   /*int get_secx(){return secx;}
   void set_secx(int sx);
   
   int get_secy(){return secy;}
   void set_secy(int sy);*/

   int get_roi(){return roi;}
   
   double get_rotangle(); //returns the current rotator angle, relative to N
   double get_rotenc(); //returns the current rotator encoder angle.
   
   double get_next_rotangle(){return next_rotangle;}
   int set_next_rotangle(double nr);
   
   
   
   //Calculates current azimuth and elevation given sidereal time, RA, and Dec.
   //Reference: J. Meeus "Astronomical Algorithms", 1991.
   int calc_az_el(bool nocalc = false);
   int calc_next_az_el();
   double get_next_az(){return next_az;}
   double get_next_el(){return next_el;}

   //Calculates current RA and DEC, given Az, El, and sidereal time.
   //Reference: J. Meeus "Astronomical Algorithms", 1991.
   int calc_ra_dec(bool nocalc = false);

   //Simulation control.
   int start_simulating();
   int stop_simulating();
   int get_simulating() {return simulating;}
   int reset_simulation();
   
   //Tracking in RA and DEC
   int start_tracking();
   int stop_tracking();
   int get_tracking() {return tracking;}

   //Guiding
   int start_guiding();
   int stop_guiding();
   int get_guiding() {return guiding;}

   //Rotator following
   int start_following();
   int stop_following();
   int get_rotmode() {return rotmode;}
   
   //Slewing
   int start_slew();
   int stop_slew();
   int offset_slew(int dir);
   int slew_az_el();
   int slew_az();
   int slew_el();
   int get_slewing() {return (slewing);}

   //Rotator moves
   int start_rotmove();
   int move_rot();
   int stop_rotmove();
   int get_rotslewing(){ return rotslewing;}
   
   
   //Interface
   int process_string(char * inpstr, FILE *fp, bool statonly);
   
   int process_command(char *response, std::string inpstr);
   
   int get_status(char *status, char *cmd_str, int statlen);
   
   int get_status(char *status, int cmd_N, int statlen);
   
   int do_command(int cmd_N, std::vector<std::string> args);
   
   std::list<std::string> lastcom;
   std::list<double> lastcom_time;
   std::list<std::string> lastresp;

   int set_tgtname(std::string tn);
   int set_cat_ra(double cra);
   int set_cat_dc(double cdc);

   int set_cat_ep(double cep);
   int set_cat_ro(double cro);
   int set_cat_rm(int crm);
};

inline
magtelescope::magtelescope()
{
   loadCommands();   //st = 0.0;
   
   for(int i=0;i<3;i++) m_dateobs[i] = 0;

   next_ra = 0.0;
   next_dec = 0.0;
   
   az = 0.0;
   el = 0.0;

   epoch = 0;
   airmass = 0;
   
   telfocus = 0;
   telx = 0;
   tely = 0;
   telh = 0;
   telv = 0;

   roi = 0;
   rotangle = 0;
   next_rotangle = 0;

   //These values set so calcs match typvals of sidereal time and zd.
   m_lat = -29.015; //degrees
   m_lon = -70.69166667; // degrees

   //Simulation flags
   simulating = 0;
   tracking = 0;
   guiding = 0;
   slewing = 0;
   rotmode = 0;
   rotslewing = 0;

   //Simulation time hacks
   simstart = 0;
   trackstart = 0;
   slewstart = 0;
   rotstart = 0;

   az_accel = 0.1;//*PI/180.;
   az_rate = 2.;//*PI/180.;
   curr_az_rate = 0.;
   el_accel = 0.1;//*PI/180.;
   el_rate = 1.0;//*PI/180.;
   curr_el_rate = 0.;

   rotrate = 6; //degrees per second.

   char logname[200];

   tgtname[0] = '\0';
   
   time_t tl;
   tl = time(0);
   strftime(logname, 200, "logs/magtel_%m%d%Y%H%M%S.log", localtime(&tl));
   tellog(0, 0, 0, "magtelescope", logname);
   TELLOG("Initialized.");

   for(int i=0;i<10;i++)
   {
      lastcom.push_back("");
      lastcom_time.push_back(-1);
      lastresp.push_back("");
   }
   pthread_mutex_init(&logmutex, 0);

   pthread_mutex_init(&comMutex, 0);
   
   srand(time(0));
}

inline
void magtelescope::loadCommands()
{
   commands.clear();
   commands[DATEOBS]= DATEOBS_N;
   commands[UT]= UT_N;
   commands[ST]= ST_N;
   commands[RA]= RA_N;
   commands[DEC]= DEC_N;
   commands[EPOCH]= EPOCH_N;
   commands[HA]= HA_N;
   commands[AIRMASS]= AIRMASS_N;
   commands[ZD]= ZD_N;
   commands[TELFOCUS]= TELFOCUS_N;
   commands[ROTANGLE]= ROTANGLE_N;
   commands[TELRA]= TELRA_N;
   commands[TELDC]= TELDC_N;
   commands[TELEP]= TELEP_N;
   commands[TELFC]= TELFC_N;
   commands[VEXSET]= VEXSET_N;
   commands[VEYSET]= VEYSET_N;
   commands[VEZSET]= VEZSET_N;
   commands[VEHSET]= VEHSET_N;
   commands[VEVSET]= VEVSET_N;
   commands[VEXENC]= VEXENC_N;
   commands[VEYENC]= VEYENC_N;
   commands[VEZENC]= VEZENC_N;
   commands[VEHENC]= VEHENC_N;
   commands[VEVENC]= VEVENC_N;
   commands[VEDATA]= VEDATA_N;
   commands[TELUT]= TELUT_N;
   commands[TELST]= TELST_N;
   commands[TELAM]= TELAM_N;
   commands[TELPA]= TELPA_N;
   commands[TELHA]= TELHA_N;
   commands[TELDM]= TELDM_N;
   commands[DMSTAT]= DMSTAT_N;
   commands[TELGUIDE]= TELGUIDE_N;
   commands[GDRMOUNTMV]= GDRMOUNTMV_N;
   commands[DATETIME]= DATETIME_N;
   commands[TELPOS]= TELPOS_N;
   commands[TELDATA]= TELDATA_N;
   commands[TELAZ]= TELAZ_N;
   commands[TELEL]= TELEL_N;
   commands[INPHA]= INPHA_N;
   commands[INPRA]= INPRA_N;
   commands[INPDC]= INPDC_N;
   commands[INPEP]= INPEP_N;
   commands[INPAZ]= INPAZ_N;
   commands[INPEL]= INPEL_N;
   commands[INPPA]= INPPA_N;
   commands[PANGLE]= PANGLE_N;
   commands[EANGLE]= EANGLE_N;
   commands[GANGLE]= GANGLE_N;
   commands[NANGLE]= NANGLE_N;
   commands[HANGLE]= HANGLE_N;
   commands[FANGLE]= FANGLE_N;
   commands[ROTOFH]= ROTOFH_N;
   commands[NROTOFF]= NROTOFF_N;
   commands[TELROI]= TELROI_N;
   commands[ROTATORE]= ROTATORE_N;
   commands[P1FILT]= P1FILT_N;
   commands[P1MASK]= P1MASK_N;
   commands[GUIDERX1]= GUIDERX1_N;
   commands[GUIDERY1]= GUIDERY1_N;
   commands[P2FILT]= P2FILT_N;
   commands[P2MASK]= P2MASK_N;
   commands[GUIDERX2]= GUIDERX2_N;
   commands[GUIDERY2]= GUIDERY2_N;
   commands[P3FILT]= P3FILT_N;
   commands[P3MASK]= P3MASK_N;
   commands[GUIDERX3]= GUIDERX3_N;
   commands[GUIDERY3]= GUIDERY3_N;
   commands[C1CUR]= C1CUR_N;
   commands[C2CUR]= C2CUR_N;
   commands[C3CUR]= C3CUR_N;
   commands[C1XY2]= C1XY2_N;
   commands[C2XY2]= C2XY2_N;
   commands[C3XY2]= C3XY2_N;
   commands[C1XY3]= C1XY3_N;
   commands[C2XY3]= C2XY3_N;
   commands[C3XY3]= C3XY3_N;
   commands[C1XY4]= C1XY4_N;
   commands[C2XY4]= C2XY4_N;
   commands[C3XY4]= C3XY4_N;
   commands[C1BOX]= C1BOX_N;
   commands[C2BOX]= C2BOX_N;
   commands[C3BOX]= C3BOX_N;
   commands[CA1]= CA1_N;
   commands[CA2]= CA2_N;
   commands[CA3]= CA3_N;
   commands[WXTEMP]= WXTEMP_N;
   commands[WXPRES]= WXPRES_N;
   commands[WXHUMID]= WXHUMID_N;
   commands[WXWIND]= WXWIND_N;
   commands[WXDIR]= WXDIR_N;
   commands[TELENV]= TELENV_N;
   commands[OFRA]= OFRA_N;
   commands[OFDC]= OFDC_N;
   commands[OFEP]= OFEP_N;
   commands[OFFP]= OFFP_N;
   commands[OFFM]= OFFM_N;
   commands[DR]= DR_N;
   commands[AEG]= AEG_N;
   commands[ZSET]= ZSET_N;
   commands[ZSTR]= ZSTR_N;
   commands[XSET]= XSET_N;
   commands[XSTR]= XSTR_N;
   commands[YSET]= YSET_N;
   commands[YSTR]= YSTR_N;
   commands[SLEW]= SLEW_N;
   commands[HALT]= HALT_N;
   commands[RO]= RO_N;
   commands[NRO]= NRO_N;
   commands[ROTMODE]= ROTMODE_N;
   commands[CATRA]= CATRA_N;
   commands[CATDC]= CATDC_N;
   commands[CATEP]= CATEP_N;
   commands[CATRO]= CATRO_N;
   commands[CATRM]= CATRM_N;
   commands[CATOBJ]= CATOBJ_N;
   commands[CATDATA]= CATDATA_N;
}

inline
void magtelescope::set_typval()
{
   dateobsNow();
   
   m_ra = st(); 
   m_dec = lat(); 
   
   next_ra = m_ra;
   next_dec = m_dec;
   
   mx::astro::calcAzEl(az, el, st()-m_ra, m_dec, m_lat);
   epoch = 2000.0;
   airmass = 1.0; 
   telfocus = 1725; 
   telx = 10;
   tely = -10;
   telh = 1.24;
   telv = -2.34;

   rotangle = 0;
   rotenc = 0;

   TELLOG("Set typical values.");
}

inline
double magtelescope::lat()
{
   return m_lat;
}

inline
double magtelescope::lon()
{
   return m_lon;
}

inline
int magtelescope::dateobsNow()
{
   time_t rawtime;
   time(&rawtime);
   struct tm * uttm = gmtime (&rawtime);

   int dobs[3];
   
   m_dateobs[0] = 1900+uttm->tm_year;
   m_dateobs[1] = uttm->tm_mon+1;
   m_dateobs[2] = 1 + uttm->tm_mday;

   return 0;
}

inline
void magtelescope::dateobs(int dobs[3])
{
   dateobsNow();
   
   for(int i=0;i<3;i++) dobs[i] = m_dateobs[i];
}

inline
double magtelescope::ut()
{
   time_t rawtime;
   time(&rawtime);
   struct tm * uttm = gmtime (&rawtime);
   
   m_ut = uttm->tm_hour + uttm->tm_min/60.0 + uttm->tm_sec/3600.0;

   return m_ut;
}

inline
double magtelescope::st()
{
   double LMST;
   int dobs[3];
   dateobs(dobs);
   
   mx::astro::getLMST(LMST, dobs[0],dobs[1], ((double) dobs[2]) +(ut()/24.0), lon());

   return LMST;
}

inline
double magtelescope::ra(bool nocalc)
{
   if(!tracking && !nocalc) calc_ra_dec();
   return m_ra;
}

inline
int magtelescope::ra(double r)
{
   m_ra = r;
   calc_az_el(true);
	
	snprintf(lmsg, lmsg_size, "Set ra: %f", r);
   TELLOG(lmsg);
   
   return 0;
}

inline
int magtelescope::dec(double d)
{
   m_dec = d;
   calc_az_el(true);
	
	snprintf(lmsg, lmsg_size, "Set dec: %f", d);
   TELLOG(lmsg);
   
   return 0;
}

inline
double magtelescope::dec(bool nocalc)
{
   if(!tracking && !nocalc) calc_ra_dec();
   return m_dec;
}

inline
void magtelescope::set_next_ra(double r)
{
   next_ra = r;
	calc_next_az_el();
   
   snprintf(lmsg, lmsg_size, "Set next ra: %f", r);
   TELLOG(lmsg);
}

inline
void magtelescope::set_off_ra(double r)
{
   off_ra = r;
   
   snprintf(lmsg, lmsg_size, "Set ra offset: %f", r);
   TELLOG(lmsg);
}

inline
double magtelescope::get_off_ra()
{
	return off_ra;
}
	
inline
void magtelescope::set_next_dec(double d)
{
   next_dec = d;
   calc_next_az_el();
   snprintf(lmsg, lmsg_size, "Set next dec: %f", d);
   TELLOG(lmsg);
}

inline
void magtelescope::set_off_dec(double d)
{
   off_dec = d;
   
   snprintf(lmsg, lmsg_size, "Set dec offset: %f", d);
   TELLOG(lmsg);
}

inline
double magtelescope::get_off_dec()
{
   return off_dec;
}

inline
void magtelescope::set_az(double a)
{
   az = a;
   calc_ra_dec(true);

   snprintf(lmsg, lmsg_size, "Set az: %f", a);
   TELLOG(lmsg);
}

inline
double magtelescope::get_az(bool nocalc)
{
   //std::cout << "get_az slewing = " << slewing << "\n";
   if(tracking && !nocalc && !slewing) 
   {
      calc_az_el();
   }
   else if(slewing && !nocalc)
   {
      slew_az_el();
   }
   return az;
}

inline
void magtelescope::set_el(double e)
{
   el = e;
   calc_ra_dec(true);

   snprintf(lmsg, lmsg_size, "Set el: %f", e);
   TELLOG(lmsg);
}

inline
double magtelescope::get_el(bool nocalc)
{
   if(tracking && !nocalc && !slewing) 
   {
      calc_az_el();
   }
   else if(slewing && !nocalc)
   {
      slew_az_el();
   }
   return el;
}

inline
int magtelescope::set_ra_dec(double nra, double ndec)
{
   m_ra = nra;
   m_dec = ndec;
   return calc_az_el(true);
}

inline
double magtelescope::get_ha(bool nocalc)
{
   double ha;
   ha = fmod(st()-ra(nocalc), 24.0);
   if(ha < 0.0) ha += 24.0;
   
   return ha;
}

inline
double magtelescope::get_next_ha()
{
   double nha;
   nha = fmod(st()-get_next_ra(), 24.0);
   if(nha < 0.0) nha += 24.0;
   return nha;
}

inline
double magtelescope::get_epoch()
{
   return epoch;
}

inline
double magtelescope::get_airmass()
{
   return 1./cos(get_zd()*3.14159/180.);
}

inline
double magtelescope::get_zd()
{
   return 90.0 - get_el();
}

inline
double magtelescope::get_pa()
{
   return mx::astro::parAngDeg(m_lat, dec(), -get_ha()/24.*360.);
   
   //return sin(get_az()*3.14159/180.)*cos(m_lat*3.14149/180.)/cos(dec()*3.14149/180.)*180./3.14159;
}

inline
double magtelescope::get_next_pa()
{
   return mx::astro::parAngDeg(m_lat, get_next_dec(), -get_next_ha()/24.*360.);
   //return sin(get_next_az()*3.14159/180.)*cos(m_lat*3.14149/180.)/cos(get_next_dec()*3.14149/180.)*180./3.14159;
}

inline
double magtelescope::get_wxtemp()
{
   return 18. + (.25 - ((double)rand())/((double)RAND_MAX)*.5);
}

inline
double magtelescope::get_wxpres()
{
	return 760. + (5. - ((double)rand())/((double)RAND_MAX)*10.);;
}

inline
double magtelescope::get_wxhumid()
{
	return 31. + (1. - ((double)rand())/((double)RAND_MAX)*2.);;
}

inline
double magtelescope::get_wxwind()
{
	return 10. + + (5. - ((double)rand())/((double)RAND_MAX)*10.);;
}

inline
double magtelescope::get_wxdir()
{
	return 90.+ (35. - ((double)rand())/((double)RAND_MAX)*70.);;
}
 
inline
void magtelescope::set_telfocus(int tf)
{
   telfocus = tf;
   
   snprintf(lmsg, lmsg_size, "Set telfocus: %i", tf);
   TELLOG(lmsg);
}

inline
void magtelescope::set_telx(int tx)
{
   telx = tx;
   
   snprintf(lmsg, lmsg_size, "Set telx: %i", tx);
   TELLOG(lmsg);
}

inline
void magtelescope::set_tely(int ty)
{
   tely = ty;
   
   snprintf(lmsg, lmsg_size, "Set tely: %i", ty);
   TELLOG(lmsg);
}

inline
void magtelescope::set_telh(double th)
{
   telh = th;
   
   snprintf(lmsg, lmsg_size, "Set telh: %f", th);
   TELLOG(lmsg);
}

inline
void magtelescope::set_telv(double tv)
{
   telv = tv;
   
   snprintf(lmsg, lmsg_size, "Set telv: %f", tv);
   TELLOG(lmsg);
}

 
/*void magtelescope::set_secx(int sx)
{
   secx = sx;
	snprintf(lmsg, lmsg_size, "Set secondary x: %i", sx);
   TELLOG(lmsg);
}*/

/*void magtelescope::set_secy(int sy)
{
   secy = sy;
   snprintf(lmsg, lmsg_size, "Set secondary y: %i", sy);
   TELLOG(lmsg);
}*/

inline
double magtelescope::get_rotangle()
{
   if(rotslewing) move_rot();
   
   if(rotmode || rotslewing)
   {
      return rotangle;
   }
   else
   {
      rotangle = rotenc - get_pa();
      return rotangle;
   }

}

inline
double magtelescope::get_rotenc()
{
   if(rotslewing) move_rot();
   
   if(rotmode || rotslewing)
   {
      rotenc = rotangle + get_pa(); 
   }

   return rotenc;
   
}

inline
int magtelescope::set_next_rotangle(double nr)
{
   next_rotangle = nr;
   return 0;
}

inline
int magtelescope::calc_az_el(bool nocalc)
{
   double ha_rad, dec_rad, lat_rad, az_rad, el_rad;

   ha_rad = get_ha(nocalc)*15.0*PI/180.0;
   dec_rad = dec(nocalc)*PI/180.0;
   lat_rad = lat()*PI/180.0;

   az_rad = atan2(sin(ha_rad), cos(ha_rad)*sin(lat_rad)-tan(dec_rad)*cos(lat_rad));

   el_rad = asin(sin(lat_rad)*sin(dec_rad) + cos(lat_rad)*cos(dec_rad)*cos(ha_rad));

   az = az_rad*180.0/PI + 180.0;
   el = el_rad*180.0/PI;

   return 0;
}

inline
int magtelescope::calc_next_az_el()
{
   double ha_rad, dec_rad, lat_rad, az_rad, el_rad;

   ha_rad = get_next_ha()*15.0*PI/180.0;
   dec_rad = next_dec*PI/180.0;
   lat_rad = lat()*PI/180.0;

   az_rad = atan2(sin(ha_rad), cos(ha_rad)*sin(lat_rad)-tan(dec_rad)*cos(lat_rad));

   el_rad = asin(sin(lat_rad)*sin(dec_rad) + cos(lat_rad)*cos(dec_rad)*cos(ha_rad));

   next_az = az_rad*180.0/PI+ 180.0;
   next_el = el_rad*180.0/PI;

   return 0;
}

inline
int magtelescope::calc_ra_dec(bool nocalc)
{
   double az_rad, lat_rad, el_rad, ha_rad, dec_rad;

   az_rad = (get_az(nocalc)-180.0)*PI/180.0;
   lat_rad = lat()*PI/180.0;
   el_rad = get_el(nocalc)*PI/180.0;

   ha_rad = atan2(sin(az_rad), cos(az_rad)*sin(lat_rad)+tan(el_rad)*cos(lat_rad));
   dec_rad = asin(sin(lat_rad)*sin(el_rad) - cos(lat_rad)*cos(el_rad)*cos(az_rad));

   m_ra = (st() - ha_rad*180.0/(15.0*PI));

   if(m_ra < 0) m_ra = m_ra + 2.*PI;

   m_dec = dec_rad*180.0/PI;

   return 0;
}

inline
int magtelescope::start_simulating()
{
   timeval stime;
   simulating = 1;

   gettimeofday(&stime, 0);

   simstart = 0;
   /*if(simstart == 0) 
   {
      simstart = stime.tv_sec+(stime.tv_usec/1000000.0);
   }*/

   if(tracking == 1)
   {
      if(trackstart == 0)
      {
         trackstart = stime.tv_sec+(stime.tv_usec/1000000.0);
      }
   }

   return 0;
}

inline
int magtelescope::stop_simulating()
{
   simulating = 0;
   return 0;
}

inline
int magtelescope::reset_simulation()
{
   simstart = 0;
   trackstart = 0;

   if(simulating) start_simulating();

   return 0;
}

inline
int magtelescope::start_tracking()
{
   timeval stime;
   tracking = 1;
   was_tracking=1;
   if(trackstart == 0) 
   {
      gettimeofday(&stime, 0);
      trackstart = stime.tv_sec+(stime.tv_usec/1000000.0);
   }
   return 0;
}

inline
int magtelescope::stop_tracking()
{
   tracking = 0;
   was_tracking = 0;
   return 0;
}

inline
int magtelescope::start_guiding()
{
   guiding = 1;
   was_guiding=1;
   return 0;
}

inline
int magtelescope::stop_guiding()
{
   guiding = 0;
   was_guiding = 0;
   return 0;
}

inline
int magtelescope::start_following()
{
   rotmode = 1;
   rotwas_following=1;
   return 0;
}

inline
int magtelescope::stop_following()
{
   rotmode = 0;
   rotwas_following=0;
   return 0;
}

inline
int magtelescope::start_slew()
{
   timeval stime;

   gettimeofday(&stime, 0);
   slewstart = stime.tv_sec+(stime.tv_usec/1000000.0);

   was_tracking = tracking;
   tracking = 0;
   slewing = 1;

   az_lastt = 0.;
   el_lastt = 0.;

   //std::cout << "starting slew\n";
   return 0;
}

inline
int magtelescope::stop_slew()
{
   tracking = was_tracking;
   slewing = 0;
   //std::cout << "slew stopped\n";
   return 0;
}

inline
int magtelescope::offset_slew(int dir)
{
   if(dir >= 0) dir = 1;
   else dir = -1;
   
   next_ra = ra(false) + dir*(off_ra/3600.)/15.;
   next_dec = dec(false) + dir*(off_dec/3600.);
   return start_slew();
}

inline
int magtelescope::slew_az_el()
{

   slew_az();
   slew_el();

   calc_ra_dec(1);

   //std::cout << next_ra << " " << ra << "\n";
   if(fabs(next_ra - m_ra) < 1e-5 && fabs(next_dec - m_dec) < 1e-5) 
   {
      slewing = 0;
      tracking = was_tracking;
   }
   //std::cout << "slewing=" << slewing << "\n";
   return 0;
}

inline
int magtelescope::slew_az()
{
   double az_rem, currt, dt, dtacc;
   timeval stime;
   static double init_az_rem;
   
   gettimeofday(&stime, 0);
   currt = stime.tv_sec+(stime.tv_usec/1000000.0);

   
   dtacc = az_rate/az_accel;
   
   calc_next_az_el();
   az_rem = fabs(next_az - az);

   
   if(az_lastt == 0)
   {
      az_lastt = currt;
      init_az_rem = az_rem;
   }
   
   dt = currt-az_lastt;
   
   //std::cout << "slew_az " << next_az << " " << az << " " << az_rem << "\n";
   //std::cout << "dtacc=" << dtacc << " dt=" << dt << "\n";
   
   curr_az_rate = fabs(curr_az_rate);
   
   if (az_rem < .5*az_accel*dtacc*dtacc && az_rem < 0.5*init_az_rem)
   {
      curr_az_rate = curr_az_rate - dt*az_accel;
      if(curr_az_rate < 0) curr_az_rate = 0.;
   }
   else if(curr_az_rate < az_rate)
   {
      curr_az_rate = curr_az_rate + dt*az_accel;
      if(curr_az_rate > az_rate) curr_az_rate = az_rate;
   }
   if((next_az-az)*curr_az_rate < 0) curr_az_rate*=-1;
   
   //std::cout << "curr_az_rate=" << curr_az_rate << " az_rate=" << az_rate << " ";
   az = az + dt*curr_az_rate;
   az_rem = fabs(next_az-az);
   //std::cout << "az_rem=" << az_rem << " dazacc=" << .5*az_accel*dtacc*dtacc << "\n";
   
   if((az_rem < .5*az_accel*dtacc*dtacc && az_rem < 0.5*init_az_rem && fabs(curr_az_rate) < 0.01) || init_az_rem < .002)
   {
      az = next_az;
   }
   
   az_lastt = currt;
   
   return 0;
}

inline
int magtelescope::slew_el()
{
   double el_rem, currt, dt, dtel;
   static double init_el_rem;
   timeval stime;

   gettimeofday(&stime, 0);
   currt = stime.tv_sec+(stime.tv_usec/1000000.0);

   
   dtel = el_rate/el_accel;
   
   calc_next_az_el();
   el_rem = fabs(next_el - el);

   if(el_lastt == 0)
   {
      el_lastt = currt;
      init_el_rem = el_rem;
   }
   
   dt = currt-el_lastt;
   curr_el_rate = fabs(curr_el_rate);
   if (el_rem < .5*el_accel*dtel*dtel && el_rem < 0.5*init_el_rem)
   {
      curr_el_rate = curr_el_rate - dt*el_accel;
      if(curr_el_rate < 0) curr_el_rate = 0.;
   }
   else if(curr_el_rate < el_rate)
   {
      curr_el_rate = curr_el_rate + dt*el_accel;
      if(curr_el_rate > el_rate) curr_el_rate = el_rate;
   }
   if((next_el-el)*curr_el_rate < 0) curr_el_rate*=-1;
      
   el = el + dt*curr_el_rate;
   el_rem = fabs(next_el-el);
   
   if((el_rem < .5*el_accel*dtel*dtel && el_rem < 0.5*init_el_rem && fabs(curr_el_rate) < 0.01) || init_el_rem < .002)
   {
      el = next_el;
   }
   //std::cout << "curr_el_rate=" << curr_el_rate << " el_rate=" << el_rate << " ";
   //std::cout << "el_rem=" << el_rem << " delacc=" << .5*el_accel*dtel*dtel << "\n";
   
   el_lastt = currt;
   
   return 0;
}

inline
int magtelescope::start_rotmove()
{
   timeval stime;

   gettimeofday(&stime, 0);
   rotstart = stime.tv_sec+(stime.tv_usec/1000000.0);

   rotwas_following = rotmode;
   
   rotmode = 0;
   rotslewing = 1;
   
   init_rotangle = rotangle;
   
   if(next_rotangle > rotangle)
   {
      if(next_rotangle - rotangle > 180.)
      {
         rotmovedir = -1;
         rotmovesz = 360.-(next_rotangle - rotangle);
      }
      else
      {
         rotmovedir = 1;
         rotmovesz = next_rotangle -rotangle;
      }
   }
   else
   {
      if(rotangle-next_rotangle > 180.)
      {
         rotmovedir = 1;
         rotmovesz = 360. - (rotangle - next_rotangle);
      }
      else
      {
         rotmovedir = -1;
         rotmovesz = rotangle-next_rotangle;
      }
   }

   //std::cout << "starting rotator move\n";
   return 0;
}

inline
int magtelescope::stop_rotmove()
{
   move_rot(); //move it one bit more.
   
   
   rotmode = rotwas_following;
   
   rotslewing = 0;
   
   
   //std::cout << "rotator move stopped\n";
   return 0;
}

inline
int magtelescope::move_rot()
{
   timeval stime;
   
   gettimeofday(&stime, 0);
   double currt = stime.tv_sec+(stime.tv_usec/1000000.0);

   double dt = currt - rotstart;

   double drotang = dt * rotrate;

   if(drotang >= rotmovesz)
   {
      rotangle = next_rotangle;
      rotmode = rotwas_following;
      rotslewing = 0;
   }
   else
   {
      rotangle = init_rotangle + rotmovedir*drotang;
      if(rotangle < 0) rotangle = rotangle + 360;
      if(rotangle >= 360) rotangle = rotangle - 360;
   }
   return 0;
}

inline
int magtelescope::process_string(char * inpstr, FILE *fp, bool statonly)
{
   char response[100];
   
   pthread_mutex_lock(&comMutex);

   
   if(inpstr[strlen(inpstr)-1]==13) inpstr[strlen(inpstr)-1] = '\0'; //from telnet

   if(strchr(inpstr, ' ') > 0)
   {
      printf ("Command received: %s.\n", inpstr);

      if(strncmp(inpstr, "tcssim", 6) == 0)
      {
         if(strlen(inpstr) < 8) return -1;

         if(strncmp(inpstr+7, "quit", 4) == 0)
         {
            exit(0);
         }
      }

      if(statonly)
      {
         std::cout << "This is the status only port.\n";
         
         pthread_mutex_unlock(&comMutex);
         
         return -1;
      }
      process_command(response, inpstr);
      
      tcs_fputs(response, fp);
   }
   else
   {
      printf ("Status request received: %s.\n", inpstr);

      get_status(response, inpstr, 100);

      tcs_fputs(response, fp);

      sprintf(response, "%s\n", response);
   }

   lastcom.push_back(inpstr);
   lastcom_time.push_back(mx::get_curr_time());
   lastcom.pop_front();
   lastcom_time.pop_front();
   lastresp.push_back(response);
   lastresp.pop_front();

   pthread_mutex_unlock(&comMutex);
   
   return 0;
}

inline
int magtelescope::process_command(char * response, std::string inpstr)
{
   std::string tmparg, comstr;
   std::vector<std::string> args;
   int lsp, nsp, cmd_N, cmdres;

   lsp = 0;
   nsp = inpstr.find(' ', lsp);
   comstr.assign(inpstr, lsp, nsp-lsp);
   while(inpstr[nsp+1] == ' ') nsp++;//drop extra spaces.
   lsp = nsp+1;

   while(lsp > -1)
   {
      nsp = inpstr.find(' ', lsp);

      if(nsp < 0) nsp = inpstr.length();//Must be the last arg.
      if(nsp-lsp > 0) //Only add an arg if it is really there
      {
         tmparg.assign(inpstr, lsp, nsp-lsp);
         args.push_back(tmparg);
      }
      if(nsp != (int)inpstr.length() && nsp > -1) //There's more to parse
      {
         while(inpstr[nsp+1] == ' ') nsp++;
         if(nsp > (int) inpstr.length()-1) nsp = -1;
         lsp = nsp + 1;
      }
      else lsp = -1; //we're done.
   }

   //execute command
   try
   {
      cmd_N = commands.at(comstr);
   }
   catch(...)
   {
      cmd_N = -1;
   }
   
   cmdres = do_command(cmd_N, args);

   sprintf(response, "%d\n", cmdres);

   //for(i=0;i<(int)args.size();i++) std::cout << args[i] << "\n";
   return 0;
}

inline
int magtelescope::get_status(char *status, char *cmd_str, int statlen)
{
   int cmd_N;

   try
   {
      cmd_N = commands.at(cmd_str);
   }
   catch(...)
   {
      cmd_N = -1;
   }
   
   return get_status(status, cmd_N, statlen);
}

inline
int magtelescope::get_status(char *status, int cmd_N, int statlen)
{
   int tmpdate[3];
   double tmparr[3];;
   
   if(cmd_N < 0)
   {
      strncpy(status, "-1", statlen);
      return -1;
   }
   
   //Adjust for synonyms
   if(cmd_N == TELUT_N) cmd_N = UT_N;
   if(cmd_N == TELST_N) cmd_N = ST_N;
   if(cmd_N == TELRA_N) cmd_N = RA_N;
   if(cmd_N == TELDC_N) cmd_N = DEC_N;
   if(cmd_N == TELEP_N) cmd_N = EPOCH_N;
   if(cmd_N == TELHA_N) cmd_N = HA_N;
   if(cmd_N == TELAM_N) cmd_N = AIRMASS_N;
   if(cmd_N == TELFC_N) cmd_N = TELFOCUS_N;
   if(cmd_N == INPEP_N) cmd_N = EPOCH_N;
   
   switch(cmd_N)
   {
      case DATEOBS_N:
      {
         if(statlen < 11)
         {
            fprintf(stderr, "Status string not long enough for m_dateobs.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         dateobs(tmpdate);
         if(mts_format_date(status, tmpdate) != 0)
         {
            fprintf(stderr, "mts_format_date returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }
      case UT_N:
      {
         if(statlen < 9)
         {
            fprintf(stderr, "Status string not long enough for ut/telut.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, ut());
         if(mts_format_time_uint(status, tmparr) != 0)
         {
            fprintf(stderr, "UT_N: mts_format_time_uint returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//UT_N
      case ST_N:
      {
         if(statlen < 9)
         {
            fprintf(stderr, "Status string not long enough for st/telst.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, st());
         if(mts_format_time_uint(status, tmparr) != 0)
         {
            fprintf(stderr, "ST_N: mts_format_time_uint returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//case ST_N:
      case RA_N:
      {
         if(statlen < 13)
         {
            fprintf(stderr, "Status string not long enough for ra/telra.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, ra());
         if(mts_format_time_doub(status, tmparr) != 0)
         {
            fprintf(stderr, "RA_N mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//case RA_N:
      case INPRA_N:
      {
         if(statlen < 13)
         {
            fprintf(stderr, "Status string not long enough for inpra.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, get_next_ra());
         if(mts_format_time_doub(status, tmparr) != 0)
         {
            fprintf(stderr, "INPRA_N mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//case INPRA_N:
      case DEC_N:
      {
         if(statlen < 13)
         {
            fprintf(stderr, "Status string not long enough for dec/teldc.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, dec());
         if(mts_format_time_doub(status, tmparr) != 0)
         {
            fprintf(stderr, "mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//case DEC_N:
      case INPDC_N:
      {
         if(statlen < 13)
         {
            fprintf(stderr, "Status string not long enough for inpdc.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, get_next_dec());
         if(mts_format_time_doub(status, tmparr) != 0)
         {
            fprintf(stderr, "mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//case INPDC_N:
      case EPOCH_N:
      {
         if(statlen < 7)
         {
            fprintf(stderr, "Status string not long enough for epoch/telep.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.2f", get_epoch());
         return 0; //success
      }//case EPOCH_N:
      case HA_N:
      {
         if(statlen < 13)
         {
            fprintf(stderr, "Status string not long enough for ha/telha.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, get_ha());
         if(mts_format_time_uint(status, tmparr) != 0)
         {
            fprintf(stderr, "mts_format_time_uint returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//case HA_N:
      case INPHA_N:
      {
         if(statlen < 13)
         {
            fprintf(stderr, "Status string not long enough for INPHA.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, get_next_ha());
         if(mts_format_time_uint(status, tmparr) != 0)
         {
            fprintf(stderr, "mts_format_time_uint returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//case INPHA_N:
      case AIRMASS_N:
      {
         if(statlen < 7)
         {
            fprintf(stderr, "Status string not long enough for airmass.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.3f", get_airmass());
         return 0; //success
      }//case AIRMASS_N:
      case ZD_N:
      {
         if(statlen < 9)
         {
            fprintf(stderr, "Status string not long enough for zd.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.6f", get_zd());
         return 0; //success
      }//case ZD_N:
      case TELPA_N:
      {
         if(statlen < 8)
         {
            fprintf(stderr, "Status string not long enough for telpa.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.4f", get_pa());
         return 0; //success
      }//case TELPA_N:
      case INPPA_N:
      {
         if(statlen < 8)
         {
            fprintf(stderr, "Status string not long enough for inppa.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.4f", get_next_pa());
         return 0; //success
      }//case INPPA_N:
      case TELFOCUS_N:
      {
         if(statlen < 7)
         {
            fprintf(stderr, "Status string not long enough for telfocus.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%06i", get_telfocus());
         return 0; //success
      }//case TELFOCUS_N:

      case VEDATA_N:
      {
         if(statlen < 36)
         {
            fprintf(stderr, "Status string not long enough for vaneset.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "  %06i %06i %06i %06i %06i %06i %.3f %.3f %.3f %.3f", get_telfocus(), get_telfocus(), get_telx(), get_telx(), get_tely(), get_tely(), get_telh(), get_telh(), get_telv(), get_telv());
         std::cout << status << "\n";
         return 0;
      }//case VEDATA_N

      case TELROI_N:
      {
         if(statlen < 2)
         {
            fprintf(stderr, "Status string not long enough for telroi.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%1d", get_roi());
         return 0; //success
      }//case TELROI_N:

      case ROTANGLE_N:
      {
         if(statlen < 9)
         {
            fprintf(stderr, "Status string not long enough for rotangle.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.4f", get_rotangle());
         return 0; //success
      }//case ROTANGLE_N:
      case NROTOFF_N:
      {
         if(statlen < 9)
         {
            fprintf(stderr, "Status string not long enough for nrotoff.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%8.4f", get_next_rotangle());
         return 0; //success
      }//case NROTOFF_N:
      case ROTATORE_N:
      {
         if(statlen < 9)
         {
            fprintf(stderr, "Status string not long enough for rotatore.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%8.4f", get_rotenc());
         return 0; //success
      }//case ROTATORE_N:
      case TELAZ_N:
      {
         if(statlen < 8)
         {
            fprintf(stderr, "Status string not long enough for telaz.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.4f", get_az());
         return 0; //success
      }//case TELAZ_N:
      case INPAZ_N:
      {
         if(statlen < 8)
         {
            fprintf(stderr, "Status string not long enough for inpaz.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.4f", get_next_az());
         return 0; //success
      }//case INPAZ_N:
      case TELEL_N:
      {
         if(statlen < 8)
         {
            fprintf(stderr, "Status string not long enough for telel.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.4f", get_el());
         return 0; //success
      }//case TELEL_N:
      case INPEL_N:
      {
         if(statlen < 8)
         {
            fprintf(stderr, "Status string not long enough for inpel.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status,statlen, "%.4f", get_next_el());
         return 0; //success
      }//case INPEL_N:
      case TELDM_N:
      {
         if(statlen < 4)
         {
            fprintf(stderr, "Status string not long enough for teldm.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%i", (int)get_az());
         return 0;
      }
      case DMSTAT_N:
      {
         if(statlen < 3)
         {
            fprintf(stderr, "Status string no long enough for dmstat.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%i", 1);
         return 0;
      }
      case WXTEMP_N:
      {
         if(statlen < 10)
         {
            fprintf(stderr, "Status string no long enough for wxtemp.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%0.2f", get_wxtemp());
         return 0;
      }
      case WXPRES_N:
      {
         if(statlen < 10)
         {
            fprintf(stderr, "Status string no long enough for wxpres.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%0.2f", get_wxpres());
         return 0;
      }
      case WXHUMID_N:
      {
         if(statlen < 10)
         {
            fprintf(stderr, "Status string no long enough for wxhumid.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%0.2f", get_wxhumid());
         return 0;
      }
      case WXWIND_N:
      {
         if(statlen < 10)
         {
            fprintf(stderr, "Status string no long enough for wxwind.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%0.2f", get_wxwind());
         return 0;
      }
      case WXDIR_N:
      {
         if(statlen < 10)
         {
            fprintf(stderr, "Status string not long enough for wxdir.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%0.2f", get_wxdir());
         return 0;
      }
      case TELENV_N:
      {
         if(statlen < 50)
         {
            fprintf(stderr, "Status string not long enough for wxdir.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         
         snprintf(status, statlen, "%0.2f %0.2f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f", get_wxtemp(), get_wxpres(), get_wxhumid(), get_wxwind(), get_wxdir(),get_wxtemp()*1.1, get_wxtemp()*.98, get_wxtemp()*.9, get_wxtemp()*1.15, get_wxtemp()-20);
         
         return 0;
      }//case TELENV_N:
      case TELGUIDE_N:
      {
         if(statlen < 3)
         {
            fprintf(stderr, "Status string not long enough for telguide.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%i%i", get_tracking(), get_guiding());
         return 0;
      }//case TELGUIDE_N
      case GDRMOUNTMV_N:
      {
         if(statlen < 4)
         {
            fprintf(stderr, "Status string not long enough for gdrmountmv.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%i%i%i", get_slewing(), 0,0);
         return 0;
      }//case GDRMOUNTMV_N
      case DATETIME_N:
      {
         char dttmp[32], uttmp[32], sttmp[32];
         if(statlen < 35)
         {
            fprintf(stderr, "Status string not long enough for gdrmountmv.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }

         dateobs(tmpdate);
         if(mts_format_date(dttmp, tmpdate) != 0)
         {
            fprintf(stderr, "mts_format_date returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }

         deg2dms(tmparr, ut());
         if(mts_format_time_uint(uttmp, tmparr) != 0)
         {
            fprintf(stderr, "UT_N: mts_format_time_uint returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }

         deg2dms(tmparr, st());
         if(mts_format_time_uint(sttmp, tmparr) != 0)
         {
            fprintf(stderr, "ST_N: mts_format_time_uint returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         
         snprintf(status, statlen, "%s %s %s", dttmp, uttmp,sttmp);
         return 0;
      }//case DATETIME_N
      case TELPOS_N:
      {
         char rastr[32], dcstr[32], hastr[32];
         if(statlen < 50)
         {
            fprintf(stderr, "Status string not long enough for telpos.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }

         deg2dms(tmparr, ra());
         if(mts_format_time_doub(rastr, tmparr) != 0)
         {
            fprintf(stderr, "TELPOS_N RA mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         


         deg2dms(tmparr, dec());
         if(mts_format_time_doub(dcstr, tmparr) != 0)
         {
            fprintf(stderr, "mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }

         deg2dms(tmparr, get_ha());
         if(mts_format_time_uint(hastr, tmparr) != 0)
         {
            fprintf(stderr, "mts_format_time_uint returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         
         snprintf(status, statlen, "%s %s %.2f %s %.3f %.4f", rastr, dcstr, get_epoch(), hastr, get_airmass(), get_rotangle());
         
         return 0;
      }//case TELPOS_N

      case TELDATA_N:
      {
         if(statlen < 50)
         {
            fprintf(stderr, "Status string not long enough for teldata.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }


         snprintf(status,statlen, "%1d %i%i %i00 %.4f %.4f %.6f %.4f %i %i", get_roi(), get_tracking(), get_guiding(), get_slewing(),get_az(), get_el(), get_zd(), get_pa(), (int)get_az(), 1);

         return 0;
         
      }//case TELDATA_N
      
      //case ROTFOLLOW_N:
      case ROTMODE_N:
      {
         if(statlen < 2)
         {
            fprintf(stderr, "Status string not long enough for rotmode.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%i", get_rotmode());
         return 0;
      }//case ROTMODE_N

      case CATRA_N:
      {
         if(statlen < 13)
         {
            fprintf(stderr, "Status string not long enough for catra.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, cat_ra);
         if(mts_format_time_doub(status, tmparr) != 0)
         {
            fprintf(stderr, "CATRA_N mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//case CATRA_N:
      case CATDC_N:
      {
         if(statlen < 13)
         {
            fprintf(stderr, "Status string not long enough for catdc.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         deg2dms(tmparr, cat_dc);
         if(mts_format_time_doub(status, tmparr) != 0)
         {
            fprintf(stderr, "mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }
         return 0; //success
      }//case CATDC_N:
      case CATEP_N:
      {
         if(statlen < 8)
         {
            fprintf(stderr, "Status string not long enough for catep.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%.2f", cat_ep);
         return 0; //success
      }//case CATEP_N:
      case CATRO_N:
      {
         if(statlen < 10)
         {
            fprintf(stderr, "Status string not long enough for catra.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         snprintf(status, statlen, "%.4f", cat_ro);
         return 0; //success
      }//case CATRO_N:

      case CATRM_N:
      {
         std::string crm = "OFF";
         if(statlen < 4)
         {
            fprintf(stderr, "Status string not long enough for catrm.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }

         if(cat_rm == 1) crm = "EQU";
         if(cat_rm == 2) crm = "GRV";
         if(cat_rm == 3) crm = "HRZ";
         
         snprintf(status, statlen, "%s", crm.c_str());
         
         return 0; //success
      }//case CATRM_N:
      case CATOBJ_N:
      {
         if(statlen < 25)
         {
            fprintf(stderr, "Status string not long enough for catobj.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }
         std::string tname = tgtname;
         tname += "(sim)";
         snprintf(status, statlen, "%s", tname.c_str());
         return 0;
      }//case CATOBJ_N
      case CATDATA_N:
      {
         char rastr[32], dcstr[32];
         if(statlen < 50)
         {
            fprintf(stderr, "Status string not long enough for catdata.\n");
            strncpy(status, "-1", statlen);
            return -2;
         }

         deg2dms(tmparr, cat_ra);
         if(mts_format_time_doub(rastr, tmparr) != 0)
         {
            fprintf(stderr, "CATRA_N mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }

         deg2dms(tmparr, cat_dc);
         if(mts_format_time_doub(dcstr, tmparr) != 0)
         {
            fprintf(stderr, "mts_format_time_doub returned error.\n");
            strncpy(status, "-1", statlen);
            return -3;
         }

         std::string crm = "OFF";
         if(cat_rm == 1) crm = "EQU";
         if(cat_rm == 2) crm = "GRV";
         if(cat_rm == 3) crm = "HRZ";
         
         std::string tname = tgtname;
         tname += "(sim)";

         
         snprintf(status, statlen, "%s %s %.2f %.4f %s %s", rastr, dcstr, cat_ep, cat_ro, crm.c_str(), tname.c_str());
         return 0;
      }//case CATDATA_N
      default:
      {
         fprintf(stderr, "Status request valid, but not implemented in magTCSsim yet.\n");
         strncpy(status, "0", statlen);
         return 1;
      }
   }
   
   
   return -10;
}//int get_status(char *status, char *cmd_str, int statlen)

inline
int magtelescope::do_command(int cmd_N, std::vector<std::string> args)
{
   if(cmd_N < 0)
   {
      return -1;
   }

   if(cmd_N == TELUT_N) cmd_N = UT_N;
   if(cmd_N == TELST_N) cmd_N = ST_N;
   if(cmd_N == TELRA_N) cmd_N = RA_N;
   if(cmd_N == TELDC_N) cmd_N = DEC_N;
   if(cmd_N == TELEP_N) cmd_N = EPOCH_N;
   if(cmd_N == TELHA_N) cmd_N = HA_N;
   if(cmd_N == TELAM_N) cmd_N = AIRMASS_N;
   if(cmd_N == TELFC_N) cmd_N = TELFOCUS_N;

   switch(cmd_N)
   {
      case RA_N:
      {
         if(args.size() != 3)
         {
            fprintf(stderr, "Incorrect number of arguments for command RA (should be 3).\n");
            return -1;
         }
         logstr = "Got ra ";
         
         snprintf(logval, 50, "%i %i %0.4f ", atoi(args[0].c_str()),atoi(args[1].c_str()), strtod(args[2].c_str(),0));
         logstr += logval;

         TELLOG(logstr.c_str());

         double nra = atof(args[0].c_str()) +atof(args[1].c_str())/60. + atof(args[2].c_str())/3600.;
         nra = nra;//*PI/180.;
         set_next_ra(nra);
         return 0;
      }
      case DEC_N:
      {
         if(args.size() != 3)
         {
            fprintf(stderr, "Incorrect number of arguments for command DEC (should be 3).\n");
            return -1;
         }
         logstr = "Got dec ";
         
         snprintf(logval, 50, "%i %i %0.4f ", atoi(args[0].c_str()),atoi(args[1].c_str()), strtod(args[2].c_str(),0));
         logstr += logval;

         TELLOG(logstr.c_str());

         double ndec = atof(args[0].c_str());
         if(ndec < 0) ndec -= atof(args[1].c_str())/60.;
         else ndec += atof(args[1].c_str())/60.;
         if(ndec < 0) ndec -=atof(args[2].c_str())/3600.;
         else ndec += atof(args[2].c_str())/3600.;

         ndec = ndec;// * PI/180.;
         set_next_dec(ndec);
         return 0;
      }
      case SLEW_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command SLEW (should be 0).\n");
            return -1;
         }
         logstr = "Got slew";
         
         TELLOG(logstr.c_str());

         start_slew();
         return 0;
      }
      case HALT_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command h (halt) (should be 0).\n");
            return -1;
         }
         logstr = "Got halt";
         
         TELLOG(logstr.c_str());

         stop_slew();
         return 0;
      }
      case OFRA_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command OFRA (should be 1).\n");
            return -1;
         }
         logstr = "Got offra ";
         
         snprintf(logval, 20, "%0.4f ", strtod(args[0].c_str(),0));
         logstr += logval;

         TELLOG(logstr.c_str());

         set_off_ra(atof(args[0].c_str()));
         return 0;
      }
      case OFDC_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command OFDC (should be 1).\n");
            return -1;
         }
         logstr = "Got ofdc ";
         
         snprintf(logval, 20, "%0.4f", strtod(args[0].c_str(),0));
         logstr += logval;

         TELLOG(logstr.c_str());

         set_off_dec(atof(args[0].c_str()));
         return 0;
      }
      case OFFP_N:
      {
        
         if(args.size() != 0)
         {
            fprintf(stderr, "Incorrect number of arguments for command OFFP (should be 0).\n");
            return -1;
         }
         logstr = "Got offp";
         
         TELLOG(logstr.c_str());

         offset_slew(1);
         return 0;
      }
      case OFFM_N:
      {
         if(args.size() != 0)
         {
            fprintf(stderr, "Incorrect number of arguments for command OFFM (should be 0).\n");
            return -1;
         }
         logstr = "Got offm";
         
         TELLOG(logstr.c_str());

         offset_slew(-1);
         return 0;
      }
      case AEG_N:
      {
         if(args.size() != 2)
         {
            fprintf(stderr, "Incorrect number of arguments for command AEG (should be 2).\n");
            return -1;
         }
         logstr = "Got aeg ";
         
         snprintf(logval, 20, "%0.2f ", strtod(args[0].c_str(),0));
         logstr += logval;
   
         snprintf(logval, 20, "%0.2f", strtod(args[1].c_str(),0));
         logstr += logval;

         TELLOG(logstr.c_str());

         set_az(get_az(false) + strtod(args[0].c_str(),0)/3600.);
         set_el(get_el(false) + strtod(args[1].c_str(),0)/3600.);
         return 0;
      }
      case ZSET_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command ZSET (should be 1).\n");
            return -1;
         }
         logstr = "Got zset ";
         
         snprintf(logval, 20, "%i ", atoi(args[0].c_str()));
         logstr += logval;
   
         TELLOG(logstr.c_str());


         set_telfocus(atoi(args[0].c_str()));
         return 0;
      }
      case ZSTR_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command ZSTR (should be 1).\n");
            return -1;
         }
         logstr = "Got zstr ";
         
         snprintf(logval, 20, "%i ", atoi(args[0].c_str()));
         logstr += logval;
   
         TELLOG(logstr.c_str());

         set_telfocus(telfocus + atoi(args[0].c_str()));
         return 0;
      }
      case XSET_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command XSET (should be 1).\n");
            return -1;
         }
         logstr = "Got xset ";
         
         snprintf(logval, 20, "%i ", atoi(args[0].c_str()));
         logstr += logval;
   
         TELLOG(logstr.c_str());

         set_telx(atoi(args[0].c_str()));
         return 0;
      }
      case XSTR_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command XSTR (should be 1).\n");
            return -1;
         }
         logstr = "Got xstr ";
         
         snprintf(logval, 20, "%i ", atoi(args[0].c_str()));
         logstr += logval;
   
         TELLOG(logstr.c_str());
         set_telx(telx + atoi(args[0].c_str()));
         return 0;
      }
      case YSET_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command YSET (should be 1).\n");
            return -1;
         }
         logstr = "Got yset ";
         
         snprintf(logval, 20, "%i ", atoi(args[0].c_str()));
         logstr += logval;
   
         TELLOG(logstr.c_str());

         set_tely(atoi(args[0].c_str()));
         return 0;
      }
      case YSTR_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command YSTR (should be 1).\n");
            return -1;
         }
         logstr = "Got ystr ";
         
         snprintf(logval, 20, "%i ", atoi(args[0].c_str()));
         logstr += logval;
   
         TELLOG(logstr.c_str());
         set_tely(tely + atoi(args[0].c_str()));
         return 0;
      }
      case RO_N:
      {
         if(args.size() != 1)
         {
            fprintf(stderr, "Incorrect number of arguments for command ro (should be 1).\n");
            return -1;
         }
         logstr = "Got ro ";
         
         snprintf(logval, 20, "%06.4f", strtod(args[0].c_str(),0));
         logstr += logval;
   
         TELLOG(logstr.c_str());

         set_next_rotangle(strtod(args[0].c_str(),0));
         start_rotmove();
         return 0;
      }
      default:
      {
         fprintf(stderr, "Command valid, but not implemented in magTCSsim yet.\n");
         return 1;
      }
   }
   
   return -10;
}//int int magtelescope::do_command(int cmd_N, std::vector<std::string> args)


inline
int magtelescope::set_tgtname(std::string tn)
{
   strncpy(tgtname, tn.c_str(), 25);
   return 0;
}

inline
int magtelescope::set_cat_ra(double cra)
{
   cat_ra = cra;
   return 0;
}

inline
int magtelescope::set_cat_dc(double cdc)
{
   cat_dc = cdc;
   return 0;
}

inline
int magtelescope::set_cat_ep(double cep)
{
   cat_ep = cep;
   return 0;
}

inline
int magtelescope::set_cat_ro(double cro)
{
   cat_ro = cro;
   return 0;
}

inline
int magtelescope::set_cat_rm(int crm)
{
   cat_rm = crm;
   return 0;
}
   


#endif //magtelescope_hpp
