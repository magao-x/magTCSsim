#ifndef magtcs_cmdstr_hpp
#define magtcs_cmdstr_hpp

#define DATEOBS "dateobs"
#define UT "ut"
#define ST "st"
#define RA "ra"
#define DEC "dec"
#define EPOCH "epoch"
#define HA "ha"
#define AIRMASS "airmass"
#define ZD "zd"
#define TELFOCUS "telfocus"
#define ROTANGLE "rotangle"
#define TELRA "telra"
#define TELDC "teldc"
#define TELEP "telep"
#define TELFC "telfc"
#define VEXSET "vexset"
#define VEYSET "veyset"
#define VEZSET "vezset"
#define VEHSET "vehset"
#define VEVSET "vevset"
#define VEXENC "vexenc"
#define VEYENC "veyenc"
#define VEZENC "vezenc"
#define VEHENC "vehenc"
#define VEVENC "vevenc"
#define VEDATA "vedata"

//#define VANESET "vaneset"
//#define VANEENC "vaneenc"

//-----------------------
#define TELUT "telut"
#define TELST "telst"
#define TELAM "telam"
#define TELPA "telpa"
#define TELHA "telha"
#define TELDM "teldm"
#define DMSTAT "dmstat"
#define TELGUIDE  "telguide"
#define GDRMOUNTMV "gdrmountmv"
#define DATETIME "datetime"
#define TELPOS   "telpos"
#define TELDATA  "teldata"

#define TELAZ "telaz"
#define TELEL "telel"
#define INPHA "inpha"
#define INPRA "inpra"
#define INPDC "inpdc"
#define INPEP "inpep"
#define INPAZ "inpaz"
#define INPEL "inpel"
#define INPPA "inppa"
#define PANGLE "pangle"
#define EANGLE "eangle"
#define GANGLE "gangle"
#define NANGLE "nangle"
#define HANGLE "hangle"
#define FANGLE "fangle"
#define ROTOFH "rotofh"
#define TELROI "telroi"
#define NROTOFF "nrotoff"
#define ROTATORE "rotatore"
#define P1FILT "p1filt"
#define P1MASK "p1mask"
#define GUIDERX1 "guiderx1"
#define GUIDERY1 "guidery1"
#define P2FILT "p2filt"
#define P2MASK "p2mask"
#define GUIDERX2 "guiderx2"
#define GUIDERY2 "guidery2"
#define P3FILT "p3filt"
#define P3MASK "p3mask"
#define GUIDERX3 "guiderx3"
#define GUIDERY3 "guidery3"
#define C1CUR "c1cur"
#define C2CUR "c2cur"
#define C3CUR "c3cur"
#define C1XY2 "c1xy2"
#define C2XY2 "c2xy2"
#define C3XY2 "c3xy2"
#define C1XY3 "c1xy3"
#define C2XY3 "c2xy3"
#define C3XY3 "c3xy3"
#define C1XY4 "c1xy4"
#define C2XY4 "c2xy4"
#define C3XY4 "c3xy4"
#define C1BOX "c1box"
#define C2BOX "c2box"
#define C3BOX "c3box"
#define CA1 "ca1"
#define CA2 "ca2"
#define CA3 "ca3"
#define WXTEMP "wxtemp"
#define WXPRES "wxpres"
#define WXHUMID "wxhumid"
#define WXWIND "wxwind"
#define WXDIR "wxwdir"
#define TELENV "telenv"

#define OFRA "ofra"
#define OFDC "ofdc"
#define OFEP "ofep"
#define OFFP "offp"
#define OFFM "offm"
#define DR   "dr"
#define AEG "aeg"
#define ZSET "zset"
#define ZSTR "zstr"
#define XSET "xset"
#define XSTR "xstr"
#define YSET "yset"
#define YSTR "ystr"
#define SLEW "slew"
#define HALT "h"
#define RO   "ro"
#define NRO  "nro"

//These are unimplemented, but should be available from 802 EDS message
//#define TELTRACK  "teltrack"

//#define ROTFOLLOW "rotfollow"

//This currently returns only whether rotator is following or not.  Should return more about the rotator.
#define ROTMODE "rotmode"

//Current target name is available in one form or another, but it might not be from TCS.  This is a placeholder
#define CATRA     "catra"
#define CATDC     "catdc"
#define CATEP     "catep"
#define CATRO     "catro"
#define CATRM     "catrm"
#define CATOBJ    "catobj"
#define CATDATA   "catdata"

/** Enumerated index of the command strings.
  * Order does not matter, should be one per command
  */
enum cmdidx { DATEOBS_N,
              UT_N,
              ST_N,
              RA_N,
              DEC_N,
              EPOCH_N,
              HA_N,
              AIRMASS_N,
              ZD_N,
              TELFOCUS_N,
              ROTANGLE_N,
              TELRA_N,
              TELDC_N,
              TELEP_N,
              TELFC_N,
              VEXSET_N,
              VEYSET_N,
              VEZSET_N,
              VEHSET_N,
              VEVSET_N,
              VEXENC_N,
              VEYENC_N,
              VEZENC_N,
              VEHENC_N,
              VEVENC_N,
              VEDATA_N,
              TELUT_N,
              TELST_N,
              TELAM_N,
              TELPA_N,
              TELHA_N,
              TELDM_N,
              DMSTAT_N,
              TELGUIDE_N,
              GDRMOUNTMV_N,
              DATETIME_N,
              TELPOS_N,
              TELDATA_N,
              TELAZ_N,
              TELEL_N,
              INPHA_N,
              INPRA_N,
              INPDC_N,
              INPEP_N,
              INPAZ_N,
              INPEL_N,
              INPPA_N,
              PANGLE_N,
              EANGLE_N,
              GANGLE_N,
              NANGLE_N,
              HANGLE_N,
              FANGLE_N,
              ROTOFH_N,
              NROTOFF_N,
              TELROI_N,
              ROTATORE_N,
              P1FILT_N,
              P1MASK_N,
              GUIDERX1_N,
              GUIDERY1_N,
              P2FILT_N,
              P2MASK_N,
              GUIDERX2_N,
              GUIDERY2_N,
              P3FILT_N,
              P3MASK_N,
              GUIDERX3_N,
              GUIDERY3_N,
              C1CUR_N,
              C2CUR_N,
              C3CUR_N,
              C1XY2_N,
              C2XY2_N,
              C3XY2_N,
              C1XY3_N,
              C2XY3_N,
              C3XY3_N,
              C1XY4_N,
              C2XY4_N,
              C3XY4_N,
              C1BOX_N,
              C2BOX_N,
              C3BOX_N,
              CA1_N,
              CA2_N,
              CA3_N,
              WXTEMP_N,
              WXPRES_N,
              WXHUMID_N,
              WXWIND_N,
              WXDIR_N,
              TELENV_N,
              OFRA_N,
              OFDC_N,
              OFEP_N,
              OFFP_N,
              OFFM_N,
              DR_N,
              AEG_N,
              ZSET_N,
              ZSTR_N,
              XSET_N,
              XSTR_N,
              YSET_N,
              YSTR_N,
              SLEW_N,
              HALT_N,
              RO_N,
              NRO_N,
              //TELTRACK_N,
              //ROTFOLLOW_N,
              ROTMODE_N,
              CATRA_N,
              CATDC_N,
              CATEP_N,
              CATRO_N,
              CATRM_N,
              CATOBJ_N,
              CATDATA_N,
              cmdidx_MAX};


#endif //magtcs_cmdstr_hpp

