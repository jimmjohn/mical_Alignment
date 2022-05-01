/*
shift in xpos_ydev and ypos_xdev only
Error in <TPad::RangeAxis>: illegal axis coordinates range: xmin=-0.500000, ymin=0.000000, xmax=1535.500000, ymax=0.000000



correction of offset is before/earlier
 align_xstr_xdev... Why const term is varying so much ?
 Fit should be within 3 to 31 & 33 to 58 separately

inefficiency_npixel
_jj0 : without any align_ystrp_xval...
_jj1 : with only align_ystrp_xval & align_xstrp_yval
_jj2 : with all four align_ystrp_xval, align_xstrp_yval, align_xstrp_xval & align_ystrp_yval
_jj3 : Again with all four but correction from _jj2 and update corrections in each iterations
_jj4 : With corrections terms from _jj3
_jj5 : Remove many efficiency plot, modified timing part and addition of correlation large multiplicity in different layers, define SHIFT0 + correction from _jj4
_jj6 : int iTimeSlopeConst = 0; no parabilic correction, raw occupancy at last iteration
_jj7 : double dx = yextloc[ij]+0.5-istr; {inittime = tmptime; tshft =ytoffset[ij][ypts[ij][yiy]]; istr = istrytime[ij] = ypts[ij][yiy];} + raw_seloccu etc in isfill loop if ( yexter[ij]<0.4 &&  xexter[ij]<0.4)
_jj8 : Without any time correction ( ?without pixel wise time correction.)
_jj9 : without CAU correction
_jj10 : istr = int(yextloc[ij]+0.5); dx = yextloc[ij]-istr;
_jj11 : without local[ij] +=align_ystr_xdev[ij][1]*otherext[ij];
_jj12 : without pixelwise correction
_jj13 : return back to jj8 with corection from jj8, but with _jj10
_jj14 : if (iiter<nmxiter-1 && ntcor==1) {xoff[iocc] += fitfx->GetParameter(1);}
        if (ntcor==0 && iiter<nmxiter-1) { align_xstr_xdev[iocc][0] =ptrx->Parameter(0);}
_jj15 : if (iiter<nmxiter-1 && ntcor==0) { xoff[iocc] += fitfx->GetParameter(1);}
        if (ntcor==1 && iiter<nmxiter-1) {  align_xstr_xdev[iocc][0] =ptrx->Parameter(0);}
_jj16 : _jj15, but if (iiter%2==0) { .... also without timing
_jj17 : use only align_xstr_ydev & anal_ystr_xdev
_jj18 : use only align_xstr_xdev & anal_ystr_ydev
_jj19 : use only align_xstr_ydev & anal_ystr_xdev, but use function cal_slope
_jj20 : all four with cal_slope
_jj21 : all four both corrections are in if (ntcor==1 && iiter<nmxiter-1)
_jj22 : all four both corrections are with four paramters and if (ntcor==1 && iiter<nmxiter-1)
_jj23 : local[ij] +=off[ij];   in GetXposInStrip, xpos_xdev etc are put outside istiming
shift in
_jj24 : no update on xoff/yoff (isCombinedOffset=false), but if (ntcor==1) align_xstr_xdev[iocc][0] +=ptrx->Parameter(0); etc
_jj25 : if (ntcor==0) align_xstr_xdev[iocc][0] +=ptrx->Parameter(0); etc
_jj26 : _jj25, but no update align_xstr_xdev[iocc][0] +=ptrx->Parameter(0); (fourr of them)
_jj27 : _jj25, but no update on align_xstr_ydev[iocc][0] +=ptrx->Parameter(0); (four of them)
_jj28 : Only three paramters in cal_slope double par3 = (par[2]-par[0]-par[1]*nstrip/4.)/(nstrip/4.); with xoff/yoff in ntcor=0 and align_ update in ntcor==1
_jj29 : _jj28 but isCombinedOffset==false & if (ntcor==0) align_xstr_xdev[iocc][0] +=ptrx->Parameter(0);
_jj30 : _jj28 but isCombinedOffset==false & if (ntcor==1) align_xstr_xdev[iocc][0] +=ptrx->Parameter(0);
_jj31 : same as _jj28, with removing dead strips
_jj32 : _jj31 with if (iiter%2==0) align_xstr_xdev[iocc][0] +=ptrx->Parameter(0);
_jj33 : _jj31 with if (iiter%2==1) align_xstr_xdev[iocc][0] +=ptrx->Parameter(0);

_jj34 : Again same as _jj31, with new align_paramters
_jj35 : without align_correction

_jj36 : With align_, but without xoff/yoff but with parabolic correction const int nmxusedtimehit=3; const int nmxusedhits=3; if (ntcor==0) {align_xstr_xdev[iocc][0] +=ptrx->Parameter(0);}
_jj37 : const int nmxusedtimehit=2; const int nmxusedhits=3;
_jj38 : const int nmxusedtimehit=2; const int nmxusedhits=2;

\\upto _jj38, align_xstr/align_ystr etc were filled within ntocr==1 option
_jj39 : Actually _jj36 with align_xstr_xdev and align_ystr_ydev correction bool whichFill = (ntcor==0) ? true : false; if (ntcor==0)
_jj40 : Actually _jj36 with align_xstr_xdev and align_ystr_ydev correction bool whichFill = (ntcor==1) ? true : false; if (ntcor==1)
_jj41 : Without ntcor=0 and and only iteration on 8th layer
_jj42 : With 5 iteration on Aug data
_jj43 : Aug 25 data    if (abs(xextloc[occulyr] - xpts[occulyr][ix])<1) { expp = xextloc[occulyr];//xpts[occulyr][ix];
_jj45 : L8 X S51 52 and 53 removed
_jj46 : double expdiff = 10000;
		  for (int ix=0; ix<ypts[occulyr].size(); ix++) {
		    if (abs(yextloc[occulyr] - ypts[occulyr][ix])<1.0) {
		      double tmpexpdiff = abs(yextloc[occulyr] - ypts[occulyr][ix]);
		      if(tmpexpdiff < expdiff) {
		      expdiff = tmpexpdiff;	expp = ypts[occulyr][ix]; break;
		      }
		    }
		  }

_jj47 : if (ypts[ij][yiy]==int(yextloc[ij]+0.5)) {
        timesy[ij] = 300.0 + 0.1*(int(event->vytdc_l[ij][iyyy[yiy]]->at(0)-event->tdc_ref_l[ij])) - (5.2*stripdelay[int(istr/8)]) +0.1*(cau_calib1[ij]+cau_calib2[ij])/4.;

_jj48 : Looked for extrapolated strips, otherwise its left or right one.
_jj49 : Without strpos_vs_time
_jj50 : Without strpos_vs_time and with only left/rigth strips, not looking for third one


_jj51 : With parabolic correction from _jj50
_jj52 : With old parabolic function and minimum timing
_jj53 : With parabolic function from _jj50 and minimum timing
_jj54 : accprng =0.5;//1.0;  for efficiency calculation
//madurai data : sim01@/home/mdridata/DaqData64
miical_jj01 : without parabolic correction. To calculate multiplicity dependent time correction
miical_jj02 : parabolic correction incorporated from jj01
miical_jj03 : correction are calculated for timing with jj02 parabolic corrections

miical_jj04 : 6000+tdc_l-tdcref_l
miical_jj05: 11500+tdc_l-tdcref_t
miical_jj06 : 6000+tdc_l_hit1-tdc_ref_l
miical_jj07 : 6000+tdc_l_hit0-tdc_ref_l , new histogram added  h_xtstrpmult_xtdchitmult[ij]
miical_jj07 : efficiency sets for 5 minutes interval

miical_jj13 : xstr_xdev.. ystr_ydev.. found not updated. These parameters are updated
miical_jj14 : without parabolic corrections for time_vs_pos

miical_jj16 : region wise study for L0


_miical_jj20_: time correction without parabolic correction
_miical_jj21_: time correction with parabolic correction from _miical_jj20_
_miical_jj22_: timeing is done for tdc-ref close to 850 ns others are not taken
_miical_jj23_: tot corrections incorporated
_miical_jj24_: tot corrections incorporated isLast = 1
_miical_noband_jj27_: no tot corrections isLast =1; The each pixel the delta t is filled to see 4 different bands
miical_noband_jj28_ : no tot corrections isLast =1; The each pixel the delta t is filled to see 4 different bands, the mean of delta t distribution is used for correction in different band except zeroth band
miical_noband_jj29_ : no tot corrections isLast =1; The each pixel the delta t is filled to see 4 different bands, the mean of delta t distribution is used for correction in different band including zeroth band
miical_noband_jj30_ :  xpos_xtdev_str filled for 4 different band regions.
miical_noband_jj31_ :  xpos_xtdev_str filled for 4 different band regions. No TOT correctiona added parabolically
miical_noband_jj32_ :  xpos_xtdev_str filled for 4 different band regions. No TOT correctiona added parabolically , No pixel wise correction
miical_noband_jj33_ : Correction of timing based on mean of delta distribution from different region using average time shift in Strip_Vs_TOT
miical_noband_jj34_ :  Pixel wise (strip vs TOT) tot corrections are included.
miical_noband_jj35_ :  Pixel wise (strip vs TOT) tot corrections are included for different multiplicity
miical_noband_jj36_ :  Pixel wise corr from functional fit fitfunctot() for whole layer
miical_noband_jj37_ :  Pixel wise corr from functional fit fitfunctot() for individual strip
miical_noband_jj38_ :  Pixel wise corr from functional fit fitfunctot2() for whole layer   parabolic 4th order
miical_noband_jj39_ :  jj38 present, xpts, ypts, chisq, ndf, rawxtime, rawytime are stored in root tuple.
miical_noband_jj40_ :  jj38 present, xpts, ypts, chisq, ndf, rawxtime, rawytime are stored in root tuple, root tuple filled without array
miical_noband_jj41_ :  position input is given from rawxtime[i]-rawxtime[0] for multiplicity 2,3,4
miical_noband_jj42_ :  position input is given from pulsewith[i]/pulsewidth[0] for multiplicity 2,3,4
miical_noband_jj43_ :  all multiplicity all timing of L8 stored in variable array.
miical_noband_jj44_ :  The xyposdev return back to 7 strips from 2 strips.
miical_noband_jj45_ : The data is run for L7 and L8 to calculated time correction for different zone of the efficiency. L7 and L8 is occulyr in timing fit
miical_noband_jj46_ : The data is run for L7 and L8 to calculated time correction for different zone of the efficiency. L7 and L8 is occulyr in timing fit only for 9.8 KV all layers are in iteration, No TOT correction
miical_noband_jj47_ : The all files of HV are run for one iteration. TOT correction included.
miical_noband_jj48_ : The all files of HV are run for one iteration. No TOT correction included.
miical_noband_jj49_ : The timing is taken strips which is having larger pulse width. No TOT correction included.
miical_noband_jj50_ : The timing is taken strips which is having larger pulse width. TOT correction included.
miical_noband_jj53_: laylast =9;
miical_noband_jj54_: laylast =9; and maxtimechisq=4.
miical_noband_jj55_: laylast =11; and maxtimechisq=2. , the datasets
miical_noband_jj56_: The pixel-wise TOT corrections are incorporated. The L7 and L8 voltage scan data analysed.
miical_noband_jj57_: pos_vs_time correction only given range M2 (-3,3),M3 (-2.5,2.5), M4 (-1.5,1.5)
miical_noband_jj58_: pos_vs_time correction only given range M2 (-5,5),M3 (-2.5,2.5), M4 (-1.5,1.5)
miical_noband_jj59_: pos_vs_tot correction only given range M2, M3, M4 (0.8,1.2)
miical_noband_jj60 : with pre-existing parabolic correction
miical_noband_jj61 : without any parabolic correction
miical_noband_jj62 : with parabolic correction from _jj61
miical_noband_jj63 : without parabolic correction ... All previous runs (jj61 and 60) timing is wrongly estimated. Consider jj63 for further analysis
miical_noband_jj64 : occulyr == 9: parabolic corrections are used other than occulyr==9
g++ -g -pthread -m64 -Wno-deprecated -I${ROOTSYS}/include -I${CLHEP_BASE_DIR}/include -o anal_rpcdata   EveTree.C StraightLineFit.cc anal_rpcdata.cc  `root-config --cflags --libs` -lMinuit -L${CLHEP_BASE_DIR}/lib -lCLHEP

g++ -g -pthread -m64 -Wno-deprecated -I${ROOTSYS}/include -o anal_rpcdata   EveTree.C StraightLineFit.cc anal_rpcdata.cc  `root-config --cflags --libs` -lMinuit


1xx. Error in extrapolation : Why so large distribution ? Also sigma_true .ne. sqrt(sigma_obs^2 - mu_ext^2)
Ans Used dist[], where we have valid timing

2. In MC : time_xreso_l11_i0->Draw("") & time_yreso_l11_i0->Draw("") is not stored properly, by
         time_xyreso_l11_i0->Draw("") does not have any problem
Ans : Problem with ubuntu@home

//xx 3. Why there is shift in xtdev in MC sample, whereas xdev is fine, initial xval[ij] was not proper 13/07/2016 :
//Ans: same as I. Used dist[], where we have valid timing
//
//4. h_xrawcorhits->Fill(-1.0, -1.0); h_xrawcorhits->Fill(ix, iy); etc are in NOT correct place.
//
////xx5. timeoffsetx[ij] & timeoffsety[ij] updated only iiter=0 and first iteration only (temporarily for all iteration in MC)
//Retrun back original one
//
//6. shift_pos, why only upto nmxiter, instead of 2*nmxiter
//
////xx7. xtime_exterr[occulyr][iiterrs]->Fill(xtexter[occulyr]); values in histogramme and in rottuple are not same
//Ans : plot and drawing were from different iteration
//
////xx8. change condition for deltatcovy[ij][jk] += dt1*dt2; (only condition on Y-side fit)
//: Done
//xx9. No bias in deltatcov etc ??
//
//Different histogrammes
//a. deltatcov[ij][jk]  -> h_deltatcov
//b. timex_correl[ij][jk][iiter+1] (wider now for ntcor=0)
//c. xtime_exterr_l11_i0->Draw(), but xtime_exterr_l11_i3->Draw() is fine ntcor=1
//   same to
//
//   Ans : //due to wrong condition
//                  int iTimeSlopeConst = ((isTimeCorOrReso) && ntcor==1) ? 0 : -1; //-1 : 0;
//
//
//
_ff : is used for time and poistion corrections
_gg : After all correction
_hh : corecting dist[ij] for time
_ii : without constraining time slope by 1/c in StraightLineFit.cc
_jj : bool isTimeCorOrReso=false; timeserrx2[iocc] =  (time_xreso[iocc][iiterrs]->GetRMS())*time_xreso[iocc][iiterrs]->GetRMS() - bias_intime_xreso2[iocc];
_kk bool isTimeCorOrReso=false; timeserrx2[iocc] =  time_xrms[iocc]*time_xrms[iocc] - bias_intime_xreso2[iocc];

kkxx : Again use constraint : bool isTimeCorOrReso=true;
kkyy : Again use constraint : bool isTimeCorOrReso=false;
kkzz : int iTimeSlopeConst = 0; also change in Straightlinefit

15th Oct : strp_xmul[ij][iiterrs]->Fill(xposinstr[ij], xhits[ij]); Moved in different place, need to return back its original position

initxtime[12][2] : Introduced for 64 strip on 27th Nov 2015
rawxtime : Measured time
rawxtime0 : +offset corection for a layer
rawxtime1 : +pathlenght correction
rawxtime2 : +Stripwise correction  //ignore
timesx    : +Pixel wise correction
xtime     : +parabolic correction
rawxtime3 : -pathlenght correction

./configure
gmake
gmake install


rm daqeventdict.cpp
rm daqeventdict.h
rootcint daqeventdict.cpp -c DaqEvent.h

g++ -g -std=gnu++11 -pthread -m64 -Wno-deprecated -I${ROOTSYS}/include -I${CLHEP_BASE_DIR}/include -o anal_rpcdata_pethu daqeventdict.cpp StraightLineFit.cc anal_rpcdata_pethu.cc  `root-config --cflags --libs` -lMinuit -L${CLHEP_BASE_DIR}/lib -lCLHEP



g++ -g -pthread -m64 -Wno-deprecated -I${ROOTSYS}/include -I${CLHEP_BASE_DIR}/include -o madurai_rpcdata madurai_rpcdata.cc  `root-config --cflags --libs` -lMinuit -L${CLHEP_BASE_DIR}/lib -lCLHEP

 StraightLineFit ytimefit(0, dist, ytime, timeserry2, yusedtime, occulyr, 4, ytcorstr, ytcorend, float(7.0));
if (jk==4) timesy[jk] -=25.0;
int ifirst=1;
 fitxy[il] = new TF1(name, fitspec, nstrip*il, nstrip*(il+0.5), 4);
 xy_occu[il]->GetXaxis()->SetRangeUser(nstrip*il, nstrip*(il+0.5));
 for (int istr=nstrip*il; istr<nstrip*(il+0.5); istr++)
*/

//#define C217STRIP       // Strip convention of C217 for time shift
//#define ONEMBYONEM      // RPC of 1m x 1m with readout in 32 strips
//#define OLDFORMAT       //old C217 formated data
//#define ONETIMEFORMAT   //Number of TDC per layer (need manual intervention to map strip and TDC)
//#define MONTECARLO      //MC simulated event
//#define FASTSIM     //Simple simulation using aperture_mc.C
//#define FULLG4SIM    //Simulation using Geant4 (C217Sim)
#define ISEFFICIENCY
#define TIMESLOPE
//#define MADURAIINTM1 //Madurai intermediate stage, where middle layers are occupied by 50%
//#define MADURAIINTM2  //Madurai intermediate stage, where 4th layer with NINO chip and 5,6,7 & 8
//      were partially filled from file INORUN_20160421_185949 (thoughOn 15/04/2016, around 12 hrs
//      SG00013 was moved to stacker to mount the NINO boards )
//#define MADURAIINTM3 //With 2mx2m trigger INORUN_20160804_153224
#define MADURAIINTM4 //With 2mx2m trigger and L1,L4,L5,L6,L7,L8,L9 RPC Daq switch only when combined IDaq and RPCdaq (RPC_INORUN_20170425_185300.ire)
#define NOISE_SIM // Noise file for simulation
//#define NEWRPCDAQ // timing data changed in new RPCDaq data in madurai (8 channel per X and Y side)
//#define CAUDATA
#define NEWRPCDAQ1 //all layers are instumented using NINO-RPCDaq and layer 11 instrumented by Anush-RPCDAQ //Data from RPC_evtraw-r231738.rre
//#define CAUDATA1
//#define SHIFT0 //for shifting 4th ninocable in layer-3 X-side,
//#define BARC_EVTBLR
//#define BARCROOT //one time stamp is there   //jim jim 2021 data
#define TIFRROOT // layer wise timestamp is present //jim jim 2018 data



#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <new>
#include<climits>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <TTimeStamp.h>
#include "TSystem.h"
#ifdef OLDFORMAT
#include "DaqEvent.h"
#endif
#ifdef NEWRPCDAQ
#include "RPCEve_old.h"
#endif
#ifdef CAUDATA
#include "CauTree_old.h"
#endif

#ifdef NEWRPCDAQ1
#ifdef TIFRROOT
#include "RPCEve.h"
#endif
#ifdef BARCROOT
#include "BARCEve.h"
#endif
#endif


#ifdef CAUDATA1
#include "CauTree.h"
#endif
#ifdef BARC_EVTBLR
#include "evtdata.h"
#include "tdcinfo.h"
#include "Evt_LinkDef.h"
#endif
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TDatime.h"

#include "TStyle.h"
#include "TAttFill.h"
#include "TPaveStats.h"
#include "TMinuit.h"
#include "TPostScript.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
//#include "EveTree.h"
//GMA161203 #include "EveTreeG4.h"
#include "StraightLineFit.h"
#include "TMath.h"
#include <vector>

struct CauInfo{
  int CENum;
  int CauTime;
  int selectline;
  double CauVal[12];
};

using namespace std;


const int plot_level=90; // Which figure want to see, more level means less plots


const double cval = 29.979; // velocity of light in cm/ns

#ifdef ONETIMEFORMAT
const int  nTDCpLayer=1; // Number of TDC per layer
#elif defined(NEWRPCDAQ) || defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
const int  nTDCpLayer=8;
const int nhitspchannel=10;
#else
const int  nTDCpLayer=2;
#endif

const int  nlayer =12; //Maximum RPC layers
const int nxybutton = 9;
const int nsplit = 10;

#ifdef ONEMBYONEM
const int  nstrip =32; // Strips in RPC chamber
const int  nusedstrip =32; // Strips connected with Amplifier/ADC in RPC chamber
const int  lastXstrip =32; // last strip in X
const int  lastYstrip =32; // last strip in Y
#else
const int  nstrip = 64; // Strips in RPC chamber
const int  nwidth = 60;
const int  nusedstrip =64; // Strips connected with Amplifier/ADC in RPC chamber
const int  lastXstrip =57;//60; // last strip in X
const int  lastYstrip =60;//63; // last strip in Y
#endif
const int nTermResSplit = nstrip/20 + 1;
bool isOnlyCom  = false;//false;//false;//true;//false;//true; //Combined fit only, no iteration, not proper efficiency //jim jim
bool isTimeCorOrReso = true; // 10th Aug 2016, must be true always //false; //true; // true: Time offset Correction and
                                    // false: Time resolution, they should come one by one
const int  firstlayer =0; //1; //Used first layer for efficiency calculation
const int  lastlayer =9; // 10; //Used last layer for efficiency calculation

const int  layfirst =0; // 1; //Used first layer in track fitting
const int  laylast =9; // 10; //Used last layer in track fitting

const int  firstXstrip =0;  //1st strip in X
const int  firstYstrip =0;  //1st strip in Y

const double  mxchisq =2.0;   //Maximum value of Normalised chi^2 (position fit);
const double  mxtimechisq=2.0; // Maximum value of Normalised chi^2 (time fit);
const double  accprng =1.0; //Acceptance region with respect tot strip centre
const double  effirng =0.25; //Additional region for efficeincy calculation
                           //Extrapolation error < 0.2 of strip width
const int xtcorstr=0;    //Starting layer for X-time correction
const int xtcorend=11;   // End layer for X-time correction
const int ytcorstr=0;    //Starting layer for Y-time correction
const int ytcorend=11;   //End layer for Y-time correction

const bool isCombinedOffset=false; //true : change xoff[iocc]/yoff[iocc], otherwise align_xstr_xdev[iocc][0] etc

const int isLast=0; // 0: also correct in last iteration, 1: no correction in last iteration //jim jim make it to 1

vector<int>deadchannel[2][nlayer];

double find_max_gap(int ix) {
  if (ix >=layfirst && ix <=laylast) return 1.5;
  if (ix <layfirst) return 1.5+0.3*(layfirst-ix);
  if (ix >layfirst) return 1.5+0.3*(ix-laylast);
  return 1.5;
}

#if defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
  //GMA 171022
const int nPosAlignPar=3;
//correction are from aug 25 26 27 jj46
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.378799, -0.0103435, 0.0294646},
					       {-0.289823, 0.005905, -0.128604},
					       {-0.0938795, -0.0015088, -0.110556},
					       {0.125222, -0.00877913, -0.208487},
					       {0.203481, -0.0029897, 0.0656587},
					       {0.0177999, -0.00396816, -0.102375},
					       {0.0179844, -0.00333615, -0.10059},
					       {0.487905, -0.00587077, 0.384449},
					       {0.0196337, 0.0083739, 0.24446},
					       {-0.178701, 0.00360379, -0.0497524},
					       {0.131686, 0.00393464, 0.238113},
					       {-0.050804, 0.000793465, -0.0253505}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.239299, 0.0107629, 0.0681905},
					       {0.242562, -0.00579628, 0.0692311},
					       {0.0944097, 0.00126593, 0.13656},
					       {-0.120545, 0.00678728, 0.085419},
					       {-0.210893, 0.00151289, -0.149695},
					       {-0.144272, 0.00472849, -0.00419812},
					       {-0.0853803, 0.00301556, -0.0325722},
					       {-0.102275, 0.000631511, -0.0904438},
					       {-0.147314, -0.00573528, -0.349159},
					       {0.327296, -0.0032148, 0.14547},
					       {-0.0968242, -0.00258652, -0.204518},
					       {-0.211821, -0.000860388, -0.27855}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/
//correction are from aug 25 26 27 jj42
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.380028, -0.0101302, 0.0342964},
					       {-0.289683, 0.00611041, -0.123519},
					       {-0.0946153, -0.00131641, -0.106738},
					       {0.123903, -0.0085883, -0.206388},
					       {0.201476, -0.00286327, 0.0674395},
					       {0.0144109, -0.00384962, -0.102246},
					       {0.0149194, -0.00315379, -0.100528},
					       {0.48219, -0.00576293, 0.382265},
					       {0.0119027, 0.00851281, 0.242229},
					       {-0.188439, 0.00372509, -0.0526404},
					       {0.120185, 0.00401643, 0.233304},
					       {-0.0627563, 0.000883728, -0.0313611}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.235746, 0.0106177, 0.0703374},
					       {0.246702, -0.00591186, 0.0737502},
					       {0.0996236, 0.00116085, 0.141489},
					       {-0.130499, 0.00669808, 0.0759722},
					       {-0.207662, 0.00142719, -0.148673},
					       {-0.143828, 0.00475575, -0.0013846},
					       {-0.0825713, 0.00299854, -0.030799},
					       {-0.0983604, 0.000548333, -0.0873971},
					       {-0.157883, -0.00578955, -0.358204},
					       {0.326391, -0.00307384, 0.147971},
					       {-0.0977867, -0.00242088, -0.203945},
					       {-0.212378, -0.000888813, -0.279605}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/
//from october 16 and jj36 conditions L8 from aug 25 jj41
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.379996, -0.0101658, 0.0320637},
					       {-0.290662, 0.00607966, -0.124038},
					       {-0.095168, -0.00115466, -0.103742},
					       {0.120866, -0.00852006, -0.206407},
					       {0.200634, -0.00284789, 0.0651593},
					       {0.0133795, -0.00368666, -0.0982465},
					       {0.0127551, -0.00300147, -0.100935},
					       {0.482941, -0.00590716, 0.377956},
					       {0.0115232, 0.00854454, 0.241911},
					       {-0.189572, 0.00384019, -0.0544424},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.239384, 0.0106047, 0.0682793},
					       {0.242856, -0.00591861, 0.0697491},
					       {0.0956601, 0.0011521, 0.138417},
					       {-0.114276, 0.00665248, 0.0901825},
					       {-0.209517, 0.00149976, -0.146055},
					       {-0.150334, 0.00460343, -0.0106465},
					       {-0.0846636, 0.00290989, -0.0310789},
					       {-0.100443, 0.000627081, -0.0882519},
					       {-0.144772, -0.00572179, -0.347843},
					       {0.32546, -0.00318339, 0.143687},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/
//from october 16 and jj36 conditions
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.379996, -0.0101658, 0.0320637},
					       {-0.290662, 0.00607966, -0.124038},
					       {-0.095168, -0.00115466, -0.103742},
					       {0.120866, -0.00852006, -0.206407},
					       {0.200634, -0.00284789, 0.0651593},
					       {0.0133795, -0.00368666, -0.0982465},
					       {0.0127551, -0.00300147, -0.100935},
					       {0.482941, -0.00590716, 0.377956},
					       {-0.229255, 0.00634179, -0.0658492},
					       {-0.189572, 0.00384019, -0.0544424},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.239384, 0.0106047, 0.0682793},
					       {0.242856, -0.00591861, 0.0697491},
					       {0.0956601, 0.0011521, 0.138417},
					       {-0.114276, 0.00665248, 0.0901825},
					       {-0.209517, 0.00149976, -0.146055},
					       {-0.150334, 0.00460343, -0.0106465},
					       {-0.0846636, 0.00290989, -0.0310789},
					       {-0.100443, 0.000627081, -0.0882519},
					       {0.0583059, -0.00368594, -0.0763008},
					       {0.32546, -0.00318339, 0.143687},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};

*/

//correction for miical
/*

double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.379996, -0.0101658, 0.0320637},
					       {0.24411, 0.00135583, 0.331437},
					       {-0.206852, -0.00241539, -0.256665},
					       {0.577799, 0.00452339, 0.623225},
					       {0.132152, 0.000202012, 0.117742},
					       {0.00508948, 0.00640389, 0.2212},
					       {-0.0995792, 0.00301801, -0.0259608},
					       {-0.075129, 0.00613611, 0.118225},
					       {0.257143, -0.0101926, -0.120563},
					       {-0.276869, 0.00803263, -0.0166791},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.239384, 0.0106047, 0.0682793},
					       {-0.135319, -0.00788478, -0.365214},
					       {0.0820874, 0.00241216, 0.156799},
					       {-0.168961, -0.00321617, -0.263009},
					       {-0.266538, 0.00390449, -0.131707},
					       {-0.128513, -0.00468823, -0.247748},
					       {0.025417, -0.00387079, -0.107385},
					       {0.0647813, -0.0040596, -0.0519331},
					       {0.00974267, 0.00953877, 0.311304},
					       {0.248968, -0.00831818, -0.0537647},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/
/*




double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.341701, 0.00232331, 0.427576},
{0.0990857, 0.00308581, 0.270982},
{0.376454, -0.00116143, 0.354883},
{0.478277, 0.00506045, 0.596225},
{0.0604978, 0.00120239, 0.0781356},
{-0.0513733, 0.00668643, 0.182061},
{-0.12846, 0.00277369, -0.0345766},
{-0.0822276, 0.00485041, 0.125525},
{0.267324, -0.0114804, -0.1115},
{-0.246483, 0.00679948, 0.0146419},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.247221, -0.00468875, -0.41369},
{0.0295191, -0.00900403, -0.258535},
{-0.349526, -0.00215533, -0.316359},
{-0.0465722, -0.00374445, -0.181353},
{-0.164872, 0.00352517, -0.0562332},
{-0.0500353, -0.0047694, -0.194131},
{0.0789477, -0.00346478, -0.0642318},
{0.0989317, -0.00374362, -0.0226601},
{0.0154978, 0.00989614, 0.336025},
{0.235398, -0.00733953, -0.0552118},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
{-0.103852, -8.21954e-05, -0.10216},
{-0.0390945, -0.00061064, -0.0339788},
{-0.00811012, 0.00190294, 0.104191},
{-0.0467445, 0.00245776, 0.0367095},
{-0.0151288, 0.00184321, 0.034453},
{-0.0323859, -0.00264712, -0.0740822},
{-0.0772644, 0.00315467, 0.00723611},
{-0.106228, 0.00321937, 0.00332138},
{-0.141061, 0.00426494, -0.0778646},
{-0.199577, 0.00860155, 0.0488217},
{-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
{0.147766, 7.88917e-05, 0.0891817},
{0.120254, 3.79443e-05, 0.0253254},
{-0.0221867, 0.00323213, -0.0412028},
{0.0377476, -0.00382866, -0.0773817},
{0.0315155, -0.00407352, -0.0725246},
{-0.0137578, -0.000103093, -0.039819},
{-0.0466412, -0.00037037, -0.0487983},
{0.03916, 0.00309259, 0.0948947},
{0.0263892, 0.00135496, 0.0742521},
{0.0918279, -0.00586882, -0.056975},
{0.147549, -0.000988079, 0.128039}};


*/

/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.328758, 0.000809782, 0.386502},
					       {0.0892291, 0.00177623, 0.237656},
					       {0.368131, -0.00224067, 0.324536},
					       {0.470926, 0.00431404, 0.569934},
					       {0.0546626, 0.000774937, 0.0566221},
					       {-0.0554943, 0.00655762, 0.166534},
					       {-0.129729, 0.00286287, -0.0437164},
					       {-0.080753, 0.00515392, 0.122512},
					       {0.271362, -0.0110227, -0.109331},
					       {-0.238964, 0.00733988, 0.0220615},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.257678, -0.00224908, -0.384689},
					       {0.0167331, -0.00684058, -0.237485},
					       {-0.362322, -0.000303646, -0.30068},
					       {-0.0595347, -0.00228428, -0.171298},
					       {-0.178324, 0.00461483, -0.0522813},
					       {-0.0640976, -0.00406399, -0.19699},
					       {0.0635583, -0.00306558, -0.0741628},
					       {0.0823044, -0.00362591, -0.0396648},
					       {-0.00217382, 0.00981419, 0.313289},
					       {0.216225, -0.00756014, -0.0830984},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.571025, -0.00390976, 0.444798},
					       {0.0916012, 0.000765553, 0.186009},
					       {0.365639, -0.00278033, 0.288237},
					       {0.470517, 0.0045195, 0.556274},
					       {0.0270887, 0.00673764, 0.230168},
					       {-0.0703424, 0.00731032, 0.188308},
					       {-0.146691, 0.00451489, -0.00648148},
					       {-0.0957742, 0.00706601, 0.174481},
					       {0.248108, -0.00877271, -0.0462803},
					       {-0.165089, -0.00069326, -0.15252},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.263288, 0.00220807, -0.171511},
					       {0.0895326, -0.00642524, -0.0839684},
					       {-0.303759, 1.44132e-06, -0.180263},
					       {-0.00992936, -0.0023968, -0.0862324},
					       {-0.178237, -0.00228524, -0.243354},
					       {-0.032866, -0.0053993, -0.184569},
					       {0.0889625, -0.00521046, -0.0959978},
					       {0.0959165, -0.0062108, -0.0947685},
					       {0.00818231, 0.00692712, 0.233862},
					       {0.226615, -0.000228641, 0.160076},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.535634, -0.00126845, 0.485571},
					       {0.060782, 0.00319624, 0.215236},
					       {0.349474, -0.000924094, 0.320466},
					       {0.456437, 0.00522901, 0.578349},
					       {0.0368425, 0.00767687, 0.260792},
					       {-0.0660684, 0.00775854, 0.197603},
					       {-0.143281, 0.00428311, -0.0109445},
					       {-0.0967556, 0.00699896, 0.160918},
					       {0.261594, -0.00935492, -0.0579318},
					       {-0.152459, -0.00146801, -0.167109},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.2417, 0.00117148, -0.197979},
					       {0.0977068, -0.00717253, -0.125698},
					       {-0.285279, -0.00126614, -0.209379},
					       {0.000135214, -0.00331508, -0.1153},
					       {-0.166418, -0.00267747, -0.256898},
					       {-0.0257564, -0.00548652, -0.184443},
					       {0.0909344, -0.00491532, -0.0912606},
					       {0.0962699, -0.00564693, -0.0833484},
					       {0.0133018, 0.0075939, 0.259944},
					       {0.223752, 0.000886112, 0.184856},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/

//corrections from 20180731-1130hrs
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.366067, 0.00710533, 0.528621},
					       {0.0332437, 0.00216177, 0.101226},
					       {0.343513, 0.00499676, 0.462708},
					       {0.511087, -0.00480403, 0.287478},
					       {-0.0444418, 0.0091628, 0.193994},
					       {0.213151, 0.00798198, 0.433},
					       {0.00543495, 0.00296659, 0.0623518},
					       {-0.161207, 0.00453709, -0.0308122},
					       {-0.0932676, 0.00376154, -0.0308502},
					       {-0.0356376, 7.2643e-05, -0.0408756},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.305639, -0.00650928, -0.546251},
					       {0.068368, -0.0061077, -0.143798},
					       {-0.32493, -0.00696254, -0.446205},
					       {-0.210888, 0.00695433, -0.0150832},
					       {-0.0702994, -0.00327929, -0.214123},
					       {0.0524762, -0.00461105, -0.104508},
					       {0.187309, -0.00260529, 0.0441212},
					       {0.0892275, -0.00170463, 0.00190568},
					       {0.129068, -0.000628512, 0.0826219},
					       {0.186677, 0.00114981, 0.121983},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/



/*
//corrections from 20180916-2038
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.585985, -0.00169607, 0.512833},
					       {0.120088, 0.00110314, 0.181508},
					       {0.547426, -0.00295906, 0.431936},
					       {0.245991, 0.00193752, 0.241715},
					       {0.00182568, 0.00864335, 0.231845},
					       {0.215798, 0.00848294, 0.460873},
					       {-0.146306, 0.00490589, -0.0310674},
					       {0.0397733, 0.00218251, 0.0898007},
					       {0.0210028, -0.000151628, -0.0464456},
					       {-0.197384, 0.00424803, -0.0769757},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.418458, 2.61958e-05, -0.436838},
					       {0.0531373, -0.00556135, -0.143379},
					       {-0.427155, 0.00091669, -0.304793},
					       {-0.10553, -0.000195184, -0.13276},
					       {-0.0559154, -0.00282129, -0.180241},
					       {0.0838446, -0.00507525, -0.0871822},
					       {0.180392, -0.00449147, -0.0215089},
					       {0.0594571, 0.000765004, 0.0505406},
					       {0.0306251, 0.00335537, 0.107452},
					       {0.303718, -0.00325192, 0.104689},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};

*/
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.881614, 0.00450472, 1.03468},
					       {0.120088, 0.00110314, 0.181508},
					       {0.543574, 0.00199643, 0.644781},
					       {0.616936, -0.00348837, 0.490372},
					       {0.0293914, 0.00880089, 0.306764},
					       {0.00817349, 0.00594236, 0.214112},
					       {-0.109646, 0.000740061, -0.0897906},
					       {-0.0968437, 0.00296114, 0.0136208},
					       {-0.0220174, 0.00659637, 0.158857},
					       {-0.0185803, -0.00237224, -0.0749279},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.270166, -0.00590999, -0.488233},
					       {0.0531373, -0.00556135, -0.143379},
					       {-0.330107, -0.00555126, -0.423603},
					       {-0.118712, 0.00378253, -0.0325714},
					       {-0.118041, -0.00423728, -0.301604},
					       {0.0818528, -0.00371114, -0.0542869},
					       {0.158509, -0.00130627, 0.0524673},
					       {0.160865, -0.000408924, 0.107198},
					       {0.132805, -0.00441858, -0.0376676},
					       {0.176432, 0.00305507, 0.16839},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};


*/
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.866167, 0.00249047, 0.995845},
					       {0.413011, -0.00509625, 0.33388},
					       {0.490783, 0.00137804, 0.578044},
					       {0.57401, -0.00368885, 0.439545},
					       {-0.00431068, 0.00911031, 0.281264},
					       {-0.0146555, 0.00613532, 0.203942},
					       {-0.114802, 0.000793162, -0.084927},
					       {-0.0915452, 0.00239189, 0.0187629},
					       {-0.0149975, 0.0064926, 0.176696},
					       {0.00147261, -0.00240046, -0.0426177},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.279774, -0.00598673, -0.509703},
					       {-0.0224384, -0.000843582, -0.0791265},
					       {-0.344824, -0.00507127, -0.429034},
					       {-0.110397, 0.00400109, -0.0208413},
					       {-0.109114, -0.00458172, -0.299551},
					       {0.0903889, -0.00420996, -0.0562883},
					       {0.166045, -0.00176419, 0.0494738},
					       {0.163433, -0.00069024, 0.104916},
					       {0.139178, -0.00475873, -0.0413109},
					       {0.184375, 0.00258219, 0.15971},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/
/*

// from RPCv4t_evtraw-20181102_174529.rre
double align_xstr_ydev[nlayer][nPosAlignPar] ={{1.04347, 0.00858077, 1.34433},
					       {0.413011, -0.00509625, 0.33388},
					       {0.490783, 0.00137804, 0.578044},
					       {0.57401, -0.00368885, 0.439545},
					       {-0.00431068, 0.00911031, 0.281264},
					       {-0.0146555, 0.00613532, 0.203942},
					       {-0.114802, 0.000793162, -0.084927},
					       {-0.0915452, 0.00239189, 0.0187629},
					       {-0.0149975, 0.0064926, 0.176696},
					       {0.00147261, -0.00240046, -0.0426177},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.210527, -0.00887566, -0.518991},
					       {-0.0224384, -0.000843582, -0.0791265},
					       {-0.344824, -0.00507127, -0.429034},
					       {-0.110397, 0.00400109, -0.0208413},
					       {-0.109114, -0.00458172, -0.299551},
					       {0.0903889, -0.00420996, -0.0562883},
					       {0.166045, -0.00176419, 0.0494738},
					       {0.163433, -0.00069024, 0.104916},
					       {0.139178, -0.00475873, -0.0413109},
					       {0.184375, 0.00258219, 0.15971},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
					       {-0.103852, -8.21954e-05, -0.10216},
					       {-0.0390945, -0.00061064, -0.0339788},
					       {-0.00811012, 0.00190294, 0.104191},
					       {-0.0467445, 0.00245776, 0.0367095},
					       {-0.0151288, 0.00184321, 0.034453},
					       {-0.0323859, -0.00264712, -0.0740822},
					       {-0.0772644, 0.00315467, 0.00723611},
					       {-0.106228, 0.00321937, 0.00332138},
					       {-0.141061, 0.00426494, -0.0778646},
					       {-0.199577, 0.00860155, 0.0488217},
					       {-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
					       {0.147766, 7.88917e-05, 0.0891817},
					       {0.120254, 3.79443e-05, 0.0253254},
					       {-0.0221867, 0.00323213, -0.0412028},
					       {0.0377476, -0.00382866, -0.0773817},
					       {0.0315155, -0.00407352, -0.0725246},
					       {-0.0137578, -0.000103093, -0.039819},
					       {-0.0466412, -0.00037037, -0.0487983},
					       {0.03916, 0.00309259, 0.0948947},
					       {0.0263892, 0.00135496, 0.0742521},
					       {0.0918279, -0.00586882, -0.056975},
					       {0.147549, -0.000988079, 0.128039}};
*/
//RPCv4t_evtraw-20181129_153631.rre
/*

double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.968782, -0.000680881, 0.974815},
{0.410075, -0.00498523, 0.309099},
{0.498086, 0.000816655, 0.556806},
{0.579883, -0.00429872, 0.420375},
{0.00459934, 0.00855027, 0.27595},
{-0.00737111, 0.00576984, 0.198652},
{-0.109568, 0.000694721, -0.0819929},
{-0.0885618, 0.00247219, 0.0272743},
{-0.01358, 0.00676457, 0.184621},
{0.00503283, -0.0020055, -0.0291174},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.384177, -0.00353378, -0.529503},
{-0.0170089, -0.000193468, -0.0556963},
{-0.331427, -0.00501745, -0.408608},
{-0.111291, 0.00422686, -0.00726831},
{-0.11229, -0.00452714, -0.29422},
{0.0897065, -0.00433023, -0.0545302},
{0.166314, -0.00221376, 0.047029},
{0.166237, -0.00103209, 0.10089},
{0.145338, -0.00532243, -0.044865},
{0.197034, 0.00201745, 0.16028},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0667488, 0.00142772, 0.0273393},
{-0.103852, -8.21954e-05, -0.10216},
{-0.0390945, -0.00061064, -0.0339788},
{-0.00811012, 0.00190294, 0.104191},
{-0.0467445, 0.00245776, 0.0367095},
{-0.0151288, 0.00184321, 0.034453},
{-0.0323859, -0.00264712, -0.0740822},
{-0.0772644, 0.00315467, 0.00723611},
{-0.106228, 0.00321937, 0.00332138},
{-0.141061, 0.00426494, -0.0778646},
{-0.199577, 0.00860155, 0.0488217},
{-0.132723, 0.00203228, 0.0491761}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.00453425, -0.00349104, -0.209847},
{0.147766, 7.88917e-05, 0.0891817},
{0.120254, 3.79443e-05, 0.0253254},
{-0.0221867, 0.00323213, -0.0412028},
{0.0377476, -0.00382866, -0.0773817},
{0.0315155, -0.00407352, -0.0725246},
{-0.0137578, -0.000103093, -0.039819},
{-0.0466412, -0.00037037, -0.0487983},
{0.03916, 0.00309259, 0.0948947},
{0.0263892, 0.00135496, 0.0742521},
{0.0918279, -0.00586882, -0.056975},
{0.147549, -0.000988079, 0.128039}};


*/

//double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.}};
// double align_ystr_xdev[nlayer][nPosAlignPar] ={{0.}};
// double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.}};
// double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.}};

/*
//corrections from miical_20190103_092942 trg L6789
 double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.974357, -0.000780252, 0.952204},
					       {0.407603, -0.00464439, 0.296278},
					       {0.502132, 0.000754752, 0.551055},
					       {0.571836, -0.00419364, 0.421177},
					       {0.0100809, 0.00822278, 0.265268},
					       {-0.00439638, 0.00578138, 0.194304},
					       {-0.106954, 0.000482919, -0.0800397},
					       {-0.090414, 0.00251342, 0.0403062},
					       {-0.0125323, 0.0069579, 0.187733},
					       {0.00879383, -0.00183703, -0.0215471},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.395435, -0.00306952, -0.519915},
					       {-0.0231225, 1.63967e-05, -0.0480968},
					       {-0.336342, -0.00470494, -0.404303},
					       {-0.115211, 0.00447614, 0.00129037},
					       {-0.110903, -0.004409, -0.287381},
					       {0.0943966, -0.00455995, -0.0525063},
					       {0.172432, -0.00231969, 0.0531122},
					       {0.172625, -0.00134351, 0.106022},
					       {0.156824, -0.0055534, -0.0481289},
					       {0.204824, 0.00186769, 0.162358},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
*/
// double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0532247, 0.000850579, 0.0727853},
// 					       {-0.13343, 0.00131589, -0.0544304},
// 					       {-0.0521189, 0.000566702, -0.0168752},
// 					       {0.0186944, 0.00103776, 0.0380581},
// 					       {-0.0383812, 0.00188794, 0.0359932},
// 					       {-0.0351564, 0.00298708, 0.0775197},
// 					       {-0.0365641, -0.00225282, -0.0719991},
// 					       {-0.0681306, 0.00199715, 0.00162731},
// 					       {-0.146001, 0.0066665, 0.0757529},
// 					       {-0.133744, 0.0036995, -0.109285},
// 					       {-0.199577, 0.00860155, 0.0488217},
// 					       {-0.132723, 0.00203228, 0.0491761}};
// double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.0495451, -0.00837708, -0.241049},
// 					       {0.149521, 0.000446819, 0.080566},
// 					       {0.0953866, -0.000301029, 0.131398},
// 					       {0.0466111, -0.00392555, -0.109068},
// 					       {-0.0348838, 0.00243936, -0.00635519},
// 					       {-0.0315368, 0.00217542, -0.0295113},
// 					       {-0.0170854, -4.14456e-05, -0.0260342},
// 					       {-0.00423381, -0.00467558, -0.0859551},
// 					       {0.0201621, 0.00530587, 0.111098},
// 					       {0.0260645, 0.00224716, 0.0685436},
// 					       {0.0918279, -0.00586882, -0.056975},
// 					       {0.147549, -0.000988079, 0.128039}};

//corrections from RPCv4t_evtraw-20181226_192924_miical_jj14_3.root
/*




double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.893067, -0.001879, 0.900119},
{0.562203, -0.00614427, 0.439102},
{0.600602, -0.000126263, 0.651504},
{0.568652, -0.00510885, 0.410989},
{-0.00312035, 0.00782, 0.258904},
{-0.0288475, 0.00535717, 0.170647},
{-0.14004, 0.000312532, -0.113553},
{-0.148819, 0.00241532, -0.0323629},
{0.0373388, 0.00676629, 0.246137},
{0.0321182, -0.00204897, 0.00499599},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.371309, -0.00241894, -0.481752},
{-0.158913, 0.000610099, -0.175186},
{-0.389387, -0.00418929, -0.449643},
{-0.0828907, 0.00487446, 0.0366989},
{-0.104944, -0.00410707, -0.278469},
{0.131809, -0.00417544, -0.0118498},
{0.164401, -0.00209133, 0.0363426},
{0.192258, -0.00121329, 0.122072},
{0.159859, -0.00544231, -0.0387917},
{0.183074, 0.0015277, 0.133429},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0}};




*/

double xpos_xtot_corr[3][nlayer-2][2] = {{
    {-0.426208 , 0.425733},
    {-0.465268 , 0.460681},
    {-0.964307 , 0.954272},
    {-0.382796 , 0.379024},
    {-0.518372 , 0.516541},
    {-0.480678 , 0.476589},
    {-0.666214 , 0.665285},
    {-0.404797 , 0.403529},
    {-0.64769 , 0.645748},
    {-0.355035 , 0.350462}
},

{
  {-1.18402 , 1.19949},
  {-1.39774 , 1.38315},
  {-2.23684 , 2.22553},
  {-1.27644 , 1.27282},
  {-1.41497 , 1.41218},
  {-1.17814 , 1.16478},
  {-1.34673 , 1.33694},
  {-1.36057 , 1.34678},
  {-1.46523 , 1.46342},
{-1.24492 , 1.24425}
},

{
  {-0.769098 , 0.75687},
  {-0.800793 , 0.780863},
  {-4.60018 , 4.6002},
  {-0.801697 , 0.780636},
  {-0.842961 , 0.836272},
  {-0.968286 , 0.933522},
  {-0.931996 , 0.934615},
  {-0.766274 , 0.744822},
  {-0.855745 , 0.844288},
{-0.927782 , 0.922621}
}};
double ypos_ytot_corr[3][nlayer-2][2] = {{
  {-0.662704 ,  0.660021},
  {-0.644938 ,  0.640753},
  {-1.11468 ,  1.11002},
  {-0.454414 ,  0.452745},
  {-0.681086 ,  0.682133},
  {-0.723326 ,  0.721711},
  {-0.726125 ,  0.721957},
  {-0.757095 ,  0.751317},
  {-0.673318 ,  0.668337},
  {-0.832125 ,  0.830638},
},

{
  {-1.01717 , 1.00661},
  {-1.00765 , 1.00276},
  {-2.3526 , 2.3449},
  {-1.17229 , 1.16379},
  {-1.0902 , 1.07228},
  {-0.978586 , 0.97538},
  {-1.21438 , 1.20635},
  {-1.11001 , 1.10223},
  {-1.18666 , 1.17626},
{-0.935049 , 0.934758}
},

{
  {-0.968616 , 0.941825},
    {-0.957148 , 0.935234},
      {-3.83931 , 3.80122},
  {-0.786205 , 0.758219},
  {-1.01094 , 0.987804},
  {-0.94055 , 0.934427},
  {-1.05563 , 1.0423},
  {-1.06539 , 1.0551},
  {-0.925681 , 0.89153},
{-0.96418 , 0.9493}
}};


double xpos_xtime_corr[3][nlayer-2][2] = {
  {
    {0.00202779 , -0.0812821},
    {-0.00138013 , -0.06553},
    {-0.000517067 , -0.0524363},
    {-0.0023242 , -0.0669679},
    {-0.00410599 , -0.0758268},
    {-0.00227463 , -0.0647343},
    {0.000668721 , -0.0721951},
    {-0.00111671 , -0.0650791},
    {-0.00113084 , -0.0723436},
    {-0.00155261 , -0.0803685}
  },
  {
    {0.0107042 , -0.154162},
    {0.00571975 , -0.127032},
    {-0.0182512 , -0.103338},
    {-0.00241358 , -0.155436},
    {-0.0125573 , -0.139422},
    {-0.00714964 , -0.130629},
    {-0.00396996 , -0.138947},
    {-0.0115882 , -0.16071},
    {0.00350838 , -0.134407},
    {0.0165607 , -0.146794}
  },
  {
    {-0.0253024 , -0.0762112},
    {-0.0144984 , -0.0819222},
    {-0.0100224 , -0.119589},
    {-0.0221553 , -0.0565737},
    {-0.0127251 , -0.0911266},
    {-0.0314252 , -0.0874673},
    {-0.0143739 , -0.0762083},
    {-0.0290836 , -0.0941933},
    {-0.0179193 , -0.0619914},
    {0.0117704 , -0.0946706}
  }};
double ypos_ytime_corr[3][nlayer-2][2] = {
  {
    {0.000308781 ,  -0.0767619},
    {0.00182891 ,  -0.0807462},
    {-0.00132333 ,  -0.0609557},
    {-0.00117311 ,  -0.0726113},
    {0.00362737 ,  -0.0819477},
    {0.000423071 ,  -0.0806853},
    {0.000759258 ,  -0.0814152},
    {-0.00190833 ,  -0.0793292},
    {-0.00279785 ,  -0.075269},
    {7.09376e-05 ,  -0.0885898}
  },
  {
    {-0.00535329 , -0.0942072},
    {0.00436107 , -0.0968358},
    {-0.00132676 , -0.107664},
    {0.0020546 , -0.145895},
    {0.00332287 , -0.110053},
    {0.000408995 , -0.0955467},
    {0.00277593 , -0.109928},
    {0.000531758 , -0.0998414},
    {-0.00498136 , -0.106872},
    {0.00255755 , -0.0886}
  },
  {
    {0.000478619 , -0.0779253},
    {-0.00299864 , -0.0889225},
    {-0.01256 , -0.0727393},
    {-0.0162453 , -0.0902802},
    {-0.00384284 , -0.0857089},
    {-0.00727766 , -0.0808863},
    {0.00775395 , -0.0801694},
    {0.00609954 , -0.0843701},
    {-0.0219408 , -0.0770351},
    {-0.000505559 , -0.0772903}
  }};




//corrections from mamta with x vs xdev, y vs ydev , x vs ydev and y vs xdev
/*


double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.902501, -0.00147741, 0.87399},
{0.558146, -0.00570852, 0.433343},
{0.606213, 0.000158929, 0.64626},
{0.564069, -0.00464128, 0.408817},
{-0.000233436, 0.00772455, 0.249388},
{-0.0261012, 0.0054246, 0.167092},
{-0.139968, 0.000253643, -0.115752},
{-0.146435, 0.002353, -0.0215648},
{0.0411409, 0.00673785, 0.242736},
{0.0373509, -0.0019646, 0.00899838},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.373262, -0.00211835, -0.477948},
{-0.16048, 0.000755804, -0.170731},
{-0.386394, -0.00401218, -0.441992},
{-0.0789207, 0.00500542, 0.0454921},
{-0.102793, -0.00413745, -0.275781},
{0.133916, -0.00423813, -0.0107951},
{0.170557, -0.00222154, 0.050112},
{0.197908, -0.00140254, 0.126709},
{0.168129, -0.00567122, -0.0421984},
{0.191933, 0.00159874, 0.141716},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.000239892, 0.000149222, -0.00969302},
{-0.0112579, 0.000485861, 0.0321615},
{0.000244042, -0.000120565, -0.00435596},
{0.00352752, 0.000169026, -0.0187873},
{-0.00712494, 0.000557995, 0.0160631},
{-0.0188293, 0.00138308, 0.0373846},
{0.040021, -0.00409907, -0.0541931},
{-0.0021807, -6.3288e-05, 0.00531951},
{-0.0527747, 0.00435077, 0.0984586},
{0.0053602, 0.00107776, -0.0385897},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.0877337, -0.00663059, -0.0775218},
{-0.0201337, 0.00148209, 0.0136654},
{-0.0318381, 0.000796397, 0.0916057},
{0.0353231, -0.00305892, -0.041516},
{-0.0423195, 0.00303366, 0.0444738},
{-0.023816, 0.00236852, 0.0179043},
{0.00174094, -8.25233e-06, 0.0125138},
{0.0446024, -0.00506314, -0.035822},
{-0.0372589, 0.00425799, 0.0347173},
{-0.00677809, 0.000968722, 0.00417498},
{0, 0, 0},
{0, 0, 0}};

*/

/*


//corrections from 9.8 KV L7,L8
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.902501, -0.00147741, 0.87399},
{0.555227, -0.00561301, 0.433036},
{0.607684, 0.000936813, 0.652025},
{0.560886, -0.00450048, 0.398619},
{-0.00230799, 0.00801645, 0.24742},
{-0.0263059, 0.00542467, 0.171042},
{-0.143405, 0.000417019, -0.122065},
{-0.142163, 0.00261052, -0.0234088},
{0.0434194, 0.00669095, 0.246633},
{0.047474, -0.00205612, 0.0157011},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.373262, -0.00211835, -0.477948},
{-0.159821, 0.000969969, -0.16404},
{-0.395012, -0.00375159, -0.440838},
{-0.0847033, 0.00525925, 0.0477489},
{-0.10607, -0.00394671, -0.274309},
{0.132832, -0.00406575, -0.0115469},
{0.177561, -0.00208875, 0.0524914},
{0.201155, -0.00127788, 0.128932},
{0.166041, -0.00569167, -0.0431354},
{0.199221, 0.00139396, 0.150199},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.000239892, 0.000149222, -0.00969302},
{-0.00914768, 0.000346942, 0.0356797},
{0.00400889, -0.000263189, -0.0133175},
{-0.0214527, 0.00185324, 0.00524513},
{0.000288077, -0.000136783, 0.00805378},
{-0.0134257, 0.00093161, 0.037025},
{0.0397634, -0.00438078, -0.0579788},
{-0.0126058, 0.000685565, 0.0180906},
{-0.0720373, 0.00603004, 0.129647},
{-0.0052696, 0.00169392, -0.022001},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.0877337, -0.00663059, -0.0775218},
{-0.0159146, 0.00110394, 0.00708923},
{-0.0373827, 0.00080425, 0.0644961},
{0.0320611, -0.00245734, -0.0308718},
{-0.0395594, 0.00282949, 0.0305381},
{-0.0210847, 0.00197196, 0.031554},
{0.00464634, 0.00042358, 0.0154522},
{0.0380908, -0.00448772, -0.0282866},
{-0.0462668, 0.00480505, 0.0547828},
{-0.0164713, 0.00167181, 0.0111886},
{0, 0, 0},
{0, 0, 0}};

*/

//corrections for RPC196_20200106_221800_9p9kv
/*


double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.841608, -0.00273157, 0.77174},
{0.484219, -0.00671078, 0.322615},
{0.525742, -0.000355885, 0.541112},
{0.456919, -0.00187056, 0.382976},
{0.217562, 0.00350088, 0.339777},
{-0.0549222, 0.00391517, 0.102113},
{-0.202141, 0.01924, 0.386629},
{-0.111026, -0.00385111, -0.189619},
{0.147128, -0.00154816, 0.0885014},
{-0.167732, 0.0108518, 0.101019},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.358378, -0.000388772, -0.40067},
{-0.14628, 0.00213272, -0.0970309},
{-0.38874, -0.00320549, -0.39087},
{-0.0385159, 0.00218853, 0.0112762},
{-0.121277, -0.00139105, -0.238883},
{0.10325, -0.00318219, 0.0028195},
{0.323539, -0.0129538, -0.182183},
{0.0418548, 0.00458161, 0.175495},
{-0.0422884, 0.00207051, 0.00524224},
{0.576053, -0.0145785, 0.173368},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0805092, -0.000866891, 0.0388358},
{0.0517821, -0.00107588, 0.0529561},
{0.0528449, -0.00189702, 0.00926318},
{-0.0817636, -3.69921e-05, -0.102498},
{0.00245831, 0.00488397, 0.15636},
{-0.00458486, -0.00122981, 0.00509552},
{-0.132489, 0.00105854, -0.0652025},
{-0.0678916, -0.00274581, -0.0918829},
{-0.0504317, 0.00268228, 0.0838753},
{0.302074, 0.00418903, 0.441678},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.0350803, -0.00724877, -0.166699},
{-0.0853027, 0.00228787, -0.0751938},
{-0.117084, 0.0014801, 0.0232025},
{0.0284745, -0.000945407, -0.0239552},
{0.34772, 0.00168454, 0.340045},
{-0.184394, 0.00371255, -0.109268},
{0.157784, 0.000685123, 0.114986},
{-0.0900028, -0.0023703, -0.120984},
{-0.0586435, 0.00712272, 0.0870417},
{-0.186131, 0.00754871, 0.00144318},
{0, 0, 0},
{0, 0, 0}};

*/
//RPC200 HV 9.9 KV L9
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.841608, -0.00273157, 0.77174},
{0.484219, -0.00671078, 0.322615},
{0.525742, -0.000355885, 0.541112},
{0.456919, -0.00187056, 0.382976},
{0.217562, 0.00350088, 0.339777},
{-0.0549222, 0.00391517, 0.102113},
{-0.202141, 0.01924, 0.386629},
{-0.111026, -0.00385111, -0.189619},
{0.147128, -0.00154816, 0.0885014},
{-0.180514, 0.00171295, -0.117173},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.358378, -0.000388772, -0.40067},
{-0.14628, 0.00213272, -0.0970309},
{-0.38874, -0.00320549, -0.39087},
{-0.0385159, 0.00218853, 0.0112762},
{-0.121277, -0.00139105, -0.238883},
{0.10325, -0.00318219, 0.0028195},
{0.323539, -0.0129538, -0.182183},
{0.0418548, 0.00458161, 0.175495},
{-0.0422884, 0.00207051, 0.00524224},
{0.44134, -0.00302878, 0.288393},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.0890824, -0.00179052, 0.00992128},
{0.05667, -0.00162227, 0.0333382},
{0.0603312, -0.00325843, 0.0095202},
{-0.081684, 5.59124e-06, -0.108244},
{-0.00720543, 0.00574159, 0.169142},
{-0.0152517, -9.33265e-05, 0.0208606},
{-0.148171, 0.00263185, -0.0375312},
{-0.0837487, -0.000866755, -0.0546602},
{-0.066307, 0.00474911, 0.127657},
{0.296443, 0.000708385, 0.35578},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.0178726, -0.00812473, -0.207566},
{-0.097703, 0.00150554, -0.107974},
{-0.123736, 0.000963833, -0.00258509},
{0.0298586, -0.00127168, -0.0288285},
{0.348022, 0.00184898, 0.349012},
{-0.18004, 0.00434755, -0.0876657},
{0.168537, 0.00162943, 0.148692},
{-0.0777029, -0.000953664, -0.0725906},
{-0.0406596, 0.00886145, 0.149566},
{-0.438533, 0.00676755, -0.251487},
{0, 0, 0},
{0, 0, 0}};

*/

/*
//jim jim from 20180603
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.690997, -0.00106525, 0.674549},
{0.437359, 6.47988e-05, 0.525401},
{0.637909, -0.00341259, 0.548055},
{0.492247, 0.00399662, 0.583821},
{-0.0126313, 0.00624234, 0.182609},
{-0.0326514, 0.00793064, 0.232742},
{-0.179706, 0.00520673, -0.0329627},
{-0.180646, 0.00762177, 0.0953763},
{0.255662, -0.00783975, -0.0331669},
{-0.00901388, 1.06885e-06, 0.0093732},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.332649, 0.000351126, -0.369147},
{-0.0832198, -0.00493023, -0.278443},
{-0.407805, 0.00101358, -0.303642},
{-0.00691401, -0.00201738, -0.117994},
{-0.0888124, -0.00143886, -0.164162},
{0.130637, -0.00544776, -0.0682027},
{0.195816, -0.00547611, -0.0268042},
{0.228186, -0.00658754, -0.00555108},
{0.0781458, 0.0065722, 0.248743},
{0.246079, -0.000215654, 0.135166},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.108475, 7.5295e-05, 0.0145186},
{-0.0212315, 0.000749986, -0.0742516},
{0.000493299, 0.00102128, -0.0250829},
{0.0239766, -0.000485094, -0.0735852},
{0.0636577, 0.00193423, 0.0934559},
{-0.122299, 0.00259041, -0.0540506},
{-0.0353896, -0.00209887, -0.106426},
{-0.0602833, 0.00212635, 0.017982},
{-0.0418195, 0.00754663, 0.169371},
{0.062212, 0.00368759, 0.124685},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.146005, -0.00728153, -0.344423},
{-0.0663927, 0.00166895, -0.0463195},
{-0.0492019, -0.000448945, 0.0840993},
{0.106776, -0.00296874, 0.0139298},
{-0.0799122, 0.00322292, 0.0123639},
{-0.0138785, 0.00325518, 0.0535273},
{0.0299057, 0.00109996, 0.0602917},
{0.0600953, -0.00321176, 0.0241669},
{-0.0178657, 0.00684659, 0.122235},
{-0.187819, 0.0034329, -0.111824},
{0, 0, 0},
{0, 0, 0}};
*/



//jim jim put from 20181226
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.893067, -0.001879, 0.900119},
{0.562203, -0.00614427, 0.439102},
{0.600602, -0.000126263, 0.651504},
{0.568652, -0.00510885, 0.410989},
{-0.00312035, 0.00782, 0.258904},
{-0.0288475, 0.00535717, 0.170647},
{-0.14004, 0.000312532, -0.113553},
{-0.148819, 0.00241532, -0.0323629},
{0.0373388, 0.00676629, 0.246137},
{0.0321182, -0.00204897, 0.00499599},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.371309, -0.00241894, -0.481752},
{-0.158913, 0.000610099, -0.175186},
{-0.389387, -0.00418929, -0.449643},
{-0.0828907, 0.00487446, 0.0366989},
{-0.104944, -0.00410707, -0.278469},
{0.131809, -0.00417544, -0.0118498},
{0.164401, -0.00209133, 0.0363426},
{0.192258, -0.00121329, 0.122072},
{0.159859, -0.00544231, -0.0387917},
{0.183074, 0.0015277, 0.133429},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0}};
*/

/*
//from jim jim 20181226 jj125 no buttons - all layers tot used. itimeslopeconst=1
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.90007, -0.00131721, 0.866598},
					       {0.558201, -0.00558347, 0.432049},
					       {0.609629, 0.000430343, 0.650158},
					       {0.562167, -0.00455759, 0.404205},
					       {-0.00581387, 0.00781396, 0.242992},
					       {-0.0180289, 0.00530058, 0.173928},
					       {-0.140532, 0.000201023, -0.118965},
					       {-0.146207, 0.00232111, -0.0225976},
					       {0.0440297, 0.00671974, 0.243564},
					       {0.0434118, -0.00201553, 0.0127627},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.373564, -0.00167021, -0.468228},
					       {-0.164041, 0.00108908, -0.164921},
					       {-0.394197, -0.00355925, -0.443333},
					       {-0.0831234, 0.00528951, 0.0484727},
					       {-0.103043, -0.00400829, -0.271215},
					       {0.135136, -0.00420663, -0.00827945},
					       {0.179543, -0.00225904, 0.059352},
					       {0.205114, -0.00140197, 0.132499},
					       {0.170852, -0.00576313, -0.0438699},
					       {0.202955, 0.00151961, 0.146974},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{-0.000420087, 7.65079e-05, -0.02938},
					       {-0.00921644, 0.000230943, 0.029509},
					       {0.0124561, -0.000454053, -0.00686321},
					       {0.0035589, 0.000215914, -0.0163051},
					       {-0.00704848, 0.000581142, 0.0196084},
					       {-0.0208716, 0.00171166, 0.0442355},
					       {0.0344102, -0.00395186, -0.054273},
					       {-0.0111212, 0.000564866, 0.0155978},
					       {-0.0606953, 0.00539653, 0.120142},
					       {-0.00665835, 0.00197361, -0.0280764},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.0974053, -0.00744922, -0.0848548},
					       {-0.0196936, 0.00132236, 0.00865583},
					       {-0.0384626, 0.000742708, 0.0816781},
					       {0.0358451, -0.00301725, -0.0404708},
					       {-0.0402897, 0.00318422, 0.0502752},
					       {-0.032797, 0.00255593, 0.0216704},
					       {0.00187542, 0.000318639, 0.0168422},
					       {0.0411863, -0.00457516, -0.027963},
					       {-0.0463964, 0.00505752, 0.0494659},
					       {-0.0178738, 0.00186994, 0.0129207},
					       {0, 0, 0},
					       {0, 0, 0}};

*/


//from jim jim 20181226 jj160 no buttons - all layers tot used. itimeslopeconst=1 --The error is modified as in the paper --- to test chisq prob
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.896717, -0.00124373, 0.860312},
					       {0.558158, -0.00552902, 0.429407},
					       {0.615939, 0.000484834, 0.654082},
					       {0.560477, -0.004517, 0.400324},
					       {-0.0072233, 0.00783906, 0.239408},
					       {-0.0131662, 0.00531898, 0.176655},
					       {-0.141062, 0.000211421, -0.121643},
					       {-0.144646, 0.00232788, -0.023237},
					       {0.0469097, 0.00672655, 0.244182},
					       {0.0492414, -0.0020106, 0.0162993},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.371748, -0.00142931, -0.460206},
					       {-0.165426, 0.00130044, -0.160972},
					       {-0.399246, -0.00337552, -0.443807},
					       {-0.0823312, 0.00543334, 0.0527669},
					       {-0.102045, -0.00390248, -0.267877},
					       {0.136592, -0.00413981, -0.00565882},
					       {0.187406, -0.00223193, 0.0671337},
					       {0.21114, -0.00140901, 0.137363},
					       {0.173892, -0.00579837, -0.0428568},
					       {0.213107, 0.00145969, 0.15437},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{-0.00133191, -0.000194871, -0.0377101},
					       {-0.00683172, 3.98178e-05, 0.027078},
					       {0.018189, -0.000550896, -0.00306874},
					       {0.00323974, 0.000227607, -0.0149807},
					       {-0.00802073, 0.000731913, 0.0242258},
					       {-0.0228318, 0.00199716, 0.0519345},
					       {0.0255597, -0.00352606, -0.0494424},
					       {-0.0186109, 0.00112422, 0.0255631},
					       {-0.0653504, 0.00606788, 0.135933},
					       {-0.0185197, 0.00274004, -0.0168129},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.102463, -0.00781714, -0.0891241},
					       {-0.0185152, 0.00105877, 0.0034941},
					       {-0.0442381, 0.000580965, 0.0720852},
					       {0.0375989, -0.00307331, -0.0394045},
					       {-0.0394436, 0.00325661, 0.0538554},
					       {-0.0389845, 0.00276792, 0.0217285},
					       {0.000397539, 0.000659379, 0.0249789},
					       {0.036935, -0.00410152, -0.0192734},
					       {-0.0525508, 0.00563551, 0.0590823},
					       {-0.0271728, 0.00254281, 0.021943},
					       {0, 0, 0},
					       {0, 0, 0}};


/*
//from jim jim 20181226 jj179 no buttons - all layers tot used. itimeslopeconst=1 --The error is modified as in the paper --- with chisq/ndf<2 criteria and ndf>=8 during alignment
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.896621, -0.00154837, 0.852423},
					       {0.559962, -0.00558859, 0.42948},
					       {0.616502, 0.000946342, 0.657815},
					       {0.562121, -0.00466892, 0.393126},
					       {-0.00723512, 0.00788305, 0.235675},
					       {-0.00672936, 0.005095, 0.178807},
					       {-0.140543, 0.000140783, -0.126776},
					       {-0.144723, 0.00236797, -0.0249279},
					       {0.0514948, 0.00685604, 0.247355},
					       {0.0539245, -0.00202286, 0.0169482},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.37278, -0.00110373, -0.455792},
					       {-0.168139, 0.00164738, -0.15657},
					       {-0.406809, -0.00341061, -0.448309},
					       {-0.0882572, 0.00574304, 0.0540729},
					       {-0.101651, -0.00373491, -0.261349},
					       {0.138204, -0.00415607, 0.00112797},
					       {0.196745, -0.00209948, 0.0806689},
					       {0.216039, -0.00134319, 0.146172},
					       {0.173297, -0.00564528, -0.0398397},
					       {0.221159, 0.00155851, 0.163112},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{-0.00787539, -0.00018496, -0.0387445},
					       {-0.00889042, 0.00016463, 0.0267919},
					       {0.0223126, -0.00044456, 0.00389219},
					       {-0.00215747, 0.000804982, -0.00307584},
					       {-0.00641841, 0.000668335, 0.0228673},
					       {-0.0221224, 0.00213218, 0.0556144},
					       {0.0194032, -0.00348286, -0.0517149},
					       {-0.0256502, 0.00158342, 0.038039},
					       {-0.0725273, 0.00698692, 0.156404},
					       {-0.0288495, 0.00333603, -0.00778559},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.108027, -0.00820645, -0.0916636},
					       {-0.0205011, 0.000913993, 6.42615e-05},
					       {-0.0522734, 0.000597786, 0.0664851},
					       {0.0384323, -0.00291638, -0.0397018},
					       {-0.0365371, 0.0031379, 0.0516219},
					       {-0.0412208, 0.00252773, 0.0315351},
					       {0.00269691, 0.000783985, 0.0268023},
					       {0.0340636, -0.00367381, -0.0135515},
					       {-0.0617262, 0.00619661, 0.0731881},
					       {-0.0338005, 0.00307091, 0.0288784},
					       {0, 0, 0},
					       {0, 0, 0}};
*/

/*
//from jim jim 20181226 jj113 no buttons
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.902712, -0.00147695, 0.875125},
{0.558546, -0.00571671, 0.435669},
{0.603269, 0.000304661, 0.646736},
{0.564058, -0.00464411, 0.40889},
{-0.00423481, 0.00776017, 0.247087},
{-0.0227705, 0.00527567, 0.17169},
{-0.140127, 0.000205112, -0.116145},
{-0.147717, 0.00231663, -0.021954},
{0.0410866, 0.0067172, 0.242893},
{0.0374556, -0.00200579, 0.00905995},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.374061, -0.00207337, -0.477451},
{-0.161527, 0.000800657, -0.170654},
{-0.388025, -0.00399897, -0.442562},
{-0.0834365, 0.0051149, 0.0430592},
{-0.1038, -0.00413164, -0.27514},
{0.133685, -0.00428195, -0.0110767},
{0.17151, -0.00226432, 0.0517283},
{0.198629, -0.00136704, 0.128137},
{0.167325, -0.00568238, -0.0439927},
{0.192265, 0.00163006, 0.140585},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{-0.000909629, 0.000179711, -0.0109235},
{-0.0116404, 0.000500213, 0.0308979},
{0.000313436, -0.000122903, -0.00491335},
{0.00395315, 0.000295329, -0.0180457},
{-0.005741, 0.000455047, 0.0142622},
{-0.0185686, 0.00143349, 0.0359402},
{0.0436332, -0.00438709, -0.0599592},
{-0.00297027, -5.07286e-05, 0.00534493},
{-0.0553429, 0.00463719, 0.103297},
{0.00594501, 0.00107335, -0.0400967},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0.0869118, -0.00656456, -0.0772511},
{-0.0215733, 0.00154277, 0.0143791},
{-0.0324116, 0.000918096, 0.0926473},
{0.0340515, -0.00293068, -0.0411811},
{-0.0408681, 0.00310311, 0.0466672},
{-0.0262495, 0.00230541, 0.0213953},
{0.00356639, -4.24595e-05, 0.00824002},
{0.045673, -0.00507943, -0.0370847},
{-0.0396657, 0.00441926, 0.0394379},
{-0.00817177, 0.00113285, 0.00330116},
{0, 0, 0},
{0, 0, 0}};
*/



/*
//jim jim put from 22092021 with button
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.737792, -0.00268072, 0.657905},
{0.516886, 0.00135763, 0.564767},
{0.576195, 0.00287897, 0.662223},
{0.447839, 0.0014215, 0.473933},
{0.105688, 0.00604315, 0.294208},
{0.128035, -0.00193063, 0.0892975},
{-0.0653755, -0.00448814, -0.202864},
{-0.0823494, -0.000278279, -0.0534066},
{0.0393771, 0.00145758, 0.0780956},
{-0.018435, -0.001612, -0.0333551},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.307544, -0.00232256, -0.366507},
{-0.145742, -0.00894343, -0.293867},
{-0.332086, -0.000666661, -0.430066},
{0.0525157, -0.00173617, -0.0134424},
{-0.0961493, -0.00419976, -0.291284},
{0.21091, 0.000565457, 0.189671},
{0.0885927, 0.00257498, 0.121435},
{0.16256, 0.00131363, 0.190193},
{0.0600231, -0.000368272, 0.0263365},
{0.404079, -0.000816889, 0.315299},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.133567, -0.00185708, 0.119415},
{-0.175305, -0.00210439, -0.132205},
{0.0525021, 0.00311956, 0.207946},
{-0.0996339, 0.00227975, -0.0247649},
{-0.108562, 0.00842646, 0.157939},
{0.331876, 0.00805335, 0.60983},
{-0.0429877, -0.00139668, -0.0232516},
{-0.0880586, 0.00264704, 0.051433},
{-0.0910401, 0.00747536, 0.196956},
{-0.27478, 0.00272147, -0.199738},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.12003, -0.00634206, -0.321027},
{0.0472891, 0.00129125, 0.219649},
{-0.0495913, 0.0029576, 0.0329466},
{-0.0665951, 0.000463857, -0.0703637},
{0.10514, 0.00209141, 0.140357},
{0.252335, 0.00126073, 0.384551},
{-0.0218036, 0.00495255, 0.130818},
{-0.0931542, 0.00117188, 0.00275089},
{-0.175988, 0.0104345, 0.10283},
{-0.196709, 0.00646535, -0.00365794},
{0, 0, 0},
{0, 0, 0}};
*/

/*
//jim jim put from 20210922 removing buttons and using tot in all layers jj130
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.754389, -0.00261143, 0.676387},
					       {0.515946, 0.00142996, 0.565728},
					       {0.578406, 0.00294917, 0.666249},
					       {0.45505, 0.00148663, 0.482867},
					       {0.103108, 0.00610883, 0.293287},
					       {0.127392, -0.00186669, 0.0901638},
					       {-0.0625184, -0.00442668, -0.198636},
					       {-0.0778195, -0.00021861, -0.0475905},
					       {0.0332435, 0.00151814, 0.0731239},
					       {-0.01642, -0.00155056, -0.0300885},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.305393, -0.0025164, -0.371716},
					       {-0.145988, -0.00911047, -0.300876},
					       {-0.332931, -0.000816311, -0.437067},
					       {0.0608652, -0.00186604, -0.0105756},
					       {-0.0940731, -0.00430113, -0.293937},
					       {0.212385, 0.000484631, 0.187039},
					       {0.0949541, 0.00251254, 0.124228},
					       {0.16828, 0.00126632, 0.192802},
					       {0.0688853, -0.000402766, 0.0324274},
					       {0.436525, -0.000843306, 0.345238},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.117472, -0.00108493, 0.13088},
					       {-0.187144, -0.00131171, -0.116279},
					       {0.0427692, 0.00392943, 0.225925},
					       {-0.117167, 0.00311388, -0.0143637},
					       {-0.118114, 0.0092896, 0.176349},
					       {0.323847, 0.00894009, 0.629972},
					       {-0.0538114, -0.00047941, -0.00543203},
					       {-0.0964959, 0.00358343, 0.071726},
					       {-0.100631, 0.00842718, 0.215963},
					       {-0.305948, 0.00370198, -0.202482},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.138911, -0.00619543, -0.326737},
					       {0.0449861, 0.00149454, 0.231733},
					       {-0.0561889, 0.00323259, 0.0426243},
					       {-0.0793661, 0.000818974, -0.0646355},
					       {0.10071, 0.00253349, 0.15678},
					       {0.244639, 0.0017733, 0.399708},
					       {-0.0339828, 0.00553101, 0.143717},
					       {-0.10812, 0.00180633, 0.0146328},
					       {-0.181375, 0.011101, 0.125594},
					       {-0.211207, 0.00717577, 0.0113109},
					       {0, 0, 0},
					       {0, 0, 0}};
*/

/*
//jim jim put from 20210922 removing buttons
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.746065, -0.00264467, 0.667236},
					       {0.516387, 0.00139477, 0.565317},
					       {0.577304, 0.00291452, 0.664299},
					       {0.451415, 0.00145521, 0.478412},
					       {0.104396, 0.00607555, 0.293758},
					       {0.127686, -0.00189803, 0.0897062},
					       {-0.0639649, -0.00445755, -0.2008},
					       {-0.0800692, -0.000249187, -0.050541},
					       {0.0363273, 0.00148705, 0.0755551},
					       {-0.0173983, -0.0015812, -0.0317495},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.306217, -0.00243565, -0.369274},
					       {-0.145695, -0.00903793, -0.297542},
					       {-0.332365, -0.000750977, -0.433699},
					       {0.0567685, -0.00180648, -0.0121065},
					       {-0.0950796, -0.00425124, -0.292654},
					       {0.211649, 0.000525869, 0.188347},
					       {0.0917223, 0.00254874, 0.122845},
					       {0.165346, 0.00129638, 0.191542},
					       {0.0643592, -0.000378775, 0.0294335},
					       {0.420243, -0.000821606, 0.330431},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.125539, -0.00147938, 0.12518},
					       {-0.181191, -0.00171576, -0.124179},
					       {0.0476569, 0.0035206, 0.216973},
					       {-0.108378, 0.00269409, -0.0195308},
					       {-0.113353, 0.00886035, 0.167118},
					       {0.327832, 0.00849802, 0.619904},
					       {-0.0484102, -0.000937738, -0.0143371},
					       {-0.0923166, 0.00311906, 0.0615808},
					       {-0.0958426, 0.00795724, 0.206412},
					       {-0.290384, 0.00322419, -0.201473},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.12931, -0.00628259, -0.324055},
					       {0.0462936, 0.00137982, 0.225576},
					       {-0.0527917, 0.00308442, 0.0376917},
					       {-0.0729039, 0.000635283, -0.0675297},
					       {0.102964, 0.00230737, 0.14856},
					       {0.248501, 0.00151587, 0.392194},
					       {-0.0279205, 0.00524421, 0.137361},
					       {-0.100714, 0.00149282, 0.00879873},
					       {-0.178789, 0.0107718, 0.11435},
					       {-0.204089, 0.00682539, 0.00394881},
					       {0, 0, 0},
					       {0, 0, 0}};

*/

//20210616 nmnhits 4 layer
/*
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.757487, -0.00343897, 0.657642},
{0.536228, 0.000365619, 0.565123},
{0.579303, 0.00231325, 0.65604},
{0.470304, -0.000734685, 0.430358},
{0.110529, 0.00556911, 0.294584},
{0.094095, -0.00265843, 0.0249921},
{-0.138294, -0.000210474, -0.139331},
{-0.0518632, -0.0015529, -0.0691958},
{0.0438181, 0.00153169, 0.0807514},
{-0.0374635, -0.000610989, -0.0287152},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.308748, -0.00014927, -0.315708},
{-0.142124, -0.00734138, -0.242085},
{-0.317736, 0.000662253, -0.391974},
{0.00513323, 0.000554557, 0.0217904},
{-0.093389, -0.00397791, -0.276653},
{0.163356, 0.0010695, 0.156878},
{0.085353, -0.00126534, 0.00818382},
{0.144314, 0.00200796, 0.210931},
{0.0357161, -0.000902928, -0.00498587},
{0.369362, -0.00180376, 0.252157},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.207383, -0.00287568, 0.118054},
{-0.116426, -0.00318945, -0.137939},
{0.10844, 0.0013424, 0.181991},
{-0.0465966, 3.78816e-05, -0.0569626},
{-0.051825, 0.00633406, 0.136257},
{0.0175829, 0.00694718, 0.305164},
{-0.150772, -0.00463718, -0.200239},
{0.0397613, -0.00124189, 0.106823},
{-0.0673569, 0.00415318, 0.141375},
{-0.0712683, -0.000156747, -0.0638069},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.109119, -0.00786714, -0.352822},
{0.0605249, 8.97977e-05, 0.182537},
{-0.0229331, 0.00181854, 0.00520725},
{0.0257394, -0.000830636, -0.0150502},
{0.112741, 0.00114588, 0.123059},
{-0.137558, 0.0046049, 0.00907458},
{0.0155312, 0.00309437, 0.0964144},
{0.0532545, -0.000976362, 0.0709323},
{-0.200369, 0.00809374, -0.00517366},
{-0.155654, 0.00480042, -0.0399305},
{0, 0, 0},
{0, 0, 0}};

*/
/*
//20210616 nmnhits 5 layer older jj79 with buttons jim jim
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.742252, -0.00215205, 0.664685},
{0.521975, 0.00148942, 0.569442},
{0.571699, 0.00330133, 0.663737},
{0.471552, 0.000139776, 0.441985},
{0.0996403, 0.00622951, 0.292578},
{0.0911996, -0.00278027, 0.011689},
{-0.135885, -1.84512e-05, -0.137365},
{-0.0469259, -0.0014692, -0.0673176},
{0.0455501, 0.00147479, 0.0770551},
{-0.0276557, -0.00073509, -0.0261322},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.289957, -0.00112752, -0.333051},
{-0.131042, -0.00811686, -0.261037},
{-0.307932, -2.21377e-05, -0.403592},
{0.0169934, 0.000120823, 0.0162288},
{-0.0852206, -0.00427646, -0.279299},
{0.164752, 0.00175756, 0.180825},
{0.0671057, -0.00109489, -0.00546781},
{0.144095, 0.00227005, 0.214272},
{0.0314972, -0.000475399, 0.00443122},
{0.375543, -0.00124537, 0.277515},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.192658, -0.00253631, 0.133657},
{-0.125528, -0.00279854, -0.114326},
{0.0948222, 0.00208507, 0.209127},
{-0.058722, 0.000624062, -0.0439148},
{-0.0636077, 0.00709237, 0.152312},
{0.00498836, 0.00833976, 0.314427},
{-0.14224, -0.00300392, -0.176849},
{0.0306129, 0.000350428, 0.11186},
{-0.0733529, 0.00538422, 0.154848},
{-0.0891262, 0.00108843, -0.0620256},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.115168, -0.00743057, -0.339593},
{0.056061, 0.000420674, 0.199022},
{-0.0329112, 0.00237916, 0.0157534},
{0.0260961, -0.00142234, -0.0156406},
{0.109441, 0.00201657, 0.134747},
{-0.133953, 0.00460737, 0.0324201},
{0.00806556, 0.00337675, 0.107735},
{0.0445589, -0.000749053, 0.0818816},
{-0.203437, 0.00826582, 0.0107982},
{-0.16507, 0.00490477, -0.0304384},
{0, 0, 0},
{0, 0, 0}};
*/

/*
//20210616 after removing buttons jj145   jim jim
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.734287, -0.000994931, 0.679548},
					       {0.514582, 0.00233569, 0.576436},
					       {0.570942, 0.0038502, 0.673084},
					       {0.47747, 0.000586943, 0.456015},
					       {0.0934983, 0.00680692, 0.296945},
					       {0.0780449, -0.00274446, 0.00173839},
					       {-0.134813, 0.000227042, -0.13095},
					       {-0.0471177, -0.00126317, -0.0608562},
					       {0.0442888, 0.00169541, 0.0776732},
					       {-0.0224388, -0.000569712, -0.0190136},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.274369, -0.00181582, -0.336287},
					       {-0.123298, -0.00879729, -0.269881},
					       {-0.300757, -0.000450844, -0.408105},
					       {0.024585, -0.000125785, 0.0144533},
					       {-0.0766508, -0.00443493, -0.276907},
					       {0.172848, 0.00219693, 0.198794},
					       {0.0467467, -0.000963447, -0.0240462},
					       {0.145059, 0.00244194, 0.22006},
					       {0.0298219, -0.000208752, 0.00928991},
					       {0.39024, -0.000994353, 0.29912},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.173592, -0.00111896, 0.131434},
					       {-0.139524, -0.00167138, -0.0926014},
					       {0.0795386, 0.00302148, 0.229698},
					       {-0.0702011, 0.00141098, -0.0321397},
					       {-0.0841537, 0.00851902, 0.174765},
					       {-0.0108142, 0.0103187, 0.314852},
					       {-0.128935, -0.00164224, -0.144848},
					       {0.0209451, 0.00172742, 0.124669},
					       {-0.0812458, 0.00688365, 0.171011},
					       {-0.111979, 0.00274722, -0.0641049},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.123341, -0.0061434, -0.328072},
					       {0.0610676, 0.000135512, 0.217324},
					       {-0.0414983, 0.00287884, 0.0290969},
					       {0.01817, -0.00145255, -0.0130714},
					       {0.105676, 0.00265743, 0.151208},
					       {-0.129883, 0.0052788, 0.0593734},
					       {0.00032655, 0.00370933, 0.120819},
					       {0.0366124, -0.000393848, 0.0944181},
					       {-0.208275, 0.00864205, 0.0272388},
					       {-0.176467, 0.0051558, -0.0206991},
					       {0, 0, 0},
					       {0, 0, 0}};
*/

/*
//20210616 after removing buttons jj169   jim jim applting the modified correction to timeserrsq2
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.732783, -0.000671442, 0.683513},
					       {0.513696, 0.00255301, 0.58042},
					       {0.572136, 0.00411452, 0.679944},
					       {0.481783, 0.000803741, 0.464276},
					       {0.0925049, 0.00698606, 0.299846},
					       {0.0678788, -0.00255032, -0.00633496},
					       {-0.133667, 0.000411944, -0.127137},
					       {-0.0465302, -0.0011234, -0.0576148},
					       {0.0431032, 0.00184789, 0.0796582},
					       {-0.0211292, -0.000414395, -0.0151361},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.267122, -0.00181587, -0.332985},
					       {-0.12184, -0.00883709, -0.272182},
					       {-0.298054, -0.000480591, -0.40812},
					       {0.0267296, -0.000160599, 0.015073},
					       {-0.0723455, -0.00444849, -0.274325},
					       {0.180556, 0.00233947, 0.210276},
					       {0.0363235, -0.000979917, -0.0335147},
					       {0.145954, 0.00246182, 0.221399},
					       {0.0301579, -0.000188701, 0.0110401},
					       {0.39823, -0.000913332, 0.309638},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.161277, -0.000512585, 0.132493},
					       {-0.145795, -0.00102742, -0.0857965},
					       {0.0725737, 0.00365035, 0.23397},
					       {-0.0767459, 0.00199663, -0.0264207},
					       {-0.0928223, 0.00916564, 0.177546},
					       {-0.020836, 0.0109495, 0.316337},
					       {-0.122732, -0.000935242, -0.126536},
					       {0.0160346, 0.00244622, 0.132268},
					       {-0.0845528, 0.00766328, 0.17859},
					       {-0.124024, 0.00353575, -0.0633055},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.119462, -0.00599511, -0.320947},
					       {0.0637856, 0.000159676, 0.225145},
					       {-0.043234, 0.003019, 0.0342901},
					       {0.0132343, -0.00128929, -0.010021},
					       {0.104654, 0.00285654, 0.159854},
					       {-0.126034, 0.00567866, 0.074269},
					       {-0.00350681, 0.00392761, 0.126669},
					       {0.0327376, -0.00015504, 0.100971},
					       {-0.211346, 0.0089179, 0.034461},
					       {-0.182859, 0.00544411, -0.0166654},
					       {0, 0, 0},
					       {0, 0, 0}};
*/

/*
//20210616 after removing buttons jj156   jim jim  no layer 5 in fit
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.733406, -0.000745347, 0.684822},
					       {0.514343, 0.00255764, 0.582172},
					       {0.57345, 0.00408913, 0.68067},
					       {0.482242, 0.000685737, 0.463317},
					       {0.0916618, 0.00704014, 0.300633},
					       {0.0666359, -0.00256626, -0.00683824},
					       {-0.133046, 0.000338157, -0.12639},
					       {-0.0459367, -0.00117765, -0.0571446},
					       {0.0424425, 0.00184425, 0.0800366},
					       {-0.0213997, -0.000438174, -0.0128063},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.267045, -0.00188437, -0.333898},
					       {-0.121303, -0.00892284, -0.272787},
					       {-0.297184, -0.000568249, -0.407665},
					       {0.0273274, -0.00020929, 0.0143329},
					       {-0.0736827, -0.00437555, -0.276632},
					       {0.186657, 0.00234928, 0.213278},
					       {0.0381931, -0.000933915, -0.0327115},
					       {0.144887, 0.00248239, 0.217255},
					       {0.0307149, -0.000105332, 0.00833897},
					       {0.400843, -0.00110789, 0.306318},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.162078, -0.000566899, 0.132529},
					       {-0.144319, -0.00111573, -0.0874882},
					       {0.0727532, 0.00367481, 0.232278},
					       {-0.0781614, 0.00207423, -0.0256978},
					       {-0.0915775, 0.00904136, 0.177746},
					       {-0.0240291, 0.0108015, 0.30336},
					       {-0.124687, -0.00119146, -0.125635},
					       {0.017068, 0.00229219, 0.134156},
					       {-0.0853808, 0.0075365, 0.178738},
					       {-0.121868, 0.00330076, -0.0659095},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.118532, -0.00591427, -0.324072},
					       {0.0612517, 0.000476988, 0.223211},
					       {-0.0442532, 0.00312496, 0.0327317},
					       {0.0139666, -0.0012983, -0.00996658},
					       {0.105931, 0.00276949, 0.159873},
					       {-0.1222, 0.00529528, 0.0747056},
					       {-0.00389333, 0.00379859, 0.128972},
					       {0.0317252, -0.000241463, 0.103624},
					       {-0.209548, 0.00862089, 0.0378411},
					       {-0.181738, 0.00526808, -0.0154664},
					       {0, 0, 0},
					       {0, 0, 0}};

*/

/*
//20210528 jj69 jim jim with button jj69
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.741326, -0.00271395, 0.658527},
{0.5174, 0.00112854, 0.561506},
{0.566494, 0.00304164, 0.655278},
{0.453584, 0.000527763, 0.44572},
{0.105824, 0.00583631, 0.294614},
{0.0476212, -0.000817665, 0.0284958},
{-0.134929, -0.000218034, -0.13944},
{-0.0432256, -0.00187107, -0.0713886},
{0.0440649, 0.00132937, 0.0745992},
{-0.0288298, -0.000948241, -0.0295263},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.292372, -0.00131762, -0.322725},
{-0.136898, -0.00789873, -0.254417},
{-0.30155, 0.000271031, -0.38957},
{0.02182, -0.000481617, 0.00337144},
{-0.0958482, -0.00409799, -0.283167},
{0.175426, 0.00224493, 0.202688},
{0.0813561, -0.000986332, 0.0085345},
{0.14981, 0.00232657, 0.222489},
{0.0300992, -0.000418857, 0.00225947},
{0.36889, -0.00121893, 0.267755},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.199196, -0.00283387, 0.121683},
{-0.123919, -0.00297523, -0.120065},
{0.0872819, 0.00149688, 0.187021},
{-0.0292871, -1.15676e-05, -0.0268187},
{-0.0615449, 0.00638985, 0.151255},
{0.0905792, 0.00756859, 0.393321},
{-0.166302, -0.00444434, -0.198722},
{0.0100476, -0.00120994, 0.0969872},
{-0.091702, 0.00446864, 0.139148},
{-0.103316, 0.000197559, -0.0782111},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.133826, -0.00737101, -0.360152},
{0.044098, 0.0005362, 0.187435},
{-0.0369324, 0.00204828, 0.00581743},
{0.0527229, -0.00146354, 0.0133141},
{0.0993133, 0.00199267, 0.130129},
{-0.190161, 0.00449231, -0.0280819},
{0.0103012, 0.00340671, 0.104818},
{0.0531437, -0.000902304, 0.0846551},
{-0.190939, 0.00821574, 0.0117984},
{-0.149253, 0.00491587, -0.024693},
{0, 0, 0},
{0, 0, 0}};
*/

/*
//20210528 jj138 jim jim without button jj138
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.732751, -0.00220113, 0.662984},
					       {0.508599, 0.00151241, 0.561901},
					       {0.562, 0.00334368, 0.657259},
					       {0.451123, 0.000751345, 0.44575},
					       {0.0988991, 0.00612039, 0.292025},
					       {0.0451608, -0.000968574, 0.0236454},
					       {-0.135093, -0.000170833, -0.138361},
					       {-0.0432834, -0.00190646, -0.0685893},
					       {0.0426515, 0.00133419, 0.0727586},
					       {-0.0243046, -0.000986726, -0.0268082},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.28253, -0.0019418, -0.330757},
					       {-0.132766, -0.00842403, -0.265862},
					       {-0.286245, -0.000193007, -0.385715},
					       {0.0228744, -0.00074339, -0.0041864},
					       {-0.0947066, -0.00420312, -0.288231},
					       {0.179138, 0.00265944, 0.217749},
					       {0.0841937, -0.000913688, 0.013637},
					       {0.150672, 0.00241543, 0.222866},
					       {0.0281541, -0.000226497, 0.00696773},
					       {0.372269, -0.00093993, 0.280239},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.190106, -0.00264739, 0.127404},
					       {-0.129471, -0.00269503, -0.106826},
					       {0.0711234, 0.00126883, 0.189006},
					       {-0.030673, 0.000102743, -0.0168343},
					       {-0.0700034, 0.00700442, 0.169556},
					       {0.0812664, 0.00860042, 0.385558},
					       {-0.172059, -0.0039396, -0.19799},
					       {0.00210812, 1.22277e-05, 0.0968575},
					       {-0.0966147, 0.00519645, 0.146991},
					       {-0.115844, 0.00100318, -0.0769351},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.141284, -0.00709274, -0.354738},
					       {0.0388608, 0.000829856, 0.197223},
					       {-0.0424528, 0.00221509, 0.0125764},
					       {0.0516127, -0.00150503, 0.0174032},
					       {0.0964011, 0.00244588, 0.134944},
					       {-0.187165, 0.00457617, -0.0162421},
					       {0.00615766, 0.00350579, 0.111542},
					       {0.0483712, -0.000720275, 0.0908526},
					       {-0.191985, 0.00838863, 0.0224544},
					       {-0.154626, 0.00503014, -0.020174},
					       {0, 0, 0},
					       {0, 0, 0}};

*/

/*
//20210528  jim jim without button jj166 the timeserrsq2 is modified by subtracting from fitted sigma
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.727743, -0.00191896, 0.663637},
					       {0.501845, 0.00176644, 0.560315},
					       {0.558465, 0.00356313, 0.658291},
					       {0.448138, 0.000939455, 0.446752},
					       {0.0930389, 0.00627464, 0.289536},
					       {0.0411767, -0.000959419, 0.0197719},
					       {-0.135612, -7.92945e-05, -0.136949},
					       {-0.0433003, -0.00184188, -0.0671973},
					       {0.0405766, 0.00137943, 0.0717297},
					       {-0.0209551, -0.00095128, -0.022666},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.276036, -0.00229582, -0.335683},
					       {-0.13151, -0.00871485, -0.274211},
					       {-0.27284, -0.000436448, -0.380213},
					       {0.0228802, -0.000932406, -0.0102737},
					       {-0.0939467, -0.00432737, -0.29152},
					       {0.185059, 0.00284768, 0.229449},
					       {0.0879058, -0.000911642, 0.0172271},
					       {0.15127, 0.00246972, 0.225079},
					       {0.0275462, -0.000113489, 0.00933455},
					       {0.378043, -0.000805407, 0.290077},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.180323, -0.00245153, 0.134796},
					       {-0.134208, -0.0024599, -0.0945094},
					       {0.054114, 0.001545, 0.188428},
					       {-0.0344821, 0.00043917, -0.00499123},
					       {-0.0748491, 0.00742672, 0.179299},
					       {0.070293, 0.0094272, 0.385781},
					       {-0.180812, -0.00332606, -0.194066},
					       {-0.00395142, 0.000709477, 0.102847},
					       {-0.101307, 0.0059527, 0.153418},
					       {-0.126653, 0.00181626, -0.0768775},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.147501, -0.00686729, -0.34948},
					       {0.0354352, 0.00103814, 0.204828},
					       {-0.0480749, 0.00241808, 0.0178115},
					       {0.0464386, -0.00130059, 0.0229881},
					       {0.0951744, 0.0026463, 0.144334},
					       {-0.184559, 0.00470618, -0.00390678},
					       {0.00183468, 0.00368764, 0.117537},
					       {0.0444902, -0.000551448, 0.0971481},
					       {-0.192951, 0.00854079, 0.0314545},
					       {-0.160252, 0.00517196, -0.0159308},
					       {0, 0, 0},
					       {0, 0, 0}};
*/


/*
//20210528 jj138 jim jim without button jj153   removing layer 5 from fit
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.727595, -0.0018988, 0.662545},
					       {0.499938, 0.00189511, 0.558729},
					       {0.557952, 0.0035714, 0.657354},
					       {0.448757, 0.000892577, 0.447254},
					       {0.0929927, 0.00635654, 0.290115},
					       {0.0406052, -0.000956772, 0.0199043},
					       {-0.134983, -0.000107744, -0.135708},
					       {-0.0428631, -0.00187415, -0.0670114},
					       {0.0403069, 0.00136852, 0.0714768},
					       {-0.0200769, -0.00102401, -0.0207191},
					       {0.117565, 0.00409309, 0.235423},
					       {-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.274944, -0.0023209, -0.334017},
					       {-0.130932, -0.00882677, -0.274378},
					       {-0.271598, -0.000548379, -0.378358},
					       {0.0225269, -0.000963731, -0.0110375},
					       {-0.0948013, -0.00428717, -0.292519},
					       {0.188874, 0.00288661, 0.231131},
					       {0.0870628, -0.00084847, 0.0152292},
					       {0.151972, 0.0024911, 0.22311},
					       {0.0279692, -3.37325e-05, 0.0066156},
					       {0.381802, -0.000922202, 0.289651},
					       {-0.0991423, -0.00246185, -0.205385},
					       {-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.180112, -0.00243994, 0.131192},
					       {-0.132811, -0.00250005, -0.0974476},
					       {0.0542204, 0.00156509, 0.185888},
					       {-0.0353579, 0.000561183, -0.00341702},
					       {-0.0733351, 0.00727003, 0.17833},
					       {0.0704254, 0.00919455, 0.376577},
					       {-0.179999, -0.00349625, -0.190002},
					       {-0.00344506, 0.000485753, 0.105964},
					       {-0.101258, 0.0057773, 0.155988},
					       {-0.125168, 0.00152611, -0.078966},
					       {0, 0, 0},
					       {0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.145967, -0.0068761, -0.352897},
					       {0.0360971, 0.00111922, 0.203397},
					       {-0.0478782, 0.00248203, 0.0170628},
					       {0.0456846, -0.00123398, 0.0217968},
					       {0.0961041, 0.00257188, 0.143088},
					       {-0.180897, 0.00440966, -0.00453133},
					       {0.00100274, 0.00364556, 0.119555},
					       {0.0433757, -0.00058388, 0.0997984},
					       {-0.191412, 0.00834395, 0.034687},
					       {-0.159748, 0.00510232, -0.0145886},
					       {0, 0, 0},
					       {0, 0, 0}};

*/

/*

double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.741326, -0.00271395, 0.658527},
{0.519379, 0.00175251, 0.57743},
{0.574443, 0.0034494, 0.674904},
{0.424608, 0.00187362, 0.46174},
{0.0993575, 0.00611151, 0.296818},
{0.0476212, -0.000817665, 0.0284958},
{-0.0725687, -0.00438123, -0.20636},
{-0.0874069, -0.000293351, -0.0627047},
{0.0449238, 0.00142292, 0.0861419},
{-0.0139342, -0.00174619, -0.0288256},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.292372, -0.00131762, -0.322725},
{-0.149787, -0.00812024, -0.281346},
{-0.302695, 0.000228822, -0.401001},
{0.0266158, -0.00174497, -0.0418445},
{-0.10312, -0.00386254, -0.290604},
{0.175426, 0.00224493, 0.202688},
{0.0739978, 0.00248977, 0.104241},
{0.174278, 0.00123536, 0.202617},
{0.0516817, -0.000388137, 0.0148537},
{0.345308, -0.000488581, 0.258937},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.199196, -0.00283387, 0.121683},
{-0.149395, -0.00280634, -0.121285},
{0.0534708, 0.00223099, 0.176893},
{-0.0348492, 0.000162304, -0.0135246},
{-0.0819426, 0.00680706, 0.156598},
{0.0905792, 0.00756859, 0.393321},
{-0.00982351, -0.00372233, -0.0264795},
{-0.0855712, -0.00029862, 0.0233469},
{-0.0736123, 0.00491329, 0.177837},
{-0.220895, 0.000927985, -0.175228},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.133826, -0.00737101, -0.360152},
{0.0578988, 0.00110751, 0.228556},
{-0.0282191, 0.00251906, 0.0405966},
{-0.0261895, 7.5693e-05, -0.0324365},
{0.12075, 0.00239904, 0.174291},
{-0.190161, 0.00449231, -0.0280819},
{0.0270188, 0.00351793, 0.14115},
{-0.0312681, -0.000901634, 0.0136154},
{-0.123457, 0.00835408, 0.0912457},
{-0.143438, 0.00455422, -0.0105844},
{0, 0, 0},
{0, 0, 0}};

*/

/*
23092021 for HV scan in layer 8
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.741326, -0.00271395, 0.658527},
{0.51705, 0.00178695, 0.577579},
{0.573441, 0.00350261, 0.675771},
{0.425079, 0.00199995, 0.463453},
{0.0983987, 0.00612751, 0.295803},
{0.0476212, -0.000817665, 0.0284958},
{-0.07299, -0.00427292, -0.205925},
{-0.0866838, -0.000306032, -0.0612786},
{0.0437987, 0.00143382, 0.0855642},
{-0.0127517, -0.00178029, -0.027483},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.292372, -0.00131762, -0.322725},
{-0.148898, -0.00806973, -0.283258},
{-0.300208, 0.000165204, -0.401538},
{0.0277785, -0.0017601, -0.0395632},
{-0.101956, -0.00384203, -0.290783},
{0.175426, 0.00224493, 0.202688},
{0.074852, 0.00250436, 0.105671},
{0.176443, 0.00118985, 0.200424},
{0.0513201, -0.000384089, 0.0158112},
{0.350542, -0.000510498, 0.264976},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0.199196, -0.00283387, 0.121683},
{-0.151634, -0.00258046, -0.122223},
{0.048711, 0.00202087, 0.182638},
{-0.0342793, 0.000314123, -0.01015},
{-0.083899, 0.00693662, 0.16039},
{0.0905792, 0.00756859, 0.393321},
{-0.0117481, -0.00361361, -0.0220824},
{-0.089055, 0.00022514, 0.0221171},
{-0.0744704, 0.00511372, 0.182186},
{-0.225478, 0.000959651, -0.184713},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{-0.133826, -0.00737101, -0.360152},
{0.054139, 0.00153737, 0.232164},
{-0.0346353, 0.00282173, 0.0432851},
{-0.0177714, -0.000747274, -0.0297125},
{0.117961, 0.00286829, 0.174233},
{-0.190161, 0.00449231, -0.0280819},
{0.0254732, 0.00362664, 0.143336},
{-0.0332224, -0.000827094, 0.0152752},
{-0.122998, 0.00840213, 0.0946516},
{-0.144547, 0.00452604, -0.00998103},
{0, 0, 0},
{0, 0, 0}};
*/
#endif

#ifdef ONEMBYONEM
const char* labels[nstrip]={"S00", "S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10",
			    "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20",
			    "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", "S30",
			    "S31"};
#else
const char* labels[nstrip]={"S00", "S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10",
			    "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20",
			    "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", "S30",
			    "S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", "S40",
			    "S41", "S42", "S43", "S44", "S45", "S46", "S47", "S48", "S49", "S50",
			    "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59", "S40",
			    "S61", "S62", "S63"};

#endif
const char* xlabels[2*nlayer]={"X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11",
			       "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10", "Y11"};



//#ifdef ONEMBYONEM
//const double layerzpos[nlayer]={0.0, 16.25, 32.05, 48.3, 64.2, 80.3, 96.4, 112.1, 128.2, 144.4, 160.3, 176.05};
//#else
//const double layerzpos[nlayer]={0.0, 17.80,34.40 ,51.10 , 68.10, 85.10,101.50 , 118.30,135.30 ,151.9 ,169.0, 185.50};
const double layerzpos[nlayer] = {0.0,10.4,20.2,30.2,40.7,51.0/*101.9*/,61.1,71.3,81.3,91.9,101.9,111.9};//9.6,19.2,28.8,38.4,48.0,57.6,67.2,76.8,86.4,96.0,105.6};   //jim jim
//#endif
const double pival=acos(-1);
const int nmxhits=4;
const int nmxusedhits=3;

const int nmxtimehit=4;
const int nmxusedtimehit=3;

//from
/*

double xposerrsq[nmxhits][nlayer] = {{0.081097, 0.0649439, 0.0548458, 0.0446785, 0.0382087, 0.0567823, 0.0405045, 0.0393467, 0.0440577, 0.0606275, 0.0499431, 0.0706403},
{0.0519754, 0.0576385, 0.0408247, 0.059 8575, 0.0482522, 0.0378107, 0.067654, 0.065411, 0.0464823, 0.0571661, 0.04932, 0.0645767},
{0.163662, 0.11067, 0.0814183, 0.0855799, 0.0813797, 0.0982121, 0.0867583, 0.0765227, 0.090884, 0.0989566, 0.0741876, 0.0968527},
{0.95571, 0.781556, 0.69278, 0.679704, 0.74384, 0.844302, 0.461556, 0.598369, 0.838896, 0.800568, 0.642708, 0.855978}};
double yposerrsq[nmxhits][nlayer] = {{0.058812, 0.0461559, 0.0468089, 0.0557433, 0.0430914, 0.0417022, 0.0545324, 0.0358047, 0.0440019, 0.0468247, 0.0413938, 0.0529024},
{0.0459345, 0.0331053, 0.0553509, 0.0482619, 0.0315608, 0.0294313, 0.0549353, 0.0402538, 0.0248514, 0.0291569, 0.0389157, 0.040075},
{0.102929, 0.115662, 0.0836209, 0.0737237, 0.117718, 0.106202, 0.078392, 0.0847044, 0.101968, 0.114648, 0.0813697, 0.102919},
{0.659345, 0.802475, 0.634551, 0.622705, 0.745852, 0.722659, 0.66842, 0.683721, 0.819542, 0.754843, 0.724981, 0.901681}};
*/
//double timeserrx2[nlayer] = {0.732514, 0.25, 0.913301, 0.948905, 0.986499, 0.921516, 0.972974, 0.25, 0.460946, 0.601381, 0.348236, 0.41084};

//double timeserry2[nlayer] = {0.862293, 0.347105, 1.41956, 1.59067, 0.25, 1.58684, 1.36139, 0.32777, 0.434341, 1.57892, 0.25, 0.550872};


Double_t gausX(Double_t* x, Double_t* par){
  return par[0]*(TMath::Gaus(x[0], par[1], par[2], kFALSE)); //kTRUE));
}

const int nmxchn = 200;
int nchannel =0;
double m_data[nmxchn];
double m_xpos[nmxchn];

void fcnsg(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t flag) {

  double fval=0;
  double x[2];
  for (int ij=0; ij<nchannel; ij++) {
    x[0] = m_xpos[ij];
    fval += pow( (m_data[ij] - gausX(x, par)), 2.) / TMath::Max(1., m_data[ij]);
  }

  f = fval;
}


Double_t fitspec(Double_t* x, Double_t* par) {
  //  int nx=int (x[0]/nstrip);
  double yy = x[0]; // - nx*nstrip;

  double yval = par[0] + par[1] + par[2]*sin(pival*yy/nusedstrip);
  //  double yval = par[2]*pow(sin(pival*yy/32), par[3]);
  //  double yval = par[2]*sin(pival*yy/32);
  return yval;
}

//bias in time resolution due to uncertainties in other layers
double bias_intime_xreso2[nlayer]={0.404097, 0.294606, 0.250296, 0.168577, 0.128642, 0.114688, 0.120965, 0.138543, 0.158922, 0.251805, 0.316103, 0.46363};
double bias_intime_yreso2[nlayer]={0.508908, 0.390572, 0.342497, 0.227587, 0.167038, 0.142859, 0.140968, 0.173794, 0.230109, 0.359566, 0.454661, 0.635356};

double bias_inpos_xreso2[nlayer]={0.0240879, 0.0176354, 0.0150597, 0.0098424, 0.00719886, 0.00583462, 0.00570846, 0.00690567, 0.00941585, 0.0144282, 0.0174868, 0.0240973};
double bias_inpos_yreso2[nlayer]={0.0242652, 0.0176931, 0.0146815, 0.0100788, 0.00727195, 0.00586189, 0.00567934, 0.00681836, 0.00907941, 0.0137591, 0.0167822, 0.0235308};

//Time shift correction based on postion of position of track in a strip
double parabolic(double x, double* par) {
  double yy = par[0]+par[1]*abs(x)+par[2]*x*x;
  //  cout<<" yyy "<< par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<x <<" "<<fabs(x)<<" "<<abs(x)<<" "<<yy<<endl;
  return yy;
}

double fitfunctot(double x, double* par) {
  double funcval;
  if(x>0.&& x<=21.) {
     return par[0]*pow(cos(par[1]*x),par[2]);
  } else if(x>21.) {
     return par[3]*TMath::Gaus(x,par[4],par[5]) + par[6]*TMath::Gaus(x,par[7],par[8]);
  }
}

double fitfunctot2(double x, double* par) {
    int arrindx;
    if(x>=0. && x<20) {
       arrindx = 0;
       return  par[arrindx]+par[arrindx+1]*pow(x,1.)+par[arrindx+2]*pow(x,2.)+par[arrindx+2]*pow(x,3.)+par[arrindx+2]*pow(x,4.);
    } else if(x>=20 && x<34.) {
      arrindx = 5;
       return  par[arrindx]+par[arrindx+1]*pow(x,1.)+par[arrindx+2]*pow(x,2.)+par[arrindx+3]*pow(x,3.)+par[arrindx+4]*pow(x,4.);
    }else if(x>=34 && x<38.) {
      arrindx = 10;
      return  par[arrindx]+par[arrindx+1]*pow(x,1.)+par[arrindx+2]*pow(x,2.)+par[arrindx+3]*pow(x,3.)+par[arrindx+4]*pow(x,4.);
    }else if(x>=38) {
      arrindx = 15;
      return  par[arrindx]+par[arrindx+1]*pow(x,1.)+par[arrindx+2]*pow(x,2.)+par[arrindx+3]*pow(x,3.)+par[arrindx+4]*pow(x,4.);
    }
}
//double strpos_vs_time[2*nmxtimehit][nlayer][3] ={0};
//from _jj50
/*



double strpos_vs_time[2*nmxtimehit][nlayer][3] ={
{{-0.123113, 0.18606, 0.813269},
{-0.120944, 0.103326, 0.760101},
{-0.0970921, 0.110158, 1.29783},
{-0.0531101, 0.104501, 1.51129},
{-0.137217, -0.00784696, 1.89188},
{-0.127418, -0.0632092, 1.83653},
{-0.159292, -0.00146545, 1.80509},
{-0.181328, -0.23882, 2.13574},
{-0.161664, -0.0574664, 1.90206},
{-0.209488, 0.0931923, 1.27755},
{-0.216192, 0.222433, 0.83137},
{-0.261329, 0.272662, 0.460641}},

{{-0.24998, 0.606532, 2.13545},
{-0.228501, 0.24806, 3.06557},
{-0.309287, -0.111194, 4.03333},
{-0.171296, -0.551247, 4.20086},
{-0.34884, -0.400734, 4.75968},
{-0.343191, -0.32021, 4.64909},
{-0.399293, -0.162118, 4.54819},
{-0.458792, 0.409844, 3.5494},
{-0.373437, -0.183276, 4.31114},
{-0.558112, 1.68538, 1.25675},
{-0.362985, 1.05872, 1.45137},
{-0.196687, 0.971843, 1.58006}},

{{-1.0718, -1.15025, 5.68805},
{-0.885998, -1.3294, 5.46625},
{-0.74276, -1.23182, 4.74245},
{-0.491446, -1.03858, 4.30676},
{-1.06557, -0.967623, 4.45857},
{-0.935407, -1.47415, 6.24446},
{-1.0148, -1.39736, 5.54038},
{-1.30535, -1.03262, 5.09559},
{-0.955787, -0.737653, 3.78017},
{-1.27183, -0.888609, 4.85362},
{-1.564, -0.368864, 4.1142},
{-1.10468, 0.231317, 4.08977}},

{{-0.755741, 0.250088, 0.563088},
{-0.302553, 0.0579232, 1.1176},
{-0.373607, 0.40424, 0.579971},
{-0.450426, 0.34948, 0.25685},
{-0.598923, -0.371326, 1.80526},
{-0.679298, 0.0239827, 0.96946},
{-0.735798, 0.49941, 0.954074},
{-1.20538, 1.53639, 0.198066},
{-0.82023, 0.262081, 1.77157},
{-1.09687, 2.4573, -2.53943},
{-0.931619, -0.29405, 4.01591},
{-0.928323, -0.265313, 2.57414}}};


*/

//from jj12 ?
/*



double strpos_vs_time[2*nmxtimehit][nlayer][3] = { // 0:1x, 1:2x, 3:3x, 4:4x, 5:1y, 6:2y, 7:3y & 8:4y
{{-0.0674297, 0.170439, 0.409651},
{-0.0611661, 0.0843116, 0.36818},
{-0.0608992, 0.0720272, 1.43811},
{-0.029488, 0.0230799, 1.55109},
{-0.0852735, -0.0772598, 2.04451},
{-0.0777609, -0.037026, 1.67501},
{-0.102506, -0.00912536, 1.81528},
{-0.107011, -0.188384, 2.03677},
{-0.075814, 0.00978978, 1.6483},
{-0.117701, 0.150895, 1.03712},
{-0.0937163, 0.208578, 0.657211},
{-0.0600049, 0.276257, 0.350935}},

{{-0.133698, 0.641704, -0.425988},
{-0.0289257, 0.783015, -0.632344},
{-0.264497, 1.28543, -0.401567},
{-0.145635, 1.00634, -0.319277},
{-0.369328, 1.34384, 0.144246},
{-0.360744, 1.73536, -0.870279},
{-0.380285, 1.82961, -0.986573},
{-0.440047, 1.88143, -0.55643},
{-0.297206, 1.37435, -0.272649},
{-0.444448, 2.7363, -2.40169},
{-0.279145, 1.79739, -1.4482},
{-0.130919, 1.3315, -0.949211}},

{{-1.08337, -0.0554113, 0.222879},
{-0.800603, -0.0311802, 0.116837},
{-0.666544, 0.0290114, 0.0521687},
{-0.45883, 0.0974538, 0.248945},
{-0.910158, 0.0556651, -0.101543},
{-0.76806, -0.219106, 1.33766},
{-0.862272, -0.275776, 0.775626},
{-1.19479, 0.0320598, 0.544137},
{-0.803924, 0.300296, -0.411876},
{-1.11993, 0.43575, -0.25818},
{-1.44454, 0.24637, 0.299847},
{-1.01937, 0.980263, -0.31354}},

{{-1.0757, 0.132411, -0.451637},
{-0.617441, -0.0612786, 0.367728},
{-0.739528, 1.30432, -2.11228},
{-0.283779, -0.181128, -0.751285},
{-0.827418, -0.418606, 0.417702},
{-0.842352, 0.996867, -2.20627},
{-0.871919, -0.0912848, 0.596201},
{-1.2719, 0.817326, 0.114753},
{-0.982607, 1.48201, -1.43956},
{-1.186, 0.636628, 0.46263},
{-1.61993, 1.47871, 0.918357},
{-0.73131, 0.668055, -0.529129}},

{{-0.0797712, 0.24008, 0.345694},
{-0.00939127, 0.0561136, 0.560917},
{-0.0531844, 0.199102, 1.52665},
{-0.00633464, -0.108056, 1.63114},
{-0.104791, 0.348658, 1.63834},
{-0.106914, 0.306611, 1.58464},
{-0.12487, 0.23325, 2.08458},
{-0.142368, 0.430819, 1.45945},
{-0.0903414, 0.381746, 1.02316},
{-0.099931, 0.337873, 1.4436},
{-0.0779324, 0.317406, 0.950911},
{0.0150185, 0.42247, 0.510248}},

{{-0.114299, 0.66488, -0.464315},
{-0.0861321, 0.632696, -0.4353},
{-0.253552, 1.11586, -0.257322},
{-0.13858, 1.057, -0.743504},
{-0.427893, 1.59802, -0.289627},
{-0.360461, 1.54384, -0.555168},
{-0.462724, 1.70897, -0.43152},
{-0.395999, 1.53944, -0.576794},
{-0.396334, 1.89273, -1.16307},
{-0.513241, 2.21021, -1.48615},
{-0.268089, 1.36834, -1.00644},
{-0.270102, 1.18027, -0.601756}},

{{-0.920801, -0.180719, 0.485889},
{-0.735903, 0.296633, -0.622558},
{-0.737952, 0.201275, 0.048042},
{-0.376173, 0.0944521, 0.131297},
{-1.03441, 0.096647, 1.02852},
{-0.808437, 0.107032, 0.323696},
{-0.882437, -0.0465669, 0.226607},
{-1.0028, -0.511033, 1.26195},
{-0.891362, -0.0159531, 0.54913},
{-1.36, 0.288394, 0.444698},
{-1.38049, 0.280667, -0.0720997},
{-0.922892, 0.233491, 0.192945}},

{{-1.04844, -0.782402, 1.6579},
{-0.702953, -0.639421, 0.825105},
{-0.631701, -0.408758, 0.945416},
{-0.435467, 0.358731, -1.13661},
{-0.862956, 0.264863, -0.510273},
{-0.872514, 0.207667, 0.157545},
{-0.951074, -0.282902, 0.735075},
{-1.2476, 0.114487, 0.710311},
{-1.11347, 0.955351, -0.0895417},
{-1.31315, -0.00438113, 0.820052},
{-1.44422, -0.647576, 1.99953},
{-1.2383, -0.0450062, 0.930323}}};

*/

//from miical_jj01
 /*
double strpos_vs_time[2*nmxtimehit][nlayer][3] = {
  {{-0.0448697, 0.221767, 0.226477},
   {-0.0722215, 0.174617, 0.74613},
   {-0.0750741, 0.0865269, 1.34962},
   {-0.123104, 0.181259, 0.925995},
   {-0.110588, 0.0382029, 1.59096},
   {-0.117586, -0.0279002, 1.56479},
   {-0.112573, -0.135599, 1.92539},
   {-0.11571, -0.0581722, 1.63511},
   {-0.118406, 0.300869, 0.687014},
   {-0.0779043, 0.0602137, 0.351175},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.133009, 0.680239, -0.416976},
   {-0.316719, 1.18976, -0.384429},
   {-0.620234, 1.82779, -0.551679},
   {-0.243464, 1.41737, -0.682817},
   {-0.340312, 1.18567, 0.22009},
   {-0.356822, 1.44935, -0.124987},
   {-0.271557, 0.980932, 0.425623},
   {-0.669687, 2.97021, -1.95365},
   {-0.297411, 1.14224, -0.387868},
   {-0.0326468, 0.458599, -0.0780037},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-1.53303, 0.31033, -0.420305},
   {-1.70447, 0.484875, -0.395405},
   {-1.72451, 0.282975, 2.74359},
   {-1.39227, 0.171867, 0.404542},
   {-1.49777, 0.363632, 0.183081},
   {-1.63348, 0.671309, 0.0261612},
   {-0.920801, 0.469904, 2.41209},
   {-1.47359, -0.0900912, 1.00886},
   {-1.39096, 1.41558, -0.139558},
   {-1.20289, 0.128169, -0.414735},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-2.11053, 1.75078, -2.32119},
   {-2.28844, 2.75009, -4.47848},
   {-1.69882, 1.14195, -1.17817},
   {-2.03384, 2.39807, -3.40608},
   {-1.73616, 0.133325, 0.487029},
   {-1.81686, 0.650243, 0.647758},
   {-1.19721, 0.0612355, 0.161366},
   {-1.76036, 0.640521, -0.730455},
   {-1.45764, -1.37109, 2.56396},
   {-1.46249, -0.887245, 1.78768},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.0183923, 0.296538, 0.0880122},
   {-0.0341522, 0.273962, 0.914879},
   {-0.192201, 0.16818, 1.43918},
   {-0.096375, 0.2471, 0.727869},
   {-0.0864231, 0.244726, 1.20165},
   {-0.102493, 0.275905, 1.31913},
   {-0.122539, 0.306003, 1.77444},
   {-0.117719, 0.26303, 1.92134},
   {-0.113409, 0.401621, 1.11629},
   {-0.0563375, 0.307794, 0.242431},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.21022, 0.593959, -0.192597},
   {-0.564823, 1.45015, -0.506134},
   {-0.119329, 1.09772, -0.481995},
   {-0.32993, 1.50685, -0.995183},
   {-0.468978, 1.61526, -0.639436},
   {-0.465197, 1.66152, -0.595886},
   {-0.419754, 1.37808, -0.0438821},
   {-0.766782, 2.68497, -1.57906},
   {-0.670304, 2.42984, -1.81317},
   {-0.304649, 0.979221, -0.616283},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-1.61049, 0.0770846, -0.0323675},
   {-2.00837, 0.644391, 0.178477},
   {-0.814113, 0.507532, 0.496697},
   {-1.71484, 0.347924, 0.516586},
   {-1.59565, 0.0264247, 0.946618},
   {-1.66754, 0.451234, 0.410325},
   {-1.07325, -0.110207, 0.520546},
   {-1.67006, 0.398388, 0.375434},
   {-1.7198, 0.58365, 0.217978},
   {-1.54951, 0.300994, -0.163665},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-1.78673, -0.772979, 1.76565},
   {-2.32276, 1.50929, -1.75524},
   {-1.38737, 0.064522, 2.46021},
   {-1.85209, 0.926195, -1.16473},
   {-1.97952, 1.09447, -0.906532},
   {-1.91389, 0.529291, -0.274773},
   {-1.5401, 0.610972, -0.696445},
   {-2.02348, 1.5105, -1.92721},
   {-2.07397, 0.490818, 0.498676},
   {-1.64799, -0.239083, 0.199888},
   {-0.1, 1, 2},
   {-0.1, 1, 2}}};
*/
// double xtimetot[nlayer][nband][2] = {{{1.76533,-0.552454},{0.571606,-0.302464},{-1.95791,0.274374},{-3.41057,0.661373}},{{1.23059,-0.655436},{1.22125,-0.28678},{0.202306,0.249155},{-2.54329,0.757466}},{{0,0},{0.,0.},{0.0},{0.,0.}},{{0.746421,-0.479509},{1.47341,-0.405422},{0.382099,0.0154258},{1.32067,-0.317953},},{{1.50034,-0.599209},{1.08593,-0.36698},{-0.0966534,0.0788617},{-1.81049,-0.66188},},{{1.27059,-0.655707},{1.01384,-0.37474},{-2.42051,-1.12939},{-4.81466,-0.535106}},{{2.56887,-0.598093},{1.84521,-0.167284},{-2.55686,-0.6868},{-1.72675,-0.0141445}},{{2.24691,-0.629948},{1.77863,-0.0601707},{-3.02019,-0.409441},{-3.80217,0.303567}},{{2.61799,-0.64872},{1.82559,-0.0761361},{-2.87257,-0.459503},{-4.51307,-0.289495}},{{2.09832,-0.638843},{1.41888,-0.0966174},{-2.68732,-0.479283},{-3.85233,-1.59173}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}};
// double ytimetot[nlayer][nband][2] = {{{0.516211,-0.697908},{0.262034,-0.06971},{-3.84632,-0.248941},{-3.70488,0.486211}},{{0.106243,-0.744873},{0.161662,-0.213244},{-2.58703,0.225818},{-2.66455,-0.3829}},{{0,0},{0,0},{0.0,0.0},{0.,0.0}},{{-0.294211,-0.522063},{0.266624,-0.376946},{-0.24816,0.0561476},{-0.609705,-0.37394}},{{0.101809,-0.628463},{0.0288145,-0.142214},{-2.76032,-0.428504},{-3.19056,-0.410711}},{{0.213043,-0.740487},{0.0135635,-0.0465408},{-4.46874,-0.530817},{-3.73315,-0.504672}},{{1.38159,-0.616003},{1.06601,-0.0495229},{-4.18792,-0.182215},{-2.01447,-0.340381}},{{1.24541,-0.627022},{0.469229,-0.0462086},{-4.48306,-0.251014},{-2.31764,-0.231737}},{{0.367391,-0.703517},{0.274061,-0.124717},{-3.96247,-0.356709},{-2.94563,-0.833886}},{{0.542087,-0.699603},{0.367609,0.0642358},{-4.57251,-0.330824},{-3.24307,-0.261537}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}};
const int nband = 4;

/*
double xtimetot[nlayer][nband][2] = {{{0.365397,-0.310688},{-0.26868,-0.161901},{-1.03894,-0.0541665},{-0.785052,0.153307}},{{0.300876,-0.410243},{-0.119811,-0.254944},{-0.759406,-0.182417},{-0.886772,-0.188281}},{{0,0},{0,0},{0,0},{0,0}},{{0.517374,-0.338227},{0.306718,-0.256946},{-0.500103,-0.0269511},{-0.480245,-0.0932887}},{{0.444027,-0.359051},{-0.109477,-0.255522},{-1.18663,-0.214307},{-0.673146,-0.13991}},{{0.228645,-0.366254},{-0.292977,-0.26158},{-1.43846,-0.00876863},{-0.00517089,0.366287}},{{0.763131,-0.365001},{0.0474142,-0.199107},{-1.08331,-0.0218876},{-0.592935,0.0888662}},{{0.678863,-0.415359},{-0.218354,-0.136743},{-1.10169,0.122598},{0.0711608,0.31746}},{{0.823585,-0.406179},{-0.132634,-0.136448},{-1.00343,0.052781},{-0.54974,-0.101239}},{{0.539525,-0.377016},{-0.204816,-0.160904},{-1.08158,-0.223523},{-1.64717,-0.404158}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}};

double ytimetot[nlayer][nband][2] = {{{0.289201,-0.345363},{-0.429917,-0.0652669},{-0.703001,-0.0391581},{-0.488214,0.485395}},{{0.194759,-0.376226},{-0.485959,-0.156851},{-0.778061,-0.068465},{-0.68208,-0.0509229}},{{0,0},{0,0},{0,0},{0,0}},{{0.599287,-0.311039},{-0.0249163,-0.221486},{-0.561424,-0.0151613},{-0.590085,-0.137404}},{{0.219232,-0.328065},{-0.408118,-0.0933984},{-0.721495,-0.140636},{-0.629362,0.0112672}},{{0.13649,-0.36608},{-0.741674,-0.0967892},{-1.26243,-0.0256001},{-1.97611,-0.392186}},{{0.659622,-0.359134},{-0.247194,-0.048673},{-0.548854,0.0304604},{-0.475263,-0.0773765}},{{0.522077,-0.335725},{-0.493541,-0.0376901},{-0.896324,-0.00136582},{-0.999284,-0.00446134}},{{0.270561,-0.357937},{-0.428979,-0.100519},{-0.707015,-0.144448},{-1.18937,-0.170653}},{{0.176224,-0.293975},{-0.45881,0.0159427},{-0.823799,0.0255578},{-1.20601,-0.113068}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}};

double xtimetot[nlayer][nband][2] = {{{0.366244,-0.310225},{-0.502852,-0.153031},{-1.06727,-0.0951597},{-0.658589,0.0635426}},{{0.302713,-0.410353},{-0.489515,-0.231509},{-0.926802,-0.0985605},{-0.566086,0.10384}},{{0,0},{0,0},{0,0},{0,0}},{{0.513234,-0.339876},{-0.0462709,-0.240169},{-0.494497,-0.0332639},{-0.60165,-0.190226}},{{0.443948,-0.359068},{-0.355798,-0.176911},{-0.879014,-0.16635},{-0.817229,-0.143853}},{{0.228838,-0.366212},{-0.696687,-0.266555},{-1.4831,-0.0122707},{-1.04602,-0.167509}},{{0.763128,-0.365003},{-0.207145,-0.148067},{-0.868259,-0.0136088},{-0.683817,-0.157035}},{{0.678796,-0.415342},{-0.423094,-0.161398},{-1.15659,0.113201},{-0.487298,-0.290046}},{{0.823246,-0.406084},{-0.330122,-0.129678},{-0.857057,-0.00583076},{-0.651664,-0.0894564}},{{0.540561,-0.377181},{-0.492987,-0.215436},{-1.21145,-0.222705},{-1.83303,-0.156608}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}};

double ytimetot[nlayer][nband][2] = {{{0.28847,-0.345392},{-0.584202,-0.0777139},{-0.865935,-0.0988868},{1.66529,1.46085}},{{0.195081,-0.376722},{-0.768721,-0.176786},{-0.958011,-0.156977},{-0.75665,-0.0948914}},{{0,0},{0,0},{0,0},{0,0}},{{0.597136,-0.310819},{-0.329471,-0.218211},{-0.659956,-0.0345034},{-0.880229,-0.203577}},{{0.219145,-0.328254},{-0.474578,-0.0805863},{-0.742617,-0.151914},{-0.416449,0.138562}},{{0.135968,-0.366271},{-0.897169,-0.118322},{-1.26962,-0.0586333},{-1.48949,-0.0241953}},{{0.659623,-0.359138},{-0.355462,-0.0786227},{-0.624252,-0.110947},{-0.664454,-0.100817}},{{0.52205,-0.335713},{-0.500857,-0.0508215},{-0.913591,-0.0152297},{-1.32065,-0.104231}},{{0.270316,-0.357818},{-0.614426,-0.115816},{-0.996073,-0.211601},{-1.06147,0.159103}},{{0.176591,-0.294069},{-0.493511,-0.0189245},{-0.850462,0.014108},{-1.66724,-0.0869364}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}};
*/
double xtimetot[nlayer][nband][2] = {{{0.198879,-0.331369},{-0.374162,-0.222662},{-0.747121,-0.18361},{-0.652328,0.120102}},{{0.309327,-0.409159},{-0.248061,-0.268369},{-0.49324,-0.299011},{-0.603298,-0.34103}},{{0,0},{0,0},{0,0},{0,0}},{{0.759908,-0.371742},{0.200017,-0.258287},{-0.504526,-0.0272435},{-0.343499,-0.1329}},{{0.448545,-0.369064},{-0.187752,-0.219339},{-0.874387,-0.259469},{-0.369771,-0.0159672}},{{0.14826,-0.359291},{-0.425297,-0.284637},{-1.30153,-0.0832935},{-0.882025,0.842308}},{{0.265746,-0.340314},{-0.283573,-0.207375},{-0.780636,-0.266771},{-0.283926,-0.0668603}},{{0.152231,-0.386965},{-0.396871,-0.203134},{-0.252621,-0.347905},{-0.310807,0.644625}},{{0.688264,-0.381714},{-0.0956937,-0.196624},{-0.745015,-0.0523803},{0.100868,0.0262983}},{{0.305827,-0.373754},{-0.360989,-0.200887},{-0.811521,-0.141851},{-0.751244,0.0410243}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}};

double ytimetot[nlayer][nband][2] = {{{0.0117474,-0.318014},{-0.562433,-0.0777283},{-0.74839,-0.089796},{-1.34379,0.0638467}},{{0.149341,-0.372949},{-0.710717,-0.14695},{-1.03407,-0.186731},{-0.64239,-0.0539065}},{{0,0},{0,0},{0,0},{0,0}},{{0.770103,-0.334476},{-0.403822,-0.230864},{-0.648029,-0.0451244},{-0.979854,-0.212098}},{{0.099382,-0.313936},{-0.57617,-0.109873},{-1.08358,-0.199357},{-1.12528,-0.216484}},{{0.0473206,-0.358023},{-0.89312,-0.101548},{-1.39471,-0.107553},{-1.16716,0.0121764}},{{0.114513,-0.324047},{-0.423507,-0.0882019},{-0.911162,-0.264624},{-0.0532677,0.0540513}},{{0.121935,-0.301594},{-0.683338,-0.100612},{-1.19258,-0.159369},{-0.904079,0.0863255}},{{-0.278587,-0.334505},{-0.821701,-0.139846},{-1.70838,-0.313187},{-0.947613,-0.124709}},{{-0.142663,-0.240049},{-0.456982,-0.00560049},{-0.769313,0.0195419},{-0.914953,0.336414}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}};

 double xt_ystr_gap[nlayer][nband+1][2] = {{{5.,5.},{25.,21.5},{31,26},{45,39},{60,60}},{{5,5},{25,20.5},{31,26},{45,39},{60,60}},{{20,20},{50,40},{60,60},{0,0},{0,0}},{{5,5},{25,20},{30,25},{45,39},{60,60}},{{5,5},{25,20},{31,26},{45,39},{60,60}},{{5,5},{25,20},{31,26},{45,39},{60,60}},{{5,5},{26,21},{32,27},{45,39},{60,60}},{{5,5},{26,21},{33,28},{45,39},{60,60}},{{5,5},{26,21},{33,28},{45,39},{60,60}},{{5,5},{26,21},{32,27},{45,39},{60,60}}};

  double yt_xstr_gap[nlayer][nband+1][2] =  {{{5,5},{25,20.5},{31,26},{43,37},{60,60}},{{5,5},{25,20},{30,25},{43,37},{60,60}},{{20,20},{50,40},{60,60},{0,0},{0,0}},{{5,5},{24,19},{29,24},{43,37},{60,60}},{{5,5},{25,19.5},{30,25},{43,37},{60,60}},{{5,5},{25,20.5},{31,26},{43,37},{60,60}},{{5,5},{26,21},{31,26},{43,37},{60,60}},{{5,5},{26,21},{31,26},{43,37},{60,60}},{{5,5},{25,20},{30,25},{43,37},{60,60}},{{5,5},{26,21},{30,25},{43,37},{60,60}}};
  double range[nband][2] = {{0,35},{25,35},{25,35},{0,35}};
  double xmean[5] = {0,0,25,0,0};
  double xtslopepar[nlayer][nband][2] = {{{16.441,0},{32.6427,-0.291593},{32.5757,0},{44.4174,-0.143605}},{{15.9542,0},{31.7869,-0.268337},{32.7905,0},{46.3708,-0.186282}},{{7.80481e-317,0},{6.95272e-310,0},{0,0},{0,0}},{{15.2415,0},{32.0825,-0.28194},{32.3217,0},{49.7497,-0.269447}},{{15.9589,0},{32.1746,-0.289861},{32.133,0},{45.8541,-0.154933}},{{16.1817,0},{32.6635,-0.296235},{32.8795,0},{45.345,-0.161264}},{{17.3394,0},{33.4887,-0.281652},{33.6761,0},{48.0628,-0.241309}},{{17.2554,0},{33.7785,-0.299409},{33.5443,0},{45.3982,-0.101133}},{{16.2518,0},{33.0899,-0.308926},{32.3748,0},{47.0624,-0.209903}},{{16.5708,0},{33.568,-0.307407},{31.8,0},{44.8725,-0.148043}}};
  double ytslopepar[nlayer][nband][2] = {{{16.7573,0},{34.3203,-0.316507},{33.405,0},{48.8086,-0.233624}},{{16.0926,0},{34.0818,-0.310135},{34.3026,0},{48.2521,-0.183208}},{{3.95253e-323,0},{0,0},{0,0},{0,0}},{{15.3963,0},{33.581,-0.291263},{33.5519,0},{52.6204,-0.312474}},{{16.3255,0},{34.1487,-0.315798},{33.9568,0},{49.7053,-0.242714}},{{16.1881,0},{32.604,-0.262157},{34.1649,0},{46.9406,-0.151591}},{{17.5134,0},{34.8198,-0.29878},{35.1536,0},{51.7693,-0.319671}},{{17.1635,0},{34.5355,-0.285974},{35.754,0},{51.0666,-0.323619}},{{17.5429,0},{35.084,-0.301724},{35.6096,0},{49.2182,-0.22165}},{{17.0539,0},{34.3719,-0.293689},{33.4808,0},{49.0093,-0.203212}}};

double xtmeanall[nlayer][nband] = {{0.1481,-0.5391,-1.187,-1.368},{0.2521,-0.5071,-0.9432,-1.228},{0.0309,-0.4786,0,0},{0.6957,-0.05258,-0.4541,-0.7248},{0.3781,-0.2893,-0.6648,-0.8027},{0.2023,-0.6489,-1.135,-1.535},{0.3062,-0.3311,-0.6902,-0.8851},{0.2641,-0.4935,-0.8641,-1.458},{0.3168,-0.3782,-0.7093,-0.9213},{0.2486,-0.4856,-1.286,-1.598},{0.,0.,0.,0.},{0.,0.,0.,0.}};
double ytmeanall[nlayer][nband] = {{0.09848,-0.6983,-1.093,-1.379},{0.2131,-0.587,-0.8423,-0.7609},{0.09415,-0.4178,0,0},{0.8001,-0.07915,-0.5239,-0.6189},{0.2664,-0.385,-0.6573,-1.04},{0.08318,-0.7917,-1.081,-1.326},{0.2589,-0.4131,-0.759,-1.223},{0.157,-0.5588,-0.9543,-1.245},{0.2359,-0.5133,-0.8163,-1.176},{0.1099,-0.5902,-1.058,-1.115},{0.,0.,0.,0.},{0.,0.,0.,0.}};
const int xnband[nlayer-2] = {4,6,5,4,6,5,6,5,6,4};
 const int npoint = 4;
const int nbandmax = 6;
  double xbandpoints[nlayer-2][nbandmax][npoint][2] = {{{{0,5},{60,5},{0,13},{60,13}},{{0,13},{60,13},{0,18},{60,18}},{{28,18},{34,18},{28,50},{34,50}},{{0,25},{28,22},{0,32},{28,26}},{{100,200},{100,200},{100,200},{100,200}},{{100,200},{100,200},{100,200}}},
						     {{{0,5},{60,5},{0,14},{60,14}},{{0,14},{60,14},{0,18},{60,18}},{{0,25},{46,19},{0,32},{46,22}},{{28,28},{42,28},{28,34},{42,34}},{{46,18},{52,18},{46,50},{52,50}},{{10,38},{40,38},{10,50},{40,50}}},
						     {{{0,24},{60,24},{0,34},{60,34}},{{0,34},{40,34,},{0,40},{40,40}},{{40,34},{60,34},{40,46},{60,40}},{{24,34},{40,40},{24,50},{40,46}},{{25,50},{40,50},{25,60},{40,60}},{{100,200},{100,200},{100,200},{100,200}}},
						     {{{0,5},{60,5},{0,14},{60,14}},{{0,14},{60,14},{0,33},{60,18}},{{6,35},{60,23},{6,35},{60,35}},{{0,42},{60,35},{0,52},{60,35}},{{100,200},{100,200},{100,200},{100,200}},{{100,200},{100,200},{100,200},{100,200}}},
						     {{{0,5},{60,5},{0,14},{60,14}},{{0,13},{40,13},{0,18},{40,18}},{{40,14},{60,14},{40,24},{60,18}},{{0,24},{40,18},{0,32},{40,24}},{{0,32},{50,25},{0,34},{50,34}},{{8,36},{32,36},{8,48},{32,48}}},
						     {{{0,5},{60,5},{0,14},{60,14}},{{0,14},{60,14},{0,18},{60,18}},{{0,25},{32,40},{0,32},{60,18}},{{36,28},{54,28},{36,34},{54,34}},{{6,37},{30,37},{6,46},{30,46}},{{100,200},{100,200},{100,200},{100,200}}},
						     {{{0,5},{60,5},{0,16},{60,16}},{{0,16},{36,16},{0,20},{36,20}},{{36,16},{60,16},{36,26},{60,20}},{{0,26},{36,20},{0,34},{36,26}},{{26,32},{60,24},{26,34},{60,34}},{{0,40},{40,40},{0,50},{40,50}}},
						     {{{0,5},{60,5},{0,15},{60,15}},{{0,15},{60,15},{0,19},{60,19}},{{0,23},{42,22},{0,36},{42,25}},{{0,36},{44,26},{0,36},{44,36}},{{16,38},{34,38},{16,46},{34,46}},{{100,200},{100,200},{100,200},{100,200}}},
						     {{{0,5},{60,5},{0,14},{60,14}},{{0,14},{60,14},{0,18},{60,18}},{{0,18},{33,18},{0,32},{33,23}},{{7,30},{20,28},{7,34},{20,34}},{{25,31},{60,23},{25,34},{60,30}},{{0,38},{40,38},{0,50},{40,50}}},
						     {{{0,5},{60,5},{0,14},{60,14}},{{0,14},{60,14},{0,19},{60,19}},{{0,26},{34,19},{0,34},{34,26}},{{28,28},{35,28},{28,44},{34,44}},{{100,200},{100,200},{100,200},{100,200}},{{100,200},{100,200},{100,200},{100,200}}}};

double xmeancorr[nlayer][nbandmax+1] = {{0.931244,0.103037,-0.307279,-0.378884,0,0,-0.566627},{0.827963,0.0993,-0.404698,-0.621681,-0.32046,-0.586387,-0.52164},{0.612368,-0.0755286,-0.107195,-0.357222,-0.791561,0,-0.600595},{0.955053,0.0477133,-0.436503,-0.506457,0,0,-0.271568},{0.76658,0.109549,-0.0683385,-0.244929,-0.499176,-0.352549,-0.392241},{0.888299,0.0881195,-0.701871,-0.702718,-0.981107,0,-0.556023},{0.728268,0.0489261,-0.0925014,-0.254262,-0.429948,-0.497811,-0.493082},{0.893734,0.118455,-0.433483,-0.764769,-0.789434,0,-0.380757},{0.920741,0.159386,-0.249202,-0.375583,-0.57254,-0.600814,-0.444042},{0.713651,0.0708401,-0.397161,-0.552502,0,0,-0.482355},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};

double ymeancorr[nlayer][nbandmax+1] = {{1.6827,0.179803,-0.38733,-0.621067,0,0,-0.725449},{1.36715,0.139773,-0.548837,-0.669828,-0.475832,-0.576993,-0.631161},{1.37468,-0.113416,-0.149322,-0.523872,-1.11956,0,-0.481297},{1.60887,0.0922841,-0.605919,-0.570078,0,0,-0.385238},{1.32067,0.183681,-0.086599,-0.354886,-0.524024,-0.509753,-0.488306},{1.43181,0.095885,-0.719635,-0.892848,-1.10951,0,-0.637609},{1.12857,0.0651868,-0.115265,-0.348407,-0.548684,-0.538714,-0.544087},{1.31028,0.142514,-0.496405,-0.920515,-0.949592,0,-0.485883},{1.57617,0.220375,-0.375749,-0.491975,-0.63029,-0.745514,-0.56238},{1.2839,0.0712272,-0.47084,-0.655361,0,0,-0.61156},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};



//20190610 fit using fitfunctot()
double xtcortotpar[nlayer][9]= {
  {17.0228,0.0509286,0.252653,14.2893,26.672,19.5558,5.33896,48.1513,5.78161},
  {18.1004,0.0452115,0.534634,13.9392,20.0666,23.0193,8.84306,55.9982,13.2114},
  {15,0.0174533,2,67.4348,-10.3135,14.8257,39.885,-495.361,370.224},
  {18.1502,0.0186777,2.96323,14.8485,24.7274,24.5304,6.0625,54.194,9.35027},
  {18.3969,0.0114734,9.18753,13.4616,19.1523,16.0811,12.3791,53.1952,14.8191},
  {18.1645,0.0143197,6.19901,38.7845,-150.095,115.767,7.97104,75.7441,29.0994},
  {17.8531,0.0562872,0.233079,14.415,24.6366,22.354,7.62142,54.8335,10.2138},
  {18.3095,0.0252605,1.77234,14.4808,27.0696,29.5511,3.40569,50.0216,5.32278},
  {18.1116,0.0391428,0.618671,14.4113,26.6726,40.8137,2.0829,49.7607,4.8693},
  {18.005,0.0457898,0.456457,14.4406,26.0558,20.9799,5.18411,50.3587,6.97156},
  {15,0.0174533,2,10,25,25,10,40,25},
  {15,0.0174533,2,10,25,25,10,40,25}};
  double ytcortotpar[nlayer][9]= {
    {17.9089,0.0485677,0.467184,1.96827,23.9787,4.58935,13.3813,37.4347,31.7976},
    {17.8583,0.0505397,0.413165,14.0673,24.3307,16.4164,9.25621,49.5231,8.48608},
    {15,0.0174533,2,4.14811,28.387,6.34851,13.3208,166.258,9678.59},
    {18.1941,0.0161499,4.19061,14.0982,22.4283,13.9552,11.9991,51.0598,11.7996},
    {17.9504,0.025599,1.73167,2.42381,22.0601,6.40228,14.2629,45.0879,37.4907},
    {18.158,0.00793206,21.2196,14.1258,9.8455,125.844,62.9597,72.5105,7.3392},
    {17.5199,0.0638656,0.154293,14.5126,26.2956,18.3263,7.20928,48.9875,6.81813},
    {17.5796,0.0578023,0.224409,8.82702,15.2484,16.394,12.7051,54.1325,26.487},
    {18.2638,0.0339195,1.06434,3.73408,20.4287,6.91681,13.8495,41.3838,25.9137},
    {17.5676,0.0465048,0.416859,11.7766,10.7977,17.3114,11.1733,44.8292,16.6698},
    {15,0.0174533,2,10,25,25,10,40,25},
    {15,0.0174533,2,10,25,25,10,40,25}};

//fit using fitfunctot2()
double xtcortotpar2[nlayer][20]= {
{17.1333,-0.0566824,-1.02641e-05,-143645,0.000405047,123.715,-15.5812,0.818858,-0.0188413,0.000159989,34.2135,1.8985,-0.218658,0.00603071,-5.19687e-05,-1375.02,125.513,-4.26788,0.0647309,-0.000369593},
{18.0417,-0.0720737,-1.84228e-05,-93127.3,0.000405047,95.0769,-13.2455,0.803926,-0.0214337,0.000211639,63.9269,1.52272,-0.240815,0.0057085,-3.54779e-05,-1371.71,126.011,-4.27892,0.0642958,-0.000360694},
{-43710.2,-2.80099,0.368543,-0.0207454,0.000405047,149.55,-15.912,0.794373,-0.0192529,0.000183515,35.4957,1.94419,-0.217749,0.00602892,-5.32142e-05,-1344.29,125.116,-4.28278,0.0646046,-0.000362624},
{18.557,-0.121168,-1.04907e-05,-52880,0.000405047,109.051,-14.3776,0.81297,-0.020179,0.000185401,-11.4523,2.58883,-0.18222,0.00651049,-8.03922e-05,-1369.66,125.875,-4.27807,0.0643831,-0.000361983},
{18.567,-0.107685,-1.4607e-05,-67057.4,0.000405047,99.5557,-13.6125,0.808176,-0.021058,0.000202968,93.9423,1.08656,-0.264241,0.00539791,-1.71812e-05,-1373.91,126.138,-4.27938,0.064244,-0.000360096},
{18.0384,-0.0791739,-1.75604e-05,-111924,0.000405047,86.638,-12.5897,0.798597,-0.0219692,0.000221429,39.6216,1.90111,-0.223958,0.00588351,-4.69936e-05,-1518.41,132.581,-4.27929,0.0610696,-0.000325122},
{17.4871,-0.011154,-1.79115e-05,-61508.8,0.000405047,100.067,-13.6912,0.809903,-0.0210034,0.00020148,-16.0769,2.65303,-0.178668,0.0065585,-8.3189e-05,-1374.51,126.047,-4.27791,0.06432,-0.000361357},
{17.6357,-0.00535185,-2.19833e-05,-75000.9,0.000405047,95.1778,-13.3774,0.810084,-0.0213258,0.000206439,15.8498,2.20189,-0.203263,0.00622568,-6.43872e-05,-1368.45,125.449,-4.27204,0.064672,-0.000366988},
{17.2789,0.0222844,-2.16557e-05,-95916.9,0.000405047,107.623,-14.3244,0.813127,-0.0202269,0.000186158,-22.568,2.73487,-0.173879,0.00662854,-8.66179e-05,-1375.98,125.923,-4.27469,0.0644174,-0.000363385},
{17.848,-0.0548057,-1.64326e-05,-96228.8,0.000405047,122.535,-15.5423,0.820242,-0.0188378,0.00015866,35.2986,1.89877,-0.219067,0.00601941,-5.19629e-05,-1362.97,125.346,-4.27391,0.064682,-0.000366713},
{24.321,-2.80099,0.368543,-0.0207454,0.000405047,123.87,-15.5838,0.818598,-0.0188464,0.00016036,34.1937,1.89521,-0.21874,0.00603115,-5.1813e-05,-1373.78,125.499,-4.26844,0.0647254,-0.000369364},
{24.321,-2.80099,0.368543,-0.0207454,0.000405047,123.87,-15.5838,0.818598,-0.0188464,0.00016036,34.1937,1.89521,-0.21874,0.00603115,-5.1813e-05,-1373.78,125.499,-4.26844,0.0647254,-0.000369364}};

double ytcortotpar2[nlayer][20]= {
{17.7733,-0.0633388,-1.93783e-05,-132931,0.000405047,94.4132,-13.5367,0.819941,-0.0212732,0.000200109,28.2708,1.87954,-0.214024,0.00617286,-5.55342e-05,-1634.69,138.007,-4.29488,0.0588901,-0.000300418},
{17.8728,-0.074639,-1.8423e-05,-124739,0.000405047,101.422,-13.9173,0.813197,-0.0206642,0.000192761,118.895,0.676628,-0.28501,0.00513673,-3.84894e-07,-1354.56,124.87,-4.27153,0.0649496,-0.000370441},
{-75507.6,-2.80099,0.368543,-0.0207454,0.000405047,101.213,-15.0859,0.845876,-0.0184743,0.000131394,37.8637,1.92553,-0.219166,0.00600414,-5.24411e-05,-1342.17,125.123,-4.28419,0.064573,-0.00036186},
{18.5697,-0.122641,-1.12138e-05,-73922.1,0.000405047,122.729,-15.5285,0.82055,-0.0188405,0.000158126,35.1753,1.9045,-0.218809,0.00602503,-5.19636e-05,-1358.57,125.297,-4.27598,0.0646641,-0.000365622},
{18.2133,-0.111082,-1.33109e-05,-173772,0.000405047,104.926,-14.1124,0.808147,-0.0201566,0.000184863,30.2554,1.91456,-0.21623,0.00608766,-5.30765e-05,-1570.36,135.043,-4.29007,0.0602309,-0.000315436},
{18.0273,-0.0829786,-1.81094e-05,-52368.1,0.000405047,102.685,-13.9861,0.812455,-0.0205984,0.000192398,0.163574,2.34453,-0.192862,0.00640175,-7.12168e-05,-1358.12,125.648,-4.28213,0.0643462,-0.000359628},
{17.2607,-0.0102866,-1.70836e-05,-103223,0.000405047,96.3385,-13.3947,0.807078,-0.0212593,0.000206631,35.0181,1.98195,-0.220132,0.00592786,-5.0011e-05,-1221.42,118.708,-4.26427,0.0678702,-0.000403719},
{17.2034,-0.00934248,-1.84231e-05,-60489.4,0.000405047,97.9028,-13.6177,0.811969,-0.0210436,0.000200358,66.8597,1.43759,-0.244026,0.00568649,-3.21442e-05,-1364.7,125.903,-4.28196,0.0642972,-0.000359841},
{18.3742,-0.102958,-1.67084e-05,-151819,0.000405047,108.944,-14.4724,0.814671,-0.0200255,0.000181333,122.236,0.62233,-0.287844,0.00510143,2.1248e-06,-1367.72,125.91,-4.27913,0.0643209,-0.000361085},
{17.8115,-0.0914744,-1.20312e-05,-191148,0.000405047,139.334,-16.8103,0.820315,-0.0171828,0.000129301,24.5602,2.08568,-0.212645,0.00605527,-5.54133e-05,-1620.89,137.29,-4.28951,0.0591472,-0.000304307},
{24.321,-2.80099,0.368543,-0.0207454,0.000405047,123.87,-15.5838,0.818598,-0.0188464,0.00016036,34.1937,1.89521,-0.21874,0.00603115,-5.1813e-05,-1373.78,125.499,-4.26844,0.0647254,-0.000369364},
{24.321,-2.80099,0.368543,-0.0207454,0.000405047,123.87,-15.5838,0.818598,-0.0188464,0.00016036,34.1937,1.89521,-0.21874,0.00603115,-5.1813e-05,-1373.78,125.499,-4.26844,0.0647254,-0.000369364}};


//from 20181226 //jim see this for 2018 data
double strpos_vs_time[2*nmxtimehit][nlayer][3]= {
{{0.0301133, 0.173192, 0.180998},
{-0.00952232, 0.214831, 0.647275},
{-0.0839172, 0.117452, 1.44721},
{-0.0169559, 0.159143, 0.820714},
{-0.0279889, 0.103518, 1.15376},
{-0.0680663, 0.0355852, 1.38724},
{-0.0505348, 0.0220958, 1.52765},
{-0.0366391, -0.102215, 1.55193},
{-0.030703, 0.0825523, 0.932763},
{-0.019575, 0.113086, 0.191272},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.223559, 0.594342, -0.0615138},
{-0.32414, 1.10071, -0.472413},
{-0.703874, 2.4027, -0.924179},
{-0.33013, 1.20633, -0.306647},
{-0.30099, 0.888727, 0.504324},
{-0.510944, 1.43281, 0.191222},
{-0.376526, 1.01101, 0.589447},
{-0.71818, 2.66586, -1.53398},
{-0.333247, 1.20845, -0.521417},
{-0.21402, 0.698503, -0.249437},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.29686, -0.0334464, 0.114018},
{-1.32659, 0.377523, -0.492212},
{-1.54595, 0.616253, 1.54542},
{-1.3644, 0.331534, -0.0902528},
{-1.13723, -0.03302, 0.430462},
{-1.86492, 0.369171, 0.435679},
{-1.26221, 0.168658, 1.94025},
{-1.62089, 0.247891, 0.329338},
{-1.39758, 0.0769055, 0.185586},
{-1.24035, 0.0471445, 0.0951117},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.7257, 1.34979, -2.66653},
{-1.769, 1.40425, -1.90233},
{-1.36299, 0.771066, 0.967465},
{-1.55748, -0.0338025, 0.549903},
{-1.39179, 0.331608, -0.268987},
{-1.90979, 0.677146, -0.870713},
{-1.48158, -0.0701869, 0.868786},
{-1.84802, -0.485351, 1.60988},
{-1.54901, 0.399841, -0.523386},
{-1.64576, 0.244088, 0.227151},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{0.0194192, 0.302901, 0.342601},
{-0.0254677, 0.469326, 0.61227},
{-0.118976, 0.213262, 1.53799},
{-0.0511441, 0.142864, 1.08353},
{-0.0542953, 0.355048, 1.27015},
{-0.0868421, 0.336026, 1.53954},
{-0.0734036, 0.261845, 1.81035},
{-0.0887113, 0.305155, 2.23019},
{-0.069478, 0.299402, 1.24791},
{-0.00807098, 0.274482, 0.415779},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.308957, 0.904941, -0.401894},
{-0.463173, 1.36705, -0.438891},
{-0.495573, 2.07782, -0.777077},
{-0.35309, 1.30112, -0.278221},
{-0.384928, 1.18199, 0.218563},
{-0.530397, 1.45613, 0.197214},
{-0.533705, 1.59568, 0.113584},
{-0.772638, 2.4454, -0.749781},
{-0.489644, 1.82011, -0.989846},
{-0.289538, 0.843905, -0.332132},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.58035, 0.581259, -0.420185},
{-1.63079, 0.0020396, 0.488355},
{-1.27329, 0.598568, 0.584731},
{-1.43202, 0.0967583, 0.829624},
{-1.33176, 0.182537, 0.545858},
{-1.90425, 0.807962, -0.0423462},
{-1.36712, 0.241394, 0.486282},
{-1.72418, 0.370829, 0.804862},
{-1.45699, 0.156131, 0.652529},
{-1.49491, 0.102651, 0.16661},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.80505, 0.473295, -0.561165},
{-2.02013, 0.950587, -1.05662},
{-1.35315, 0.756275, 0.00397401},
{-1.53042, 0.758478, -0.918424},
{-1.54903, -0.0816109, 0.6537},
{-1.90495, -0.462841, 1.06018},
{-1.62926, 0.303324, 0.310038},
{-1.98481, 1.51705, -2.14954},
{-1.62989, -0.0230485, 0.688019},
{-1.63011, -0.117912, 0.495087},
{-0.1, 1, 2},
{-0.1, 1, 2}}};





//rpc 196 (L9) HV 9.9 kV
/*





double strpos_vs_time[2*nmxtimehit][nlayer][3]= {
{{-0.024813, 0.103942, 0.405951},
{-0.0615963, 0.228577, 0.755601},
{0.0107021, 0.08793, 1.68551},
{-0.0858886, 0.0483545, 1.28225},
{-0.0975513, 0.0678466, 1.33934},
{-0.136087, 0.109142, 1.51309},
{-0.0903227, 0.0671507, 1.40504},
{-0.0928458, -0.0195053, 1.48384},
{-0.0867691, 0.151662, 0.569281},
{-0.0634191, 0.210318, 0.00347566},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.103042, 0.707875, -0.31348},
{-0.272647, 1.23583, -0.616469},
{-0.727225, 2.89996, -1.4571},
{-0.370881, 1.9127, -1.29516},
{-0.415003, 1.84039, -0.938632},
{-0.646062, 2.09298, -0.410501},
{-0.679022, 2.05477, -0.714463},
{-0.768165, 2.72513, -1.46072},
{-0.297796, 1.25793, -0.482567},
{-0.221918, 0.372983, 0.0518777},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.885624, -0.0513916, -0.0594131},
{-1.26372, -0.022945, 0.432721},
{-1.57012, 0.432434, 1.74575},
{-1.35653, 0.0795492, 0.419293},
{-1.52174, 0.332944, 0.334091},
{-1.91115, 0.140807, 1.21597},
{-1.83803, -0.0465923, 1.48333},
{-1.76326, 0.382983, 0.190098},
{-1.62747, 0.244055, 0.211982},
{-1.88131, -0.210955, 0.827474},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.46257, 0.381236, -1.11183},
{-1.71989, 1.35704, -1.79692},
{-1.49305, 0.842176, 0.582097},
{-1.57259, 0.102313, 0.297015},
{-1.65703, 0.350623, 0.0326045},
{-2.18718, 2.01134, -3.04987},
{-1.81073, -1.47192, 3.02426},
{-2.00697, 0.609581, -0.869049},
{-1.83596, 0.484977, -0.691684},
{-2.02916, -3.11552, 7.36227},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.0047398, 0.343696, 0.271653},
{-0.019589, 0.388349, 0.922749},
{-0.106389, 0.102845, 1.96308},
{-0.0900993, 0.0892048, 1.34571},
{-0.111545, 0.162665, 1.55963},
{-0.117605, 0.374565, 1.63268},
{-0.115235, 0.177293, 1.77504},
{-0.0922093, 0.327424, 1.70837},
{-0.0738832, 0.443556, 0.698847},
{-0.0277417, 0.209864, 0.140987},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.254306, 0.702771, -0.120218},
{-0.468255, 1.37996, -0.573594},
{-0.736037, 2.56837, -1.08896},
{-0.404066, 1.85362, -1.14386},
{-0.56454, 2.04362, -1.0109},
{-0.740184, 2.09282, -0.492297},
{-0.874996, 2.64379, -1.35948},
{-0.785011, 2.42694, -1.00934},
{-0.48219, 1.65884, -0.883886},
{-0.358676, 0.603862, -0.398758},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.50216, 0.526046, -0.633382},
{-1.5919, 0.243748, 0.109221},
{-1.64072, 0.686116, 1.00181},
{-1.5746, 0.252301, 0.788362},
{-1.69625, 0.296772, 0.698091},
{-2.08064, 0.397223, 1.06381},
{-2.27544, 1.06054, 0.531966},
{-1.98618, 0.572692, 0.579524},
{-1.72922, 0.353061, 0.180112},
{-2.23296, 0.0202918, 0.439145},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.92091, 1.67357, -2.67366},
{-1.91255, 0.050568, 0.354469},
{-1.54658, 1.72699, -1.73822},
{-1.66401, 0.344944, -0.254318},
{-1.7351, -0.342942, 1.13719},
{-2.34697, 0.530928, -0.458358},
{-2.37634, 0.887588, -0.1734},
{-2.12516, -0.333302, 1.12593},
{-2.00344, -0.177138, 0.987953},
{-2.32491, 0.779882, -0.205638},
{-0.1, 1, 2},
{-0.1, 1, 2}},
};


*/

//20210614 layer 5 at top position
/*
double strpos_vs_time[2*nmxtimehit][nlayer][3]= {
{{-0.0694753, 0.249005, 0.0962191},
{-0.106682, 0.197162, 0.540986},
{-0.0237624, -0.180051, 1.90404},
{-0.137892, 0.149238, 1.1891},
{-0.152609, 0.0335269, 1.25188},
{-0.0762019, 0.106238, 0.271251},
{-0.129098, 0.200511, 1.41201},
{-0.141513, 0.108769, 1.10961},
{-0.167681, 0.236882, 0.937598},
{-0.10784, 0.154339, 0.588031},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{0.0574771, 0.455925, -0.194198},
{-0.0889343, 1.10756, -0.660657},
{-0.607655, 2.68616, -1.70001},
{-0.322184, 1.96891, -1.25603},
{-0.17875, 1.53385, -0.822418},
{-0.0541689, 0.344349, 0.276223},
{-0.689507, 2.23869, -0.863572},
{-0.376301, 2.3363, -1.67267},
{-0.393702, 2.35701, -1.81874},
{-0.283185, 1.23089, -0.463042},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.430817, 0.0793913, -0.202517},
{-1.2364, 0.00177095, 0.307448},
{-1.55732, 0.172024, 0.744367},
{-1.45457, 0.249587, 0.823089},
{-1.18943, 0.306447, 0.0346651},
{-2.0009, -0.0328066, 0.0357835},
{-1.98162, 0.211273, 0.926739},
{-1.20765, 0.214585, -0.0506553},
{-1.40242, 0.102023, 0.617472},
{-1.40796, 0.294745, -0.486604},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.19504, 0.188806, -0.240826},
{-1.57059, -0.247368, 1.05832},
{-1.914, 0.442999, -0.201863},
{-1.85655, 0.0991492, 0.987859},
{-1.53209, 0.21479, 0.268317},
{-2.6515, 6.50356, -8.99532},
{-2.39236, 2.07015, -3.30686},
{-1.71977, 1.03617, -2.21113},
{-1.70539, -0.396404, 1.12108},
{-2.12691, 1.43237, -1.74559},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.0311051, 0.283405, 0.241715},
{-0.0997738, 0.165202, 0.951622},
{-0.116859, 0.178968, 1.95888},
{-0.143133, 0.128171, 1.49086},
{-0.140216, 0.143797, 1.50709},
{-0.101563, 0.307029, 0.175919},
{-0.124547, 0.378391, 1.47692},
{-0.151211, 0.291247, 1.62645},
{-0.163832, 0.375247, 1.30196},
{-0.0758934, 0.308914, 0.777237},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.0589989, 0.520785, -0.227418},
{-0.163674, 1.04252, -0.453349},
{-0.808646, 2.55192, -1.19591},
{-0.42858, 2.04184, -1.11842},
{-0.351499, 1.50686, -0.39771},
{-0.224452, 0.9074, -0.608845},
{-0.765228, 2.06242, -0.382055},
{-0.48862, 1.92399, -0.810765},
{-0.637436, 2.60659, -1.74588},
{-0.585099, 1.68618, -0.817738},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.674392, 0.149512, -0.166899},
{-1.33432, -0.000752371, 0.694189},
{-2.07281, 0.622441, 1.05227},
{-1.59407, 0.164453, 0.997109},
{-1.5093, -0.101276, 1.18926},
{-2.55282, 1.3056, -1.22282},
{-1.85067, -0.136305, 1.42491},
{-1.49353, 0.270031, 0.312402},
{-1.47903, 0.271199, 0.165127},
{-2.03218, 0.253582, 0.465873},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.22694, 0.048959, 0.29221},
{-1.85004, 0.615764, -0.653695},
{-2.30007, 1.19611, -1.14182},
{-1.88936, 0.589444, 0.305982},
{-1.87076, -0.270641, 1.84676},
{-2.97215, 5.3695, -7.46708},
{-2.1989, -0.646064, 1.94764},
{-1.81163, -0.0378112, 0.268362},
{-1.89383, 0.606995, -0.544506},
{-2.35919, -0.350748, 1.9024},
{-0.1, 1, 2},
{-0.1, 1, 2}}};
*/

/* //Just before jim modifying


double strpos_vs_time[2*nmxtimehit][nlayer][3]= { //jim jim
{{-0.060314, 0.223242, 0.115526},
{-0.0956059, 0.232567, 0.468032},
{-0.0704615, -0.208044, 1.93179},
{-0.122478, 0.195499, 0.934874},
{-0.122791, 0.06982, 1.10783},
{-0.0483741, 0.111169, 0.216602},
{-0.108357, 0.279822, 1.15076},
{-0.106302, 0.14241, 1.07882},
{-0.101599, 0.236903, 0.822841},
{-0.0779696, 0.132122, 0.590308},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{0.0623496, 0.458916, -0.194421},
{-0.0651495, 1.19205, -0.891641},
{-0.275775, 0.717207, 0.801334},
{-0.291927, 1.9166, -1.30588},
{-0.107297, 1.33577, -0.620287},
{-0.00325278, 0.333676, 0.247755},
{-0.533622, 1.38739, 0.243739},
{-0.343375, 2.24185, -1.50743},
{-0.291719, 2.21478, -1.72387},
{-0.210766, 1.10911, -0.394979},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.422287, 0.0542924, -0.183341},
{-1.18648, 0.0492624, 0.179883},
{-1.52744, 0.829058, 2.64087},
{-1.45731, 0.260798, 0.63846},
{-1.13364, 0.427437, -0.22807},
{-1.93512, -0.107336, 0.345925},
{-1.84136, 1.35912, 1.22521},
{-1.13593, 0.110628, 0.112438},
{-1.34026, 0.126589, 0.453432},
{-1.34926, 0.32542, -0.449545},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.18912, 0.370183, -0.632572},
{-1.57863, 0.556484, -0.803466},
{-1.52383, -2.34646, 4.91422},
{-0.764262, -1.10195, 2.62571},
{-1.50669, 0.483989, -0.168618},
{-2.32882, 7.26796, -10.0543},
{-2.21531, 1.94616, -3.07942},
{-1.63341, 0.909786, -2.00117},
{-1.63659, -0.386073, 0.937396},
{-2.20061, 2.23539, -3.179},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.0261495, 0.296724, 0.194045},
{-0.0961549, 0.182669, 0.884351},
{-0.1061, 0.228002, 1.78591},
{-0.130476, 0.137584, 1.35748},
{-0.115042, 0.163002, 1.39774},
{-0.0871426, 0.315616, 0.132084},
{-0.115925, 0.361524, 1.40677},
{-0.0998681, 0.318975, 1.38135},
{-0.0948436, 0.38085, 1.14036},
{-0.0463567, 0.3024, 0.635138},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.0474032, 0.487042, -0.193853},
{-0.170602, 1.15193, -0.584116},
{-0.753918, 2.41076, -1.14189},
{-0.3626, 1.82042, -0.905762},
{-0.307104, 1.44083, -0.384248},
{-0.205385, 0.929816, -0.661647},
{-0.707041, 1.78735, -0.117967},
{-0.393849, 1.70323, -0.636605},
{-0.504525, 2.2547, -1.37314},
{-0.513376, 1.4495, -0.623667},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.654653, 0.123764, -0.0752707},
{-1.32262, -0.0739905, 0.767044},
{-2.02463, 0.473912, 1.20624},
{-1.54274, 0.654261, 1.09606},
{-1.4361, -0.075547, 1.01619},
{-2.53048, 1.03376, -0.453812},
{-1.81855, -0.478361, 1.90183},
{-1.40586, 0.158123, 0.470562},
{-1.41241, 0.33155, -0.0668696},
{-1.99014, -0.000397913, 1.00189},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.16673, -0.00416365, 0.214266},
{-1.82638, 0.792902, -0.950764},
{-2.20445, 0.642546, -0.240816},
{-1.75712, 0.348231, 0.443689},
{-1.86123, 0.317524, 0.397121},
{-2.8209, -1.63735, 12.9245},
{-2.23335, 0.168096, 0.0823259},
{-1.69695, -0.293204, 0.645877},
{-1.82848, 0.711683, -0.908215},
{-2.20979, -0.582413, 2.0145},
{-0.1, 1, 2},
{-0.1, 1, 2}},
  };


*/

/*
//After jim modifying 20210922 with buttons
double strpos_vs_time[2*nmxtimehit][nlayer][3]= { //jim jim
{{0.0241002, 0.203819, 0.156117},
{-0.0496038, 0.174001, 0.66763},
{-0.0740419, 0.0174115, 1.71636},
{-0.0959757, 0.118462, 1.3227},
{-0.12143, -0.0876975, 1.76855},
{-0.0668307, 0.0188426, 1.43294},
{-0.10548, -0.0859754, 2.25875},
{-0.0834303, -0.0219795, 1.38564},
{-0.0977418, 0.121281, 0.730005},
{-0.0572004, 0.20388, 0.14916},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.0849193, 0.61549, -0.177967},
{-0.262115, 1.51674, -1.21133},
{-0.926295, 3.4342, -2.67681},
{-0.471602, 2.10991, -1.30669},
{-0.306582, 1.7907, -1.07591},
{-0.497628, 1.89786, -0.765964},
{-0.827091, 2.40183, -0.588677},
{-0.595867, 2.42288, -1.38293},
{-0.292561, 1.44946, -0.925886},
{-0.434795, 0.989224, -0.41286},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.473341, 0.156654, -0.245993},
{-1.45986, 0.482811, -0.378995},
{-1.73266, 0.159814, 0.495506},
{-1.32956, -0.095357, 0.821744},
{-1.211, 0.247635, -0.216042},
{-1.27315, 0.159703, 0.278703},
{-1.7251, -0.200836, 1.57162},
{-1.12596, 0.1955, -0.0187315},
{-1.18458, 0.198627, -0.239126},
{-1.4287, 0.240301, -0.581358},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.998187, 0.35982, -0.667077},
{-1.46804, -0.357071, 0.48977},
{-2.16955, 1.40192, -2.10823},
{-1.59335, 0.317624, -0.743295},
{-1.53521, -0.615522, 1.71505},
{-1.59508, -0.0381995, 1.36688},
{-2.16472, 1.39861, -1.34042},
{-1.41386, -1.25753, 2.68608},
{-1.44091, 0.196897, -1.11589},
{-1.67272, -2.69045, 4.7011},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{0.0791078, 0.250784, 0.449635},
{-0.0469119, 0.230735, 0.829066},
{-0.108418, 0.20993, 2.31557},
{-0.0942722, 0.0536717, 1.69814},
{-0.0809028, 0.0603304, 1.90532},
{-0.138519, 0.134922, 2.07089},
{-0.105316, 0.316547, 2.17205},
{-0.0782661, 0.234058, 2.03159},
{-0.0850383, 0.401674, 0.868347},
{-0.0319437, 0.36005, 0.324946},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.240336, 0.796463, -0.240223},
{-0.364143, 1.55878, -0.916299},
{-1.11264, 3.10887, -1.28127},
{-0.606933, 2.2984, -1.00428},
{-0.588009, 1.90636, -0.273446},
{-0.55225, 1.74824, -0.207721},
{-1.01563, 2.83266, -0.951055},
{-0.677876, 2.42772, -1.13578},
{-0.509247, 1.80929, -0.987427},
{-0.521422, 1.26583, -0.635209},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-0.792573, 0.152562, 0.0857114},
{-1.55259, 0.147977, 0.489447},
{-2.25056, 0.699833, 1.40078},
{-1.64044, 0.233563, 0.877477},
{-1.6703, 0.28306, 0.911145},
{-1.48948, -0.0283451, 0.979202},
{-1.88906, -0.0436606, 1.54477},
{-1.44096, 0.121886, 0.691478},
{-1.4061, 0.110522, 0.282977},
{-2.01762, -0.0290871, 0.551294},
{-0.1, 1, 2},
{-0.1, 1, 2}},

{{-1.00981, -0.203586, 0.766933},
{-1.74449, 0.37487, -0.332377},
{-2.18107, 0.061734, 1.36904},
{-1.63273, -0.627522, 2.46345},
{-1.69014, 0.905881, -1.12403},
 {-1.46942, -0.365314, 0.655367},
{-2.19462, 0.915752, -0.464729},
{-1.77652, 1.56786, -2.03768},
{-1.60149, 0.258702, 0.330725},
{-2.42455, 2.05731, -3.4655},
{-0.1, 1, 2},
{-0.1, 1, 2}}};
*/

/*
//After removing buttons 20210922 jim jim
double strpos_vs_time[2*nmxtimehit][nlayer][3]= {
  {{0.0283716, 0.23694, 0.109458},
   {-0.0391088, 0.183605, 0.646273},
   {-0.0616182, 0.0470743, 1.68405},
   {-0.0788627, 0.106218, 1.34156},
   {-0.105831, -0.0628332, 1.70798},
   {-0.0412041, 0.00796176, 1.47027},
   {-0.08094, -0.122703, 2.32804},
   {-0.0551267, -0.0398852, 1.42957},
   {-0.0680604, 0.111178, 0.743814},
   {-0.025431, 0.212903, 0.140969},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.071762, 0.654081, -0.292784},
   {-0.237872, 1.49508, -1.16844},
   {-0.924031, 3.50166, -2.76408},
   {-0.45713, 2.16856, -1.42647},
   {-0.285877, 1.82365, -1.10279},
   {-0.470156, 1.92467, -0.826766},
   {-0.816962, 2.39076, -0.468805},
   {-0.575677, 2.42653, -1.34953},
   {-0.266309, 1.49215, -0.982326},
   {-0.389988, 0.992132, -0.44443},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.447162, 0.204593, -0.370515},
   {-1.42066, 0.707193, -0.876186},
   {-1.73029, 0.452863, -0.00899947},
   {-1.26211, -0.372793, 1.28876},
   {-1.19018, 0.350612, -0.343072},
   {-1.21211, -0.0280769, 0.550029},
   {-1.70413, -0.283893, 1.69383},
   {-1.07101, 0.183423, -0.162738},
   {-1.12947, 0.139893, -0.0653493},
   {-1.27707, 0.176559, -0.77979},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.991939, 0.346866, -0.604155},
   {-1.41544, -0.250802, 0.0714012},
   {-2.12894, 1.21268, -1.67601},
   {-1.56211, 0.252376, -0.29215},
   {-1.46744, -0.547965, 1.51761},
   {-1.64323, 0.924415, -0.633367},
   {-1.99641, 0.0223809, 1.06844},
   {-1.43777, -1.25217, 3.1043},
   {-1.34445, 0.0761109, -1.28926},
   {-1.57173, -3.25447, 5.4505},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{0.0868812, 0.218399, 0.525979},
   {-0.0274652, 0.188276, 0.902464},
   {-0.0863952, 0.176992, 2.34296},
   {-0.0750694, 0.0342211, 1.74855},
   {-0.0610966, 0.0663597, 1.88417},
   {-0.108195, 0.092996, 2.1466},
   {-0.0794371, 0.285108, 2.24945},
   {-0.0458681, 0.210881, 2.06561},
   {-0.0546162, 0.417907, 0.819773},
   {0.00457328, 0.359317, 0.328092},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.227117, 0.811955, -0.26948},
   {-0.350761, 1.60214, -0.986853},
   {-1.11453, 3.25542, -1.50507},
   {-0.586994, 2.34276, -1.08974},
   {-0.564299, 1.952, -0.379796},
   {-0.527473, 1.73689, -0.153963},
   {-0.993695, 2.81848, -0.868748},
   {-0.658519, 2.50284, -1.23021},
   {-0.476457, 1.8256, -0.987368},
   {-0.480805, 1.34373, -0.768588},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.753799, 0.10559, 0.170346},
   {-1.52632, 0.310451, 0.237999},
   {-2.20937, 0.64011, 1.43153},
   {-1.60365, 0.255249, 0.84515},
   {-1.6358, 0.258867, 0.9878},
   {-1.48323, 0.246597, 0.425577},
   {-1.84786, -0.0255499, 1.4592},
   {-1.38113, 0.0467761, 0.788908},
   {-1.35751, 0.139108, 0.219833},
   {-1.98368, 0.414587, -0.27138},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-1.00503, -0.0865787, 0.662209},
   {-1.70522, 0.913322, -1.36611},
   {-2.24505, 0.789926, 0.233422},
   {-1.5683, -1.22502, 3.47842},
   {-1.63092, 0.455343, -0.0781352},
   {-1.47378, 0.412449, -1.0925},
   {-2.15467, 0.760409, -0.248523},
   {-1.72803, 1.41527, -1.74827},
   {-1.51241, -0.0670172, 0.75439},
   {-2.38843, 2.47071, -3.68155},
   {-0.1, 1, 2},
   {-0.1, 1, 2}}

};

*/

/*
// 20210528 data with layer 9 on top all 39 ohm resistors jim jim
double strpos_vs_time[2*nmxtimehit][nlayer][3]= {
{{-0.00270629, -0.0649766, 0.0869805},
 {-0.0058035, -0.0633149, 0.0929842},
 {-0.0227812, 0.298467, -0.578908},
 {-0.0213641, -0.00574275, 0.0243322},
 {-0.010217, 0.0580722, -0.0761703},
 {-0.0221749, 0.0775263, -0.0926914},
 {0.00519485, -0.102412, 0.278815},
 {0.0153744, -0.0497015, 0.0991808},
 {0.0232899, -0.00731983, 0.0229845},
 {-0.0141752, 0.0993492, -0.147479},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.021723, 0.0243535, -0.00580714},
 {-0.0427331, 0.141102, -0.153931},
 {-0.0888936, 0.108998, -0.115428},
 {-0.0188348, -0.0626928, 0.150642},
 {-0.0119128, 0.267599, -0.460338},
 {-0.0724769, 0.359499, -0.668967},
 {-0.0776955, 0.0290053, 0.0878901},
 {0.0203039, -0.168884, 0.3226},
 {0.0506757, -0.213582, 0.388639},
 {-0.0261986, -0.0304933, 0.133203},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.1328, -0.0275749, 0.156351},
 {-0.0483867, -0.0170261, 0.255773},
 {-0.15045, 0.0519303, -0.298444},
 {0.0166571, 0.70809, 0.675044},
 {0.060884, -0.625288, 1.18124},
 {0.079793, -0.140515, 0.485889},
 {-0.132667, -0.311034, 0.188228},
 {-0.0229836, 0.110262, -0.265201},
 {-0.0350799, 0.45544, -0.890373},
 {-0.0820902, -0.175579, 0.351729},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{0.0288087, -0.717877, 0.868314},
 {0.0999026, -0.420168, 0.553312},
 {-0.267528, 2.2552, -4.52044},
 {0.0707259, 1.99449, -1.40567},
 {0.122608, -0.703962, 0.836778},
 {-0.0225985, -5.57253, 8.05145},
 {0.0284424, -1.30063, 1.99973},
 {0.107185, -1.17828, 2.38777},
 {0.0374289, 0.288466, -0.508058},
 {0.142175, -1.2796, 1.7765},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{0.0261649, 0.00878101, 0.061208},
 {-0.012225, 0.0354965, -0.0901814},
 {-0.00115676, 0.0416588, -0.0645994},
 {-0.0202394, 0.0346996, 0.00968428},
 {0.00223241, -0.0157623, 0.0769124},
 {-0.0195384, 0.0447614, -0.0290576},
 {0.00190605, -0.0351909, 0.167568},
 {0.0132143, 0.0706521, -0.0179991},
 {0.0241024, -0.00416777, 0.189427},
 {-0.0196118, 0.146448, -0.0537321},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.062859, 0.0650832, 0.0191323},
 {-0.0272263, 0.007697, -0.0177765},
 {-0.0655359, -0.106713, 0.249177},
 {-0.0400718, 0.0934716, -0.105229},
 {-0.0265129, 0.151516, -0.164733},
 {-0.024202, -0.113674, 0.258106},
 {-0.0707518, -0.0341041, 0.179806},
 {-0.0333661, 0.190338, -0.205601},
 {-0.0483168, 0.355016, -0.456926},
 {-0.0790944, 0.201795, -0.159397},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.160459, 0.00680401, 0.0470617},
 {-0.0619665, 0.187014, -0.355453},
 {-0.157164, 0.438895, -0.85327},
 {-0.0337004, 0.0471162, -0.0439314},
 {0.0030825, 0.402743, -0.915745},
 {0.0173801, -0.943926, 1.54134},
 {-0.174107, 0.975808, -1.8476},
 {-0.0131928, 0.00536785, -0.00265619},
 {0.026024, -0.0766667, 0.166978},
 {-0.112818, 0.900909, -1.50322},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.150983, 0.0944222, 0.0215536},
 {0.0277063, -0.528389, 1.06222},
 {-0.278792, 0.631818, -0.825591},
 {-0.102587, 0.699219, -1.56199},
 {0.101604, -0.0406687, -0.499123},
 {0.537094, 0.855371, -11.806},
 {-0.374902, 1.56157, -2.3281},
 {-0.107747, 0.608274, -0.456393},
 {0.0414123, -0.038398, -0.136865},
 {-0.184756, 0.788524, -1.77174},
 {-0.1, 1, 2},
 {-0.1, 1, 2}}

};
*/

/*
// 20210616 data with layer 9 on top modified 15,18,22 ohm resistors jim jim something wrong in layer 5 correction
double strpos_vs_time[2*nmxtimehit][nlayer][3]= {
{{-0.0694753, 0.249005, 0.0962191},
 {-0.106682, 0.197162, 0.540986},
 {-0.0237624, -0.180051, 1.90404},
 {-0.137892, 0.149238, 1.1891},
 {-0.152609, 0.0335269, 1.25188},
 {-0.0762019, 0.106238, 0.271251},
 {-0.129098, 0.200511, 1.41201},
 {-0.141513, 0.108769, 1.10961},
 {-0.167681, 0.236882, 0.937598},
 {-0.10784, 0.154339, 0.588031},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{0.0574771, 0.455925, -0.194198},
 {-0.0889343, 1.10756, -0.660657},
 {-0.607655, 2.68616, -1.70001},
 {-0.322184, 1.96891, -1.25603},
 {-0.17875, 1.53385, -0.822418},
 {-0.0541689, 0.344349, 0.276223},
 {-0.689507, 2.23869, -0.863572},
 {-0.376301, 2.3363, -1.67267},
 {-0.393702, 2.35701, -1.81874},
 {-0.283185, 1.23089, -0.463042},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.430817, 0.0793913, -0.202517},
 {-1.2364, 0.00177095, 0.307448},
 {-1.55732, 0.172024, 0.744367},
 {-1.45457, 0.249587, 0.823089},
 {-1.18943, 0.306447, 0.0346651},
 {-2.0009, -0.0328066, 0.0357835},
 {-1.98162, 0.211273, 0.926739},
 {-1.20765, 0.214585, -0.0506553},
 {-1.40242, 0.102023, 0.617472},
 {-1.40796, 0.294745, -0.486604},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-1.19504, 0.188806, -0.240826},
 {-1.57059, -0.247368, 1.05832},
 {-1.914, 0.442999, -0.201863},
 {-1.85655, 0.0991492, 0.987859},
 {-1.53209, 0.21479, 0.268317},
 {-2.6515, 6.50356, -8.99532},
 {-2.39236, 2.07015, -3.30686},
 {-1.71977, 1.03617, -2.21113},
 {-1.70539, -0.396404, 1.12108},
 {-2.12691, 1.43237, -1.74559},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.0311051, 0.283405, 0.241715},
 {-0.0997738, 0.165202, 0.951622},
 {-0.116859, 0.178968, 1.95888},
 {-0.143133, 0.128171, 1.49086},
 {-0.140216, 0.143797, 1.50709},
 {-0.101563, 0.307029, 0.175919},
 {-0.124547, 0.378391, 1.47692},
 {-0.151211, 0.291247, 1.62645},
 {-0.163832, 0.375247, 1.30196},
 {-0.0758934, 0.308914, 0.777237},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.0589989, 0.520785, -0.227418},
 {-0.163674, 1.04252, -0.453349},
 {-0.808646, 2.55192, -1.19591},
 {-0.42858, 2.04184, -1.11842},
 {-0.351499, 1.50686, -0.39771},
 {-0.224452, 0.9074, -0.608845},
 {-0.765228, 2.06242, -0.382055},
 {-0.48862, 1.92399, -0.810765},
 {-0.637436, 2.60659, -1.74588},
 {-0.585099, 1.68618, -0.817738},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.674392, 0.149512, -0.166899},
 {-1.33432, -0.000752371, 0.694189},
 {-2.07281, 0.622441, 1.05227},
 {-1.59407, 0.164453, 0.997109},
 {-1.5093, -0.101276, 1.18926},
 {-2.55282, 1.3056, -1.22282},
 {-1.85067, -0.136305, 1.42491},
 {-1.49353, 0.270031, 0.312402},
 {-1.47903, 0.271199, 0.165127},
 {-2.03218, 0.253582, 0.465873},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-1.22694, 0.048959, 0.29221},
 {-1.85004, 0.615764, -0.653695},
 {-2.30007, 1.19611, -1.14182},
 {-1.88936, 0.589444, 0.305982},
 {-1.87076, -0.270641, 1.84676},
 {-2.97215, 5.3695, -7.46708},
 {-2.1989, -0.646064, 1.94764},
 {-1.81163, -0.0378112, 0.268362},
 {-1.89383, 0.606995, -0.544506},
 {-2.35919, -0.350748, 1.9024},
 {-0.1, 1, 2},
 {-0.1, 1, 2}}

};
*/

/*
// 20210616 data with layer 9 on top modified 15,18,22 ohm resistors jim jim jj144
double strpos_vs_time[2*nmxtimehit][nlayer][3]= {
{{-0.0806618, 0.281413, 0.0760067},
 {-0.126003, 0.259846, 0.438288},
 {-0.0938806, -0.272223, 2.16565},
 {-0.142127, 0.123312, 1.29194},
 {-0.151309, 0.00162875, 1.38503},
 {-0.080901, 0.180052, 0.183207},
 {-0.119752, 0.184072, 1.57627},
 {-0.124054, 0.0759701, 1.32031},
 {-0.163245, 0.273702, 0.989947},
 {-0.108008, 0.188302, 0.609577},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{0.0153953, 0.55041, -0.266645},
 {-0.128082, 1.38872, -1.08381},
 {-0.727889, 2.7867, -1.63574},
 {-0.362746, 2.18362, -1.5105},
 {-0.205253, 1.72287, -1.01836},
 {-0.0537331, 0.59686, -0.232766},
 {-0.733487, 2.57621, -1.32856},
 {-0.408081, 2.64223, -2.09023},
 {-0.409481, 2.59315, -2.15008},
 {-0.274078, 1.37879, -0.771757},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.461422, 0.234945, -0.447654},
 {-1.25175, 0.254569, -0.133443},
 {-1.62504, 0.0515469, 1.23736},
 {-1.47437, 0.390067, 0.400951},
 {-1.21008, 0.530649, -0.324346},
 {-2.03055, 0.0349942, 0.31229},
 {-2.03298, 0.553492, 0.521512},
 {-1.19972, 0.321005, -0.209965},
 {-1.40397, 0.495368, -0.120466},
 {-1.31763, 0.0662547, -0.100912},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-1.15208, -0.111478, 0.43971},
 {-1.41845, -0.670551, 1.50485},
 {-1.93455, 0.535635, -0.827151},
 {-1.88155, 0.401836, 0.26851},
 {-1.48602, -0.340619, 1.27841},
 {-2.47629, 5.96759, -10.3391},
 {-2.2344, 0.731253, -0.809958},
 {-1.5703, -0.26708, 0.523602},
 {-1.69988, 0.463641, -0.664264},
 {-1.78581, -0.295727, -0.021295},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.0428639, 0.21683, 0.372328},
 {-0.123596, 0.222431, 0.760835},
 {-0.127161, 0.166869, 1.99891},
 {-0.146029, 0.0805903, 1.6121},
 {-0.131821, 0.107431, 1.59336},
 {-0.0950991, 0.293858, 0.248522},
 {-0.118911, 0.410393, 1.558},
 {-0.129754, 0.252939, 1.8608},
 {-0.143487, 0.303799, 1.50239},
 {-0.0748659, 0.355612, 0.760985},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.11348, 0.534821, -0.194716},
 {-0.228022, 1.36797, -0.899927},
 {-0.897038, 2.95443, -1.66496},
 {-0.467745, 2.20315, -1.23506},
 {-0.388332, 1.70861, -0.560852},
 {-0.238157, 0.976908, -0.708774},
 {-0.781906, 2.01569, -0.113116},
 {-0.541491, 2.33098, -1.27941},
 {-0.675305, 2.87593, -2.02687},
 {-0.590292, 1.88866, -1.19545},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-0.64828, 0.0392538, -0.037283},
 {-1.35915, 0.311922, 0.152854},
 {-2.08203, 0.906964, 0.53969},
 {-1.64481, 0.577711, 0.268608},
 {-1.56554, 0.367473, 0.470653},
 {-2.4151, -0.059317, 1.39474},
 {-1.89857, 0.370572, 0.656341},
 {-1.46547, 0.20485, 0.557135},
 {-1.4152, 0.0923303, 0.419576},
 {-2.00754, 0.308416, 0.56101},
 {-0.1, 1, 2},
 {-0.1, 1, 2}},

{{-1.3052, 0.492897, -0.535473},
 {-1.75904, 0.20583, 0.0424339},
 {-2.30428, 0.972711, -0.621546},
 {-1.90939, 0.870508, -0.321482},
 {-1.86785, 0.466113, -0.185377},
 {-2.53907, -0.547731, 2.75272},
 {-2.28732, 0.171257, 0.53104},
 {-1.91186, 0.605763, -0.377298},
 {-1.77496, 0.171042, 0.0662591},
 {-2.28099, -0.285826, 1.21485},
 {-0.1, 1, 2},
 {-0.1, 1, 2}}

};
*/


/*
double strpos_vs_time[2*nmxtimehit][nlayer][3]= {
  {{-0.1, 1, 2},
   {-0.00742301, 0.184462, 0.131238},
   {-0.0726225, 0.19087, 0.901143},
   {-0.0826275, 0.107989, 1.09605},
   {-0.111848, 0.0269024, 1.29682},
   {-0.1, 1, 2},
   {-0.0856062, 0.0102598, 1.77302},
   {-0.0616103, -0.0836363, 1.44963},
   {-0.0894981, 0.229494, 0.323477},
   {-0.0236952, 0.0950921, 0.345416},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.1, 1, 2},
   {-0.0889827, 0.682834, -0.500263},
   {-0.468595, 1.34504, -0.544125},
   {-0.393341, 1.97517, -1.48003},
   {-0.232037, 1.48109, -0.815321},
   {-0.1, 1, 2},
   {-0.677158, 1.89591, -0.377848},
   {-0.525811, 1.96111, -0.821341},
   {-0.257523, 1.47692, -1.03861},
   {-0.39786, 0.617511, 0.213837},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.1, 1, 2},
   {-1.36196, 0.211256, -0.191627},
   {-1.76212, 0.224115, 0.773333},
   {-1.39156, 0.273621, 0.139286},
   {-1.32924, 0.0592473, 0.459536},
   {-0.1, 1, 2},
   {-1.80897, -0.339814, 1.94393},
   {-1.22102, 0.239915, 0.171449},
   {-1.22055, 0.34978, -0.290603},
   {-1.44083, -0.252155, 0.0582215},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.1, 1, 2},
   {-1.40524, -1.06702, 1.67554},
   {-1.99329, 0.60474, -0.699506},
   {-1.70249, 0.993747, -1.47345},
   {-1.61717, 0.25576, 0.00827773},
   {-0.1, 1, 2},
   {-1.99803, 0.293025, -0.324337},
   {-1.68301, 0.330479, 0.200454},
   {-1.50223, 0.54449, -0.705361},
   {-1.8255, 0.103049, -2.21487},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.1, 1, 2},
   {0.00265726, 0.16601, 0.381162},
   {-0.0746955, 0.360941, 1.30056},
   {-0.0871353, 0.161841, 1.28421},
   {-0.0679057, 0.0866404, 1.66127},
   {-0.1, 1, 2},
   {-0.0597036, 0.212156, 2.01078},
   {-0.0690801, 0.447488, 1.2893},
   {-0.0643342, 0.396644, 0.825977},
   {-0.00406406, 0.382384, 0.226233},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.1, 1, 2},
   {-0.173856, 0.767059, -0.402311},
   {-0.845798, 2.15616, -0.977794},
   {-0.526448, 2.18964, -1.3461},
   {-0.532106, 1.61888, -0.29008},
   {-0.1, 1, 2},
   {-0.858846, 2.02028, -0.134518},
   {-0.627734, 2.301, -1.14392},
   {-0.487517, 1.70202, -0.846552},
   {-0.51528, 1.14884, -0.468393},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.1, 1, 2},
   {-1.53309, 0.223198, -0.0376532},
   {-2.23402, 0.624626, 0.654162},
   {-1.66224, 0.438846, 0.413713},
   {-1.85663, 0.664562, 0.0642409},
   {-0.1, 1, 2},
   {-2.02302, 0.342375, 0.679817},
   {-1.5079, 0.242389, 0.500221},
   {-1.54026, 0.268444, 0.422453},
   {-2.07277, 0.0356366, 0.439713},
   {-0.1, 1, 2},
   {-0.1, 1, 2}},

  {{-0.1, 1, 2},
   {-1.70596, 0.239288, -0.168292},
   {-2.20324, 0.39034, 0.52843},
   {-1.74516, 1.40578, -1.55704},
   {-1.85196, 0.729668, -0.712391},
   {-0.1, 1, 2},
   {-2.16977, 1.19327, -1.49381},
   {-1.60373, -0.142435, 0.821837},
   {-1.58868, 0.175762, -0.0210929},
   {-2.20285, -0.342131, 0.84556},
   {-0.1, 1, 2},
   {-0.1, 1, 2}}};
*/


//Button positon for 20181226 jim jim
double button_xpos[10][9] = {{6.15, 12.5, 19, 25.4, 31.6, 38.1, 44.5, 50.9,1000}, {6.35, 12.75, 19.2, 25.6, 32, 38.4, 44.8, 51.2,1000}, {6.1, 12.4, 18.75, 25.1, 31.6, 37.9, 44.4, 50.7,1000}, {6.35, 12.7, 19.1, 25.5, 31.9, 38.3, 44.7, 51.13,1000}, {6.25, 12.65, 19, 25.4, 31.75, 38.13, 44.5, 50.9,1000}, {3.4, 9.7, 16.05, 22.4, 28.62, 34.9, 41.3, 47.6, 53.85}, {3.8, 10.9, 18, 25.12, 32.36, 39.5, 46.75, 54,1000}, {3.4, 9.85, 16.1, 22.4, 28.75, 35.1, 41.4, 47.65, 53.9},{6.5, 12.8, 19.125, 25.5, 31.85, 38.25, 44.6, 51,1000}, {3.6, 9.9, 16.2, 22.1, 28.8, 35.1, 41.4, 47.85, 54}};
double button_ypos[10][9] = {{6.4, 13.1, 20, 26.8, 33.6, 40.4, 47.2, 54,1000}, {6.1, 12.8, 19.8, 26.6, 33.35, 40.1, 46.9, 53.65,1000}, {6.4, 13.1, 19.9, 26.9, 33.75, 40.4, 47.1, 53.9,1000}, {6.3, 13.1, 19.85, 26.6, 33.6, 40.4, 47.1, 53.85,1000}, {6.15, 12.9, 19.7, 26.55, 33.4, 40.12, 46.9, 53.6,1000}, {3.4, 10.1, 16.6, 23.38, 30, 36.7, 43.4, 50, 56.61}, {3.8, 11.25, 18.7, 26.2, 33.65, 41.1, 48.6, 56,1000}, {3.4, 10.1, 16.6, 23.4, 30, 36.6, 43.4, 50, 56.6}, {6.3, 13.1, 19.9, 26.6, 33.5, 40.3, 47.1, 53.9,1000}, {3.6, 10.1, 16.9, 23.4, 30.4, 36.9, 43.6, 50.1, 56.8}};


/*
//Button position from 20210922 jim jim
double button_xpos[10][9] = {
  {6.2,12.6,18.9,25.3,31.7,38,44.55,50.85,1000},
  {5.8,12.1,18.4,24.85,31.3,37.7,44.1,50.5,1000},
  {3.15,9.38,15.73,22.1,28.4,34.65,41,47.38,53.4},
  {6.4,12.75,19.15,25.6,31.9,38.3,44.7,51.1,1000},
  {6.3,12.5,19.1,25.5,31.9,38.4,44.7,51.15,1000},
  {3.5,9.65,16,22.4,28.75,35.125,41.4,47.8,54.2},
  {3.8,10.1,16.375,22.55,28.875,35.2,41.5,47.75,54.1},
  {3.5,9.85,16.125,22.4,28.7,35,41.35,47.625,53.9},
  {6.325,12.75,19.1,25.5,31.875,38.35,44.6,51.125,1000},
  {3.625,10,16.3,22.5,28.75,35.1,41.325,47.65,54}};
double button_ypos[10][9] = {
  {6.6,13.35,19.95,26.85,33.55,40.3,47.1,53.9,1000},
  {6.5,13.2,20.1,26.8,33.65,40.5,47.4,54.15,1000},
  {3.8,10.4,17,23.8,30.4,37.25,43.8,50.5,57.1},
  {6.3,13.125,19.9,26.65,33.6,40.35,47.05,53.8,1000},
  {6.4,13.125,19.9,26.65,33.55,40.4,47.1,53.9,1000},
  {3.8,10.3,16.8,23.6,30.3,36.9,43.5,50.3,56.7},
  {3.15,9.7,16.25,23,29.6,36.25,43,49.6,56.125},
  {3.4,10.1,16.7,23.375,30,36.7,43.35,49.9,56.6},
  {6.1,12.825,19.75,26.625,33.4,40.2,47,53.8,1000},
  {3.325,10,16.725,23.4,30.25,36.9,43.6,50.15,56.825}};
*/

double biasinxtime[nlayer]={0};
double biasinytime[nlayer]={0};

const double stripwidth = 3.0; // cm
const double stripgap = 0.2; // cm
const double layergap = 16.0; // cm
const double fact = stripwidth; //layergap;
const double sigmaz = 0.1; //in cm as the smallest division in scale
//  const double sigmapos = 0.8; //in cm as 2.8/sqrt(12). its a uniform distribution
const double sigmarpc = 1.5; //in ns RPC time resolution
const int nmnhits =  6;//3;//4;//9;//5;//4; //9; //Minimum number of layers required for position and time measurements (must be more than 2); // jim jim modifiy here
const int nmnentry = 10; // Statistics requred in a strip for time correction and use in next iteration

bool isTiming = true;//true;//true;//true; // true //Don't do timing fit etc
bool useallpos=true; //true; //While use all strip irrespecive of aligned or not
bool usealltime=false; //While use timing of all strip irrespecive of aligned or not
bool timefit;

const float xyPosDev=1000.0; // seven sigma 2.0; //maximum deviation of points from fit line

#ifdef ONEMBYONEM
const int trigly1 =2; //0;
const int trigly2 =4;
const int trigly3 =7;
const int trigly4 =9; //11;
#else
const int trigly1 =6;//1;//1;//4;//4;//1;//4;//1;//8;//1;//0;//1;//0;//1;//2;
const int trigly2 =7;//3;//3;//5;//5;//2;//5;//2;//9;//2;//9;//1;//2;//1;//2;//5;
const int trigly3 =8;//5;//5;//6;//6;//9;//6;//3;//10;//9;//10;//2;//9;//2;//9;//7;
const int trigly4 =9;//7;//7;//10;//7;//4;//10;//11;//3;//10;//3;//10;//9;

#endif

const int nDelayPar=6;
double timecorpar[2][nlayer][nstrip][nDelayPar]={0.0};

Double_t polfunc(Double_t* x, Double_t* par) {
  return par[0]+
    par[1]*x[0]+
    par[2]*x[0]*x[0]+
    par[3]*x[0]*x[0]*x[0]+
    par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t cal_slope(Double_t* x, Double_t* par) {
  if (x[0]<nstrip/2. -0.5) {
    return par[0] + par[1]*(x[0] - nstrip/4. +0.5);
  } else {
    double par3 = (par[2]-par[0]-par[1]*nstrip/4.)/(nstrip/4.);
    return par[2] + par3*(x[0] - 3*nstrip/4. +0.5);
  }
}

Double_t cal_slope2(Double_t x, Double_t* par) {
  if (x<nstrip/2. -0.5) {
    return par[0] + par[1]*(x - nstrip/4. +0.5);
  } else {
    double par3 = (par[2]-par[0]-par[1]*nstrip/4.)/(nstrip/4.);
    return par[2] + par3*(x - 3*nstrip/4. +0.5);
  }
}

double path_correction_time(int ixy, int il, int ist, double xx) {
  double xval[2];
  double par[nDelayPar];
  xval[0] = xx;

  for (int ij=0; ij<nDelayPar; ij++) {
    par[ij] = timecorpar[ixy][il][ist][ij];
  }

  double cor= polfunc(xval, par) - par[nDelayPar-1]; //compared wrt central one

  //  xval[0] = 16.0;
  //  if (gRandom->Uniform()<1.e-3) cout<<"ixy "<<ixy<<" "<<il<<" "<<ist<<" "<<xx<<" "<< cor<<" "<<par[0]<<" "<< timecorpar[ixy][il][ist][5]<<" "<<polfunc(xval, par)<<endl;


  if (cor > 5.0) cor =  5.0;
  if (cor < -5.0) cor = -5.0;

  return cor;
}

void GetXposInStrip(int ixy, double* ext, double* otherext, double* off, double* pos, double* local) {
  for (int ij=0; ij<nlayer; ij++) {
    local[ij] = ext[ij]; //+off[ij];
#if defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
    if (ixy==0) {
      local[ij] +=cal_slope2(otherext[ij], &align_ystr_xdev[ij][0]);
      local[ij] +=cal_slope2(local[ij], &align_xstr_xdev[ij][0]);
    } else {
      local[ij] +=cal_slope2(otherext[ij], &align_xstr_ydev[ij][0]);
      local[ij] +=cal_slope2(local[ij], &align_ystr_ydev[ij][0]);
    }
#endif
    local[ij] +=off[ij];
    int istr = int(local[ij]);
    if (local[ij]<0.0) {
      pos[ij] = local[ij] - istr + 0.5;
    } else {
      pos[ij] = local[ij] - istr - 0.5;
    }

    local[ij] -= 0.5;

    //    cout <<"ij "<<ij<<" "<< ext[ij]<<" "<<off[ij]<<" "<<pos[ij]<<endl;
  }
}

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -1;
  for (int ix=1; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return 1000;
}

// double Median(const TH1D* h1) {
//   int n = h1->GetXaxis()->GetNbins();
//   std::vector<double>  x(n);
//   h1->GetXaxis()->GetCenter( &x[0] );
//   const double* y = h1->GetArray();
//   // exclude underflow/overflows from bin content array y
//   return TMath::Median(n, &x[0], &y[1]);
// }

double get_average(TH2F* histx) { //GMA
  int ncnt=0;
  double sum=0;
  for (int ij=2; ij<(lastXstrip/2)-1; ij++) {
    for (int jk=2; jk<(lastYstrip/2)-1; jk++) {
      ncnt++;
      sum += histx->GetBinContent(ij, jk);
    }
  }
  cout<<sum/max(1, ncnt)<<endl;
  return sum/max(1, ncnt);

}
double get_average1(TH2F* histx) { //GMA
  int ncnt=0;
  double sum=0;
  for (int ij=1; ij<lastXstrip-1; ij++) {
    for (int jk=1; jk<lastYstrip-1; jk++) {
      ncnt++;
      sum += histx->GetBinContent(ij, jk);
    }
  }
  cout<<sum/max(1, ncnt)<<endl;
  return sum/max(1, ncnt);

}


int main(int argc, char **argv) {
  int hvset = stoi(argv[2]);
  int obslay = stoi(argv[3]);

  if (nmnhits <=2) { //GMA
    cout<<"nmnhits must be more than two, but value is ,"<<nmnhits <<",Change this value and rerun the code"<<endl;
    return 0;
  }
  //gSystem->Load("libTree");

  // ----------------- Don't Change Anything Starting from Here ------------       ********************   ---------------------------
  TFile* fileEff_L9 = new TFile("RPC155_CorX_L8_HV_9p4kV.root","read");//RPC196_TrigEff_L9_HV_9p6kV.root","read");
  TH2F* efficiency_x_l9 = (TH2F*)fileEff_L9->Get("Efficiency_X_L8");
  //fileEff_L9->Close();

  TFile* fileToT_Pixel = new TFile("fileout_noband_20181226_jj161.root","read");//fileout_noband_20181226_jj180.root","read");//fileout_noband_20210616_jj146.root","read");//fileout_noband_20181226_jj161.root","read");//fileout_noband_20210616_jj157.root","read");//fileout_noband_20210528_jj154.root","read");//fileout_noband_20210616_jj146.root","read");//fileout_noband_20210528_jj139.root","read");//fileout_noband_20210528_172700.root","read");//fileout_noband_20181226_jj114.root","read");//fileout_noband_20210922_jj131.root","read");//fileout_noband_20181226_jj128.root","read");//fileout_noband_20210922_jj118,jj114.root","read");//fileout_noband_20210922_jj110.root","read");//fileout_noband_20210922.root","read"); //jim jim

    TH2F* xstr_vs_ytot_corr[nlayer];
    TH2F* ystr_vs_xtot_corr[nlayer];

    //TH2F* xstr_vs_ytot_corr[nlayer][nstrip/20];
    //TH2F* ystr_vs_xtot_corr[nlayer][nstrip/20];

  char titname[200];
  for(int ij=0;ij<nlayer-2;ij++) {


    sprintf(titname,"xstr_ytwidth_ratio_l%i", ij);
    xstr_vs_ytot_corr[ij] = (TH2F*)fileToT_Pixel->Get(titname);
    sprintf(titname,"ystr_xtwidth_ratio_l%i", ij);
    ystr_vs_xtot_corr[ij] = (TH2F*)fileToT_Pixel->Get(titname);

    /*
    for(int jk=0;jk<nstrip/20;jk++) {
      sprintf(titname,"xstr_ytwidth_ratio_l%i_str%i_to_%i", ij, jk*20,(jk+1)*20);
      xstr_vs_ytot_corr[ij][jk] = (TH2F*)fileToT_Pixel->Get(titname);
      sprintf(titname,"ystr_xtwidth_ratio_l%i_str%i_to_%i", ij, jk*20,(jk+1)*20);
      ystr_vs_xtot_corr[ij][jk] = (TH2F*)fileToT_Pixel->Get(titname);
    }
    */

  }

  for (int ixy=0; ixy<2; ixy++) {
    for (int il=0; il<nlayer; il++) {
      deadchannel[ixy][il].clear();
    }
  }

  static unsigned int mypow_2[32];
  for (int ij=0; ij<32; ij++) {
    mypow_2[ij] = pow(2, ij);
  }

  cout<<int(-.5)<<" "<<int(-1.5)<<" "<<int(0.5)<<" "<<int(1.5)<<endl;

  gStyle->SetPalette(1,0);
  gStyle->SetFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatStyle(1001);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetTitleYOffset(0.05);

  gStyle->SetStatFont(22);        // Times New Roman
  gStyle->SetTextFont(22);        // Times New Roman
  gStyle->SetTitleFont(22,"XYZ"); // Times New Roman
  gStyle->SetLabelFont(22,"XYZ"); // Times New Roman
  gStyle->SetLabelSize(0.07, "XYZ"); // Times New Roman
  gStyle->SetNdivisions(606, "XYZ");
  gStyle->SetPaintTextFormat("6.4f");

  //  gStyle->SetStatFontSize(.015);

  //  gStyle->SetPadGridX(3);
  //  gStyle->SetPadGridY(3);
  //  gStyle->SetGridStyle(3);

  gStyle->SetOptTitle(0);
  gStyle->SetFuncWidth(1);
  gStyle->SetFuncColor(2);
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(101);
  gStyle->SetOptLogy(0);
  gStyle->SetStatW(.18);
  gStyle->SetStatH(.08);
  gStyle->SetPadTopMargin(.001); //0.09
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadLeftMargin(0.001);
  gStyle->SetPadRightMargin(0.001);



  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.12);
  latex.SetTextFont(42);
  latex.SetTextAlign(1); //(31); // align right

#ifdef MONTECARLO
TH2F* tmp_sel_theta_phi[10];
TH2F* sel_theta_phi[10]; //Added to store theta phi distributions for different sets, mainly to obtain acceptance (+selection) efficiency to compare data and MC
TH2F* tmp_sel_theta_phi_8[10][8];
TH2F* sel_theta_phi_8[10][8];
TH2F* sel_theta_phi_16[10][16];
 TH2F* tmp_sel_theta_phi_16[10][16];
#ifdef FASTSIM
  float thgen, phgen, xpgen, ypgen;
  float xtrue[nlayer];
  float ytrue[nlayer];
  float xposa[nlayer];
  float yposa[nlayer];
  int   xystrp[nlayer];
  int   xstrp[nlayer][3];
  int   ystrp[nlayer][3];
  float xtimea[nlayer];
  float ytimea[nlayer];

  //TH2F* tmp_sel_theta_phi[10];
  //TH2F* sel_theta_phi[10]; //Added to store theta phi distributions for different sets, mainly to obtain acceptance (+selection) efficiency to compare data and MC

#else //FULLG4SIM

  static const unsigned int ngenmx=50;
  static const unsigned int nlayermx = 12;
  static const unsigned int nsimhtmx= 1000;
  /*
  unsigned   irun;                // Run number of these events
  unsigned   ievt;                //Event number
  unsigned   ngent;
  float ievt_wt;                //*GMa
  int           intxn_id;               //*GMa


  int   pidin[ngenmx];    //PID of incident particle
  float momin[ngenmx];    //Energy of incident particle
  float thein[ngenmx];    //Initial polar angle of incident particle
  float phiin[ngenmx];     //Initial azimuthal angle of incident particle
  float posxin[ngenmx];   //Initial X-position
  float posyin[ngenmx];     //Initial Y-position
  float poszin[ngenmx];     //Initial Z-position

  // For Simulation output
  unsigned int nsimhtx;  // Number of Hits in X side
  unsigned int nsimhty;  // Number of Hits in Y side

  unsigned int nlayert; // Number of Layers
  unsigned int nstript; // Number of Layers
  int digiXtime[nlayermx]; // Time for each layer. in multiples of 100 ps. for each layer. Time stamp using the earliest signal. X side.
  int digiYtime[nlayermx]; // Time for each layer. in multiples of 100 ps. for each layer. Time stamp using the earliest signal. Y side.
  unsigned long int xdata[nlayermx]; // 32 strips data bit wise for 12 layers. X side.
  unsigned long int ydata[nlayermx]; // 32 strips data bit wise for 12 layers. X side.
  */

  UInt_t    irun;                // Run number of these events
  UInt_t    ievt;                //Event number
  UInt_t    ngent;
  Float_t   ievt_wt;                //*GMa
  Int_t      intxn_id;               //*GMa


  Int_t   pidin[ngenmx];    //PID of incident particle
  Float_t  momin[ngenmx];    //Energy of incident particle
  Float_t  thein[ngenmx];    //Initial polar angle of incident particle
  Float_t  thegen[ngenmx];    //Initial polar angle of incident particle
  Float_t  phiin[ngenmx];     //Initial azimuthal angle of incident particle
  Float_t  posxin[ngenmx];   //Initial X-position
  Float_t  posyin[ngenmx];     //Initial Y-position
  Float_t  poszin[ngenmx];     //Initial Z-position

  // For Simulation output
  UInt_t nsimhtx;  // Number of Hits in X side
  UInt_t nsimhty;  // Number of Hits in Y side

  Int_t nlayert; // Number of Layers
  Int_t nstript; // Number of Layers
  Int_t digiXtime[nlayermx]; // Time for each layer. in multiples of 100 ps. for each layer. Time stamp using the earliest signal. X side.
  Int_t digiYtime[nlayermx]; // Time for each layer. in multiples of 100 ps. for each layer. Time stamp using the earliest signal. Y side.
  ULong64_t xdata[nlayermx]; // 32 strips data bit wise for 12 layers. X side.
  ULong64_t ydata[nlayermx]; // 32 strips data bit wise for 12 layers. X side.
  UInt_t triggerinfoX;
  UInt_t triggerinfoY;  // Number of Hits in X side
 // UInt_t nsimhty;

  // unsigned int xzlay[nsimhtmx];
  // unsigned int xstrip[nsimhtmx];
  // unsigned int yzlay[nsimhtmx];
  // unsigned int ystrip[nsimhtmx];
#endif
#endif

  int EvtIndx,yyindx,xxindx,BID;

  bool posinuse[8][nlayer][nstrip]; //Reject strips in timing, which are not aligned
  bool timeinuse[8][nlayer][nstrip]; //Reject strips in timing, which are not aligned

  for ( int ij=0; ij<8; ij++) {
    for ( int jk=0; jk<nlayer; jk++) {
      for ( int kl=0; kl<nstrip; kl++) {
	posinuse[ij][jk][kl] = timeinuse[ij][jk][kl]=true;
      }
    }
  }
  // To calcualte average X-time and Y-time to synchronise them.
  int ntotxtime=0;
  double totxtime=0;
  int ntotytime=0;
  double totytime=0;
  //  double ytimeshift=0; //included in time_corr_2479_150126.txt

  //  int isequence[nlayer]={0,11,1,10,2,9,3,8,4,7,5,6};
  //  int isequence[nlayer]={6,5,7,4,8,3,9,2,10,1,11,0};
  int isequence[nlayer]={0,1,2,3,4,5,6,7,8,9,10,11};
  int nentrymx=-1;
  int ntotal = 0;
  int nseltot=0;
  //int hvset;

  UInt_t nsec=0;
  int isalign=0;

  char outfile[100];
  char outfilx[100];

  char name[100];
  char title[100];
  char infile[200];
  char datafile[100];
  char rootfiles[100];
  // cout <<"Give the input file name"<<endl;;
  // cin>> rootfiles;
  //  sprintf(rootfiles, "test_2479.log");
  //  sprintf(rootfiles, "test_mc.log");
  if (isOnlyCom) {
    isalign = 0;
  } else {
    //    cout <<"Do you want to align ? yes/no 1/0" <<endl;;
    //    cout <<"for alignment, it can be any number greater than zero"<<endl;
    //    cin>> isalign;
    isalign =3;
  }

  //# of iteration where all layers are included
  // For the time being it is not implemented, but can be used
  // Can be use in other way, first few iteration with combined+individual, then
  // Only combined
  //Total # of iteration all layers + individual layers
//#ifdef MONTECARLO
 // const int nmxiter = 1;
//#else
  const int nmxiter = (isalign>0) ? 5 : 1;// 0;//1; //Less than 12, otherwise change the pad in canvas //jim jim make it 1
//#endif
  const int nlayerit = nlayer; // (isalign >0) ? nlayer : 1;

  int ievt, nhits;
  double xslope, xinters, yslope, yinters, timexslope, timeyslope, timex2slope, timey2slope;
  float txslop, tyslop, timexinters, timeyinters, errtimexslope,errtimeyslope,errtimexinters, errtimeyinters;
  double xchi2, ychi2, xt0chi2, yt0chi2;
  int nxstrip, nystrip,Nx,Ny,nxtime,nytime, ntxyla,tmpnxtime,tmpnytime;
  double zen;
  double xtexter[nlayer];
  double ytexter[nlayer];
  double tmpxtexter;
  double tmpytexter;
  double tmpxexter;
  double tmpyexter;
  double tmpxdev;
  double tmpydev;
  int tmpNx;
  int tmpNy;
  double tmpxchis2;
  double tmpychis2;
  int tmpxpts[nstrip];
  int tmpypts[nstrip];
  double tmpxtime;
  double tmpytime;
  double tmpxtext;
  double tmpytext;
  double tmpxext;
  double tmpyext;
  double tmpxtdev;
  double tmpytdev;
  int tmpxmult;
  int tmpymult;
  double xtpulsewidth;
  double ytpulsewidth;
  double xtpulsewidth_mul[nstrip];
  double ytpulsewidth_mul[nstrip];
  double rawtmpxtime[nstrip];
  double rawtmpytime[nstrip];
  int layoccu;
  double Xpos[nlayer], Xdev[nlayer];
  bool Xusedpos[nlayer];
  double Ypos[nlayer], Ydev[nlayer];
  bool Yusedpos[nlayer];//=new float[nlayer];

  double Xpos1[nlayer], Xdev1[nlayer];
  bool Xusedpos1[nlayer];
  double Ypos1[nlayer], Ydev1[nlayer];
  bool Yusedpos1[nlayer];//=new float[nlayer];

  double Xpos2[nlayer], Xdev2[nlayer];
  bool Xusedpos2[nlayer];
  double Ypos2[nlayer], Ydev2[nlayer];
  bool Yusedpos2[nlayer];//=new float[nlayer];

  double xslope1, xinters1, yslope1, yinters1;
  double xslope2, xinters2, yslope2, yinters2;
  double xchi21, ychi21,xchi22, ychi22;
  int Nx1,Ny1,Nx2,Ny2;
  int xbina;
  bool isfid_l9_eff;

  //int len = strlen(rootfiles);
  int len = strlen(argv[1]);
  //strncpy(outfilx, rootfiles, len-4);
  strncpy(outfilx, argv[1], len-4);
  outfilx[len-4]='\0';
  sprintf(outfilx, "%s_miical_jj1_", outfilx); //jim jim change
  len = strlen(outfilx);
  outfilx[len]='\0';

  //#if ndefined(MONTECARLO) && defined(NOISE_SIM)
#if defined(NOISE_SIM)
  //#ifndef MONTECARLO
  //  int ntot_sim=0;
  int ntot_uncsim=0;
  sprintf(outfile, "%snoise_%i.root", outfilx, isalign);
  TFile* NoisefileOut = new TFile(outfile, "recreate");

  TH1F* n_multall[2];

  TH1F* layer_hitsall[2][nlayer];
  TH1F* layer_multall[2][nlayer];
  TH1F* layer_timeall[2][nlayer];

  TH2F* total_vs_indmulall[2][nlayer];

  TH2F* mult_correlall[2][nlayer];
  TH2F* hist_correlall[2][nlayer];

  //After removal of muon hits
  TH1F* n_mult[2];
  TH1F* layer_hits[2][nlayer];
  TH1F* layer_mult[2][nlayer];
  TH1F* layer_time[2][nlayer];

  TH2F* total_vs_indmul[2][nlayer];

  TH2F* mult_correl[2][nlayer];
  TH2F* hist_correl[2][nlayer];

  for (int ixy=0; ixy<2; ixy++) {
    sprintf(title, "%s", (ixy==0) ? "x" : "y");
    sprintf(name, "%sn_mult", title);
    n_mult[ixy] = new TH1F(name, name, nlayer*nstrip+1, -0.5, nlayer*nstrip+0.5);
    sprintf(name, "%sn_multall", title);
    n_multall[ixy] = new TH1F(name, name, nlayer*nstrip+1, -0.5, nlayer*nstrip+0.5);

    for (int ij = 0; ij<nlayer; ij++) {
      sprintf(name, "%slayer_hitsall_%i", title, ij);
      layer_hitsall[ixy][ij] = new TH1F(name, name, nstrip, -0.5, nstrip-0.5);

      sprintf(name, "%slayer_multall_%i", title, ij);
      layer_multall[ixy][ij] = new TH1F(name, name, nstrip+1, -0.5, nstrip+0.5);

      sprintf(name, "%slayer_timeall_%i", title, ij);
      layer_timeall[ixy][ij] = new TH1F(name, name, 200, 0.0, 200.);

      sprintf(name, "%shist_correlall_l%i", title, ij);
      hist_correlall[ixy][ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(name, "%stotal_vs_indmulall_l%i", title, ij);
      total_vs_indmulall[ixy][ij] = new TH2F(name, name, nlayer*nstrip+1, -0.5, nlayer*nstrip+0.5, nstrip+1, -0.5, nstrip+0.5);

      sprintf(name, "%smult_correlall_l%i", title, ij);
      mult_correlall[ixy][ij] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5, nstrip+1, -0.5, nstrip+0.5);

      sprintf(name, "%slayer_hits_%i", title, ij);
      layer_hits[ixy][ij] = new TH1F(name, name, nstrip, -0.5, nstrip-0.5);

      sprintf(name, "%slayer_mult_%i", title, ij);
      layer_mult[ixy][ij] = new TH1F(name, name, nstrip+1, -0.5, nstrip+0.5);

      sprintf(name, "%slayer_time_%i", title, ij);
      layer_time[ixy][ij] = new TH1F(name, name, 200, 0.0, 200.);

      sprintf(name, "%shist_correl_l%i", title, ij);
      hist_correl[ixy][ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(name, "%stotal_vs_indmul_l%i", title, ij);
      total_vs_indmul[ixy][ij] = new TH2F(name, name, nlayer*nstrip+1, -0.5, nlayer*nstrip+0.5, nstrip+1, -0.5, nstrip+0.5);

      sprintf(name, "%smult_correl_l%i", title, ij);
      mult_correl[ixy][ij] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5, nstrip+1, -0.5, nstrip+0.5);
    }
  }
  //#endif
#endif

  sprintf(outfile, "%s%i.root", outfilx, isalign);
  TFile* fileOut = new TFile(outfile, "recreate");

  //  sprintf(outfile, "cor_%s%i.root", outfilx, isalign);
  // TFile* filecorOut = new TFile(outfile, "recreate");
  Int_t EvNum;
  UShort_t XYhitall[nlayer];
  // UShort_t Yhitall[nlayer];
  ULong64_t XYhitcor;
  ULong64_t XYhitcor2; //good strips to correlate
  // ULong_t Yhitcor[12];



  sprintf(outfile, "%s%i.ps", outfilx, isalign);
  TPostScript ps(outfile,111);
  ps.Range(20,30); //ps.Range(10,20);

  sprintf(outfile, "%s%i.txt", outfilx, isalign);
  ofstream file_out(outfile);

  sprintf(outfile, "%s%i.h", outfilx, isalign);
  ofstream file_out_h(outfile);

  sprintf(outfile, "%s%i_str.txt", outfilx, isalign);
  ofstream file_outstr(outfile);

  TTree* T2 = new TTree("T2", "store"); //("Tree Name","Tree Title")

  //  T2->Branch("Evt",&ievt,"ievt/I");
  T2->Branch("nhits", &nhits, "nhits/I");
  T2->Branch("xslope", &xslope, "xslope/D");
  T2->Branch("xinters", &xinters, "xinters/D");
  T2->Branch("yslope", &yslope, "yslope/D");
  T2->Branch("yinters", &yinters, "yinters/D");

  T2->Branch("xchi2", &xchi2, "xchi2/D");
  T2->Branch("ychi2", &ychi2, "ychi2/D");
  T2->Branch("zen",&zen,"zen/D");


  TTree* T3 = new TTree("T3","store_T3");
  T3->Branch("EvNum",&EvNum,"EvNum/I");
  T3->Branch("XYhitall",XYhitall,"XYhitall[12]/s");
  T3->Branch("XYhitcor",&XYhitcor,"XYhitcor/l");
  T3->Branch("XYhitcor2",&XYhitcor2,"XYhitcor2/l");

  TTree* T4 = new TTree("T4","store_T4");
  T4->Branch("ievt",&ievt,"ievt/I");
  T4->Branch("layoccu",&layoccu,"layoccu/I");

  T4->Branch("tmpxext",&tmpxext,"tmpxext/D");
  T4->Branch("tmpyext",&tmpyext,"tmpyext/D");
  T4->Branch("tmpxdev",&tmpxdev,"tmpxdev/D");
  T4->Branch("tmpydev",&tmpydev,"tmpydev/D");
  T4->Branch("tmpxmult",&tmpxmult,"tmpxmult/I");
  T4->Branch("tmpymult",&tmpymult,"tmpymult/I");
  T4->Branch("tmpxexter",&tmpxexter,"tmpxexter/D");
  T4->Branch("tmpyexter",&tmpyexter,"tmpyexter/D");
  T4->Branch("tmpNx", &tmpNx, "tmpNx/I");
  T4->Branch("tmpNy", &tmpNy, "tmpNy/I");
  T4->Branch("tmpxchis2", &tmpxchis2, "tmpxchis2/D");
  T4->Branch("tmpychis2", &tmpychis2, "tmpychis2/D");
  T4->Branch("tmpxpts",tmpxpts,"tmpxpts[tmpxmult]/I");
  T4->Branch("tmpypts",tmpypts,"tmpypts[tmpymult]/I");
  T2->Branch("xslope", &xslope, "xslope/D");
  T4->Branch("xinters", &xinters, "xinters/D");
  T4->Branch("yslope", &yslope, "yslope/D");
  T4->Branch("yinters", &yinters, "yinters/D");
  T4->Branch("xchi2", &xchi2, "xchi2/D");
  T4->Branch("ychi2", &ychi2, "ychi2/D");
  T4->Branch("zen",&zen,"zen/D");
  if (isTiming) {
    T4->Branch("tmpnxtime", &tmpnxtime, "tmpnxtime/I");
    T4->Branch("tmpnytime", &tmpnytime, "tmpnytime/I");
    T4->Branch("tmpxtdev",&tmpxtdev,"tmpxtdev/D");
    T4->Branch("tmpytdev",&tmpytdev,"tmpytdev/D");
    T4->Branch("xtpulsewidth",&xtpulsewidth,"xtpulsewidth/D");
    T4->Branch("ytpulsewidth",&ytpulsewidth,"ytpulsewidth/D");
    T4->Branch("xtpulsewidth_mul",xtpulsewidth_mul,"xtpulsewidth_mul[4]/D");
    T4->Branch("ytpulsewidth_mul",ytpulsewidth_mul,"ytpulsewidth_mul[4]/D");
    T4->Branch("rawtmpxtime",rawtmpxtime,"rawtmpxtime[tmpxmult]/D");
    T4->Branch("rawtmpytime",rawtmpytime,"rawtmpytime[tmpymult]/D");
    T4->Branch("tmpxtime",&tmpxtime,"tmpxtime/D");
    T4->Branch("tmpytime",&tmpytime,"tmpytime/D");
    T4->Branch("tmpxtext",&tmpxtext,"tmpxtext/D");
    T4->Branch("tmpytext",&tmpytext,"tmpytext/D");


    T2->Branch("txslop", &txslop, "txslop/F");
    T2->Branch("tyslop", &tyslop, "tyslop/F");
    T2->Branch("xt0chi2", &xt0chi2, "xt0chi2/D");
    T2->Branch("yt0chi2", &yt0chi2, "yt0chi2/D");
    T2->Branch("nxtime", &nxtime, "nxtime/I");
    T2->Branch("nytime", &nytime, "nytime/I");
    T2->Branch("ntxyla", &ntxyla, "ntxyla/I");
    T2->Branch("timexslope",&timexslope,"timexslope/D");
    T2->Branch("timeyslope",&timeyslope,"timeyslope/D");
    //    T2->Branch("timexinters",&timexinters,"timexinters/F");
    //    T2->Branch("timeyinters",&timeyinters,"timeyinters/F");
    T2->Branch("errtimexslope",&errtimexslope,"errtimexslope/F");
    //    T2->Branch("errtimeyslope",&errtimeyslope,"errtimeyslope/F");
    //    T2->Branch("errtimexinters",&errtimexinters,"errtimexinters/F");
    //    T2->Branch("errtimeyinters",&errtimeyinters,"errtimeyinters/F");
    // T2->Branch("tmpxtdev",tmpxtdev,"tmpxtdev[12]/D");
    //T2->Branch("tmpytdev",tmpytdev,"tmpytdev[12]/D");
    // T2->Branch("xtpulsewidth",xtpulsewidth,"xtpulsewidth[12]/D");
    // T2->Branch("ytpulsewidth",ytpulsewidth,"ytpulsewidth[12]/D");
    T2->Branch("xtexter", xtexter, "xtexter[12]/D");
    T2->Branch("ytexter", ytexter, "ytexter[12]/D");
  }






#ifdef MONTECARLO
#ifdef FASTSIM

  T2->Branch("thgen", &thgen, "thgen/F");
  T2->Branch("phgen", &phgen, "phgen/F");
  T2->Branch("xpgen", &xpgen, "xpgen/F");
  T2->Branch("ypgen", &ypgen, "ypgen/F");

  T2->Branch("xystrp", xystrp, "xystrp[12]/I");
  T2->Branch("xstrp", xstrp, "xstrp[12][3]/I");
  T2->Branch("ystrp", ystrp, "ystrp[12][3]/I");
  //  T2->Branch("xposa",  xposa,  "xposa[12]/F");
  //  T2->Branch("yposa",  yposa,  "yposa[12]/F");
  //  T2->Branch("xtrue",  xtrue,  "xtrue[12]/F");
  //  T2->Branch("ytrue",  ytrue,  "ytrue[12]/F");
  T2->Branch("xtimea", xtimea, "xtimea[12]/F");
  T2->Branch("ytimea", ytimea, "ytimea[12]/F");

#else  //FULLG4SIM
  T2->Branch("irun",&irun,"irun/i"); //VALGRIND
  T2->Branch("ievt",&ievt,"ievt/i");

  T2->Branch("ngent",&ngent,"ngent/i");
  T2->Branch("pidin",pidin,"pidin[ngent]/I");
  T2->Branch("ievt_wt",&ievt_wt,"ievt_wt/F");
  T2->Branch("intxn_id",&intxn_id,"intxn_id/I");
  T2->Branch("momin",momin,"momin[ngent]/F");
  T2->Branch("thein",thein,"thein[ngent]/F");
  T2->Branch("phiin",phiin,"phiin[ngent]/F");
  T2->Branch("posxin",posxin,"posxin[ngent]/F");
  T2->Branch("posyin",posyin,"posyin[ngent]/F");
  T2->Branch("poszin",poszin,"poszin[ngent]/F");
  T2->Branch("nsimhtx", &nsimhtx, "nsimhtx/i");
  T2->Branch("nsimhty", &nsimhty, "nsimhty/i");
  T2->Branch("nlayert", &nlayert, "nlayert/I");
  T2->Branch("xdata",xdata,"xdata[nlayert]/L");
  T2->Branch("ydata",ydata,"ydata[nlayert]/L");
  T2->Branch("triggerinfoX",&triggerinfoX,"triggerinfoX/i");
  T2->Branch("triggerinfoY",&triggerinfoY,"triggerinfoY/i");

  if (isTiming) {
    T2->Branch("digiXtime", digiXtime, "digiXtime[nlayert]/I");
    T2->Branch("digiYtime", digiYtime, "digiYtime[nlayert]/I");
  }
#endif
#endif



  int narray =90; //75;//60;
  const int nselcrit=16;//6;
  TH1F* costhe[nselcrit];
  costhe[0] = new TH1F("costhe_all", "costhe_all", narray, 0., narray);//filling as dn/d(theta)
  costhe[1] = new TH1F("costhe_trig", "costhe_trig", narray, 0., narray);
  costhe[2] = new TH1F("costhe_selecx", "costhe_selecx", narray, 0., narray);
  costhe[3] = new TH1F("costhe_selecy", "costhe_selecy", narray, 0., narray);
  costhe[4] = new TH1F("costhe_accep", "costhe_accep", narray, 0., narray);
  costhe[5] = new TH1F("costhe_accep_H", "costhe_accep_H", 10, 0.5, 1.0);//10Nov Honda also give 0 to 60 deg in same bin filling as dn/d(costheta)

  // To find the multiple scattering in Layer 5 due to lead block
  const int npixel_theta = 16;
  TH2F* pixel_diff_theta[16];;
  TH1F* theta_diff_1 = new TH1F("theta_diff_1","theta_diff_1",100.,0.,10.);
  TH1F* theta_diff_2 = new TH1F("theta_diff_2","theta_diff_2",100.,0.,10.);

 for(int ij=0;ij<npixel_theta;ij++) {
    sprintf(name,"diff_theta_L5_i%i",ij);
     pixel_diff_theta[ij] = new TH2F(name,name,stripwidth*nstrip/2,0.,96.,stripwidth*nstrip/2,0.,96.);
    }


 TH1F* pixel_scatang1[nstrip/2][nstrip/2] = {{0}};
 TH1F* pixel_scatang[nstrip/2][nstrip/2] = {{0}};
 double val_scatang[nstrip/2][nstrip/2]={{0.}};
 for(int ij =0;ij< nstrip/2 ;ij++) {
   for(int jk=0;jk<nstrip/2;jk++) {
     sprintf(name,"pixel_scatang1_x%i_y%i",ij,jk);
     pixel_scatang1[ij][jk] = new TH1F(name,name,100,0.,10.);
      sprintf(name,"pixel_scatang_x%i_y%i",ij,jk);
       pixel_scatang[ij][jk] = new TH1F(name,name,100,0.,10.);
   }
 }

 TH2F* pixel_scatmean1= new TH2F("scatang_mean1","scatang_mean1",32,-0.5,31.5,32,-0.5,31.5);
 TH2F* pixel_scatmean= new TH2F("scatang_mean","scatang_mean",32,-0.5,31.5,32,-0.5,31.5);
  // 9th July 2016 : Added to test "n" for different criteria
  const int ntrkselcrit=140;
  TH1F* zenithang[nselcrit][ntrkselcrit];
  TH1F* azimuthang[nselcrit][ntrkselcrit];
  TH1F* zenithang_azimuth[nselcrit][ntrkselcrit];
  TH1F* zenithang_azimuth_8[nselcrit][ntrkselcrit];
  TH2F* sel_theta_phi_reco[nselcrit][ntrkselcrit];
  for (int ij=0; ij<nselcrit; ij++) {
    for (int jk=0; jk<ntrkselcrit; jk++) {
      sprintf(name, "sel_theta_phi_reco_%i_%i", ij, jk);
      sel_theta_phi_reco[ij][jk] = new TH2F(name,name,360,-180.,180.,2*narray, 0., narray);
      sprintf(name, "zenithang_%i_%i", ij, jk);
      zenithang[ij][jk] = new TH1F(name, name, 2*narray, 0., narray);
      sprintf(name, "azimuthang_%i_%i", ij, jk);
      azimuthang[ij][jk] = new TH1F(name, name, 36, 0.0, 360.0);
      sprintf(name, "zenithang_azimuth_%i_%i", ij, jk);
      zenithang_azimuth[ij][jk] = new TH1F(name, name, 2*narray, 0., narray);
      sprintf(name, "zenithang_azimuth_8_%i_%i", ij, jk);
      zenithang_azimuth_8[ij][jk] = new TH1F(name, name, 2*narray, 0., narray);

    }
  }

  TH1F* phiang[6];

  phiang[0] = new TH1F("phiang_all", "phiang_all", 36, -pival, pival);
  phiang[1] = new TH1F("phiang_trig", "phiang_trig", 36, -pival, pival);
  phiang[2] = new TH1F("phiang_selecx", "phiang_selecx", 36, -pival, pival);
  phiang[3] = new TH1F("phiang_selecy", "phiang_selecy", 36, -pival, pival);
  phiang[4] = new TH1F("phiang_accep", "phiang_accep", 36, -pival, pival);

  TH2D* strp_xmul[nlayer][2*nmxiter];
  TH2D* strp_ymul[nlayer][2*nmxiter];

  TH2D* strp_xmulsim[nlayer][2*nmxiter];
  TH2D* strp_ymulsim[nlayer][2*nmxiter];
  TH2D* strp_xmulsim_split[nlayer][nsplit][2*nmxiter];
  TH2D* strp_ymulsim_split[nlayer][nsplit][2*nmxiter];

#ifdef ISEFFICIENCY

  TH2D* defefficiency_uncx[nlayer];
  TH2D* defefficiency_uncy[nlayer];

  TH2D* deftriggereffi_x[nlayer];
  TH2D* deftriggereffi_y[nlayer];

  TH1D* efficiency_xallpixel[nlayer][2*nmxiter];
  TH1D* efficiency_yallpixel[nlayer][2*nmxiter];

  TH1D* efficiency_xpixel[nlayer][2*nmxiter];
  TH1D* efficiency_ypixel[nlayer][2*nmxiter];

  TH2D* inefficiency_uncx[nlayer][2*nmxiter];
  TH2D* inefficiency_uncy[nlayer][2*nmxiter];

  TH2D* inefficiency_corx[nlayer][2*nmxiter];
  TH2D* inefficiency_cory[nlayer][2*nmxiter];

  TH2D* inefficiencytrue_corx[nlayer][2*nmxiter];
  TH2D* inefficiencytrue_cory[nlayer][2*nmxiter];


  TH2D* inefficiency_xt[nlayer][2*nmxiter];
  TH2D* inefficiency_yt[nlayer][2*nmxiter];

  TH2D* total_xt[nlayer][2*nmxiter];
  TH2D* total_yt[nlayer][2*nmxiter];

  TH2D* triggereffi_x[nlayer][2*nmxiter];
  TH2D* triggereffi_y[nlayer][2*nmxiter];

  TH2D* fine_triggereffi_x[nlayer][2*nmxiter];
  TH2D* fine_triggereffi_y[nlayer][2*nmxiter];


  TH2D* triggereffi_xevt[nlayer][2*nmxiter];
  TH2D* triggereffi_yevt[nlayer][2*nmxiter];

  TH2D* totalentry[nlayer][2*nmxiter];
  TH2D* fine_totalentry[nlayer][2*nmxiter];

  TH2D* difefficiency_uncx[nlayer][2*nmxiter];
  TH2D* difefficiency_uncy[nlayer][2*nmxiter];

  TH2D* difefficiency_xt[nlayer][2*nmxiter];
  TH2D* difefficiency_yt[nlayer][2*nmxiter];

  TH2D* diftriggereffi_x[nlayer][2*nmxiter];
  TH2D* diftriggereffi_y[nlayer][2*nmxiter];

  TH2D* diftriggereffi_xevt[nlayer][2*nmxiter];
  TH2D* diftriggereffi_yevt[nlayer][2*nmxiter];

  //corr efficiencies for Monte Carlo

#endif

  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<2*nmxiter; jk++) {
      sprintf(title, "strp_xmul_l%i_i%i", ij, jk);
      strp_xmul[ij][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);

      sprintf(title, "strp_ymul_l%i_i%i", ij, jk);
      strp_ymul[ij][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);

      sprintf(title, "strp_xmulsim_l%i_i%i", ij, jk);
      strp_xmulsim[ij][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);

      sprintf(title, "strp_ymulsim_l%i_i%i", ij, jk);
      strp_ymulsim[ij][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);
      for(int um=0;um<nsplit;um++) {
	sprintf(title, "strp_xmulsim_l%i_splt%i_i%i", ij, um,jk);
      strp_xmulsim_split[ij][um][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);

      sprintf(title, "strp_ymulsim_l%i_splt%i_i%i", ij,um, jk);
      strp_ymulsim_split[ij][um][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);
      }

    }



#ifdef ISEFFICIENCY

    sprintf(title, "defefficiency_uncx_l%i", ij);
    defefficiency_uncx[ij]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "defefficiency_uncy_l%i", ij);
    defefficiency_uncy[ij]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "deftriggereffi_x_l%i", ij);
    deftriggereffi_x[ij]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "deftriggereffi_y_l%i", ij);
    deftriggereffi_y[ij]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    for (int jk=0; jk<2*nmxiter; jk++) {
      sprintf(title, "totalentry_l%i_i%i", ij, jk);
      totalentry[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "fine_totalentry_l%i_i%i", ij, jk);
      fine_totalentry[ij][jk]=new TH2D(title, title, 4*nstrip, -0.5, nstrip-0.5, 4*nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncx_l%i_i%i", ij, jk);
      inefficiency_uncx[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncy_l%i_i%i", ij, jk);
      inefficiency_uncy[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_corx_l%i_i%i", ij, jk);
      inefficiency_corx[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_cory_l%i_i%i", ij, jk);
      inefficiency_cory[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiencytrue_corx_l%i_i%i", ij, jk);
      inefficiencytrue_corx[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiencytrue_cory_l%i_i%i", ij, jk);
      inefficiencytrue_cory[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "efficiency_xpixel_l%i_i%i", ij, jk);
      efficiency_xpixel[ij][jk]=new TH1D(title, title, 60, -0.5, 0.5);

      sprintf(title, "efficiency_ypixel_l%i_i%i", ij, jk);
      efficiency_ypixel[ij][jk]=new TH1D(title, title, 60, -0.5, 0.5);

      sprintf(title, "efficiency_xallpixel_l%i_i%i", ij, jk);
      efficiency_xallpixel[ij][jk]=new TH1D(title, title, 60, -0.5, 0.5);

      sprintf(title, "efficiency_yallpixel_l%i_i%i", ij, jk);
      efficiency_yallpixel[ij][jk]=new TH1D(title, title, 60, -0.5, 0.5);

      sprintf(title, "inefficiency_xt_l%i_i%i", ij, jk);
      inefficiency_xt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_yt_l%i_i%i", ij, jk);
      inefficiency_yt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "total_xt_l%i_i%i", ij, jk);
      total_xt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "total_yt_l%i_i%i", ij, jk);
      total_yt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_x_l%i_i%i", ij, jk);
      triggereffi_x[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_y_l%i_i%i", ij, jk);
      triggereffi_y[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "fine_triggereffi_x_l%i_i%i", ij, jk);
      fine_triggereffi_x[ij][jk]=new TH2D(title, title, 4*nstrip, -0.5, nstrip-0.5, 4*nstrip, -0.5, nstrip-0.5);

      sprintf(title, "fine_triggereffi_y_l%i_i%i", ij, jk);
      fine_triggereffi_y[ij][jk]=new TH2D(title, title, 4*nstrip, -0.5, nstrip-0.5, 4*nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_xevt_l%i_i%i", ij, jk);
      triggereffi_xevt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_yevt_l%i_i%i", ij, jk);
      triggereffi_yevt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "difefficiency_uncx_l%i_i%i", ij, jk);
      difefficiency_uncx[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "difefficiency_uncy_l%i_i%i", ij, jk);
      difefficiency_uncy[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "difefficiency_xt_l%i_i%i", ij, jk);
      difefficiency_xt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "difefficiency_yt_l%i_i%i", ij, jk);
      difefficiency_yt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "diftriggereffi_x_l%i_i%i", ij, jk);
      diftriggereffi_x[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "diftriggereffi_y_l%i_i%i", ij, jk);
      diftriggereffi_y[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "diftriggereffi_xevt_l%i_i%i", ij, jk);
      diftriggereffi_xevt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "diftriggereffi_yevt_l%i_i%i", ij, jk);
      diftriggereffi_yevt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    }
#endif
  }

  TH2F* time_layer = new TH2F("time_layer", "time_layer", 2*nTDCpLayer*nlayer+1, -0.5, 2*nTDCpLayer*nlayer+0.5, 140, -0.5, 1399.5);

  TH1F* rawtimex = new TH1F("time_rawx", "time_rawx",1800,0.,1800.);
  TH1F* rawtimey = new TH1F("time_rawy", "time_rawy",1800,0.,1800.);
  TH1F* xlay_timediff[11];
  TH1F* ylay_timediff[11];
  TH1F* xlay_timediff_l1[6];
  TH1F* xlay_timediff_l4[5];
  TH1F* ylay_timediff_l1[6];
  TH1F* ylay_timediff_l4[5];
  for(int ij=0;ij<nlayer-6;ij++) {
     sprintf(name,"xlayer_timediff_l1_%i%i",1,ij+4);
     xlay_timediff_l1[ij] = new TH1F(name,name,500,-49.5,49.5);
     sprintf(name,"ylayer_timediff_l1_%i%i",1,ij+4);
     ylay_timediff_l1[ij] = new TH1F(name,name,500,-49.5,49.5);
  }
 for(int ij=0;ij<nlayer-7;ij++) {
     sprintf(name,"xlayer_timediff_l4_%i%i",4,ij+5);
     xlay_timediff_l4[ij] = new TH1F(name,name,500,-49.5,49.5);
     sprintf(name,"ylayer_timediff_l4_%i%i",4,ij+5);
     ylay_timediff_l4[ij] = new TH1F(name,name,500,-49.5,49.5);
  }

  for(int ij=0;ij<nlayer-1;ij++){
    sprintf(name,"xlayer_timediff_%i%i",ij,ij+1);
    xlay_timediff[ij] = new TH1F(name,name,500,-49.5,49.5);
    sprintf(name,"ylayer_timediff_%i%i",ij,ij+1);
    ylay_timediff[ij] = new TH1F(name,name,500,-49.5,49.5);
  }
  TH2F* time_layerstrip = new TH2F("time_layerstrip", "time_layerstrip", 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5, 8000, -0.5, 7999.5);
  TH2F* time_layerstripcpy = new TH2F("time_layerstripcpy", "time_layerstripcpy", 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5, 20000, -0.5, 19999.5);
  TH1F* h_time_layerprx[nlayer];
  TH1F* h_time_layerpry[nlayer];
  TH1F* h_time_layerprx1[nlayer];
  TH1F* h_time_layerpry1[nlayer];
  // = new TH1F("time_layer_1d","time_layer_1d",8000,-0.5,7999.5);

  TH1F* xtimecpy[nlayer][3];
  TH1F* ytimecpy[nlayer][3];
  TH1F* xytimecpy[nlayer][3];
  TH1F* xytimecpy_nTermSplit[nlayer][nTermResSplit][3];
  for(int ij=0;ij<nlayer;ij++) {

    for(int jk=0;jk<3;jk++) {
      sprintf(title,"xtime_l%i_c%i",ij,jk);
      xtimecpy[ij][jk] = new TH1F(title, title, 800,800,1200);
      sprintf(title,"ytime_l%i_c%i",ij,jk);
      ytimecpy[ij][jk] = new TH1F(title, title, 800,800,1200);
      sprintf(title,"xytime_l%i_c%i",ij,jk);
      xytimecpy[ij][jk] = new TH1F(title, title, 800,-20,20);
      for(int nT=0;nT<nTermResSplit;nT++) {
	sprintf(title,"xytime_l%i_str%i_to_%i_c%i",ij, nT*20,(nT+1)*20,jk);
	xytimecpy_nTermSplit[ij][nT][jk] = new TH1F(title, title, 800,-20,20);
    }
    }
  }
  TH2F* xlayer_timeoccu[nlayer];
  TH2F* ylayer_timeoccu[nlayer];
  TH2F* xylayer_timeoccu[nlayer];
  for(int ij=0;ij<nlayer;ij++) {
    sprintf(name,"xlayer_timeoccu_%i",ij);
    xlayer_timeoccu[ij] = new TH2F(name,name,nstrip,-0.5,nstrip-0.5,nstrip,-0.5,nstrip-0.5);
    sprintf(name,"ylayer_timeoccu_%i",ij);
    ylayer_timeoccu[ij] = new TH2F(name,name,nstrip,-0.5,nstrip-0.5,nstrip,-0.5,nstrip-0.5);
    sprintf(name,"xylayer_timeoccu_%i",ij);
    xylayer_timeoccu[ij] = new TH2F(name,name,nstrip,-0.5,nstrip-0.5,nstrip,-0.5,nstrip-0.5);
    sprintf(name,"xtime_layer_1d_%i",ij);
    h_time_layerprx[ij] = new TH1F(name,name,8000,-0.5,7999.5);
    sprintf(name,"ytime_layer_1d_%i",ij);
    h_time_layerpry[ij] = new TH1F(name,name,8000,-0.5,7999.5);
    sprintf(name,"xtime_layer1_1d_%i",ij);
    h_time_layerprx1[ij] = new TH1F(name,name,25000,-0.5,249999.5);
    sprintf(name,"ytime_layer1_1d_%i",ij);
    h_time_layerpry1[ij] = new TH1F(name,name,25000,-0.5,249999.5);
  }
  TH1F* xlayer_alloccu=new TH1F("xlayer_alloccu","xlayer_alloccu",nlayer*nstrip,-0.5, nlayer*nstrip-0.5);
  TH1F* ylayer_alloccu=new TH1F("ylayer_alloccu","ylayer_alloccu",nlayer*nstrip,-0.5, nlayer*nstrip-0.5);

  TH1F* xlayer_alloccusel=new TH1F("xlayer_alloccusel","xlayer_alloccusel",nlayer*nstrip,-0.5, nlayer*nstrip-0.5);
  TH1F* ylayer_alloccusel=new TH1F("ylayer_alloccusel","ylayer_alloccusel",nlayer*nstrip,-0.5, nlayer*nstrip-0.5);

  int nbin = 400000;

  TH1F* trigrate = new TH1F("trigrate", "Trigger rate (Hz)", 300, -0.5, 299.5);
  TH1F* muposxrate = new TH1F("muposxrate", "Muon (xpos) rate (Hz)", 300, -0.5, 299.5);
  TH1F* muposyrate = new TH1F("muposyrate", "Muon (ypos) rate (Hz)", 300, -0.5, 299.5);
  TH1F* mutimexrate = new TH1F("mutimexrate", "Muon (xtime) rate (Hz)", 300, -0.5, 299.5);
  TH1F* mutimeyrate = new TH1F("mutimeyrate", "Muon (ytime) rate (Hz)", 300, -0.5, 299.5);

  TH1F* trigratex = new TH1F("trigratex", "Trigger ratex (Hz)", nbin, -0.5, nbin-0.5);
  TH1F* muposxratex = new TH1F("muposxratex", "Muon (xpos) rate (Hz)", nbin, -0.5, nbin-0.5);
  TH1F* muposyratex = new TH1F("muposyratex", "Muon (ypos) rate (Hz)", nbin, -0.5, nbin-0.5);
  TH1F* mutimexratex = new TH1F("mutimexratex", "Muon (xtime) rate (Hz)", nbin, -0.5, nbin-0.5);
  TH1F* mutimeyratex = new TH1F("mutimeyratex", "Muon (ytime) rate (Hz)", nbin, -0.5, nbin-0.5);


  TH1F* xlayer_occu[nlayer];
  TH1F* ylayer_occu[nlayer];

  TH1F* xlayer_seloccu[nlayer];
  TH1F* ylayer_seloccu[nlayer];

  TH1F* xlayer_sel2occu[nlayer];
  TH1F* ylayer_sel2occu[nlayer];

  TH1F* xlayer_noiseoccu[nlayer];
  TH1F* ylayer_noiseoccu[nlayer];

  TH2F* raw_occu[nlayer];
  TH2F* raw_seloccu[nlayer];
  TH2F* raw_noiseoccu[nlayer];

  for (int ij=0; ij<nlayer; ij++) {
    sprintf(title, "xlayer_occu_l%i", ij);
    xlayer_occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "ylayer_occu_l%i", ij);
    ylayer_occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "raw_occu_l%i", ij);
    raw_occu[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "raw_seloccu_l%i", ij);
    raw_seloccu[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "raw_noiseoccu_l%i", ij);
    raw_noiseoccu[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);


     sprintf(title, "xlayer_seloccu_l%i", ij);
    xlayer_seloccu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "ylayer_seloccu_l%i", ij);
    ylayer_seloccu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "xlayer_sel2occu_l%i", ij);
    xlayer_sel2occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "ylayer_sel2occu_l%i", ij);
    ylayer_sel2occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "xlayer_noiseoccu_l%i", ij);
    xlayer_noiseoccu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "ylayer_noiseoccu_l%i", ij);
    ylayer_noiseoccu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
  }

  TH1F* xlayer_mult[nlayer];
  TH1F* ylayer_mult[nlayer];

  TH1F* xlayer_allmult[nlayer];
  TH1F* ylayer_allmult[nlayer];

  TH1F* xlayer_allmumult[nlayer];
  TH1F* ylayer_allmumult[nlayer];

  TH1F* xlayer_allmutimemult[nlayer];
  TH1F* ylayer_allmutimemult[nlayer];

  const int nEnv=4; //Multiplicity in different env
  const char* envname[nEnv] = {"BothXY", "onlyX", "OnlyY", "neither"};
  const int nMuselec=8;
  const char* muselname[nMuselec] = {"Xmu", "Ymu", "XYmu", "X!Ymu", "!X!ymu", "!Xmu", "!Ymu"};
  TH2F* mult2d_selecxy[nlayer][nMuselec][nEnv];

  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<nMuselec; jk++) {
      for (int kl=0; kl<nEnv; kl++) {
	sprintf(name, "mult2d_selecxy_l%i_m%i_e%i", ij, jk, kl);
	sprintf(title, "mult2d_selecxy_l%i %s %s", ij, muselname[jk], envname[kl]);
	mult2d_selecxy[ij][jk][kl] = new TH2F(name, title, nstrip+1, -0.5, nstrip+0.5, nstrip+1, -0.5, nstrip+0.5);
      }
    }


    sprintf(title, "xlayer_mult_l%i", ij);
    xlayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "ylayer_mult_l%i", ij);
    ylayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "xlayer_allmult_l%i", ij);
    xlayer_allmult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "ylayer_allmult_l%i", ij);
    ylayer_allmult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "xlayer_allmumult_l%i", ij);
    xlayer_allmumult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "ylayer_allmumult_l%i", ij);
    ylayer_allmumult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "xlayer_allmutimemult_l%i", ij);
    xlayer_allmutimemult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "ylayer_allmutimemult_l%i", ij);
    ylayer_allmutimemult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

  }

  TH1F* xstrip_mult = new TH1F("xstrip_mult", "xstrip_mult", 6*nlayer, -0.5,  6*nlayer-0.5);
  TH1F* ystrip_mult = new TH1F("ystrip_mult", "ystrip_mult", 6*nlayer, -0.5,  6*nlayer-0.5);

  TH2F* dir_cxchi = new TH2F("dir_cxchi", "dir_cxchi", 100, 0., 100., 200, -20., 20.);
  TH2F* dir_cychi = new TH2F("dir_cychi", "dir_cychi", 100, 0., 100., 200, -20., 20.);

  TH2F* dir_cxy = new TH2F("dir_cxy", "dir_cxy", 100, -2., 2., 100, -2., 2.);

  TH1F* timex_shift[nlayer][nmxiter+1];
  TH1F* timey_shift[nlayer][nmxiter+1];

  TH2F* timex_2dshift[nlayer][nmxiter+1];
  TH2F* timey_2dshift[nlayer][nmxiter+1];

  //  TH1F* time_xstshift[nlayer][nstrip][2*nmxiter];
  //  TH1F* time_ystshift[nlayer][nstrip][2*nmxiter];

  TH2F* timex_fy[nlayer];
  TH2F* timey_fx[nlayer];

  TH2F* timex_fy2[nlayer];
  TH2F* timey_fx2[nlayer];

  TProfile* timex_pfy[nlayer];
  TProfile* timey_pfx[nlayer];

  TProfile* timex_pfy2[nlayer];
  TProfile* timey_pfx2[nlayer];
  if (isTiming) {
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<=nmxiter; jk++) {
	sprintf(title, "timex_shift_l%i_i%i", ij, jk);
	timex_shift[ij][jk] = new TH1F(title, title, 120, -30.0, 30.0);

	sprintf(title, "timey_shift_l%i_i%i", ij, jk);
	timey_shift[ij][jk] = new TH1F(title, title, 120, -30.0, 30.0);

	sprintf(title, "timex_2dshift_l%i_i%i", ij, jk);
	timex_2dshift[ij][jk] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 120, -30.0, 30.0);

	sprintf(title, "timey_2dshift_l%i_i%i", ij, jk);
	timey_2dshift[ij][jk] = new TH2F(title, title,  nstrip, -0.5, nstrip-0.5, 120, -30.0, 30.0);
      }


      sprintf(title, "timex_fy_l%i", ij);
      timex_fy[ij] = new TH2F(title, title, nstrip*nstrip, -0.5, nstrip*nstrip-0.5, 120, -12., 12.);
      sprintf(title, "timey_fx_l%i", ij);
      timey_fx[ij] = new TH2F(title, title, nstrip*nstrip, -0.5, nstrip*nstrip-0.5, 120, -12., 12.);
      sprintf(title, "timex_fy2_l%i", ij);
      timex_fy2[ij] = new TH2F(title, title, nstrip*nstrip, -0.5, nstrip*nstrip-0.5, 120, -12., 12.);
      sprintf(title, "timey_fx2_l%i", ij);
      timey_fx2[ij] = new TH2F(title, title, nstrip*nstrip, -0.5, nstrip*nstrip-0.5, 120, -12., 12.);

      sprintf(title, "timex_pfy_l%i", ij);
      timex_pfy[ij] = new TProfile(title, title, nstrip*(nstrip+1), -0.5, nstrip*(nstrip+1)-0.5, -12., 12.);
      sprintf(title, "timey_pfx_l%i", ij);
      timey_pfx[ij] = new TProfile(title, title, nstrip*(nstrip+1), -0.5, nstrip*(nstrip+1)-0.5, -12., 12.);
      sprintf(title, "timex_pfy2_l%i", ij);
      timex_pfy2[ij] = new TProfile(title, title, nstrip*(nstrip+1), -0.5, nstrip*(nstrip+1)-0.5, -12., 12.);
      sprintf(title, "timey_pfx2_l%i", ij);
      timey_pfx2[ij] = new TProfile(title, title, nstrip*(nstrip+1), -0.5, nstrip*(nstrip+1)-0.5, -12., 12.);
    }
  }

  TH1F* xlayer_reso[nlayer][2*nmxiter];
  TH1F* ylayer_reso[nlayer][2*nmxiter];

  TH1F* xlayer_reso_nTermSplit[nlayer][nTermResSplit][2*nmxiter];
  TH1F* ylayer_reso_nTermSplit[nlayer][nTermResSplit][2*nmxiter];

  TH1F* xlayer_reso_split[nlayer][nsplit][2*nmxiter]; //umesh
  TH1F* ylayer_reso_split[nlayer][nsplit][2*nmxiter];

  TH2F* xylayer_reso[nlayer][2*nmxiter];
  bool is_xyreso[nlayer][2*nmxiter];

  TH1F* xlayer_reso_mul[nlayer][2*nmxiter][nmxhits];
  TH1F* ylayer_reso_mul[nlayer][2*nmxiter][nmxhits];
  TH1F* xlayer_reso_mul_split[nlayer][nsplit][2*nmxiter][nmxhits];
  TH1F* ylayer_reso_mul_split[nlayer][nsplit][2*nmxiter][nmxhits];

  TH1F* xlayer_exterr[nlayer][2*nmxiter];
  TH1F* ylayer_exterr[nlayer][2*nmxiter];

  TH1F* xlayer_exterr_nTermSplit[nlayer][nTermResSplit][2*nmxiter];
  TH1F* ylayer_exterr_nTermSplit[nlayer][nTermResSplit][2*nmxiter];

  // TH2F* xstrip_xdev[nlayer][2*nmxiter];//
  // TH2F* ystrip_ydev[nlayer][2*nmxiter];//

  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<2*nmxiter; jk++) {
      double xposmx=(jk<nmxiter) ? 5.9 : 6.0;
      sprintf(title, "xlayer_reso_l%i_i%i", ij, jk);
      xlayer_reso[ij][jk]=new TH1F(title, title, 150, -xposmx, xposmx);
      sprintf(title, "ylayer_reso_l%i_i%i", ij, jk);
      ylayer_reso[ij][jk]=new TH1F(title, title, 150, -xposmx, xposmx);
      for(int nT=0;nT<nTermResSplit;nT++) {
	//	double xposmx=(jk<nmxiter) ? 5.9 : 6.0;
	sprintf(title, "xlayer_reso_l%i_str%i_to_%i_i%i", ij, nT*20,(nT+1)*20, jk);
	xlayer_reso_nTermSplit[ij][nT][jk]=new TH1F(title, title, 150, -xposmx, xposmx);
	sprintf(title, "ylayer_reso_l%i_str%i_to_%i_i%i", ij, nT*20,(nT+1)*20, jk);
	ylayer_reso_nTermSplit[ij][nT][jk]=new TH1F(title, title, 150, -xposmx, xposmx);
	sprintf(title, "xlayer_exterr_l%i_str%i_to_%i_i%i", ij, nT*20,(nT+1)*20, jk);
	xlayer_exterr_nTermSplit[ij][nT][jk]=new TH1F(title, title, 120, 0.0, 1.0);
	sprintf(title, "ylayer_exterr_l%i_str%i_to_%i_i%i", ij, nT*20,(nT+1)*20, jk);
	ylayer_exterr_nTermSplit[ij][nT][jk]=new TH1F(title, title, 120, 0.0, 1.0);
      }

      for(int um = 0;um<nsplit;um++) {
	sprintf(title, "xlayer_reso_l%i_splt%i_i%i", ij,um, jk);
	xlayer_reso_split[ij][um][jk]=new TH1F(title, title, 150, -xposmx, xposmx);
	sprintf(title, "ylayer_reso_l%i_splt%i_i%i", ij,um, jk);
	ylayer_reso_split[ij][um][jk]=new TH1F(title, title, 150, -xposmx, xposmx);
      }
      sprintf(title, "xylayer_reso_l%i_i%i", ij, jk);
      xylayer_reso[ij][jk]=new TH2F(title, title, 90, -xposmx, xposmx, 90, -xposmx, xposmx);

      sprintf(title, "xlayer_exterr_l%i_i%i", ij, jk);
      xlayer_exterr[ij][jk]=new TH1F(title, title, 120, 0.0, 1.0);
      sprintf(title, "ylayer_exterr_l%i_i%i", ij, jk);
      ylayer_exterr[ij][jk]=new TH1F(title, title, 120, 0.0, 1.0);


      //   sprintf(name, "xstrip_xdev_l%i_i%i", ij,jk);
      // xstrip_xdev[ij][jk] = new TH2F(name,name,nstrip,-0.5,nstrip-0.5,1000,-64.0,64.0);
      // sprintf(name, "ystrip_ydev_l%i_i%i", ij,jk);
      // ystrip_ydev[ij][jk] = new TH2F(name,name,nstrip,-0.5,nstrip-0.5,1000 ,-64.0,64.0);

      for (int kl=0; kl<nmxhits; kl++) {
	sprintf(title, "xlayer_reso_l%i_i%i_mul%i", ij, jk, kl+1);
	xlayer_reso_mul[ij][jk][kl]=new TH1F(title, title, 150, -xposmx, xposmx);
	sprintf(title, "ylayer_reso_l%i_i%i_mul%i", ij, jk, kl+1);
	ylayer_reso_mul[ij][jk][kl]=new TH1F(title, title, 150, -xposmx, xposmx);

      }

      for (int kl=0; kl<nmxhits; kl++) {
	for(int um=0;um<nsplit;um++) {
	  sprintf(title, "xlayer_reso_l%i_splt%i_i%i_mul%i", ij,um, jk, kl+1);
	  xlayer_reso_mul_split[ij][um][jk][kl]=new TH1F(title, title, 150, -xposmx, xposmx);
	  sprintf(title, "ylayer_reso_l%i_splt%i_i%i_mul%i", ij,um, jk, kl+1);
	  ylayer_reso_mul_split[ij][um][jk][kl]=new TH1F(title, title, 150, -xposmx, xposmx);
	}
      }
    }
  }

  const int nstr_posmx=nmxhits; //3; //This could be same as nmxhits;
  TH2F* xstr_xdev[nlayer][nmxiter][nstr_posmx];
  TH2F* ystr_xdev[nlayer][nmxiter][nstr_posmx];
  TH2F* xstr_ydev[nlayer][nmxiter][nstr_posmx];
  TH2F* ystr_ydev[nlayer][nmxiter][nstr_posmx];
  TH2F* xstr_xdev_occu[nlayer][nmxiter][nstr_posmx];
  TH2F* ystr_ydev_occu[nlayer][nmxiter][nstr_posmx];
  TProfile* prof_xstr_xdev[nlayer][nmxiter];
  TProfile* prof_ystr_xdev[nlayer][nmxiter];
  TProfile* prof_xstr_ydev[nlayer][nmxiter];
  TProfile* prof_ystr_ydev[nlayer][nmxiter];

  TH2F* xstr_xtdev[nlayer][nmxiter];
  TH2F* ystr_xtdev[nlayer][nmxiter];
  TH2F* xstr_ytdev[nlayer][nmxiter];
  TH2F* ystr_ytdev[nlayer][nmxiter];




  TH2F* nxystr_xtdev[nlayer][nmxiter+6];
  TH2F* xystr_xtdev[nlayer][nmxiter+6]; // + rawytime, rawytime1, timesy & ytimewo, rawytime2 & ryawtime3

  TH2F* nxystr_ytdev[nlayer][nmxiter+6];
  TH2F* xystr_ytdev[nlayer][nmxiter+6];

  TH2F* nxypos_xytdev[8];
  for (int ij=0; ij<8; ij++) {
    sprintf(name, "nxypos_xytdev_%i", ij);
    nxypos_xytdev[ij] = new TH2F(name, name, 2*(nstrip+2), -1.5, nstrip+0.5, 2*(nstrip+2), -1.5, nstrip+0.5);
  }
  TH2F* xtdev_ytdev[4];
  for (int ij=0; ij<4; ij++) {
    sprintf(name, "xtdev_ytdev_%i", ij);
    xtdev_ytdev[ij] = new TH2F(name, name, 120, -30., 30., 120, -30., 30.);
  }



  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<nmxiter; jk++) {
      for (int kl=0; kl<nstr_posmx; kl++) {
	sprintf(name, "xstr_xdev_l%i_i%i_mul%i", ij, jk, kl+1);
	xstr_xdev[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
	sprintf(name, "ystr_xdev_l%i_i%i_mul%i", ij, jk, kl+1);
	ystr_xdev[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
	sprintf(name, "xstr_xdev_occu_l%i_i%i_mul%i", ij, jk, kl+1);
	xstr_xdev_occu[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
	sprintf(name, "ystr_ydev_occu_l%i_i%i_mul%i", ij, jk, kl+1);
	ystr_ydev_occu[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
	sprintf(name, "xstr_ydev_l%i_i%i_mul%i", ij, jk, kl+1);
	xstr_ydev[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
	sprintf(name, "ystr_ydev_l%i_i%i_mul%i", ij, jk, kl+1);
	ystr_ydev[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
      }

      sprintf(name, "prof_xstr_xdev_l%i_i%i", ij, jk);
      prof_xstr_xdev[ij][jk] = new TProfile(name, name, nstrip, -0.5, nstrip-0.5, -1.5, 1.5);
      sprintf(name, "prof_ystr_xdev_l%i_i%i", ij, jk);
      prof_ystr_xdev[ij][jk] = new TProfile(name, name, nstrip, -0.5, nstrip-0.5, -1.5, 1.5);
      sprintf(name, "prof_xstr_ydev_l%i_i%i", ij, jk);
      prof_xstr_ydev[ij][jk] = new TProfile(name, name, nstrip, -0.5, nstrip-0.5, -1.5, 1.5);
      sprintf(name, "prof_ystr_ydev_l%i_i%i", ij, jk);
      prof_ystr_ydev[ij][jk] = new TProfile(name, name, nstrip, -0.5, nstrip-0.5, -1.5, 1.5);
    }
  }

  if (isTiming) {
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<nmxiter; jk++) {
	sprintf(name, "xstr_xtdev_l%i_i%i", ij, jk);
	xstr_xtdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 120, -18., 18.);
	sprintf(name, "ystr_xtdev_l%i_i%i", ij, jk);
	ystr_xtdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 120, -18., 18.);
	sprintf(name, "xstr_ytdev_l%i_i%i", ij, jk);
	xstr_ytdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 120, -18., 18.);
	sprintf(name, "ystr_ytdev_l%i_i%i", ij, jk);
	ystr_ytdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 120, -18., 18.);
      }

      for (int jk=0; jk<nmxiter+6; jk++) {
	sprintf(name, "nxystr_xtdev_l%i_i%i", ij, jk);
	nxystr_xtdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

	sprintf(name, "nxystr_ytdev_l%i_i%i", ij, jk);
	nxystr_ytdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

	sprintf(name, "xystr_xtdev_l%i_i%i", ij, jk);
	xystr_xtdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

	sprintf(name, "xystr_ytdev_l%i_i%i", ij, jk);
	xystr_ytdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      }
    }
  }

  TH2F* shift_pos;
  TH2F* rms_pos;
  if (isalign>0) {
    shift_pos = new TH2F("shift_pos", "shift_pos", 2*nlayer, -0.5, 2*nlayer-0.5, nmxiter, -0.5, nmxiter-0.5);
    rms_pos = new TH2F("rms_pos", "rms_pos", 2*nlayer, -0.5, 2*nlayer-0.5, nmxiter, -0.5, nmxiter-0.5);
  }

  TH1F* passed_strip[nmxiter];
  for (int ij=0; ij<nmxiter; ij++) {
    sprintf(title, "passed_strip_i%i", ij);
    passed_strip[ij] = new TH1F(title, title, 4*(nlayer+2), -0.5, 4*(nlayer+2)-0.5);
  }

  TH1F* h_chisqx = new TH1F("chisqx", "chisqx", 120, 0.0, 150.0);
  TH1F* h_reduchisqx = new TH1F("reduchisqx", "reduced chisqx", 90, 0.0, 30.0/*45.0*/);
  TH1F* h_chisqy = new TH1F("chi2y", "chisqy", 120, 0.0, 150.0);
  TH1F* h_reduchisqy = new TH1F("reduchisqy", "reduced chisqy", 90, 0.0, 30.0/*45.0*/);
  TH1F* h_xndf = new TH1F("xndf", "xndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_yndf = new TH1F("yndf", "yndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_xprob = new TH1F("xprob", "xprob", 120, 0.0, 1.0);
  TH1F* h_yprob = new TH1F("yprob", "yprob", 120, 0.0, 1.0);

  TH1F* h_tmp1xndf = new TH1F("tmp1xndf", "tmp1xndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_tmp2xndf = new TH1F("tmp2xndf", "tmp2xndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_tmp3xndf = new TH1F("tmp3xndf", "tmp3xndf", nlayer, 0.5, nlayer+0.5);


  TH1F* h_tchisqx = new TH1F("tchisqx", "tchisqx", 120, 0.0, 90.0);
  TH1F* h_treduchisqx = new TH1F("treduchisqx", "reduced tchisqx", 90, 0.0, 30.0);
  TH1F* h_tchisqy = new TH1F("tchi2y", "tchisqy", 120, 0.0, 90.0);
  TH1F* h_treduchisqy = new TH1F("treduchisqy", "reduced tchisqy", 90, 0.0, 30.0);
  TH1F* h_txndf = new TH1F("txndf", "txndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_tyndf = new TH1F("tyndf", "tyndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_xtprob = new TH1F("xtprob", "xtprob", 120, 0.0, 1.0);
  TH1F* h_ytprob = new TH1F("ytprob", "ytprob", 120, 0.0, 1.0);

  TH2F* xtprob_vs_tchisqx = new TH2F("xtprob_vs_tchisqx", "xtprob_vs_tchisqx", 500, 0.0, 1.0, 500, 0.0, 90.0);
  TH2F* ytprob_vs_tchisqy = new TH2F("ytprob_vs_tchisqy", "ytprob_vs_tchisqy", 500, 0.0, 1.0, 500, 0.0, 90.0);

  TH2F* xtprob_vs_treduchisqx = new TH2F("xtprob_vs_treduchisqx", "xtprob_vs_treduchisqx", 500, 0.0, 1.0, 500, 0.0, 30.0);
  TH2F* ytprob_vs_treduchisqy = new TH2F("ytprob_vs_treduchisqy", "ytprob_vs_treduchisqy", 500, 0.0, 1.0, 500, 0.0, 30.0);

  TH2F* xtprob_vs_ndf = new TH2F("xtprob_vs_ndf", "xtprob_vs_ndf", 500, 0.0, 1.0, 10, 0.5, 10.5);
  TH2F* ytprob_vs_ndf = new TH2F("ytprob_vs_ndf", "ytprob_vs_ndf", 500, 0.0, 1.0, 10, 0.5, 10.5);

  TH2F* tchisqx_vs_ndf = new TH2F("tchisqx_vs_ndf", "tchisqx_vs_ndf", 500, 0.0, 90.0, 10, 0.5, 10.5);
  TH2F* tchisqy_vs_ndf = new TH2F("tchisqy_vs_ndf", "tchisqy_vs_ndf", 500, 0.0, 90.0, 10, 0.5, 10.5);

  const int nprob=6;
  double probs[nprob+1]={3.5, 5.5, 7.5, 8.5, 9.5, 10.5, 12.5};
  //  double probs[nprob+1]={2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 12.5};

  TH1F* h_xnprob[nprob];
  TH1F* h_ynprob[nprob];
  TH1F* h_xtnprob[nprob];
  TH1F* h_ytnprob[nprob];

  for (int ij=0; ij<nprob; ij++) {
    sprintf(name, "h_xnprob_%i", ij);
    sprintf(title, "h_xnprob_%i_%i", int(probs[ij]+1), int(probs[ij+1]));
    h_xnprob[ij] = new TH1F(name, title, 120, 0.0, 1.0);

    sprintf(name, "h_ynprob_%i", ij);
    sprintf(title, "h_ynprob_%i_%i", int(probs[ij]+1), int(probs[ij+1]));
    h_ynprob[ij] = new TH1F(name, title, 120, 0.0, 1.0);

    sprintf(name, "h_xtnprob_%i", ij);
    sprintf(title, "h_xtnprob_%i_%i", int(probs[ij]+1), int(probs[ij+1]));
    h_xtnprob[ij] = new TH1F(name, title, 120, 0.0, 1.0);

    sprintf(name, "h_ytnprob_%i", ij);
    sprintf(title, "h_ytnprob_%i_%i", int(probs[ij]+1), int(probs[ij+1]));
    h_ytnprob[ij] = new TH1F(name, title, 120, 0.0, 1.0);
  }

  TH2F* correction_xtime[nlayer];
  TH2F* correction_ytime[nlayer];

  TH2F* fitted_rms_xtime[nlayer];
  TH2F* fitted_rms_ytime[nlayer];

  TH1F* muon_xmul_eq[nlayer][2*nmxiter][nlayer][4];//layer number, iiteration, ndf== , 3 conditoons ---- 0-No Condition, 1-Only Position, 2-Only Time, 3- Both Position and Time
  TH1F* muon_ymul_eq[nlayer][2*nmxiter][nlayer][4];

  TH1F* muon_xmul_gr[nlayer][2*nmxiter][nlayer][4];//layer number, iiteration, ndf== , 3 conditoons ---- 0-No Condition, 1-Only Position, 2-Only Time, 3- Both Position and Time
  TH1F* muon_ymul_gr[nlayer][2*nmxiter][nlayer][4];


  TH1F* time_xreso[nlayer][2*nmxiter];
  TH1F* time_yreso[nlayer][2*nmxiter];

  TH1F* time_xreso_chi2_ndf[nlayer][2*nmxiter][13];//0-chisq/ndf<infinity, 1-chisq/ndf<10, 2-chisq/ndf<5, 3-chisq/ndf<2, 4-chisq/ndf<1.5, 5-chisq/ndf<1
  TH1F* time_yreso_chi2_ndf[nlayer][2*nmxiter][13];//13------<infinity,20,10,9,8,7,6,5,4,3,2,1.5,1

  TH1F* timeshift_xreso[nlayer-1][nlayer];//layer0-layer1, layer1-layer2, layer2-layer3, layer3-layer4,--------xside
  TH1F* timeshift_yreso[nlayer-1][nlayer];//layer0-layer1, layer1-layer2, layer2-layer3, layer3-layer4,--------yside
  TH1F* timeshift_xyreso[nlayer-1][nlayer];//(layerx0+layery0)/2 - (layerx1+layery1)/2,   etc......................


  TH2F* time_xtdev_ytdev[nlayer][2*nmxiter];

  TH2F* time_xreso_vs_ndf[nlayer][2*nmxiter];
  TH2F* time_yreso_vs_ndf[nlayer][2*nmxiter];

  TH1F* time_xreso_nTermSplit[nlayer][nTermResSplit][2*nmxiter];
  TH1F* time_yreso_nTermSplit[nlayer][nTermResSplit][2*nmxiter];

  TH1F* time_xreso_split[nlayer][nsplit][2*nmxiter];
  TH1F* time_yreso_split[nlayer][nsplit][2*nmxiter];

  TH1F* tmptime_xreso[nlayer][2*nmxiter];
  TH1F* tmptime_yreso[nlayer][2*nmxiter];

  TH2F* time_xyreso[nlayer][2*nmxiter];
  bool istime_xyreso[nlayer][2*nmxiter];

  bool xposEffPass[nlayer][2*nmxiter];
  bool yposEffPass[nlayer][2*nmxiter];

  TH1F* time_xstrreso[nlayer][nstrip][2*nmxiter];
  TH1F* time_ystrreso[nlayer][nstrip][2*nmxiter];

  TH2F* time_mean_reso = new TH2F("time_mean_reso", "time_mean_reso", 2*nlayer, -0.5, 2*nlayer-0.5, 2*nmxiter, -0.5, 2*nmxiter-0.5);
  TH2F* time_rms_reso = new TH2F("time_rms_reso", "time_rms_reso", 2*nlayer, -0.5, 2*nlayer-0.5, 2*nmxiter, -0.5, 2*nmxiter-0.5);
  TH2F* time_corrms_reso = new TH2F("time_corrms_reso", "time_corrms_reso", 2*nlayer, -0.5, 2*nlayer-0.5, 2*nmxiter, -0.5, 2*nmxiter-0.5);
  TH2F* time_exterr_reso = new TH2F("time_exterr_reso", "time_exterr_reso", 2*nlayer, -0.5, 2*nlayer-0.5, 2*nmxiter, -0.5, 2*nmxiter-0.5);



  TH1F* dir_cxlay[nlayerit+1][nmxiter][nprob];
  TH1F* dir_cylay[nlayerit+1][nmxiter][nprob];

  TH1F* dir_cx[nlayerit+1][nmxiter];
  TH1F* dir_cy[nlayerit+1][nmxiter];

  TH1F* dir_cx2[nlayerit+1][nmxiter];
  TH1F* dir_cy2[nlayerit+1][nmxiter];

  TH1F* dir_c0x[nlayerit+1][nmxiter];
  TH1F* dir_c0y[nlayerit+1][nmxiter];

  TH1F* time_entry[nmxiter];
  TH1F* time_underflow[nmxiter];
  TH1F* time_overflow[nmxiter];

  TH1F* time_offset[nmxiter];
  TH1F* shift_time_mnft[nmxiter];
  TH1F* statmean_time[nmxiter];
  TH1F* statrms_time[nmxiter];
  TH1F* statskew_time[nmxiter];
  TH1F* statkurt_time[nmxiter];

  TH1F* rms_time[nmxiter];
  TH1F* rms_timeused[nmxiter];

  TH1F* time_offsetx[nmxiter];
  TH1F* shift_time_mnftx[nmxiter];
  TH1F* statmean_timex[nmxiter];
  TH1F* statrms_timex[nmxiter];
  TH1F* statskew_timex[nmxiter];
  TH1F* statkurt_timex[nmxiter];

  TH1F* rms_timex[nmxiter];
  TH1F* rms_timeusedx[nmxiter];

  TH1F* time_offsety[nmxiter];
  TH1F* shift_time_mnfty[nmxiter];
  TH1F* statmean_timey[nmxiter];
  TH1F* statrms_timey[nmxiter];
  TH1F* statskew_timey[nmxiter];
  TH1F* statkurt_timey[nmxiter];

  TH1F* rms_timey[nmxiter];
  TH1F* rms_timeusedy[nmxiter];

  TH1F* time_mulxreso[nlayer][2*nmxiter][nmxtimehit];
  TH1F* time_mulyreso[nlayer][2*nmxiter][nmxtimehit];
  const int nbandregion = 4;
  TH1F* time_mulxreso_band[nlayer][2*nmxiter][nmxtimehit][nbandregion];
  TH1F* time_mulyreso_band[nlayer][2*nmxiter][nmxtimehit][nbandregion];

  TH1F* timex_mulxyreso_band[nlayer][2*nmxiter][nstrip/2][nstrip/2][nbandregion];
  TH1F* timey_mulxyreso_band[nlayer][2*nmxiter][nstrip/2][nstrip/2][nbandregion];

  TH1F* time_mulxreso_split[nlayer][nsplit][2*nmxiter][nmxtimehit];
  TH1F* time_mulyreso_split[nlayer][nsplit][2*nmxiter][nmxtimehit];

  TH1F* xtime_exterr[nlayer][2*nmxiter];
  TH1F* ytime_exterr[nlayer][2*nmxiter];

  TH1F* xtime_exterr_chi2_ndf[nlayer][2*nmxiter][13];
  TH1F* ytime_exterr_chi2_ndf[nlayer][2*nmxiter][13];

  TH2F* xtime_exterr_vs_ndf[nlayer][2*nmxiter];
  TH2F* ytime_exterr_vs_ndf[nlayer][2*nmxiter];

  TH1F* xtime_exterr_nTermSplit[nlayer][nTermResSplit][2*nmxiter];
  TH1F* ytime_exterr_nTermSplit[nlayer][nTermResSplit][2*nmxiter];

  TH2F* rawhits_corr_xymul[nlayer];
  TH2F* rawhits_xlay_corr_mul[nlayer][nlayer];
  TH2F* rawhits_ylay_corr_mul[nlayer][nlayer];

  for (int ij=0; ij<nlayer; ij++) {
    sprintf(name, "rawhits_corr_xymul_l%i", ij);
    rawhits_corr_xymul[ij] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5,  nstrip+1, -0.5, nstrip+0.5);
    for (int jk=ij+1; jk<nlayer; jk++) {
      sprintf(name, "rawhits_xlay_corr_mul_l%i_l%i", ij, jk);
      rawhits_xlay_corr_mul[ij][jk] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5,  nstrip+1, -0.5, nstrip+0.5);

      sprintf(name, "rawhits_ylay_corr_mul_l%i_l%i", ij, jk);
      rawhits_ylay_corr_mul[ij][jk] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5,  nstrip+1, -0.5, nstrip+0.5);
    }
  }

  if (isTiming) {
    if (isalign>0) {
      for (int ij=0; ij<nlayer; ij++) {
	sprintf(title, "correction_xtime_l%i", ij);
	correction_xtime[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nmxiter, -0.5, nmxiter-0.5);

	sprintf(title, "correction_ytime_l%i", ij);
	correction_ytime[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nmxiter, -0.5, nmxiter-0.5);

	sprintf(title, "fitted_rms_xtime_l%i", ij);
	fitted_rms_xtime[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nmxiter, -0.5, nmxiter-0.5);

	sprintf(title, "fitted_rms_ytime_l%i", ij);
	fitted_rms_ytime[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nmxiter, -0.5, nmxiter-0.5);
      }
    }


    for (int ij=0; ij<nlayer-1; ij++) {
      double trange2 =  10.;// 35.0 : 35.0; // 25. : 20.; //7.5 : 10.;
      for(int jk=0; jk<nlayer; jk++) {
	sprintf(title, "timeshift_xreso_l%i_l%i", ij, jk);
	timeshift_xreso[ij][jk] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);
	sprintf(title, "timeshift_yreso_l%i_l%i", ij, jk);
	timeshift_yreso[ij][jk] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);
	sprintf(title, "timeshift_xyreso_l%i_l%i", ij, jk);
	timeshift_xyreso[ij][jk] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);
      }
    }

    for (int kl=0; kl<2*nmxiter; kl++ ){
      double trange = (kl<nmxiter) ? 10. : 25.;// 35.0 : 35.0; // 25. : 20.;//7.5 : 10.;
      double trange2 = (kl<nmxiter) ? 10. : 25.;// 35.0 : 35.0; // 25. : 20.; //7.5 : 10.;
      for (int ij=0; ij<nlayer; ij++) {

	for(int eq=0; eq<nlayer; eq++) {
	  for(int cond=0; cond<4; cond++) {
	    sprintf(title, "muon_xmul_l%i_i%i_eq%i_cond%i", ij, kl, eq, cond);
	    muon_xmul_eq[ij][kl][eq][cond] = new TH1F(title, title, 64, -0.5, 63.5); //150, -80.0, 70.0); //120, -trange2, trange2);
	    sprintf(title, "muon_ymul_l%i_i%i_eq%i_cond%i", ij, kl, eq, cond);
	    muon_ymul_eq[ij][kl][eq][cond] = new TH1F(title, title, 64, -0.5, 63.5); //150, -80.0, 70.0); //120, -trange2, trange2);

  	    sprintf(title, "muon_xmul_l%i_i%i_gr%i_cond%i", ij, kl, eq, cond);
	    muon_xmul_gr[ij][kl][eq][cond] = new TH1F(title, title, 64, -0.5, 63.5); //150, -80.0, 70.0); //120, -trange2, trange2);
	    sprintf(title, "muon_ymul_l%i_i%i_gr%i_cond%i", ij, kl, eq, cond);
	    muon_ymul_gr[ij][kl][eq][cond] = new TH1F(title, title, 64, -0.5, 63.5); //150, -80.0, 70.0); //120, -trange2, trange2);

	  }
	}

	sprintf(title, "time_xreso_l%i_i%i", ij, kl);
	time_xreso[ij][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);
	sprintf(title, "time_yreso_l%i_i%i", ij, kl);
	time_yreso[ij][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);

	for(int mn=0; mn<13; mn++) {
	  sprintf(title, "time_xreso_l%i_i%i_c%i", ij, kl, mn);
	  time_xreso_chi2_ndf[ij][kl][mn] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);
	  sprintf(title, "time_yreso_l%i_i%i_c%i", ij, kl, mn);
	  time_yreso_chi2_ndf[ij][kl][mn] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);
	}


	sprintf(title, "time_xtdev_ytdev_l%i_i%i", ij, kl);
	time_xtdev_ytdev[ij][kl] = new TH2F(title, title, 120, -trange2, trange2, 120, -trange2, trange);


	sprintf(title, "time_xreso_vs_ndf_l%i_i%i", ij, kl);
	time_xreso_vs_ndf[ij][kl] = new TH2F(title, title, 120, -trange2, trange2, 10, 0.5, 10.5);
	sprintf(title, "time_yreso_vs_ndf_l%i_i%i", ij, kl);
	time_yreso_vs_ndf[ij][kl] = new TH2F(title, title, 120, -trange2, trange2, 10, 0.5, 10.5);

	for(int um=0;um<nsplit;um++) {
	sprintf(title, "time_xreso_l%i_splt%i_i%i", ij,um, kl);
	time_xreso_split[ij][um][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);

	sprintf(title, "time_yreso_l%i_splt%i_i%i", ij, um,kl);
	time_yreso_split[ij][um][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);
	}

	for(int um=0;um<nTermResSplit;um++) {
	  sprintf(title, "time_xreso_l%i_str%i_to_%i_i%i", ij,um*20,(um+1)*20,kl);
	time_xreso_nTermSplit[ij][um][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);

	sprintf(title, "time_yreso_l%i_str%i_to_%i_i%i", ij, um*20,(um+1)*20,kl);
	time_yreso_nTermSplit[ij][um][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);
	}


	sprintf(title, "tmptime_xreso_l%i_i%i", ij, kl);
	tmptime_xreso[ij][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);

	sprintf(title, "tmptime_yreso_l%i_i%i", ij, kl);
	tmptime_yreso[ij][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);

	sprintf(title, "time_xyreso_l%i_i%i", ij, kl);
	time_xyreso[ij][kl] = new TH2F(title, title, 120, -trange2-1., trange2-1., 120, -trange2-1., trange2-1.);

	for (int lm=0; lm<nmxtimehit; lm++) {
	  sprintf(title, "time_mulxreso_l%i_i%i_m%i", ij, kl, lm+1);
	  time_mulxreso[ij][kl][lm] = new TH1F(title, title, 120, -trange2, trange2);

	  sprintf(title, "time_mulyreso_l%i_i%i_m%i", ij, kl, lm+1);
	  time_mulyreso[ij][kl][lm] = new TH1F(title, title, 120, -trange2, trange2);
	  for(int nbd =0;nbd<nbandregion;nbd++) {
	    sprintf(title, "time_mulxreso_l%i_i%i_m%i_nbd%i", ij, kl, lm+1,nbd);
	  time_mulxreso_band[ij][kl][lm][nbd] = new TH1F(title, title, 120, -trange2, trange2);

	  sprintf(title, "time_mulyreso_l%i_i%i_m%i_nbd%i", ij, kl, lm+1,nbd);
	  time_mulyreso_band[ij][kl][lm][nbd] = new TH1F(title, title, 120, -trange2, trange2);
	  }
	}

	  for(int ix=0;ix<nstrip/2;ix++) {
	    for(int iy=0;iy<nstrip/2;iy++) {
	      for(int nbd =0;nbd<nbandregion;nbd++) {
	      sprintf(title, "timex_mulxyreso_l%i_i%i_x%i_y%i_nbd%i", ij, kl, ix,iy,nbd);
	      timex_mulxyreso_band[ij][kl][ix][iy][nbd] = new TH1F(title, title, 120, -trange2, trange2);

	      sprintf(title, "timey_mulxyreso_l%i_i%i_x%i_y%i_nbd%i", ij, kl, ix,iy,nbd);
	      timey_mulxyreso_band[ij][kl][ix][iy][nbd] = new TH1F(title, title, 120, -trange2, trange2);
	    }
	  }
	}

	for (int lm=0; lm<nmxtimehit; lm++) {
	  for(int um=0;um<nsplit;um++) {
	    sprintf(title, "time_mulxreso_l%i_splt%i_i%i_m%i", ij,um, kl, lm+1);
	  time_mulxreso_split[ij][um][kl][lm] = new TH1F(title, title, 120, -trange2, trange2);

	  sprintf(title, "time_mulyreso_l%i_splt%i_i%i_m%i", ij,um, kl, lm+1);
	  time_mulyreso_split[ij][um][kl][lm] = new TH1F(title, title, 120, -trange2, trange2);
	  }
	  }


	sprintf(title, "xtime_exterr_l%i_i%i", ij, kl);
	xtime_exterr[ij][kl]=new TH1F(title, title, 120, 0.0, 2.1);
	sprintf(title, "ytime_exterr_l%i_i%i", ij, kl);
	ytime_exterr[ij][kl]=new TH1F(title, title, 120, 0.0, 2.1);

	for(int mn=0; mn<13; mn++) {
	  sprintf(title, "xtime_exterr_l%i_i%i_c%i", ij, kl, mn);
	  xtime_exterr_chi2_ndf[ij][kl][mn]=new TH1F(title, title, 120, 0.0, 2.1);
	  sprintf(title, "ytime_exterr_l%i_i%i_c%i", ij, kl, mn);
	  ytime_exterr_chi2_ndf[ij][kl][mn]=new TH1F(title, title, 120, 0.0, 2.1);
	}

	sprintf(title, "xtime_exterr_vs_ndf_l%i_i%i", ij, kl);
	xtime_exterr_vs_ndf[ij][kl]=new TH2F(title, title, 120, 0.0, 2.1, 10, 0.5, 10.5);
	sprintf(title, "ytime_exterr_vs_ndf_l%i_i%i", ij, kl);
	ytime_exterr_vs_ndf[ij][kl]=new TH2F(title, title, 120, 0.0, 2.1, 10, 0.5, 10.5);
	for(int um=0;um<nTermResSplit;um++) {
	  sprintf(title, "xtime_exterr_l%i_str%i_to_%i_i%i", ij,um*20,(um+1)*20, kl);
	  xtime_exterr_nTermSplit[ij][um][kl]=new TH1F(title, title, 120, 0.0, 2.1);
	  sprintf(title, "ytime_exterr_l%i_str%i_to_%i_i%i", ij,um*20,(um+1)*20, kl);
	  ytime_exterr_nTermSplit[ij][um][kl]=new TH1F(title, title, 120, 0.0, 2.1);
	}
	for (int jk=0; jk<nstrip; jk++) {
	  sprintf(title, "time_xstrreso_l%i_s%i_i%i", ij, jk, kl);
	  time_xstrreso[ij][jk][kl] = new TH1F(title, title, 120, -trange, trange);

	  sprintf(title, "time_ystrreso_l%i_s%i_i%i", ij, jk, kl);
	  time_ystrreso[ij][jk][kl] = new TH1F(title, title, 120, -trange, trange);

	  //	  sprintf(title, "time_xstshift_l%i_s%i_i%i", ij, jk, kl);
	  //	  time_xstshift[ij][jk][kl] = new TH1F(title, title, 180, -trange, trange);

	  //	  sprintf(title, "time_ystshift_l%i_s%i_i%i", ij, jk, kl);
	  //	  time_ystshift[ij][jk][kl] = new TH1F(title, title, 180, -trange, trange);
	}
      }
    }
    for (int kl=0; kl<nmxiter; kl++) {
      for (int ij=0; ij<=nlayerit; ij++) {
	for (int lm=0; lm<nprob; lm++) {
	  sprintf(name, "dir_cxlay_l%i_i%i_n%i", ij, kl, lm);
	  sprintf(title, "dir_cxlay_l%i_i%i_n[%i-%i]", ij, kl, int(probs[lm]+1), int(probs[lm+1]));
	  dir_cxlay[ij][kl][lm] = new TH1F(name, title, 180, -3.5, 5.5);
	  sprintf(name, "dir_cylay_l%i_i%i_n%i", ij, kl, lm);
	  sprintf(title, "dir_cylay_l%i_i%i_n[%i-%i]", ij, kl, int(probs[lm]+1), int(probs[lm+1]));

	  dir_cylay[ij][kl][lm] = new TH1F(name, title, 180, -3.5, 5.5);
	}

	sprintf(title, "dir_cx_l%i_i%i", ij, kl);
	dir_cx[ij][kl] = new TH1F(title, title, 180, -3.5, 5.5);

	sprintf(title, "dir_cy_l%i_i%i", ij, kl);
	dir_cy[ij][kl] = new TH1F(title, title, 180, -3.5, 5.5);

	sprintf(title, "dir_cx2_l%i_i%i", ij, kl);
	dir_cx2[ij][kl] = new TH1F(title, title, 180, -3.5, 5.5);

	sprintf(title, "dir_cy2_l%i_i%i", ij, kl);
	dir_cy2[ij][kl] = new TH1F(title, title, 180, -3.5, 5.5);

	sprintf(title, "dir_c0x_l%i_i%i", ij, kl);
	dir_c0x[ij][kl] = new TH1F(title, title, 120, 70., 130.);

	sprintf(title, "dir_c0y_l%i_i%i", ij, kl);
	dir_c0y[ij][kl] = new TH1F(title, title, 120, 70., 130.);
      }
    }

    if (isalign>0) {
      for (int kl=0; kl<nmxiter; kl++) {
	sprintf(title, "time_entry_i%i", kl);
	time_entry[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "time_underflow_i%i", kl);
	time_underflow[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "time_overflow_i%i", kl);
	time_overflow[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "time_offset_i%i", kl);
	time_offset[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "shift_time_mnft_i%i", kl);
	shift_time_mnft[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "statmean_time_i%i", kl);
	statmean_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "statrms_time_i%i", kl);
	statrms_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "statskew_time_i%i", kl);
	statskew_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "statkurt_time_i%i", kl);
	statkurt_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "rms_time_i%i", kl);
	rms_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	sprintf(title, "rms_timeused_i%i", kl);
	rms_timeused[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);

	/////////////////////////////////////////////
	sprintf(title, "time_offsetx_i%i", kl);
	time_offsetx[kl] = new TH1F(title, title, 60, -7.0, 8.0);

	sprintf(title, "shift_time_mnftx_i%i", kl);
	shift_time_mnftx[kl] = new TH1F(title, title, 60, -0.3, 0.3);

	sprintf(title, "statmean_timex_i%i", kl);
	statmean_timex[kl] = new TH1F(title, title, 60, -1.0, 1.0);

	sprintf(title, "statrms_timex_i%i", kl);
	statrms_timex[kl] = new TH1F(title, title, 60, 0.6, 3.0);

	sprintf(title, "statskew_timex_i%i", kl);
	statskew_timex[kl] = new TH1F(title, title, 60, -1.0, 3.0);

	sprintf(title, "statkurt_timex_i%i", kl);
	statkurt_timex[kl] = new TH1F(title, title, 60, -2.0, 6.0);

	sprintf(title, "rms_timex_i%i", kl);
	rms_timex[kl] = new TH1F(title, title, 60, 0.6, 2.7);

	sprintf(title, "rms_timeusedx_i%i", kl);
	rms_timeusedx[kl] = new TH1F(title, title, 60, 0.6, 2.4);

	/////////////////////////////////////////////
	sprintf(title, "time_offsety_i%i", kl);
	time_offsety[kl] = new TH1F(title, title, 60, -7.0, 8.0);

	sprintf(title, "shift_time_mnfty_i%i", kl);
	shift_time_mnfty[kl] = new TH1F(title, title, 60, -0.3, 0.3);

	sprintf(title, "statmean_timey_i%i", kl);
	statmean_timey[kl] = new TH1F(title, title, 60, -1.0, 1.0);

	sprintf(title, "statrms_timey_i%i", kl);
	statrms_timey[kl] = new TH1F(title, title, 60, 0.6, 3.0);

	sprintf(title, "statskew_timey_i%i", kl);
	statskew_timey[kl] = new TH1F(title, title, 60, -1.0, 3.0);

	sprintf(title, "statkurt_timey_i%i", kl);
	statkurt_timey[kl] = new TH1F(title, title, 60, -2.0, 6.0);

	sprintf(title, "rms_timey_i%i", kl);
	rms_timey[kl] = new TH1F(title, title, 60, 0.6, 2.7);

	sprintf(title, "rms_timeusedy_i%i", kl);
	rms_timeusedy[kl] = new TH1F(title, title, 60, 0.5, 2.4);
      }
    }
  } //  if (isTiming)
  TH2F* h_xcorhits = new TH2F("xcorhits","xcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_ycorhits = new TH2F("ycorhits","ycorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_xycorhits = new TH2F("xycorhits","xycorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  TH2F* h_xtcorhits = new TH2F("xtcorhits","xtcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_ytcorhits = new TH2F("ytcorhits","ytcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_xytcorhits = new TH2F("xytcorhits","xytcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  TH2F* h_xrawcorhits = new TH2F("xrawcorhits","xrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_yrawcorhits = new TH2F("yrawcorhits","yrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_xyrawcorhits = new TH2F("xyrawcorhits","xyrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  TH2F* h_xtrawcorhits = new TH2F("xtrawcorhits","xtrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_ytrawcorhits = new TH2F("ytrawcorhits","ytrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_xytrawcorhits = new TH2F("xytrawcorhits","xytrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  TH2F* h_xcorstrips[nlayer];
  TH2F* h_ycorstrips[nlayer];
  TH2F* h_xycorstrips[nlayer];
  TH2F* h_yxcorstrips[nlayer];

  TH2F* h_corr_layer_mult = new TH2F("corr_layer_mult", "corr_layer_mult", 2*nlayer, -0.5, 2*nlayer-0.5, 2*nlayer, -0.5, 2*nlayer-0.5);
  TH2F* h_raw_xcorstrips[nlayer];
  TH2F* h_raw_ycorstrips[nlayer];
  TH2F* h_raw_xystrpnhits[nlayer];
  TH2F* h_raw_yxstrpnhits[nlayer];
  TH2F* h_raw_xstrpnhits[nlayer];
  TH2F* h_raw_ystrpnhits[nlayer];

  TH2F* h_xtcorstrips[nlayer];
  TH2F* h_ytcorstrips[nlayer];
  TH2F* h_xytcorstrips[nlayer];
  TH2F* h_yxtcorstrips[nlayer];
  TH2F* h_xmucorstrips[nlayer];
  TH2F* h_ymucorstrips[nlayer];
  TH2F* h_xymucorstrips[nlayer];
  TH2F* h_yxmucorstrips[nlayer];
  TH2F* h_xmucornhits[nlayer];
  TH2F* h_ymucornhits[nlayer];
  TH2F* h_xymucornhits[nlayer];
  TH2F* h_yxmucornhits[nlayer];

  TH2F* h_xtstrpmult_xtdchitmult[nlayer];
  TH2F* h_ytstrpmult_ytdchitmult[nlayer];
  TH2F* h_xtstrpmult_ytdchitmult[nlayer];
  TH2F* h_ytstrpmult_xtdchitmult[nlayer];
  TH2F* h_corr_tdc_ref = new TH2F("corr_tdc_ref","corr_tdc_ref",nlayer,-0.5,nlayer-0.5,nlayer,-0.5,nlayer-0.5);
  TH1F* h_tdc_ref = new TH1F("tdc_ref","tdc_ref",nlayer+1,-0.5,nlayer+0.5);



  for (int ij=0; ij<nlayer; ij++) {
    sprintf(name, "xcorstrips_l%i", ij);
    h_xcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "ycorstrips_l%i", ij);
    h_ycorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "xycorstrips_l%i", ij);
    h_xycorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "yxcorstrips_l%i", ij);
    h_yxcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);


    sprintf(name, "xmucorstrips_l%i", ij);
    h_xmucorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "ymucorstrips_l%i", ij);
    h_ymucorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "xymucorstrips_l%i", ij);
    h_xymucorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "yxmucorstrips_l%i", ij);
    h_yxmucorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name,"xmucornhits_l%i",ij);
    h_xmucornhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name,"ymucornhits_l%i",ij);
    h_ymucornhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name,"xymucornhits_l%i",ij);
    h_xymucornhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name,"yxmucornhits_l%i",ij);
    h_yxmucornhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "raw_xcorstrips_l%i", ij);
    h_raw_xcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "raw_ycorstrips_l%i", ij);
    h_raw_ycorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "xstrpnhits_l%i", ij);
    h_raw_xstrpnhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "ystrpnhits_l%i", ij);
    h_raw_ystrpnhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "xystrpnhits_l%i", ij);
    h_raw_xystrpnhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "yxstrpnhits_l%i", ij);
    h_raw_yxstrpnhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "xtstrpmult_xtdchitmult_l%i",ij);
    h_xtstrpmult_xtdchitmult[ij] =  new TH2F(name,name,nstrip,-0.5,nstrip-0.5, nTDCpLayer*nhitspchannel,-0.5,nTDCpLayer*nhitspchannel-0.5);
    sprintf(name, "ytstrpmult_ytdchitmult_l%i",ij);
    h_ytstrpmult_ytdchitmult[ij] =  new TH2F(name,name,nstrip,-0.5,nstrip-0.5, nTDCpLayer*nhitspchannel,-0.5,nTDCpLayer*nhitspchannel-0.5);
    sprintf(name, "xtstrpmult_ytdchitmult_l%i",ij);
    h_xtstrpmult_ytdchitmult[ij] =  new TH2F(name,name,nstrip,-0.5,nstrip-0.5, nTDCpLayer*nhitspchannel,-0.5,nTDCpLayer*nhitspchannel-0.5);
    sprintf(name, "ytstrpmult_xtdchitmult_l%i",ij);
    h_ytstrpmult_xtdchitmult[ij] =  new TH2F(name,name,nstrip,-0.5,nstrip-0.5, nTDCpLayer*nhitspchannel,-0.5,nTDCpLayer*nhitspchannel-0.5);


    sprintf(name, "xtcorstrips_l%i", ij);
    h_xtcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "ytcorstrips_l%i", ij);
    h_ytcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "xytcorstrips_l%i", ij);
    h_xytcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "yxtcorstrips_l%i", ij);
    h_yxtcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

  }

  TH1F* time_xraw[nlayer][nstrip];
  TH1F* time_yraw[nlayer][nstrip];
  TH2F* time_xraw2d = new TH2F("time_xraw2d", "time_xraw2d", nlayer, -0.5, nlayer-0.5, nstrip, -0.5, nstrip-0.5);
  TH2F* time_yraw2d = new TH2F("time_yraw2d", "time_yraw2d", nlayer, -0.5, nlayer-0.5, nstrip, -0.5, nstrip-0.5);

  //  double xt_slope_cor[nlayer][nstrip][nstrip]={0};
  //  double yt_slope_cor[nlayer][nstrip][nstrip]={0};

#ifdef TIMESLOPE
  TH2F* time_xslope[nlayer][nTDCpLayer];
  TH2F* time_yslope[nlayer][nTDCpLayer];

  TH2F* time_xslope_pr[nlayer][nstrip][nmxiter];
  TH2F* time_yslope_pr[nlayer][nstrip][nmxiter];

  TH2F* time_X_L7_corr[nlayer];
  TH2F* time_X_corr[nlayer][nTDCpLayer];
   TH2F* time_Y_corr[nlayer][nTDCpLayer];
  for(int ij=0;ij<nlayer;ij++) {
      sprintf(name,"time_corr_L7_L%i",ij);
      time_X_L7_corr[ij] = new TH2F(name,name,nstrip,-0.5,nstrip-0.5,150,250,400);
      for(int jk=0;jk<nTDCpLayer;jk++) {
	sprintf(name,"time_corr_X_L%i_tdc%i",ij,jk);
	time_X_corr[ij][jk] = new TH2F(name,name,nstrip,-0.5,nstrip-0.5,350,250,600);
	sprintf(name,"time_corr_Y_L%i_tdc%i",ij,jk);
	time_Y_corr[ij][jk] = new TH2F(name,name,nstrip,-0.5,nstrip-0.5,350,250,600);
      }
  }

#endif
  for (int ij=0; ij <nlayer; ij++) {
    for (int jk=0; jk<nstrip; jk++) {
      sprintf(title, "time_xraw_l%i_s%i", ij, jk);
      time_xraw[ij][jk] = new TH1F(title, title, 120, 750., 1350.);
      sprintf(title, "time_yraw_l%i_s%i", ij, jk);
      time_yraw[ij][jk] = new TH1F(title, title, 120, 750., 1350.);
#ifdef TIMESLOPE
      for (int kl=0; kl<nmxiter; kl++) {

	sprintf(title, "time_xslope_pr_l%i_s%i_i%i", ij, jk, kl);
	time_xslope_pr[ij][jk][kl] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 90, -6.0, 6.0);
	sprintf(title, "time_yslope_pr_l%i_s%i_i%i", ij, jk, kl);
	time_yslope_pr[ij][jk][kl] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 90, -6.0, 6.0);
      }
#endif

    }
#ifdef TIMESLOPE
    for (int jk=0; jk<nTDCpLayer; jk++) {
      sprintf(title, "time_xslope_l%i_tdc%i", ij, jk);
      time_xslope[ij][jk] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 200, 250., 450.);
      sprintf(title, "time_yslope_l%i_tdc%i", ij, jk);
      time_yslope[ij][jk] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 200, 250., 450.);
    }
#endif
  }
  TH1F* timex_correl[nlayer][nlayer][nmxiter+1];
  TH1F* timey_correl[nlayer][nlayer][nmxiter+1];

  TH2F* time_both_cormean[nmxiter+1];
  TH2F* time_both_corrms[nmxiter+1];

  const int npixel=16;
  TH1F* timexy_correl[nlayer][npixel+1][nmxiter+1];

  TH2F* timexy_cormean[npixel+1];
  TH2F* timexy_corrms[npixel+1];

  const int ntimecor=2; //how many correlation do we want
  TProfile* indtimexy_prof[nlayer][ntimecor];
  TH1F* indtimexy_correl[nlayer][nstrip][nstrip][ntimecor];
  TH2F* indtimexy_cormean[nlayer][ntimecor];
  TH2F* indtimexy_corrms[nlayer][ntimecor];
  TH2F* indtimexy_fitmean[nlayer][ntimecor];
  TH2F* indtimexy_fitrms[nlayer][ntimecor];

  if (isTiming) {
    for (int ij=0; ij<nlayer-1; ij++) {
      for (int jk=ij+1; jk<nlayer; jk++) {
	for (int kl=0; kl<nmxiter+1; kl++) {
	  sprintf(title, "timex_correl_l%i_l%i_i%i", ij, jk, kl);
	  timex_correl[ij][jk][kl] = new TH1F(title, title, 120, -8., 8.);

	  sprintf(title, "timey_correl_l%i_l%i_i%i", ij, jk, kl);
	  timey_correl[ij][jk][kl] = new TH1F(title, title, 120, -8., 8.);
	}
      }
    }

    for (int kl=0; kl<nmxiter+1; kl++) {
      sprintf(title, "time_both_cormean_i%i", kl);
      time_both_cormean[kl] = new TH2F(title, title, nlayer+1, -0.5, nlayer+0.5, nlayer, -0.5, nlayer-0.5);

      sprintf(title, "time_both_corrms_i%i", kl);
      time_both_corrms[kl] = new TH2F(title, title, nlayer+1, -0.5, nlayer+0.5, nlayer, -0.5, nlayer-0.5);
    }

    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<npixel+1; jk++) {
	for (int kl=0; kl<nmxiter+1; kl++) {
	  sprintf(title, "timexy_correl_l%i_pixel%i_i%i", ij, jk, kl);
	  timexy_correl[ij][jk][kl] = new TH1F(title, title, 120, -8., 8.);
	}
      }
    }

    for (int kl=0; kl<npixel+1; kl++) {
      sprintf(title, "timexy_cormean_pixel%i", kl);
      timexy_cormean[kl] = new TH2F(title, title, nlayer, -0.5, nlayer-0.5, nmxiter+1, -0.5, nmxiter+0.5);
      sprintf(title, "timexy_corrms_pixel%i", kl);
      timexy_corrms[kl] = new TH2F(title, title, nlayer, -0.5, nlayer-0.5, nmxiter+1, -0.5, nmxiter+0.5);
    }


    for (int ij=0; ij<nlayer; ij++) {
      for (int lm=0; lm<ntimecor; lm++) {
	sprintf(title, "indtimexy_prof_l%i_%i", ij, lm);
	indtimexy_prof[ij][lm] = new TProfile(title, title, 2*nstrip, -0.5, 2*nstrip-0.5, -5.0, 15.0);

	sprintf(title, "indtimexy_cormean_l%i_%i", ij, lm);
	indtimexy_cormean[ij][lm] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
	sprintf(title, "indtimexy_corrms_l%i_%i", ij, lm);
	indtimexy_corrms[ij][lm] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

	sprintf(title, "indtimexy_fitmean_l%i_%i", ij, lm);
	indtimexy_fitmean[ij][lm] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
	sprintf(title, "indtimexy_fitrms_l%i_%i", ij, lm);
	indtimexy_fitrms[ij][lm] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

	for (int jk=0; jk<nstrip; jk++) {
	  for (int kl=0; kl<nstrip; kl++) {
	    sprintf(name, "indtimexy_correl_l%i_x%i_y%i_%i", ij, jk, kl, lm);
	    sprintf(title, "indtimexy_correl_l%i_x%i_y%i_%i", ij, jk+8, kl+8, lm);
	    indtimexy_correl[ij][jk][kl][lm] = new TH1F(name, title, 120, -5.0, 15.0);
	  }
	}
      }
    }
  } // if (isTiming)

#ifdef OLDFORMAT
  DaqEvent *event = new DaqEvent();
#endif
#ifdef BARC_EVTBLR
  EvtData *event = new EvtData();

#endif

  int xhits[nlayer],yhits[nlayer];   //number of hits after noise rejection for position fit
  int xallhits[nlayer][nTDCpLayer],yallhits[nlayer][nTDCpLayer]; // raw : for one hit, return strip number otherwise -10-multiplicity

  bool passxtime[nlayer], passytime[nlayer]; // passed through timing criteria in all hitted strips in x/ystr_x/ytdev
  bool passxtmy[nlayer], passytmx[nlayer]; //for X/Y strips, use boundary of Y/X extrapolation
  int istrxtime[nlayer], istrytime[nlayer]; //strip with early timing

  float errxco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};
  float erryco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};

  double xxerr[nlayer], yyerr[nlayer];
  for ( int ix=0; ix<nlayer; ix++) {xxerr[ix] = yyerr[ix] = errxco[ix]*errxco[ix];}

  double xrms[nlayer]={0};
  double yrms[nlayer]={0};

  double deltaposxcov[nlayer][nlayer];
  int deltaposxCount[nlayer][nlayer];

  double deltaposycov[nlayer][nlayer];
  int deltaposyCount[nlayer][nlayer];

  double deltatcov[nlayer][nlayer];
  int deltatCount[nlayer][nlayer];

  double deltatcov2[nlayer][nlayer];
  int deltatCount2[nlayer][nlayer];

  double deltatcovy[nlayer][nlayer];
  int deltatCounty[nlayer][nlayer];

  TH2F* h_deltaposxcov = new TH2F("deltaposxcov", "deltaposxcov", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_deltaposycov = new TH2F("deltapUosycov", "deltaposycov", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_deltatcov = new TH2F("deltatcov", "deltatcov", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_deltatcov2 = new TH2F("deltatcov2", "deltatcov2", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_deltatcovy = new TH2F("deltatcovy", "deltatcovy", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  //Pethu layerwise corr effic

  double posxcoreff[nlayer][nlayer];
  double posxcoreffCount;

  double posycoreff[nlayer][nlayer];
  double posycoreffCount;

  TH2F* h_posxcoreff = new TH2F("posxcoreff", "posxcoreff", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posycoreff = new TH2F("posycoreff", "posycoreff", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posxcoreff_def = new TH2F("posxcoreff_def", "posxcoreff_def", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posycoreff_def = new TH2F("posycoreff_def", "posycoreff_def", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posxcoreff_dif = new TH2F("posxcoreff_dif", "posxcoreff_dif", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posycoreff_dif = new TH2F("posycoreff_dif", "posycoreff_dif", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  for(int ij=0;ij<nlayer;ij++){
    for(int jk=0;jk<nlayer;jk++){
      deltatcov[ij][jk]= deltatcov2[ij][jk]=deltatcovy[ij][jk]=0.0;
      deltatCount[ij][jk]=deltatCount2[ij][jk]=deltatCounty[ij][jk]=0;
      deltaposxcov[ij][jk]= deltaposycov[ij][jk]=0.0;
      deltaposxCount[ij][jk]=deltaposyCount[ij][jk]=0;
      posxcoreff[ij][jk]= posycoreff[ij][jk]=0.0;


    }
  }
posxcoreffCount=0;posycoreffCount=0;
  //correlation of extrapolated position and time delay

  TH2F* xpos_xdev[nlayer][nstr_posmx];
  TH2F* xpos_ydev[nlayer][nstr_posmx];
  TH2F* ypos_xdev[nlayer][nstr_posmx];
  TH2F* ypos_ydev[nlayer][nstr_posmx];
  TH2F* xpos_xtdev[nlayer][nmxtimehit];
  TH2F* xpos_ytdev[nlayer][nmxtimehit];
  TH2F* ypos_xtdev[nlayer][nmxtimehit];
  TH2F* ypos_ytdev[nlayer][nmxtimehit];

  TH2F* xpos_xtdev_str[nlayer][nmxtimehit]; //Without strip position correction
  TH2F* ypos_ytdev_str[nlayer][nmxtimehit];

  TH2F* xpos_xtdev_glb[nlayer][nmxtimehit]; //Without strip position correction
  TH2F* ypos_ytdev_glb[nlayer][nmxtimehit];

  TH2F* xpos_ytdev_str[nlayer][nmxtimehit]; //Without strip position correction
  TH2F* ypos_xtdev_str[nlayer][nmxtimehit];

  TH2F* xpos_ytdev_glb[nlayer][nmxtimehit]; //Without strip position correction
  TH2F* ypos_xtdev_glb[nlayer][nmxtimehit];

  TH2F* xpos_xtdev_band_str[nlayer][nband][nmxtimehit]; //Without strip position correction
  TH2F* ypos_ytdev_band_str[nlayer][nband][nmxtimehit];

  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<nstr_posmx; jk++) {
      sprintf(name, "xpos_xdev_l%i_mul%i", ij, jk+1);
      xpos_xdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -2.0, 2.0);

      sprintf(name, "xpos_ydev_l%i_mul%i", ij, jk+1);
      xpos_ydev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -2.0, 2.0);

      sprintf(name, "ypos_xdev_l%i_mul%i", ij, jk+1);
      ypos_xdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -2.0, 2.0);

      sprintf(name, "ypos_ydev_l%i_mul%i", ij, jk+1);
      ypos_ydev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -2.0, 2.0);
    }
    for(int kl=0;kl<nband;kl++) {
      for(int jk=0;jk<nmxtimehit;jk++) {
	sprintf(name, "xpos_xtdev_str_l%i_nbd%i_mul%i", ij, kl,jk+1);
	xpos_xtdev_band_str[ij][kl][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

	sprintf(name, "ypos_ytdev_str_l%i_nbd%i_mul%i", ij,kl,jk+1);
	ypos_ytdev_band_str[ij][kl][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);
      }
    }
    for (int jk=0; jk<nmxtimehit; jk++) {
      sprintf(name, "xpos_xtdev_l%i_mul%i", ij, jk+1);
      xpos_xtdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "xpos_ytdev_l%i_mul%i", ij, jk+1);
      xpos_ytdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "ypos_xtdev_l%i_mul%i", ij, jk+1);
      ypos_xtdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "ypos_ytdev_l%i_mul%i", ij, jk+1);
      ypos_ytdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "xpos_xtdev_str_l%i_mul%i", ij, jk+1);
      xpos_xtdev_str[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "ypos_ytdev_str_l%i_mul%i", ij, jk+1);
      ypos_ytdev_str[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "xpos_xtdev_glb_l%i_mul%i", ij, jk+1);
      xpos_xtdev_glb[ij][jk] = new TH2F(name, name, nstrip+2, -1.5, nstrip+0.5, 60, -15.0, 15.0);

      sprintf(name, "ypos_ytdev_glb_l%i_mul%i", ij, jk+1);
      ypos_ytdev_glb[ij][jk] = new TH2F(name, name, nstrip+2, -1.5, nstrip+0.5, 60, -15.0, 15.0);

      sprintf(name, "xpos_ytdev_str_l%i_mul%i", ij, jk+1);
      xpos_ytdev_str[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "ypos_xtdev_str_l%i_mul%i", ij, jk+1);
      ypos_xtdev_str[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "xpos_ytdev_glb_l%i_mul%i", ij, jk+1);
      xpos_ytdev_glb[ij][jk] = new TH2F(name, name, nstrip+2, -1.5, nstrip+0.5, 60, -15.0, 15.0);

      sprintf(name, "ypos_xtdev_glb_l%i_mul%i", ij, jk+1);
      ypos_xtdev_glb[ij][jk] = new TH2F(name, name, nstrip+2, -1.5, nstrip+0.5, 60, -15.0, 15.0);
    }
  }

  TH1F* pos_xslope[nlayerit+1][nmxiter];
  TH1F* pos_yslope[nlayerit+1][nmxiter];

  TH1F* pos_theta[nlayerit+1][nmxiter];
  TH1F* pos_phi[nlayerit+1][nmxiter];

  for (int ij=0; ij<nlayerit+1; ij++) {
    for (int jk=0; jk<nmxiter; jk++) {
      sprintf(name, "pos_xslope_l%i_i%i", ij, jk);
      pos_xslope[ij][jk] = new TH1F(name, name, 120, -3.6, 3.6);

      sprintf(name, "pos_yslope_l%i_i%i", ij, jk);
      pos_yslope[ij][jk] = new TH1F(name, name, 120, -3.6, 3.6);

      sprintf(name, "pos_theta_l%i_i%i", ij, jk);
      pos_theta[ij][jk] = new TH1F(name, name, 180, 0.0, 90.0);

      sprintf(name, "pos_phi_l%i_i%i", ij, jk);
      pos_phi[ij][jk] = new TH1F(name, name, 120, -180.0, 180.0);

    }
  }

  TH2F* widthx_timex[nlayer][nmxiter]; //Not yet implemented in data. Width of pulse hegight vs time delay
  TH2F* widthy_timey[nlayer][nmxiter];
  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<nmxiter; jk++) {
      sprintf(name, "widthx_timex_l%i_i%i", ij, jk);
      widthx_timex[ij][jk] = new TH2F(name, name, 60., -20., 10., 90, 100., 1000.);
      sprintf(name, "widthy_timey_l%i_i%i", ij, jk);
      widthy_timey[ij][jk] = new TH2F(name, name, 60., -20., 10., 90, 100., 1000.);
    }
  }

  TH2F* sel_theta_phi_g4 = new TH2F("sel_theta_phi_g4", "sel_theta_phi_g4", 180,-180.,180.0, 90, 0.0, 90.0);
  TH2F* sel_theta_phi_g4_trig = new TH2F("sel_theta_phi_g4_trig", "sel_theta_phi_g4_trig", 180,-180.,180.0, 90, 0.0, 90.0);
  TH2F* sel_theta_mom = new TH2F("sel_theta_mom","sel_theta_mom",240,0.,60.,1000,0.150,3.);
  TH2F* sel_phi_mom = new TH2F("sel_phi_mom","sel_phi_mom",180,-180.,180.0,1000,0.150,3.);
  TH2F* sel_theta_rec_gen = new TH2F("sel_theta_rec_gen","sel_theta_rec_gen",90,0.,90.0,90,0.,90.0);

 const int nCount=250000; //0; //Number of events //GMA
 const int ntCount = 300;//1800;//3600; //Time is second
  const int nsetmx=2000; //Maximum number of set
   const int ntsetmx=1500;
  int iset=0;
  int isetold=-1;
  int itset = 0;
  int itsetold = -1;

  TH2F* strp_count_set[nsetmx]={0};
  TH1F* strp_xmult_set[nlayer][nsetmx]={0};
  TH1F* strp_ymult_set[nlayer][nsetmx]={0};
  TH2F* raw_occu_set[nlayer][nsetmx]={0};

  TH1F* xlayer_reso_set[nlayer][nsetmx]={0};
  TH1F* ylayer_reso_set[nlayer][nsetmx]={0};
  TH1F* time_xreso_set[nlayer][nsetmx]={0};
  TH1F* time_yreso_set[nlayer][nsetmx]={0};
  TH2F* inefficiency_corx_set[nlayer][nsetmx]={0};
  TH2F* triggereffi_x_set[nlayer][nsetmx]={0};
  TH2F* triggereffi_y_set[nlayer][nsetmx]={0};
  TH2F* totalentry_set[nlayer][nsetmx]={0};
#ifdef MONTECARLO
  //All these will go to the above four include files
  //#include "effic_pos_time_madurai64_aug_24_25_26_eff_change.txt"
 //   #include "effic_pos_time_madurai64_aug_24_25_26.txt" //"effi_postime_150705.txt"
  // #include "effic_pos_time_madurai64_aug_19.txt"
   // #include "effic_pos_time_madurai64_combined_eff_change.txt"


  //   #include "effic_pos_time_madurai64_combined.txt"
    //  #include "effic_trig_madurai64_aug_24_25_26.txt" //"effic_trig_150705.txt"
    //  #include "time_corr_madurai64_aug_24_25_26.txt" //time_corr_150705.txt" //"time_corr_2479_150228.txt"
  //  // #include "pos_time_inuse_madurai64_july_14_15.txt" //pos_time_inuse_150705.txt" //"pos_time_inuse_150228.txt"

  //  #include "pos_time_inuse_byhand_110117.txt"
   //#include "pos_time_inuse_madurai64_aug_24_25_26.txt" //pos_time_inuse_150705.txt" //"pos_time_inuse_150228.txt"

  /*
#include "effic_pos_time_madurai64_aug252627.txt"
#include "effic_trig_madurai64_aug252627.txt"
#include "time_corr_madurai64_aug252627.txt"
#include "pos_time_inuse_madurai64_aug252627.txt"
*/
#include "pos_time_inuse_madurai64_trg4567_aug252627_jj46.txt"
#include "time_corr_madurai64_trg4567_aug252627_jj46.txt"
#include "effic_trig_madurai64_trg4567_aug252627_jj46.txt"
#include "effic_pos_time_madurai64_trg4567_aug252627_jj46.txt"


/*
  //   const double slope_path=0.140;
  double ytimeshift=0.0;
double xoff[nlayer] = {-0.067085, -0.134262, -0.0177441, 0.0899896, 0.117874, 0.0472966, 0.169018, 0.242022, 0.239679, 0.139407, 0.296545, 0.428804};
double yoff[nlayer] = {0.447182, -0.00889266, 0.14596, 0.0519742, 0.0189781, 0.211516, -0.000755971, 0.110225, 0.112833, 0.0794779, 0.171389, 0.029938};

   double timeoffsetx[nlayer]={0}; //0, -6.302, -4.438, -4.826, -4.765, -3.144, -7.204, -5.564, -4.616, 0.864, -0.241, 1.585};
   double timeoffsety[nlayer]={0}; //0, -3.465, -3.741, -1.967, -2.528, -2.700, -3.613, -0.950, -2.360, 1.519,  2.767, 16.15};
   double xtoffset[nlayer][nstrip]={{0}};
   double ytoffset[nlayer][nstrip]={{0}};
   double xtoffystr[nlayer][nstrip]={{0}};
   double ytoffxstr[nlayer][nstrip]={{0}};
   double xt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};
   double yt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};
   double ineffi_table_uncX[nlayer][nstrip][nstrip]={{{0}}};
   double ineffi_table_uncY[nlayer][nstrip][nstrip]={{{0}}};
   double effi_trig_X[nlayer][nstrip][nstrip]={{{0}}};
   double effi_trig_Y[nlayer][nstrip][nstrip]={{{0}}};

   //   double timeserrx2[nlayer] = {0.787813, 0.728151, 0.594104, 0.590319, 0.899884, 1.05821, 0.586255, 0.630633, 1.20044, 0.843687, 0.623494, 0.707891};

   //   double timeserry2[nlayer] = {0.937061, 0.899824, 0.799044, 0.797195, 1.05882, 1.13075, 1.71724, 1.46413, 1.63097, 1.16403, 0.930633, 0.947063};

   //   for (int ij=0; ij<nlayer; ij++) {
   //     timeserrx2[ij] = timeserry2[ij] = 1.0*1.0;  //GMA
   //   }
//from test_2479_jj_3.root, 18th July 2016 truerm2
   double timeserrx2[nlayer] = {1.03837,0.910886,0.753216,0.677123,0.926671,1.0769,0.595245,0.660482,1.27677,1.01035,0.804336,0.991293};
   //   double timeserrx2[nlayer] = {0.951517,0.903104,0.821113,0.814621,0.97833,  1.04022,0.798109,0.83788, 1.12377,0.965306,0.845179,0.900186};
   //3               true[0] = {0.955993,0.907748,0.829315,0.819892,0.979419, 1.041,  0.799486,0.840678,1.12759,0.97259, 0.851285,0.902587,};
   //4               true[0] = {1.02796, 0.973469,0.957819,0.909363,1.01335,  1.07755,0.890681,0.966507,1.16705,1.08586, 1.01474, 1.08993,};
   //5               true[0] = {1.02841, 0.97228, 0.953233,0.902027,1.01046,  1.07046,0.884924,0.962161,1.16122,1.08813, 1.00987, 1.0872,}

   double timeserry2[nlayer] = {1.36734,1.1458,1.00839,0.886361,1.08254,1.16014,1.72595,1.52067,1.77919,1.47349,1.28456,1.29094};
   //   double timeserry2[nlayer] = {1.02475, 0.990309,0.93482,  0.918863,1.04715,1.08065,1.30887, 1.2327, 1.31418, 1.13814, 1.00606,1.04457};
   //3               true[1] = {1.02925, 0.99668, 0.946999, 0.926776,1.04953,1.08246,1.30951, 1.23694,1.32058, 1.1524 , 1.01621,1.046,};
   //4               true[1] = {1.10615, 1.07182, 1.03763,  1.03538, 1.10543,1.14718,1.41974, 1.34301,1.36697, 1.23174, 1.20712,1.24222,};
   //5               true[1] = {1.10502, 1.07266, 1.03751,  1.03516, 1.10506,1.1478, 1.41965, 1.34465,1.36691, 1.23091, 1.21031,1.24375,}

   for (int ij=0; ij<nlayer; ij++) {
     timeserrx2[ij] = timeserrx2[ij]* timeserrx2[ij];
     timeserry2[ij] = timeserry2[ij]* timeserry2[ij];
   }*/

#elif defined(OLDFORMAT) || defined(ONEMBYONEM)
  // #include "effic_trig_150807.txt"
  //#include "effi_postime_150807.txt"
   //#include "effic_trig_160629.txt" //Using trigger on 0-1-3-4 & 7-8-10-11
   //#include "effi_postime_160629.txt" //Using trigger on 0-1-3-4 & 7-8-10-11
#include "effic_trig_160707.txt" //Using trigger on 0-1-3-4 & 7-8-10-11
#include "effi_postime_160707.txt" //Using trigger on 0-1-3-4 & 7-8-10-11

#include "pos_time_inuse_160711.txt" //using only 2479 trigger data (otherwise wider time
                                     //resolution paticularly in few strip of Layer-1 & 2 in X-side)
#include "time_corr_160711.txt"

#elif defined(NEWRPCDAQ) && defined(CAUDATA)

#include "effic_pos_time_madurai_rpcdaq1.txt"
#include "effic_trig_madurai_rpcdaq1.txt"
#include "time_corr_madurai_rpcdaq2.txt"
#include "pos_time_unused_madurai_rpcdaq2.txt"


  /*  const double slope_path=0.140;
     double ytimeshift=0.0;


     double xoff[nlayer] = {0.151723, 0.24867, 0.186928, -0.0389714, 0.238133, 0.0667214, 0.110758, 0.0844209, -0.0294068, 0.177622, -0.00548693, -0.221559};
     double yoff[nlayer] = {-0.189908, -0.00563298, -0.024175, -0.113663, 0.0830386, 0.0607469, -0.168968, 0.509867, -0.199473, -0.12742, 0.0698661, -0.0766159};



     double xposerrsq[nmxhits][nlayer] = {{0.0909937, 0.0781988, 0.0642407, 0.0705577, 127.58, 0.118255, 0.0720771, 0.0757826, 0.110619, 0.0994707, 0.0856205, 0.117994},
					  {0.0669101, 0.0391113, 0.0434733, 0.0660287, 1040.53, 0.138427, 0.0394018, 0.0210766, 0.044745, 0.0349779, 0.0326737, 0.0374608},
					  {0.109605, 0.179713, 0.109713, 0.121659, 0.02, 0.10815, 0.111044, 0.180156, 0.0925861, 0.389084, 0.322293, 0.202311},
					  {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};
     double yposerrsq[nmxhits][nlayer] = {{0.08127, 0.0826758, 0.0808017, 0.14994, 0.02, 0.2128, 0.102187, 0.0880633, 0.120059, 0.064384, 0.0872874, 0.0935829},
					  {0.0694922, 0.0532706, 0.0562377, 0.0826594, 0.02, 0.0960678, 0.0594412, 0.0457572, 0.0535982, 0.0573025, 0.0524933, 0.054638},
					  {0.129846, 0.120791, 0.184287, 0.0910893, 0.02, 0.124909, 0.155805, 0.299488, 0.160473, 0.228101, 0.268159, 0.138742},
					  {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};


     double timeoffsetx[nlayer]={0}; //0, -6.302, -4.438, -4.826, -4.765, -3.144, -7.204, -5.564, -4.616, 0.864, -0.241, 1.585};
     double timeoffsety[nlayer]={0}; //0, -3.465, -3.741, -1.967, -2.528, -2.700, -3.613, -0.950, -2.360, 1.519,  2.767, 16.15};
     double xtoffset[nlayer][nstrip]={{0}};
     double ytoffset[nlayer][nstrip]={{0}};
     double xtoffystr[nlayer][nstrip]={{0}};
     double ytoffxstr[nlayer][nstrip]={{0}};
     double xt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};
     double yt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};
     double ineffi_table_uncX[nlayer][nstrip][nstrip]={{{0}}};
     double ineffi_table_uncY[nlayer][nstrip][nstrip]={{{0}}};
     double effi_trig_X[nlayer][nstrip][nstrip]={{{0}}};
     double effi_trig_Y[nlayer][nstrip][nstrip]={{{0}}};

     double timeserrx2[nlayer] = {0.578446, 0.771683, 0.79788, 0.920692, 0.25, 1.3069, 0.94565, 1.12215, 0.852368, 1.08287, 0.665644, 0.78667};

     double timeserry2[nlayer] = {0.744046, 0.808488, 1.40582, 1.25652, 0.25, 1.62414, 1.53565, 1.45712, 1.32676, 1.04824, 1.36222, 0.740011};

     double Xlay_Corr_eff[nlayer][nlayer] = {{0.9}};
     double Ylay_Corr_eff[nlayer][nlayer] = {{0.9}};
  */


#elif defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
/* && defined(CAUDATA1)*/

  /*#include "effic_pos_time_madurai_rpcdaq1.txt"
#include "effic_trig_madurai_rpcdaq1.txt"
#include "time_corr_madurai_rpcdaq1.txt"
#include "pos_time_unused_madurai_rpcdaq1.txt"
  */
  //  #include "effic_pos_time_rpcdaq.txt"
  //  #include "effic_trig_rpcdaq.txt"
    // #include "time_corr_rpcdaq.txt"
    //   #include "pos_time_inuse_rpcdaq.txt"
  //#include "effic_pos_time_madurai_rpcdaq3.txt"
  //#include "effic_trig_madurai_rpcdaq3.txt"
  //#include "time_corr_madurai_rpcdaq3.txt"
  //#include "pos_time_inuse_madurai_rpcdaq3.txt"
  // #include "effic_pos_time_madurai_rpcdaq4.txt"
  //#include "effic_trig_madurai_rpcdaq4.txt"
  //#include "time_corr_madurai_rpcdaq4.txt"
  //#include "pos_time_inuse_madurai_rpcdaq4.txt"

  //#include "effic_pos_time_madurai_rpcdaq5.txt"
  //#include "effic_trig_madurai_rpcdaq5.txt"
  //#include "time_corr_madurai_rpcdaq5.txt"
  //#include "pos_time_inuse_madurai_rpcdaq5.txt"

  //#include "effic_pos_time_madurai_rpcdaq6.txt"
  //#include "effic_trig_madurai_rpcdaq6.txt"
  //#include "time_corr_madurai_rpcdaq6.txt"
  //#include "pos_time_inuse_madurai_rpcdaq6.txt"

  //#include "effic_pos_time_madurai64_aug252627.txt"
  //#include "effic_trig_madurai64_aug252627.txt"
  //#include "time_corr_madurai64_aug252627.txt"
  // #include "pos_time_inuse_madurai64_aug252627.txt"

    // #include "effic_pos_time_madurai64_trg4567_sept27.txt"
    //#include "pos_time_inuse_madurai64_trg4567_sept27.txt"
  //#include "effic_trig_madurai64_trg4567_sept27.txt "
  //#include "time_corr_madurai64_trg4567_sept27.txt"

  //#include "effic_pos_time_madurai64_trg4567_oct16.txt"
  //  #include "effic_trig_madurai64_trg4567_oct16.txt"
  //    //#include "time_corr_madurai64_trg4567_oct16.txt"
  //    #include "pos_time_inuse_madurai64_trg4567_oct16.txt"

  //#include "effic_pos_time_madurai64_trg4567_oct16_jj4.txt"
  //#include "effic_trig_madurai64_trg4567_oct16_jj4.txt"
  //#include "pos_time_inuse_madurai64_trg4567_oct16_jj4.txt"
  //#include "time_corr_madurai64_trg4567_oct16_jj4.txt"

// #include "effic_pos_time_madurai64_trg4567_oct16_jj8.txt"
// #include "effic_trig_madurai64_trg4567_oct16_jj8.txt"
// #include "pos_time_inuse_madurai64_trg4567_oct16_jj8.txt"
// #include "time_corr_madurai64_trg4567_oct16_jj8.txt"




//#include "effic_pos_time_madurai64_trg4567_aug25_jj41.txt"
//#include "pos_time_inuse_madurai64_trg4567_aug252627_jj42.txt"
//#include "time_corr_madurai64_trg4567_aug252627_jj42.txt"
//#include "effic_trig_madurai64_trg4567_aug252627_jj42.txt"
//#include "effic_pos_time_madurai64_trg4567_aug252627_jj42.txt"

#ifdef SHIFT0
#include "pos_time_inuse_madurai64_trg4567_aug252627_jj46.txt"
#include "time_corr_madurai64_trg4567_aug252627_jj46.txt"

#include "effic_trig_madurai64_trg4567_aug252627_jj46.txt"
#include "effic_pos_time_madurai64_trg4567_aug252627_jj46.txt"

#else
    // #include "/home/inodaq/Anal_madurai_64_09-11-17/effic_pos_time_madurai_miical_26052018.txt"
    // #include "/home/inodaq/Anal_madurai_64_09-11-17/effic_trig_madurai_miical_26052018.txt"
    // #include "/home/inodaq/Anal_madurai_64_09-11-17/time_corr_madurai_miical_26052017.txt"
 // #include "effic_trig_madurai64_trg4567_oct16_jj36.txt"
 // #include "pos_time_inuse_madurai64_trg4567_oct16_jj36.txt"
 // #include "effic_pos_time_madurai64_trg4567_oct16_jj36.txt"
 // #include "time_corr_madurai64_trg4567_oct16_jj36.txt"

#include "effic_pos_time_madurai64_miical_20180602.txt"
#include "effic_trig_madurai64_miical_20180602.txt"
  //#include "time_corr_madurai64_miical_xxx_20180602.txt"
  //#include "time_cor_madurai64_miical_june142018_trg6789.txt"
  //#include "time_corr_madurai64_miical_20180731103.txt"
  //#include "time_cor_madurai64_miical_20180916_2038_trg1357.txt"
  // #include "time_corr_madurai64_miical_20180928_trg4567.txt"
//#include "pos_time_inuse_madurai64_miical_20180916_2038_trg1357.txt"
// #include "pos_time_inuse_madurai64_miical_xxx_20180602.txt"
//#include "time_corr_madurai64_20181005_trg6789.txt"
// #include "time_cor_madurai64_miical_20181129_153631_trg6789.txt"
//  #include "time_corr_miical_20190103_092942.txt"
// #include "time_corr_miical_20181226_192924.txt"
//#include "time_corr_madurai64_miical_mamta.txt"
//#include "time_corr_miical_26dec2018_noparaboliccorr.txt"
//#include "time_corr_miical_26dec2018_withparabolcorr.txt" //used where all the button are there nmnhits>=4
 #include "time_corr_miniical_20181226_jj160.txt"    //jim jim 2018 data after removing button and with all layer tot applied
  //#include "time_corr_miniical_20181226_jj113.txt" //jim jim 2018 data removing all button but not all layer tot applied****not*****
//#include "time_corr_madurai_miniical_20181128_021221.txt"
//#include "20180603165542.txt"
#include "pos_time_inuse_madurai64_26dec2018.txt"   // jim jim 2018 data
  //#include "time_corr_miical64_10layer_phase2_20210616.txt"

  //    #include "time_corr_madurai64_10layer_20210616_jj169.txt"    //jim jim data with the changed resistor rpc on top but all the resistors are 15,18,22 ohm
  //  #include "time_corr_madurai64_10layer_20210616_jj156.txt"    //jim jim data with the changed resistor rpc on top but all the resistors are 15,18,22 ohm  layer 5 not used in fit
  //#include "time_corr_miniical_20210923_174123_jj75.txt" //jim jim disable this and include new time corr
  //   #include "time_corr_madurai64_10layer_20210922_jj109.txt"  //jim jim  used finally with button
  //#include "time_corr_madurai64_10layer_20210922_jj117.txt"  //jim jim  used finally removing button
  //         #include "time_corr_madurai64_10layer_20210922_jj130.txt"  //jim jim  used finally removing button and using tot in all layers
  //#include "time_corr_madurai64_10layer_20210528_jj166.txt"    //jim jim data with 5th layer on top and all the resistors are 39 ohm
//#include "time_corr_madurai64_rpc196_9p9kv_without_paracorr_jj63.txt"
  //#include "time_corr_madurai64_rpc196_9p9_without_paracorr_jj64.txt"
//#include "TOT_Corr_Values.txt"
  //#include "TOT_Corr_Values_mul.txt"
//#include "TOT_Corr_Values_fitfunc.txt"

/*
    double xposerrsq[nmxhits][nlayer] = {{0.113279, 0.0818007, 0.0744173, 0.0789816, 0.0751351, 0.0735476, 0.0615256, 0.063695, 0.072483, 0.117045, 0.02, 0.02},
					 {0.0667296, 0.0326339, 0.0348177, 0.0350458, 0.0298755, 0.0262117, 0.0685308, 0.02, 0.0906959, 0.0785113, 0.02, 0.02},
					 {0.206786, 0.186776, 0.250294, 0.202714, 0.168759, 0.216815, 0.181322, 0.192127, 0.232744, 0.218847, 0.02, 0.02},
					 {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};
    double yposerrsq[nmxhits][nlayer] = {{0.114141, 0.0736622, 0.0746076, 0.084452, 0.0749809, 0.0669596, 0.0565289, 0.0576602, 0.0580676, 0.0898826, 0.02, 0.02},
					 {0.0699685, 0.0354721, 0.029979, 0.0452407, 0.0414917, 0.0417897, 0.0400507, 0.0360985, 0.0305168, 0.0617391, 0.02, 0.02},
					 {0.186147, 0.204842, 0.211932, 0.236854, 0.189943, 0.192353, 0.123058, 0.192505, 0.216022, 0.217481, 0.02, 0.02},
					 {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};
    */
    /*
    double xposerrsq[nmxhits][nlayer] = {{0.106357, 0.0745162, 0.0665068, 0.0689515, 0.072208, 0.0696568, 0.0615053, 0.063785, 0.0772467, 0.117856, 0.02, 0.02},
					 {0.0868731, 0.0357328, 0.041457, 0.0479257, 0.0415874, 0.0429817, 0.0353604, 0.0255034, 0.0760438, 0.0952902, 0.02, 0.02},
					 {0.161657, 0.127963, 0.17764, 0.167898, 0.114854, 0.15253, 0.102063, 0.109315, 0.142508, 0.170694, 0.02, 0.02},
					 {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};
    double yposerrsq[nmxhits][nlayer] = {{0.0989658, 0.0653573, 0.0673575, 0.0696284, 0.0652101, 0.0664139, 0.0516923, 0.0525604, 0.0477285, 0.0804387, 0.02, 0.02},
					 {0.080124, 0.0411679, 0.0385191, 0.0444887, 0.0516316, 0.0537554, 0.0517678, 0.0475527, 0.0394105, 0.0689746, 0.02, 0.02},
					 {0.147055, 0.128808, 0.139476, 0.153457, 0.127888, 0.129357, 0.101453, 0.108971, 0.0919941, 0.153493, 0.02, 0.02},
					 {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};


    //from 20180916-2038
    double xoff[nlayer] = { 0.};
    double yoff[nlayer] = { 0.};
    double xposerrsq[nmxhits][nlayer] = {{0.108909, 0.0768391, 0.0675894, 0.0639632, 0.0612128, 0.0632426, 0.0567263, 0.0582503, 0.0681779, 0.109901, 0.02, 0.02},
    {0.077319, 0.0360595, 0.0243164, 0.029018, 0.029708, 0.025834, 0.0266059, 0.02, 0.0470693, 0.0747822, 0.02, 0.02},
    {0.160539, 0.135314, 0.180577, 0.14747, 0.123004, 0.154799, 0.107872, 0.115905, 0.13441, 0.185262, 0.02, 0.02},
    {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};
    double yposerrsq[nmxhits][nlayer] = {{0.105558, 0.0889334, 0.066422, 0.0678411, 0.0623633, 0.0596872, 0.0528411, 0.0515247, 0.0496765, 0.080447, 0.02, 0.02},
    {0.077234, 0.0416191, 0.031046, 0.0377991, 0.0439917, 0.0453914, 0.0423967, 0.038345, 0.0372904, 0.0598037, 0.02, 0.02},
    {0.157018, 0.157472, 0.159579, 0.175056, 0.14038, 0.140228, 0.103738, 0.124047, 0.122436, 0.180575, 0.02, 0.02},
    {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};


    */
     double Xlay_Corr_eff[nlayer][nlayer] = {{0.9}};
     double Ylay_Corr_eff[nlayer][nlayer] = {{0.9}};
     //const double slope_path=0.140;
     // double ytimeshift=0.0;
 // double timeoffsetx[nlayer] = {0, -3.18426, -2.84865, -3.0988, -2.1637, -0.769563, -1.30187, -0.250336, 0.383466, 0.0424222, 0.133142, 1.89628};
 //   double timeoffsety[nlayer] = {0, -3.37519, -2.8192, -3.02597, -1.58974, -0.640005, -1.01107, -0.159466, 0.451872, 0.151471, 0.00930853, 0.882566};
 //   double xtoffset[nlayer][nstrip]={{0}};
 //   double ytoffset[nlayer][nstrip]={{0}};
 //   double xtoffystr[nlayer][nstrip]={{0}};
 //   double ytoffxstr[nlayer][nstrip]={{0}};
 //   double xt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};
 //   double yt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};

      double xtw_slope_cor[nlayer][nstrip][nwidth]={{{0}}};
      double ytw_slope_cor[nlayer][nstrip][nwidth]={{{0}}};

 //    double ineffi_table_uncX[nlayer][nstrip][nstrip]={{{0}}};
 //    double ineffi_table_uncY[nlayer][nstrip][nstrip]={{{0}}};
 //    double effi_trig_X[nlayer][nstrip][nstrip]={{{0}}};
 //    double effi_trig_Y[nlayer][nstrip][nstrip]={{{0}}};

 //   double timeserrx2[nlayer] = {0.846541, 0.70648, 0.880431, 0.837493, 1.03849, 1.13667, 1.09692, 1.03072, 0.958941, 0.832274, 0.772451, 0.851193};

 //     double timeserry2[nlayer] = {0.868608, 0.953603, 0.993743, 1.08104, 1.15238, 1.2563, 1.21322, 1.1519, 1.15236, 0.959422, 0.842085, 0.961558};
#endif
  /*
 const double slope_path=0.140;
     double ytimeshift=0.0;


     double xoff[nlayer] = {0.151723, 0.24867, 0.186928, -0.0389714, 0.238133, 0.0667214, 0.110758, 0.0844209, -0.0294068, 0.177622, -0.00548693, -0.221559};
     double yoff[nlayer] = {-0.189908, -0.00563298, -0.024175, -0.113663, 0.0830386, 0.0607469, -0.168968, 0.509867, -0.199473, -0.12742, 0.0698661, -0.0766159};



     double xposerrsq[nmxhits][nlayer] = {{0.0909937, 0.0781988, 0.0642407, 0.0705577, 127.58, 0.118255, 0.0720771, 0.0757826, 0.110619, 0.0994707, 0.0856205, 0.117994},
					  {0.0669101, 0.0391113, 0.0434733, 0.0660287, 1040.53, 0.138427, 0.0394018, 0.0210766, 0.044745, 0.0349779, 0.0326737, 0.0374608},
					  {0.109605, 0.179713, 0.109713, 0.121659, 0.02, 0.10815, 0.111044, 0.180156, 0.0925861, 0.389084, 0.322293, 0.202311},
					  {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};
     double yposerrsq[nmxhits][nlayer] = {{0.08127, 0.0826758, 0.0808017, 0.14994, 0.02, 0.2128, 0.102187, 0.0880633, 0.120059, 0.064384, 0.0872874, 0.0935829},
					  {0.0694922, 0.0532706, 0.0562377, 0.0826594, 0.02, 0.0960678, 0.0594412, 0.0457572, 0.0535982, 0.0573025, 0.0524933, 0.054638},
					  {0.129846, 0.120791, 0.184287, 0.0910893, 0.02, 0.124909, 0.155805, 0.299488, 0.160473, 0.228101, 0.268159, 0.138742},
					  {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};

  */
  //  double ytimeshift = 0.0104378;
  //  const double slope_path = 0.140;
  //  double timeoffsetx[nlayer] = {0, -3.18426, -2.84865, -3.0988, -2.1637, -0.769563, -1.30187, -0.250336, 0.383466, 0.0424222, 0.133142, 1.89628};
  //  double timeoffsety[nlayer] = {0, -3.37519, -2.8192, -3.02597, -1.58974, -0.640005, -1.01107, -0.159466, 0.451872, 0.151471, 0.00930853, 0.882566};
  //  double xtoffset[nlayer][nstrip]={{0}};
  //  double ytoffset[nlayer][nstrip]={{0}};
  //  double xtoffystr[nlayer][nstrip]={{0}};
  //  double ytoffxstr[nlayer][nstrip]={{0}};
  //  double xt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};
  //  double yt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};

    // double ineffi_table_uncX[nlayer][nstrip][nstrip]={{{0}}};
    // double ineffi_table_uncY[nlayer][nstrip][nstrip]={{{0}}};
    // double effi_trig_X[nlayer][nstrip][nstrip]={{{0}}};
    // double effi_trig_Y[nlayer][nstrip][nstrip]={{{0}}};

  //  double timeserrx2[nlayer] = {0.846541, 0.70648, 0.880431, 0.837493, 1.03849, 1.13667, 1.09692, 1.03072, 0.958941, 0.832274, 0.772451, 0.851193};

  //  double timeserry2[nlayer] = {0.868608, 0.953603, 0.993743, 1.08104, 1.15238, 1.2563, 1.21322, 1.1519, 1.15236, 0.959422, 0.842085, 0.961558};


    // double Xlay_Corr_eff[nlayer][nlayer] = {{0.9}};
     // double Ylay_Corr_eff[nlayer][nlayer] = {{0.9}};



#else

   //double xtoffset[nlayer][nstrip]={{0}};
   //double ytoffset[nlayer][nstrip]={{0}};
   //double xtoffystr[nlayer][nstrip]={{0}};
   //double ytoffxstr[nlayer][nstrip]={{0}};
   //double xt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};
   //double yt_slope_cor[nlayer][nstrip][nstrip]={{{0}}};

   //double timeoffsetx[nlayer] = {0};
   //double timeoffsety[nlayer] = {0};

   //double timeserrx2[nlayer] = {0.60027, 0.566616, 0.647098, 0.658842, 0.25, 0.965348, 0.701985, 0.76533, 0.617521, 0.630508, 0.566835, 0.668583};

   //double timeserry2[nlayer] = {0.507792, 0.868581, 0.878023, 0.903707, 0.25, 1.18415, 1.12585, 1.20028, 0.826304, 0.872657, 1.09714, 0.37428};

   //for (int ij=0; ij<nlayer; ij++) {
     //timeserrx2[ij] = timeserry2[ij] = 2.0;
   //}

   //#include "effic_pos_time_madurai64_aug_05_64strip.txt" //"effi_postime_150705.txt"
   //#include "effic_trig_madurai64_aug_05_64strip.txt" //"effic_trig_150705.txt"
   //#include "time_corr_madurai64_aug_05_64strip.txt" //time_corr_150705.txt" //"time_corr_2479_150228.txt"
  //// #include "pos_time_inuse_madurai64_july_14_15.txt" //pos_time_inuse_150705.txt" //"pos_time_inuse_150228.txt"
   //#include "pos_time_inuse_madurai64_aug_05_64strip.txt" //pos_time_inuse_150705.txt" //"pos_time_inuse_150228.txt"
	//#include "effic_pos_time_madurai64_aug_24_25_26_eff_change.txt"
    //// #include "effic_pos_time_madurai64_combined_eff_change.txt"
     //#include "effic_trig_madurai64_combined.txt"
 #include "effic_pos_time_madurai64_combined.txt"
  // #include "effic_pos_time_madurai64_aug_24_25_26.txt" //"effi_postime_150705.txt"
  #include "effic_trig_madurai64_aug_24_25_26.txt" //"effic_trig_150705.txt"
   #include "time_corr_madurai64_aug_24_25_26.txt" //time_corr_150705.txt" //"time_corr_2479_150228.txt"
  // #include "pos_time_inuse_madurai64_july_14_15.txt" //pos_time_inuse_150705.txt" //"pos_time_inuse_150228.txt"
   //#include "pos_time_inuse_madurai64_aug_24_25_26.txt" //pos_t
  #include "pos_time_inuse_byhand_110117.txt"

//#include "effic_pos_time_madurai64_july_14_15.txt"
//#include "effic_trig_madurai64_july_14_15.txt"
//#include "pos_time_inuse_madurai64_july_14_15.txt"
//#include "time_corr_madurai64_july_14_15.txt"

   //   //#include "effic_trig_madurai64_160316.txt" //"effi_postime_150705.txt"
   //   //#include "effi_postime_madurai64_160316.txt" //"effic_trig_150705.txt"
   //   //#include "time_corr_madurai64_160316.txt" //time_corr_150705.txt" //"time_corr_2479_150228.txt"
   //   //#include "pos_time_inuse_madurai64_160316.txt" //pos_time_inuse_150705.txt" //"pos_time_inuse_150228.txt"
#endif

// 051116 changed ---only in MC
#ifdef MONTECARLO
#ifdef FASTSIM
for(int ij =0;ij<nmxhits;ij++) {
	for(int jk=0;jk<nlayer;jk++) {
		xposerrsq[ij][jk] *= 1.1015;//1.0358; 1.0815;
		yposerrsq[ij][jk] *= 1.1207;//1.0490; 1.105;
	}
}
#endif
#ifdef FULLG4SIM
for(int ij =0;ij<nmxhits;ij++) {
	for(int jk=0;jk<nlayer;jk++) {
	  xposerrsq[ij][jk] *= 1.;//1.647/1.556;//(1.696/1.556);//1.500;//00.5000;//1.200;//0.900;//1.1015;//1.0358; 1.0815;
	  yposerrsq[ij][jk] *= 1.;//1.64/1.472;//;//(1.705/1.472);//1.500;//0.5000;//1.200;//0.900;//1.1207;//1.0490; 1.105;
	}
}
#endif
#endif

#ifdef ISEFFICIENCY
  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<nstrip; jk++) {
      for (int kl=0; kl<nstrip; kl++) {
	defefficiency_uncx[ij]->Fill(jk, kl, ineffi_table_uncX[ij][jk][kl]);
	defefficiency_uncy[ij]->Fill(jk, kl, ineffi_table_uncY[ij][jk][kl]);
	deftriggereffi_x[ij]->Fill(jk, kl, effi_trig_X[ij][jk][kl]);
	deftriggereffi_y[ij]->Fill(jk, kl, effi_trig_Y[ij][jk][kl]);
      }
    }
  }

for(int ij=0;ij<nlayer;ij++) {
	for(int jk=0;jk<nlayer;jk++) {
	 	h_posxcoreff_def->Fill(ij,jk,Xlay_Corr_eff[ij][jk]);
	 	h_posycoreff_def->Fill(ij,jk,Ylay_Corr_eff[ij][jk]);
     }
    }
#endif
//-------------------------------------------------------------------------

  if (isalign >0) {
    for (int ij=0; ij<nlayer; ij++) {
      //      xoff[ij] =  yoff[ij] = 0;
      // timeoffsetx[ij] =  timeoffsety[ij] =0;

      for (int jk=0; jk<nstrip; jk++) {
	//	xtoffset[ij][jk] =  ytoffset[ij][jk] = 0;
	//	xtoffystr[ij][jk] = ytoffxstr[ij][jk]= 0;
	for (int kl=0; kl<nstrip; kl++) {
	  //	  xt_slope_cor[ij][jk][kl] = yt_slope_cor[ij][jk][kl] = 0.0;
	}
      }
    }
  }

  //nonlinear shift in x/y timing with the number in in y/x strip
  //  double xtoffystr[nlayer][nstrip]={0};
  //  double ytoffxstr[nlayer][nstrip]={0};


  double xtrms[nlayer][nstrip]={0};
  double ytrms[nlayer][nstrip]={0};

  double timesx[nlayer]={0};
  double timesy[nlayer]={0};

  int cau_calib1[nlayer] = {0};
  int cau_calib2[nlayer] = {0};


  double widthx[nlayer]={0};
  double widthy[nlayer]={0};

  //Position resolution in layers. 150219
  double pos_xrms[nlayer][nmxhits];
  double pos_yrms[nlayer][nmxhits];

  //Timing resolution in layers. 150112
  double time_xrms[nlayer];
  double time_yrms[nlayer];

  //  double time_tmpxrms[nlayer];
  //  double time_tmpyrms[nlayer];


  for (int ij=0; ij<nlayer; ij++) {
    time_xrms[ij] = time_yrms[ij] = 1.1; //1.1//ns
    for (int jk=0; jk<nmxhits; jk++) {
      if (jk<=1) {
	pos_xrms[ij][jk] = pos_yrms[ij][jk] =0.25;
      } else {
	pos_xrms[ij][jk] = pos_yrms[ij][jk] =0.30;
      }
    }
  }

  bool filloccu = true; // For first entries fill occuplancy plot and remove noisy channels

  double errcst, errcov, errlin;
  double dist[nlayer], xtime[nlayer],ytime[nlayer],  xtdev[nlayer],ytdev[nlayer];
  double initxtime[nlayer][nTDCpLayer], rawxtime[nlayer], rawxtime0[nlayer], rawytime0[nlayer], rawxtime1[nlayer], rawytime1[nlayer], rawxtime2[nlayer], rawytime2[nlayer];
  double initytime[nlayer][nTDCpLayer], rawytime[nlayer], rawxtime3[nlayer], rawytime3[nlayer]; //150628 : Time without path lenght correction, but with all other correction

  bool xusedtime[nlayer], yusedtime[nlayer];

  double xval, yval;
  int nxtfail, nytfail, nentry, nTotalp, nTotalt,nTotallx,nTotally;
  double xc0inters,yc0inters,DDx,errcst_tx,errcov_tx,errlin_tx,DDy,errcst_ty,errcov_ty,errlin_ty;
  int firstiter,ntcormx, iiterrs, lstr, lend,occulyr ,isfill,ntrigX,ntrigY;

  double xtext[nlayer]; //, xtexter[nlayer];
  double ytext[nlayer]; //, ytexter[nlayer];


  if(nmxiter==1) {
    sprintf(infile,"%s%i_align.root", outfilx, isalign);
    TFile *fileIn_align = new TFile(infile, "read");
    TTree *in_align= (TTree*)fileIn_align->Get("Alignment");
    in_align->SetBranchAddress("xoff",xoff);
    in_align->SetBranchAddress("yoff",yoff);
    in_align->SetBranchAddress("align_xstr_ydev",align_xstr_ydev);
    in_align->SetBranchAddress("align_ystr_xdev",align_ystr_xdev);
    in_align->SetBranchAddress("align_xstr_xdev",align_xstr_xdev);
    in_align->SetBranchAddress("align_ystr_ydev",align_ystr_ydev);
    in_align->SetBranchAddress("xposerrsq",xposerrsq);
    in_align->SetBranchAddress("yposerrsq",yposerrsq);
    in_align->SetBranchAddress("timeoffsetx",timeoffsetx);
    in_align->SetBranchAddress("timeoffsety",timeoffsety);
    in_align->SetBranchAddress("timeserrx2",timeserrx2);
    in_align->SetBranchAddress("timeserry2",timeserry2);
    in_align->SetBranchAddress("xtoffset",xtoffset);
    in_align->SetBranchAddress("ytoffset",ytoffset);
    in_align->SetBranchAddress("xtoffystr",xtoffystr);
    in_align->SetBranchAddress("ytoffxstr",ytoffxstr);
    in_align->SetBranchAddress("xt_slope_cor",xt_slope_cor);
    in_align->SetBranchAddress("yt_slope_cor",yt_slope_cor);
    fileIn_align->cd();
    in_align->GetEntry(0);    // while running for entire event file put "ij<numentries".
  }


  sprintf(outfile, "%s%i_align.root", outfilx, isalign);
  TFile* fileOut_align = new TFile(outfile, "recreate");

  TTree *align = new TTree("Alignment","Alignment");
  align->Branch("xoff", xoff, "xoff[12]/D");
  align->Branch("yoff", yoff, "yoff[12]/D");
  align->Branch("align_xstr_ydev", align_xstr_ydev, "align_xstr_ydev[12][3]/D");
  align->Branch("align_ystr_xdev", align_ystr_xdev, "align_ystr_xdev[12][3]/D");
  align->Branch("align_xstr_xdev", align_xstr_xdev, "align_xstr_xdev[12][3]/D");
  align->Branch("align_ystr_ydev", align_ystr_ydev, "align_ystr_ydev[12][3]/D");
  align->Branch("xposerrsq", xposerrsq, "xposerrsq[4][12]/D");
  align->Branch("yposerrsq", yposerrsq, "yposerrsq[4][12]/D");
  align->Branch("timeoffsetx", timeoffsetx, "timeoffsetx[12]/D");
  align->Branch("timeoffsety", timeoffsety, "timeoffsety[12]/D");
  align->Branch("timeserrx2", timeserrx2, "timeserrx2[12]/D");
  align->Branch("timeserry2", timeserry2, "timeserry2[12]/D");
  align->Branch("xtoffset", xtoffset, "xtoffset[12][64]/D");
  align->Branch("ytoffset", ytoffset, "ytoffset[12][64]/D");
  align->Branch("xtoffystr", xtoffystr, "xtoffystr[12][64]/D");
  align->Branch("ytoffxstr", ytoffxstr, "ytoffxstr[12][64]/D");
  align->Branch("xt_slope_cor", xt_slope_cor, "xt_slope_cor[12][64][64]/D");
  align->Branch("yt_slope_cor", yt_slope_cor, "yt_slope_cor[12][64][64]/D");




  //*********************************************************************************************
  //               Iteration Starts
  //*********************************************************************************************

  fileOut->cd();

  nTotalp = nTotalt = 0;
  nTotallx=nTotally=0;


  const int timebin=60; //60 second bin
  int tmptimerate=-1;
  int tmp1evetime = -1;
  double ttmplsttime = -1.;
  int tmpoldtime=-1;
  int tmpoldmuposxrate=-1;
  int tmpoldmuposyrate=-1;
  int tmpoldmutimexrate=-1;
  int tmpoldmutimeyrate=-1;

  int ntimerate = 0;  //counter for number of events in a second
  int nmuposxrate = 0;
  int nmuposyrate = 0;
  int nmutimexrate = 0;
  int nmutimeyrate = 0;

  int nalltimerate = 0; //counter for seconds
  int nallmuposxrate = 0;
  int nallmuposyrate = 0;
  int nallmutimexrate = 0;
  int nallmutimeyrate = 0;

  //Total time of data taken in the file
  double init_time=0,end_time=0,total_time=0;
  double start_time =0.;
  double time_diff = 0.;
  double tmptime_diff = 0.;
  int tmptime=0;
  double start_time_l[nlayer] = {0.};
  double tmptime_diff_l[nlayer] = {0.};
  double xtbins[101];
  for(int ij=0;ij<101;ij++) {
    double xx = -.30+0.026*ij;
    xtbins[ij] = pow(10.,xx);
  }
  TH1F* event_time_diff = new TH1F("event_time_diff","event_time_diff",100, xtbins);
  TH1F* event_time_diff_l[nlayer];
  for(int jk=0;jk<nlayer;jk++) {
    sprintf(name,"event_time_diff_l%i",jk);
    event_time_diff_l[jk] = new TH1F(name,name,20000,0.,200.);
  }
  for (int iiter=0; iiter<nmxiter; iiter++) {

    ntotxtime = totxtime = ntotytime = totytime = 0;

    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<nstrip; jk++) {
    	xtrms[ij][jk]= ytrms[ij][jk]=0.0;
    //	xtmean[ij][jk]= ytmean[ij][jk]=-1.0;
      }
    }

    firstiter=0; // before fit, calculate time shift in each layer with all hists

    //    int ntcormx = 2;
    //Check this condition properly, why should we use this
    ntcormx =(layfirst-firstlayer>0 || lastlayer-laylast>0) ? 1 : 2;
    if (isOnlyCom) ntcormx=1;
    for (int ntcor=0; ntcor< ntcormx; ntcor++) { //jim jim ntcor
      iiterrs=nmxiter*ntcor+iiter;

      lstr = max(firstlayer,0);  //Always 0
      //      lstr = max(firstlayer,2);  //GMA140130

      lend = min(lastlayer+1,nlayerit); //nlayerit=nlayer=12     //Always 12

      if (ntcor==0) { lend = lstr+1;}

      file_out <<"lay1 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<endl;

      for (int laye= lstr; laye<lend; laye++) {
	occulyr = (ntcor==1) ? isequence[abs(laye)] : nlayer;
	//	if (occulyr!=8) continue;
	//	if (occulyr!=nlayer && occulyr!=8) continue;

       	//if (ntcor==1 && (occulyr==0||occulyr==1||occulyr==2||occulyr==3||occulyr==4 || occulyr==5|| occulyr==8 ||occulyr==10||occulyr==11)) continue; //only when RPCdaq data
	//	if (ntcor==1 && (occulyr!=0)) continue;
	//	if (ntcor==1 && (occulyr ==6 || occulyr==7|| occulyr==8 || occulyr==9||occulyr==10||occulyr==11)) continue;
	//	if(ntcor==1 && occulyr!=0) continue;
	//	if(ntcor==1 && (occulyr == 10 || occulyr == 11 || occulyr == 9 || occulyr == 6 || occulyr == 5 || occulyr == 4 || occulyr == 3 || occulyr == 2 || occulyr==1 || occulyr == 0)) continue;
	//	if(ntcor==1 && occulyr!=3) continue;
	//	if(ntcor==1 && occulyr!=8) continue;

	//if(ntcor==1 && occulyr!=3) continue; //jim jim enable this
	if(ntcor==1 && (occulyr == 10 || occulyr == 11)) continue; //jim jim disable this
	//	if(ntcor==1 && (occulyr == 0 || occulyr == 5 || occulyr == 10 || occulyr == 11)) continue;
	isfill = (iiter==nmxiter-1 && ntcor==0) ? true : false; //Fill rootuple and other histogramme for the final iteration with all layers

	ifstream file_db;
	//	file_db.open(rootfiles);
	file_db.open(argv[1]);

	file_out <<"lay2 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<" "<<laye<<" "<<occulyr<<" "<<ytimeshift<<endl;
	cout <<"lay2 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<" "<<laye<<" "<<occulyr<<" "<<ytimeshift<<endl;

	int nfile=0; //Initialisation for new histogramme booking
	while(!(file_db.eof())) {
	  file_db >> datafile>>nentrymx;
	  if (strstr(datafile,"#")) continue;
	  if(file_db.eof()) break;
	  nfile++;
	   int nx4 =0,nx5=0,nx6=0,nx7=0,ny4=0,ny5=0,ny6=0,ny7=0,nxy4=0,nxy5=0,nxy6=0,nxy7=0,nxory4=0,nxory5=0,nxory6=0,nxory7=0, nxytrig=0;
	   int ntotxext[nlayer][nTDCpLayer] ={0};//
	   int ntotyext[nlayer][nTDCpLayer] ={0};//
	   int ntxtime[nlayer][nTDCpLayer] = {0};//
	   int ntytime[nlayer][nTDCpLayer] ={0};//

	   int ntxxref[nlayer][nTDCpLayer] = {0};//--
	   int ntyyref[nlayer][nTDCpLayer] = {0};//--
	   int ntxyttref[nlayer] = {0};//
	   int ntxreftdc[nlayer] = {0};//
	   int ntyreftdc[nlayer] = {0};//


	   int ntotxext2[nlayer][nTDCpLayer] ={0};//
	   int ntotyext2[nlayer][nTDCpLayer] ={0};//
	   int ntxtime2[nlayer][nTDCpLayer] = {0};//
	   int ntytime2[nlayer][nTDCpLayer] ={0};//

	   int ntxytime2[4] = {0};//
	   int ntxytime1[4] = {0};
	   int ntxytref[4] = {0};//
	   int ntxext[nlayer][6] = {0};//
	   int ntyext[nlayer][6] = {0};//
	   int ntxref[nlayer][6] = {0};//
	   int ntyref[nlayer][6] = {0};//
	   int nttxext[nlayer][6] = {0};//
	   int nttyext[nlayer][6] = {0};//



#ifdef MONTECARLO
#ifdef FASTSIM

	 sprintf(infile, "../%s", datafile); //Standalone Monte carlo
	  //	 sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_Geant4_Data/%s", datafile);
	  //	  sprintf(infile, "/home/pethuraj/Pictures/RPCStackSim27122016_BackMuon/%s", datafile);
	  //  sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_East_West_Ayymetry_Data/%s",datafile);
	  //  sprintf(infile, "/media/9CFCFF46FCFF196A/Muon_BackScatt_DATA/%s",datafile);
	   TFile *fileIn = new TFile(infile, "read");
	  if (isfill) {
	    fileOut->cd();
	    for (int ij=0; ij<10; ij++) {
	      sprintf(name, "sel_theta_phi_%i", ij);
	      tmp_sel_theta_phi[ij] = (TH2F*)fileIn->Get(name);
	      if (nfile==1) {
		sel_theta_phi[ij] = new TH2F(name,
					     tmp_sel_theta_phi[ij]->GetTitle(),
					     tmp_sel_theta_phi[ij]->GetNbinsX(),
					     tmp_sel_theta_phi[ij]->GetXaxis()->GetXmin(),
					     tmp_sel_theta_phi[ij]->GetXaxis()->GetXmax(),
					     tmp_sel_theta_phi[ij]->GetNbinsY(),
					     tmp_sel_theta_phi[ij]->GetYaxis()->GetXmin(),
					     tmp_sel_theta_phi[ij]->GetYaxis()->GetXmax());
	      }
	      sel_theta_phi[ij]->Add(tmp_sel_theta_phi[ij]);


	    }
	    fileIn->cd();
	  }

	  //#ifdef FASTSIM
	  TTree *T1= (TTree*)fileIn->Get("T1");

	  T1->SetBranchAddress("thgen", &thgen);
	  T1->SetBranchAddress("phgen", &phgen);
	  T1->SetBranchAddress("xpgen", &xpgen);
	  T1->SetBranchAddress("ypgen", &ypgen);

	  T1->SetBranchAddress("xystrp", xystrp);
	  T1->SetBranchAddress("xstrp", xstrp);
	  T1->SetBranchAddress("ystrp", ystrp);
	  //	  T1->SetBranchAddress("xposa",  xposa);
	  //	  T1->SetBranchAddress("yposa",  yposa);
	  //	  T1->SetBranchAddress("xtrue",  xtrue);
	  //	  T1->SetBranchAddress("ytrue",  ytrue);
	  T1->SetBranchAddress("xtimea", xtimea);
	  T1->SetBranchAddress("ytimea", ytimea);

	  nentry = T1->GetEntries();
	  nentry = min(nentry,nentrymx);

	  cout <<"MC file "<< datafile<<" entries "<<nentry<<" iiter= "<<iiter<<" Occuly= "<<occulyr<<endl;
	  for(int iev=0;iev<nentry;iev++) {
	    //	    for(int iev=1746; iev<1750; iev++) {

	    xslope= xinters= yslope= yinters= timexslope= timeyslope= timex2slope= timey2slope=-100.;
	    xchi2= ychi2= xt0chi2= yt0chi2=-100.;
	    nxstrip= nystrip=Nx=Ny=nxtime=nytime=ntxyla=0;

	    for (int ix=0; ix<nlayer; ix++) {
	      for (int iy=0; iy<2*nmxiter; iy++) {
		istime_xyreso[ix][iy]=false;
		is_xyreso[ix][iy]=false;
		xposEffPass[ix][iy] = false;
		yposEffPass[ix][iy] = false;
	      }
	    }

	    fileIn->cd();
	    T1->GetEntry(iev);
	    //
         //   for (int ij=0; ij<nlayer; ij++) { Xdev[ij] = 100; Xpos[ij]=0.0;}
          //  for (int ij=0; ij<nlayer; ij++) { Ydev[ij] = 100; Ypos[ij]=0.0;}

	    vector<int> xpts[nlayer];
	    vector<int> ypts[nlayer];

	    vector<int> xptsall[nlayer];
	    vector<int> yptsall[nlayer];

	    int tmpnx=0;
          int xx[nlayer][nmxusedhits]= {-100};
	      int yy[nlayer][nmxusedhits]= {-100};
	      int xmultcount[nlayer]={0};
	      int ymultcount[nlayer]={0};
	    for (int ij=0; ij<nlayer; ij++) {
	     // Xpos[ij] = Ypos[ij] = -100;
	     for(int xxx=0;xxx<nmxusedhits;xxx++) {
	       xx[ij][xxx]= -100;
	       yy[ij][xxx]= -100;
	      }

	       xmultcount[ij]=0;
	       ymultcount[ij]=0;
	      int nxmul = xystrp[ij]%10;
	      int nymul = int(xystrp[ij]/10)%10;
	      int ix = int(xystrp[ij]/1000)%100;
	      int iy = int(xystrp[ij]/100000);
	      for(int nml=0;nml<nmxusedhits;nml++) {
	       xx[ij][nml] = xstrp[ij][nml];if(xx[ij][nml]>-100) {xmultcount[ij]++;}
	       yy[ij][nml] = ystrp[ij][nml];if(yy[ij][nml]>-100) {ymultcount[ij]++;}
	       //cout<<ij<<"		"<<xx[ij][nml]<<"	"<<yy[ij][nml]<<endl;
	      }
	     // cout<<"********************************************"<<endl;
	      //cout<<"mult"<<"	"<<xmultcount<<"	"<<ymultcount<<endl;
	      if (nxmul>0) { tmpnx++;}
	      float tmpxpos=0;
	      for (int jk=0; jk<xmultcount[ij]; jk++) {
		if (xx[ij][jk]>=nstrip) {
		  cout <<"Error in MC sample, X-strip id = "<<iev<<" "<<jk<<"+"<<ix<<" "<<xystrp[ij]<<endl;
		} else {
		//   if(( ij!=5 && ij!=6 && ij!=7  && ij!=8) ||  xx[ij][jk] <29) {
		// if(ij!= 4) {
		  xptsall[ij].push_back(xx[ij][jk]);
		  xpts[ij].push_back(xx[ij][jk]);
		  xlayer_occu[ij]->Fill(xptsall[ij][jk]);
		  tmpxpos +=(xx[ij][jk]);
		//}
		}
	//}
	      }

	      tmpxpos = tmpxpos/max(1, xmultcount[ij]) + 0.5 - xoff[ij];

	      float tmpypos = 0;
	      for (int jk=0; jk<ymultcount[ij]; jk++) {
		if (yy[ij][jk]>=nstrip) {
		  cout <<"Error in MC sample, Y-strip id = "<<iev<<" "<<jk<<"+"<<iy<<" "<<xystrp[ij]<<endl;
		} else {
		// if((ij!=5 && ij!=6 && ij!=7  && ij!=8) ||  yy[ij][jk] <29) {
		//  if( ij !=4 ) {
		  yptsall[ij].push_back(yy[ij][jk]);
		  ypts[ij].push_back(yy[ij][jk]);
		  ylayer_occu[ij]->Fill(yptsall[ij][jk]);
		  tmpypos +=(yy[ij][jk]);
		//}
	//}
		}
	      }
	      tmpypos = tmpypos/max(1, ymultcount[ij]) + 0.5 - yoff[ij];


	      //	      if ( nxmul>0 && (xposa[ij]>-10 && abs(xposa[ij] -tmpxpos)>0.1)) cout <<"XXXXXXXXXX "<<iev<<" "<<ij<<" "<<xystrp[ij]<<" "<<xtrue[ij]/3.<<" "<<xposa[ij]<<" = "<<tmpxpos<<" "<<ytrue[ij]/3.<<" "<<yposa[ij]<<" = "<<tmpypos<<endl;

	      //	      if ( nymul>0 && (yposa[ij]>-10 && abs(yposa[ij] -tmpypos)>0.1)) cout <<"YYYYYYYYYY "<<iev<<" "<<ij<<" "<<xystrp[ij]<<" "<<xtrue[ij]/3.<<" "<<xposa[ij]<<" = "<<tmpxpos<<" "<<ytrue[ij]/3.<<" "<<yposa[ij]<<" = "<<tmpypos<<endl;


	      xhits[ij] = xpts[ij].size();
	      yhits[ij] = ypts[ij].size();



#ifdef MADURAIINTM1
	      //Intermediate Madurai stack
	     	      if (ij>=5 && ij<=8 && tmpxpos>nstrip/2) {
	     // if (tmpxpos>nstrip/2) {

		nxmul = 0;
	      }
	      if (ij>=5 && ij<=8 && tmpypos>nstrip/2) {
		nymul = 0;
	      }
#endif
#ifdef MADURAIINTM2
	      //Intermediate Madurai stack
	      if (ij==4) { xhits[ij] = nymul = 0;}
	      //	      if (ij>=5 && ij<=8 && tmpxpos>nstrip/2) {
	      if (tmpxpos>nstrip/2) {
		xhits[ij] = 0;
	      }
	      if (ij>=5 && ij<=8 && tmpypos>nstrip/2) {
		nymul = 0;
	      }
#endif
#ifdef MADURAIINTM3
	      //Intermediate Madurai stack
	      if (ij==4) { xhits[ij] = yhits[ij] = 0;}
	      if (ij>=5 && ij<=8 && tmpxpos>nstrip/2) {
		xhits[ij] = 0;
	      }
	      if (ij>=5 && ij<=8 && tmpypos>nstrip/2) {
		yhits[ij] = 0;
	      }
#endif


	      if (tmpxpos >lastXstrip) { xhits[ij] = 0;}
	      if (tmpypos >lastYstrip) { yhits[ij] = 0;}
	     // cout<<"xhits/yhits"<<"		"<<ij<<"	"<<tmpxpos<<"	"<<tmpypos<<"	"<<xhits[ij]<<"		"<<yhits[ij]<<endl;


/*

	      if (nxmul>0 ) {

		Xpos[ij] =  tmpxpos; // xposa[ij]; // tmpxpos
		//if(Xpos[ij]<=0) cout<<Xpos[ij]<<endl;
		if (Xpos[ij]>-100) {xxerr[ij] = xposerrsq[nxmul-1][ij];}
	      }
	      xtime[ij] = xtimea[ij];
	      //	      timeserrx2[ij] =

	      if (nymul>0 ) {
		Ypos[ij] = tmpypos; // yposa[ij]; // tmpypos
		//if(Ypos[ij]<=0) cout<<Ypos[ij]<<endl;
		if (Ypos[ij]>-100) {yyerr[ij] = yposerrsq[nymul-1][ij];}
	      }
	      ytime[ij] = ytimea[ij];
	      //	      timeserry2[ij] =
*/
	      if (nxmul>0 && nxmul<4 && xxerr[ij]<=0.001) cout <<"xij "<<iiterrs<<" "<<ij<<" "<< nxmul<<" "<<nymul<<" "<<xposerrsq[nxmul-1][ij]<<" "<<xxerr[ij]<<endl;
	      if (nymul>0 && nymul<4 && yyerr[ij]<=0.001) cout <<"yij "<<iiterrs<<" "<<ij<<" "<< nxmul<<" "<<nymul<<" "<<yposerrsq[nymul-1][ij]<<" "<<yyerr[ij]<<endl;
	    }
	    h_tmp1xndf->Fill(tmpnx);
#else	//FULLG4SIM
	    //	  sprintf(infile, "%s", datafile);
	    //          TFile *fileIn = new TFile(infile, "read");
	    //	   sprintf(infile, "/media/9CFCFF46FCFF196A/Muon_BackScatt_DATA/%s",datafile);
	    //  sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_East_West_Ayymetry_Data/%s",datafile);
	    //  sprintf(infile, "/home/pethuraj/Documents/RPCStackSim27122016/%s",datafile);
	    //  sprintf(infile, "/home/pethuraj/Documents/RPCStackSim20170906/%s",datafile);
	    sprintf(infile, "/media/9CFCFF46FCFF196A/Geant4_Data_20171214/%s", datafile);

	   // sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_Geant4_Data/%s", datafile); //media/pethuraj/9CFCFF46FCFF196A/Madurai_Geant4_Data/
     TFile *fileIn = new TFile(infile, "read");
     if (isfill) {
	    fileOut->cd();
	    for (int ij=0; ij<10; ij++) {
	      sprintf(name, "sel_theta_phi_%i", ij);
	      tmp_sel_theta_phi[ij] = (TH2F*)fileIn->Get(name);
	      if (nfile==1) {
		sel_theta_phi[ij] = new TH2F(name,
					     tmp_sel_theta_phi[ij]->GetTitle(),
					     tmp_sel_theta_phi[ij]->GetNbinsX(),
					     tmp_sel_theta_phi[ij]->GetXaxis()->GetXmin(),
					     tmp_sel_theta_phi[ij]->GetXaxis()->GetXmax(),
					     tmp_sel_theta_phi[ij]->GetNbinsY(),
					     tmp_sel_theta_phi[ij]->GetYaxis()->GetXmin(),
					     tmp_sel_theta_phi[ij]->GetYaxis()->GetXmax());
	      }
	      sel_theta_phi[ij]->Add(tmp_sel_theta_phi[ij]);

	       for(int jk=0;jk<8;jk++) {
		sprintf(name, "sel_theta_phi_%i_%i", ij,jk);
		tmp_sel_theta_phi_8[ij][jk] = (TH2F*)fileIn->Get(name);
		if (nfile==1) {
		  sel_theta_phi_8[ij][jk] = new TH2F(name,
					       tmp_sel_theta_phi_8[ij][jk]->GetTitle(),
					       tmp_sel_theta_phi_8[ij][jk]->GetNbinsX(),
					       tmp_sel_theta_phi_8[ij][jk]->GetXaxis()->GetXmin(),
					       tmp_sel_theta_phi_8[ij][jk]->GetXaxis()->GetXmax(),
					       tmp_sel_theta_phi_8[ij][jk]->GetNbinsY(),
					       tmp_sel_theta_phi_8[ij][jk]->GetYaxis()->GetXmin(),
					       tmp_sel_theta_phi_8[ij][jk]->GetYaxis()->GetXmax());
		}
		sel_theta_phi_8[ij][jk]->Add(tmp_sel_theta_phi_8[ij][jk]);

	      }

	       for(int jk=0;jk<16;jk++) {
		sprintf(name, "sel_theta_phi_16_%i_%i", ij,jk);
		tmp_sel_theta_phi_16[ij][jk] = (TH2F*)fileIn->Get(name);
		if (nfile==1) {
		  sel_theta_phi_16[ij][jk] = new TH2F(name,
					       tmp_sel_theta_phi_16[ij][jk]->GetTitle(),
					       tmp_sel_theta_phi_16[ij][jk]->GetNbinsX(),
					       tmp_sel_theta_phi_16[ij][jk]->GetXaxis()->GetXmin(),
					       tmp_sel_theta_phi_16[ij][jk]->GetXaxis()->GetXmax(),
					       tmp_sel_theta_phi_16[ij][jk]->GetNbinsY(),
					       tmp_sel_theta_phi_16[ij][jk]->GetYaxis()->GetXmin(),
					       tmp_sel_theta_phi_16[ij][jk]->GetYaxis()->GetXmax());
		}
		sel_theta_phi_16[ij][jk]->Add(tmp_sel_theta_phi_16[ij][jk]);
	       }

	    }
	    fileIn->cd();
	  }
     //  cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

	  TTree *T1= (TTree*)fileIn->Get("T1");
	 /* EveTreeG4 *eventg4 = new EveTreeG4(T1);
      eventg4->Loop();*/
	  T1->SetBranchAddress("irun", &irun);
	  T1->SetBranchAddress("ievt", &ievt);

	  T1->SetBranchAddress("ngent", &ngent);
	  T1->SetBranchAddress("pidin", pidin);
	  T1->SetBranchAddress("ievt_wt", &ievt_wt);
	  T1->SetBranchAddress("intxn_id", &intxn_id);
	  T1->SetBranchAddress("momin", momin);
	  T1->SetBranchAddress("thein", thein);
	  T1->SetBranchAddress("thegen", thegen);
	  T1->SetBranchAddress("phiin", phiin);
	  T1->SetBranchAddress("posxin", posxin);
	  T1->SetBranchAddress("posyin", posyin);
	  T1->SetBranchAddress("poszin", poszin);
	  T1->SetBranchAddress("nsimhtx", &nsimhtx);
	  T1->SetBranchAddress("nsimhty", &nsimhty);
	  T1->SetBranchAddress("nlayert", &nlayert);
	  T1->SetBranchAddress("xdata", xdata);
	  T1->SetBranchAddress("ydata", ydata);
	  T1->SetBranchAddress("triggerinfoX", &triggerinfoX);
	  T1->SetBranchAddress("triggerinfoY", &triggerinfoY);


	  if (isTiming) {
	    T1->SetBranchAddress("digiXtime", digiXtime);
	    T1->SetBranchAddress("digiYtime", digiYtime);
	  }

	//  cout<<"nsimhtx-----------------"<<"		"<<nsimhtx<<endl;
	 // fileOut->cd();


	  nentry = T1->GetEntries();
	  nentry = min(nentry,nentrymx);
	  // if (filloccu) {ntotal += nentry;}
	  cout<<"nentry "<< nentry<<endl;
	  for(int iev=0;iev<nentry;iev++) {
	    //  cout<<"nentry "<< nentry<<" "<<iev<<" "<<endl;
	    xslope= xinters= yslope= yinters= timexslope= timeyslope= timex2slope= timey2slope=-100.;
	    xchi2= ychi2= xt0chi2= yt0chi2=-100.;
	    nxstrip= nystrip=Nx=Ny=nxtime=nytime=ntxyla=0;

	    for (int ix=0; ix<nlayer; ix++) {
	      for (int iy=0; iy<2*nmxiter; iy++) {
		istime_xyreso[ix][iy]=false;
		is_xyreso[ix][iy]=false;
		xposEffPass[ix][iy] = false;
		yposEffPass[ix][iy] = false;
	      }
	    }

	    fileIn->cd();

	    //T1->Show(iev);
        T1->GetEntry(iev);
        if (isfill) {
        fileOut->cd();
        sel_theta_phi_g4->Fill(180.*(phiin[0]/pival),180.-180.*(thein[0]/pival));
        //fileIn->cd();
	   }

	//	bool trginfo = false;
	//	trginfo = (triggerinfoX==15 || triggerinfoY==15) ? 1:0;
	//	if(triggerinfoX!=15 && triggerinfoY!=15) continue;
	//	if(!trginfo) continue;
	if(triggerinfoX!=15) continue;
	ntotal++;
        if (isfill) {
        fileOut->cd();
        sel_theta_phi_g4_trig->Fill(180.*(phiin[0]/pival),180.-180.*(thein[0]/pival));
        //fileIn->cd();
	   }
        // cout<<1<<endl;
            for (int ij=0; ij<nlayer; ij++) { Xdev[ij] = 100; Xpos[ij]=0.0;}
            for (int ij=0; ij<nlayer; ij++) { Ydev[ij] = 100; Ypos[ij]=0.0;}
      // cout<<"nsimhitx------------------"<<"		"<<
	    vector<int> xpts[nlayer];
	    vector<int> ypts[nlayer];

	    vector<int> xptsall[nlayer];
	    vector<int> yptsall[nlayer];

	    for (int ij=0; ij<nlayer; ij++) {
	      for (int jk=0; jk<nstrip; jk++) {
		if ((xdata[ij]>>jk)&0x01) {
		  xptsall[ij].push_back(jk);
		}
		if ((ydata[ij]>>jk)&0x01) {
		  yptsall[ij].push_back(jk);
		}
		if (isTiming) {
		  xtime[ij] = digiXtime[ij]*0.1;
		  ytime[ij] = digiYtime[ij]*0.1;
		}
	      }
	      //  cout<<"iev "<< iev<<" "<<ij<<" "<<" "<<xptsall[ij].size()<<" "<<yptsall[ij].size()<<" "<<xtime[ij]<<" "<<ytime[ij]<<endl;
	    }



	    xslope = 0;
	    xinters = 0;
	    yslope = 0;
	    yinters = 0;
	    ntrigX = ntrigY = 0;
	    for (int iz=0; iz<nlayer; iz++) {
	      xhits[iz] = yhits[iz] =0;
	      //Don't use it here	      if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && xptsall[iz].size()>0) {ntrigX++;} //18 01 2012
	      if (filloccu) {
		xlayer_allmult[iz]->Fill(xptsall[iz].size());

	      }
	      for (int ix=0; ix<xptsall[iz].size(); ix++) {
		if (filloccu) {
		  xlayer_alloccu->Fill(nstrip*iz+xptsall[iz][ix]);
		  xlayer_occu[iz]->Fill(xptsall[iz][ix]);

		}
		bool failed=false;
		//	if (!posinuse[0][iz][xptsall[iz][ix]]) {failed=true;}
		//		if(iz==1 && xptsall[iz][ix]>=25 && xptsall[iz][ix]<=31) {failed=true;}
		//		if(iz==2 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=7)  {failed=true;}
		//		if(iz==4 && xptsall[iz][ix]>=24 && xptsall[iz][ix]<=31)  {failed=true;}
		//		if(iz==8 && (xptsall[iz][ix]==17 || xptsall[iz][ix]==30))  {failed=true;}
		//		if(iz==10 && xptsall[iz][ix]==10)  {failed=true;}
		//		if(iz==11 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=7)  {failed=true;}
		/*
		if(iz==4 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=nstrip) {failed=true;}
		if(iz==5 && xptsall[iz][ix]>=29 && xptsall[iz][ix]<=nstrip) {failed=true;}
		if(iz==6 && xptsall[iz][ix]>=29 && xptsall[iz][ix]<=nstrip) {failed=true;}
		if(iz==7 && xptsall[iz][ix]>=29 && xptsall[iz][ix]<=nstrip) {failed=true;}
		if(iz==8 && xptsall[iz][ix]>=29 && xptsall[iz][ix]<=nstrip) {failed=true;}
		*/
		if(iz==2 && (xptsall[iz][ix]==36 || xptsall[iz][ix]==52)) {failed=true;}
		if(iz==3 && xptsall[iz][ix]==16) {failed=true;}
		if(iz==9 && xptsall[iz][ix]==47) {failed=true;}
		//	if(xptsall[iz][ix]>=nstrip/2) { failed = true;}
		if (!failed){
		  xpts[iz].push_back(xptsall[iz][ix]);
		  if (filloccu) {xlayer_alloccusel->Fill(nstrip*iz+xptsall[iz][ix]);}
		}

		//		if (failed) {xpts[iz].erase(xpts[iz].begin()+ix); ix--;} else {
		//		  if (filloccu) {xlayer_alloccusel->Fill(nstrip*iz+xpts[iz][ix]);}
		//		}
	      }
	      xhits[iz] = xpts[iz].size();
	      //      cout <<"xhits[iz] "<< iz<<" "<<xhits[iz]<<endl;
	     // if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && yptsall[iz].size()>0) { ntrigY++;} //18 01 2012
	      if (filloccu) {
		ylayer_allmult[iz]->Fill(yptsall[iz].size());
	      }
	      for (int iy=0; iy<yptsall[iz].size(); iy++) {
		// cout<<"yptsall"<<"	"<<iz<<"	"<<iy<<"	"<<yptsall[iz][iy]<<endl;
		if (filloccu) {
		  ylayer_occu[iz]->Fill(yptsall[iz][iy]);
		  ylayer_alloccu->Fill(nstrip*iz+yptsall[iz][iy]);
		}
		bool failed=false;
		//	if (!posinuse[1][iz][yptsall[iz][iy]]) {failed=true;}

		//if(iz==1 && yptsall[iz][iy]>=41 && yptsall[iz][iy]<=42) {failed=true;}
	    //if(iz==5 && yptsall[iz][iy]>=14 && yptsall[iz][iy]<=15)  {failed=true;}
	    //if(iz==5 && yptsall[iz][iy]>=21 && yptsall[iz][iy]<=22)  {failed=true;}
	    //if(iz==5 && yptsall[iz][iy]>=29 && yptsall[iz][iy]<=30)  {failed=true;}
		//if(iz==4 && xptsall[iz][ix]>=24 && xptsall[iz][ix]<=31)  {failed=true;}
		/*
		if(iz==4 && yptsall[iz][iy]>=0 && yptsall[iz][iy]<=nstrip) {failed=true;}
		if(iz==5 && yptsall[iz][iy]>=31 && yptsall[iz][iy]<=nstrip) {failed=true;}
		if(iz==6 && yptsall[iz][iy]>=31 && yptsall[iz][iy]<=nstrip) {failed=true;}
		if(iz==7 && yptsall[iz][iy]>=31 && yptsall[iz][iy]<=nstrip) {failed=true;}
		if(iz==8 && yptsall[iz][iy]>=31 && yptsall[iz][iy]<=nstrip) {failed=true;}
		*/
		if(iz==2 && yptsall[iz][iy]==9) {failed=true;}
		if(iz==3 && yptsall[iz][iy]==58) {failed=true;}
		if(iz==5 && yptsall[iz][iy]==35) {failed=true;}
		//	if(yptsall[iz][iy]>=nstrip/2) { failed = true;}
		if (!failed){
		  ypts[iz].push_back(yptsall[iz][iy]);
		  if (filloccu) {ylayer_alloccusel->Fill(nstrip*iz+yptsall[iz][iy]);}
		}
		//		if (failed) {ypts[iz].erase(ypts[iz].begin()+iy); iy--;} else {
		//		  if (filloccu) {ylayer_alloccusel->Fill(nstrip*iz+ypts[iz][iy]);}
		//		}
	      }
	      yhits[iz] = ypts[iz].size();
	    }

	    if(filloccu){
	    for (int ij=0; ij<nlayer; ij++) {
	       	xlayer_mult[ij]->Fill(xpts[ij].size());
	       	ylayer_mult[ij]->Fill(ypts[ij].size());

		rawhits_corr_xymul[ij]->Fill(xptsall[ij].size(), yptsall[ij].size());

		for (int ix=0; ix<xptsall[ij].size(); ix++) {
                  for (int iy=0; iy<yptsall[ij].size(); iy++) {
		    raw_occu[ij]->Fill(xptsall[ij][ix], yptsall[ij][iy]);
		    //  raw_occu_set[ij][iset]->Fill(xptsall[ij][ix], yptsall[ij][iy]);
		  }
		}
		for (int jk=ij+1; jk<nlayer; jk++) {
		  rawhits_xlay_corr_mul[ij][jk]->Fill(xptsall[ij].size(), xptsall[jk].size());
		  rawhits_ylay_corr_mul[ij][jk]->Fill(yptsall[ij].size(), yptsall[jk].size());
		}
	      }
	  }


#endif
#else

#ifdef OLDFORMAT
	  sprintf(infile, "daq/%s", datafile);
	  TFile *fileIn = new TFile(infile, "read");
	  TTree *event_tree= (TTree*)fileIn->Get("EVE");
	  event_tree->SetBranchAddress("InoEvents",&event);

          nentry = event_tree->GetEntries();
	  nentry = min(nentry,nentrymx);
	  cout <<"Reading data file "<< datafile<<" with entries "<<nentry<<endl;
	  for(int iev=0;iev<nentry;iev++) {    //ij is event loop variable. while checking put ij<3 or a small number.
	    fileIn->cd();
	    event_tree->GetEntry(iev);    // while running for entire event file put "ij<numentries".
	    if (filloccu) {
	      if (ntotal==0) {nsec=event->t;}
	      tmptime = tmptimerate = event->t; //-nsec;
#else
#ifdef ONEMBYONEM
	      //	  sprintf(infile, "hrdc217/%s", datafile);
	  sprintf(infile, "newdaq/%s", datafile);
#else
	  //	  sprintf(infile, "madurai64/daq/%s", datafile);
	  //  sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_Data_64/%s", datafile);
	  // sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_Merged_data_64/%s", datafile);
	  //	  sprintf(infile, "/home/kcr/esha/madurai_lab/rpcdaq64/evtblr_data/%s",datafile);
	  // sprintf(infile, "/media/9CFCFF46FCFF196A/new_ire_data/%s", datafile);
	  //  sprintf(infile, "/home/pethuraj/Documents/Daq_Data64_080716/%s", datafile);
	  // sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_Data2015/DaqData/%s",datafile);
	  // sprintf(infile, "/esha/madurai_lab/madurai_data/evt_data/%s", datafile);
	  // sprintf(infile, "/home/pethuraj/Documents/BARC_EVTBLR_data/%s",datafile);
	  // sprintf(infile, "/media/9CFCFF46FCFF196A/Mi-Ical_data/%s", datafile);
	  //sprintf(infile, "./data/%s", datafile);
	  sprintf(infile, "%s",datafile);//media/INO1_mical_SSD/Pethu_Desktop/Documents/Madurai_MiiCal_20190130/data/%s", datafile);
	   // sprintf(infile, "/media/pethuraj/Pethu_SSD2/Mini_ICAL_RPC_SCAN/%s",datafile);
	  // sprintf(infile, "/media/pethuraj/Pethu_SSD2/MiICAL_Data_201908/%s",datafile);
	   //sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_New_Data_64_20082017/%s",datafile);
#endif
	 TFile *fileIn = new TFile(infile, "read");
#ifdef NEWRPCDAQ1
#ifdef TIFRROOT
 //  TFile *fileIn = new TFile(infile, "read");
	  TTree *event_tree=(TTree*)fileIn->Get("evetree");
          RPCEve *event = new RPCEve(event_tree);
          event->Loop();
#endif
#ifdef BARCROOT
	  TTree *event_tree=(TTree*)fileIn->Get("evetree");
          BARCEve *event = new BARCEve(event_tree);
          event->Loop();
#endif

#elif defined(BARC_EVTBLR)
	  // cout<<"SSSSSS"<<endl;

	  TTree *event_tree=(TTree*)fileIn->Get("evetree");
	  event_tree->SetBranchAddress("ENum", &event->ENum);
	  event_tree->SetBranchAddress("REnum", &event->REnum);
	  event_tree->SetBranchAddress("CEnum", &event->CEnum);
	   event_tree->SetBranchAddress("Evetime", &event->EveTS);
	   printf("assigned branch for enum to evetime");
	   for (int kk = 0; kk<nlayer; kk++)
	     {
	       char x_tmp[100];
	       char y_tmp[100];
	       sprintf(x_tmp, "xstriphitsL%d", kk);
	       sprintf(y_tmp, "ystriphitsL%d", kk);

	       event_tree->SetBranchAddress(x_tmp, &event->xstriphits[kk]);
	       event_tree->SetBranchAddress(y_tmp, &event->ystriphits[kk]);
	     }
	   event_tree->SetBranchAddress("tdc_ref_l", &event->tdc_ref_l);

	   //sprintf(rtdc_tmp_t, "tdc_ref_t[%d]/l", NLAYERS);
	   event_tree->SetBranchAddress("tdc_ref_t", &event->tdc_ref_t);
	   for (int i = 0; i<nlayer; i++)
	     {
	       for (int j = 0; j<nTDCpLayer; j++)
		 {
		   char xtdc_tmp_l[100], ytdc_tmp_l[100], xtdc_tmp_t[100], ytdc_tmp_t[100];
		   sprintf(xtdc_tmp_l, "xtdc_l_%d_%d", i, j);
		   sprintf(ytdc_tmp_l, "ytdc_l_%d_%d", i, j);
		   event_tree->SetBranchAddress(xtdc_tmp_l, &event->xtdc_l[i][j]);
		   event_tree->SetBranchAddress(ytdc_tmp_l, &event->ytdc_l[i][j]);

		   sprintf(xtdc_tmp_t, "xtdc_t_%d_%d", i, j);
		   sprintf(ytdc_tmp_t, "ytdc_t_%d_%d", i, j);
		   event_tree->SetBranchAddress(xtdc_tmp_t, &event->xtdc_t[i][j]);
		   event_tree->SetBranchAddress(ytdc_tmp_t, &event->ytdc_t[i][j]);
		 }
	     }

#else
	  //  TFile *fileIn = new TFile(infile, "read");
	  TTree *event_tree=(TTree*)fileIn->Get("evetree");
          EveTree *event = new EveTree(event_tree);
          event->Loop();
#endif
#ifdef CAUDATA

	  TTree *cau_tree=(TTree*)fileIn->Get("cautree");
          CauTree *cau = new CauTree(cau_tree);
          cau->Loop();
	  vector<CauInfo> CauData;
	  CauInfo cau_db_info;
	  CauData.clear();
	  for(int jkk =0;jkk<cau_tree->GetEntries();jkk++) {
	    //  int CENum;
	    // int CauTime;
	    // double CauVal[nlayer];
	    cau_tree->GetEntry(jkk);

	    cau_db_info.CENum = cau->CauNum;
	    cau_db_info.CauTime = cau->Cautime->GetSec();

	    for(int nl =0;nl<nlayer;nl++) {
	      cau_db_info.CauVal[nl]=-10000.;
	    }
	     if(cau->CauNum%2==0){
	     if(cau->cau_tdc[0]->size()) cau_db_info.CauVal[8] = cau->cau_tdc_ref1[0]->at(0);
	     if(cau->cau_tdc[1]->size()) cau_db_info.CauVal[9] = cau->cau_tdc_ref1[1]->at(0);
	    }else {
	     if(cau->cau_tdc[0]->size()) cau_db_info.CauVal[1] = cau->cau_tdc_ref1[0]->at(0);
	     if(cau->cau_tdc[8]->size()) cau_db_info.CauVal[4] = cau->cau_tdc_ref1[8]->at(0);
	     if(cau->cau_tdc[10]->size()) cau_db_info.CauVal[5] = cau->cau_tdc_ref1[10]->at(0);
	     if(cau->cau_tdc[12]->size()) cau_db_info.CauVal[6] = cau->cau_tdc_ref1[12]->at(0);
	     if(cau->cau_tdc[14]->size()) cau_db_info.CauVal[7] = cau->cau_tdc_ref1[14]->at(0);
	      // cau_db_info.CauVal[9] = cau_tdc_ref1[1]->at(0);
	    }

	     // for(int ly=0;ly<nlayer;ly++) {
	     //    cout<<"CAU Input"<<"   "<<cau->CauNum<<"    "<<ly<<"   "<<cau_db_info.CauVal[ly]<<endl;
	     //   }
	    CauData.push_back(cau_db_info);
	  }


#elif defined(CAUDATA1)
	  TTree *cau_tree=(TTree*)fileIn->Get("cautree");
          CauTree *cau = new CauTree(cau_tree);
          cau->Loop();
	  vector<CauInfo> CauData;
	  CauInfo cau_db_info;
	  CauData.clear();
	  for(int jkk =0;jkk<cau_tree->GetEntries();jkk++) {
	    //  int CENum;
	    // int CauTime;
	    // double CauVal[nlayer];
	    cau_tree->GetEntry(jkk);

	    cau_db_info.CENum = cau->CauNum;
	    cau_db_info.CauTime = cau->Cautime->GetSec();
	    //  cout<<jkk<<" "<<"cautime"<<"  "<<cau->Cautime->GetSec()<<"  "<< cau->CauNum<<endl;
	    cau_db_info.selectline = cau->select_line;
	    for(int nl =0;nl<nlayer;nl++) {
	      cau_db_info.CauVal[nl]=-1000.;
	    }
	     if(cau->select_line==0){
	     if(cau->cau_tdc_l[0]->size()) cau_db_info.CauVal[0] = cau->cau_tdc_ref1_l[0]->at(0);
	     if(cau->cau_tdc_l[2]->size()) cau_db_info.CauVal[1] = cau->cau_tdc_ref1_l[2]->at(0);
	     if(cau->cau_tdc_l[4]->size()) cau_db_info.CauVal[2] = cau->cau_tdc_ref1_l[4]->at(0);
	     if(cau->cau_tdc_l[6]->size()) cau_db_info.CauVal[3] = cau->cau_tdc_ref1_l[6]->at(0);
	     if(cau->cau_tdc_l[8]->size()) cau_db_info.CauVal[4] = cau->cau_tdc_ref1_l[8]->at(0);
	     if(cau->cau_tdc_l[10]->size()) cau_db_info.CauVal[5] = cau->cau_tdc_ref1_l[10]->at(0);
	     if(cau->cau_tdc_l[12]->size()) cau_db_info.CauVal[6] = cau->cau_tdc_ref1_l[12]->at(0);
	     if(cau->cau_tdc_l[14]->size()) cau_db_info.CauVal[7] = cau->cau_tdc_ref1_l[14]->at(0);
	    }else {
	     if(cau->cau_tdc_l[0]->size()) cau_db_info.CauVal[8] = cau->cau_tdc_ref1_l[0]->at(0);
	     if(cau->cau_tdc_l[1]->size()) cau_db_info.CauVal[9] = cau->cau_tdc_ref1_l[1]->at(0);
	     if(cau->cau_tdc_l[2]->size()) cau_db_info.CauVal[10] = cau->cau_tdc_ref1_l[2]->at(0);
	     if(cau->cau_tdc_l[3]->size()) cau_db_info.CauVal[11] = cau->cau_tdc_ref1_l[3]->at(0);
	      // cau_db_info.CauVal[9] = cau_tdc_ref1[1]->at(0);
	    }

	     // for(int ly=0;ly<nlayer;ly++) {
	     //    cout<<"CAU Input"<<"   "<<cau->CauNum<<"    "<<ly<<"   "<<cau_db_info.CauVal[ly]<<endl;
	     //   }
	    CauData.push_back(cau_db_info);
	  }

#endif



          nentry = event_tree->GetEntries();


	  //	  nentry = 2;
	  nentry = min(nentry,nentrymx);
	 // if(nentry == 0)
// 	  if (filloccu) {
// 	    for (int il=0; il<nlayer; il++) {
// 	      sprintf(name, '''');
// 	      hist[ij] = new TH1F("trig_diff[ij]", "name", nentry, -0.5, nentry-0.5);

// 	      hist[ij] = new TH1F("counterf[ij]", "name", nentry, -0.5, nentry-0.5);
// 	    }
// 	  }

	  int cauentry =1;//=1 if cau eve number starts from 1
	                  // = 0 if cau eve number starts from 2


	  // for(int ievt =1 ;ievt<nentry;ievt++) {
	  // 	    event_tree->GetEntry(ievt);
	  // 	    if(ievt == 1) { init_time = event->Evetime->GetSec();/*AsGMST()*3600.);*/ cout<<"init_time"<<"		"<<init_time<<endl; }
	  // 	    if(ievt == nentry-1) { end_time = event->Evetime->GetSec();/*AsGMST()*3600.;*/
	  // 	      cout<<"end_time"<<"		"<<end_time<<endl;
	  // 	      total_time +=  end_time -init_time ;

	  // 	      cout<<"total_time"<<"			"<<end_time -init_time <<"		"<<total_time<<endl;
	  // 	    }
	  // 	  }


	  cout <<infile<<" has "<< nentry<<" events "<<endl;
	  //  int nx4 =0,nx5=0,nx6=0,nx7=0,ny4=0,ny5=0,ny6=0,ny7=0,nxy4=0,nxy5=0,nxy6=0,nxy7=0,nxory4=0,nxory5=0,nxory6=0,nxory7=0, nxytrig=0;


	  double initevttime, lsttime ;
          for(int iev=0/*0*/;iev<nentry;iev++) {    //ij is event loop variable. while checking put ij<3 or a small number.

	    fileIn->cd();
	    event_tree->GetEntry(iev);    // while running for entire event file put "ij<numentries".
	    ievt = iev;//eventnumbxx[iev];
	    //  cout<<"YYYYYYYYYYYYYYYYYY"<<endl;
	    if (filloccu) {
#ifdef BARC_EVTBLR
	      // int tmptime = 0;

	        if (ntotal==0) {nsec=event->EveTS->GetSec(); cout<<"nsec"<<"	"<<nsec<<endl;}

	      if(iev == 1) { init_time = event->EveTS->GetSec();; }
	      if(iev == nentry-1) { end_time = event->EveTS->GetSec();
	      	total_time +=  end_time -init_time ;
	      	////total_time + =
	      	cout<<"total_time"<<"			"<<end_time -init_time <<"		"<<total_time<<endl;
	      }

	       tmptime = tmptimerate = event->EveTS->GetSec(); //-nsec;
	       if(ntotal==0) { tmp1evetime = tmptimerate;}
	      tmptime_diff = 1.e+3*event->EveTS->GetSec()+1.e-6*event->EveTS->GetNanoSec()-start_time;
	      start_time = 1.e+3*event->EveTS->GetSec()+1.e-6*event->EveTS->GetNanoSec();

#else
	      // if (ntotal==0) {nsec=event->Evetime->GetSec(); cout<<"nsec"<<"	"<<nsec<<endl;}

	      // if(iev == 1) { init_time = event->Evetime->GetSec();/*(event->Evetime->AsGMST()*3600.)*/; }
	      // if(iev == nentry-1) { end_time = event->Evetime->GetSec();/* event->Evetime->AsGMST()*3600.;*/
	      // 	total_time +=  end_time -init_time ;
	      // 	////total_time + =
	      // 	cout<<"total_time"<<"			"<<end_time -init_time <<"		"<<total_time<<endl;
	      // }

	      // int tmptime = tmptimerate = event->Evetime->GetSec(); //-nsec;
	      // tmptime_diff = 1.e+3*event->Evetime->GetSec()+1.e-6*event->Evetime->GetNanoSec()-start_time;
	      // start_time = 1.e+3*event->Evetime->GetSec()+1.e-6*event->Evetime->GetNanoSec();
	      double tmpevttime[nlayer] = {0.};
	      double tmplsttime[nlayer] = {0.};
#ifdef TIFRROOT
	      for(int jj=0;jj<nlayer;jj++) {
	      tmpevttime[jj] = 1.e+3*event->EveTS[jj]->GetSec()+1.e-6*event->EveTS[jj]->GetNanoSec();
	      tmplsttime[jj] = event->EveTS[jj]->AsGMST()+5.30;
	      }
	      int ntstamp=0;
	      for(int jkl=0;jkl<nlayer-2;jkl++) {
		if(tmpevttime[jkl]>0) {
		initevttime = tmpevttime[jkl] ;
		ntstamp++;
		break;
		}
	      }
	      for(int jkl =ntstamp;jkl<nlayer-2;jkl++) {
		//	if(jkl==2) continue;
		if(tmpevttime[jkl]<=0) continue;
		if(tmpevttime[jkl]<initevttime) {
		  initevttime = tmpevttime[jkl];
		  lsttime = tmplsttime[jkl];
		}
	      }




	      for(int kl=0;kl<nlayer-2;kl++) {
		//	if(kl==2) continue;
		tmptime_diff_l[kl] =  tmpevttime[kl]-start_time_l[kl];
		start_time_l[kl] = tmpevttime[kl];
      		event_time_diff_l[kl]->Fill(tmptime_diff_l[kl]);
	      }


	      tmptime = tmptimerate = 1.e-3*initevttime;//event->Evetime->GetSec(); //-nsec;
	      ttmplsttime = lsttime;
	      if(ntotal==0) { tmp1evetime = tmptimerate;}
	      //cout<<"tmptime"<<"		"<<iev<<"		"<<tmptime<<endl;
	      //	      tmptime_diff = event->Evetime->AsGMST()*3600.*1000000.-start_time;
	      //	      start_time = event->Evetime->AsGMST()*3600.*1000000.;
	      tmptime_diff = initevttime -start_time;
	      start_time = initevttime;



#elif defined(BARCROOT)
	      //  if (ntotal==0) {nsec=event->Evetime->GetSec(); cout<<"nsec"<<"	"<<nsec<<endl;}
	      if(ntotal==0) { tmp1evetime = tmptimerate;}
	      if(iev == 1) { init_time = event->Evetime->GetSec(); }
	      if(iev == nentry-1) { end_time = event->Evetime->GetSec();
	      	total_time +=  end_time -init_time ;
	      	////total_time + =
	      	cout<<"total_time"<<"			"<<end_time -init_time <<"		"<<total_time<<endl;
	      }

	      tmptime = tmptimerate = event->Evetime->GetSec(); //-nsec;
	      if(ntotal==0) { tmp1evetime = tmptimerate;}
	      tmptime_diff = 1.e+3*event->Evetime->GetSec()+1.e-6*event->Evetime->GetNanoSec()-start_time;
	      start_time = 1.e+3*event->Evetime->GetSec()+1.e-6*event->Evetime->GetNanoSec();

#endif

	      // tmpevttime[1] = 1.e+3*event->Evetime_1->GetSec()+1.e-6*event->Evetime_1->GetNanoSec();
	      // tmpevttime[2] = 1.e+3*event->Evetime_2->GetSec()+1.e-6*event->Evetime_2->GetNanoSec();
	      // tmpevttime[3] = 1.e+3*event->Evetime_3->GetSec()+1.e-6*event->Evetime_3->GetNanoSec();
	      // tmpevttime[4] = 1.e+3*event->Evetime_4->GetSec()+1.e-6*event->Evetime_4->GetNanoSec();
	      // tmpevttime[5] = 1.e+3*event->Evetime_5->GetSec()+1.e-6*event->Evetime_5->GetNanoSec();
	      // tmpevttime[6] = 1.e+3*event->Evetime_6->GetSec()+1.e-6*event->Evetime_6->GetNanoSec();
	      // tmpevttime[7] = 1.e+3*event->Evetime_7->GetSec()+1.e-6*event->Evetime_7->GetNanoSec();
	      // tmpevttime[8] = 1.e+3*event->Evetime_8->GetSec()+1.e-6*event->Evetime_8->GetNanoSec();
	      // tmpevttime[9] = 1.e+3*event->Evetime_9->GetSec()+1.e-6*event->Evetime_9->GetNanoSec();
	      // tmpevttime[10] = 1.e+3*event->Evetime_10->GetSec()+1.e-6*event->Evetime_10->GetNanoSec();
	      // tmpevttime[11] = 1.e+3*event->Evetime_11->GetSec()+1.e-6*event->Evetime_11->GetNanoSec();












#endif
	      event_time_diff->Fill(tmptime_diff);
	      //event_time_diff_1->Fill(tmptime_diff);

#endif
	      if (iev>0 && tmptimerate >=tmpoldtime+timebin) {
		file_out <<"totrate "<<datafile<<" "<<iev<<" "<<nalltimerate<<" "<<ntimerate<<" "<<tmptimerate<<endl;
		trigrate->Fill(ntimerate*(1.0/timebin));
		trigratex->Fill(nalltimerate, ntimerate*(1.0/timebin));
		ntimerate = 0;
		nalltimerate++;
		tmpoldtime = tmptimerate;
	      }
	      ntimerate++;

	      ntotal++;



	      iset = int(ntotal/nCount);
	      if (iset >=nsetmx) { iset = nsetmx-1;}
	      if (isetold != iset) {
		fileOut->cd();
		TDatime T0x(tmptime);
		sprintf(name, "strp_count_set_%i", iset);
		sprintf(title, "strp_count_set_%i [%s]", iset, T0x.AsString());
		strp_count_set[iset] = new TH2F(name, title, 2*nlayer, -0.5,  2*nlayer-0.5, nstrip, -0.5, nstrip-0.5);
		for (int ix=0; ix<nlayer; ix++) {
		  sprintf(name, "strp_xmult_set_l%i_%i", ix, iset);
		  sprintf(title, "strp_xmult_set_l%i_%i [%s]", ix, iset, T0x.AsString());
		  strp_xmult_set[ix][iset] = new TH1F(name, title, nstrip+1, -0.5, nstrip+0.5);

		  sprintf(name, "strp_ymult_set_l%i_%i", ix, iset);
		  sprintf(title, "strp_ymult_set_l%i_%i [%s]", ix, iset, T0x.AsString());
		  strp_ymult_set[ix][iset] = new TH1F(name, title, nstrip+1, -0.5, nstrip+0.5);

		  sprintf(name, "raw_occu_set_l%i_%i", ix, iset);
		  sprintf(title, "raw_occu_set_l%i_%i [%s]", ix, iset, T0x.AsString());
		  raw_occu_set[ix][iset] = new TH2F(name, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

		  sprintf(name, "xlayer_reso_set_l%i_i%i", ix, iset);
		  sprintf(title, "xlayer_reso_set_l%i_i%i [%s]", ix, iset, T0x.AsString());
		  xlayer_reso_set[ix][iset]=new TH1F(name, title, 60, -1.8, 1.8);

		  sprintf(name, "ylayer_reso_set_l%i_i%i", ix, iset);
		  sprintf(title, "ylayer_reso_set_l%i_i%i [%s]", ix, iset, T0x.AsString());
		  ylayer_reso_set[ix][iset]=new TH1F(name, title, 60, -1.8, 1.8);

		  sprintf(name, "time_xreso_set_l%i_i%i", ix, iset);
		  sprintf(title, "time_xreso_set_l%i_i%i [%s]", ix, iset, T0x.AsString());
		  time_xreso_set[ix][iset] = new TH1F(name, title, 60, -7.5, 7.5);

		  sprintf(name, "time_yreso_set_l%i_i%i", ix, iset);
		  sprintf(title, "time_yreso_set_l%i_i%i [%s]", ix, iset, T0x.AsString());
		  time_yreso_set[ix][iset] = new TH1F(name, title, 60, -7.5, 7.5);




		}
		isetold = iset;
		fileIn->cd();
	      }

	      itset = int((tmptimerate-tmp1evetime)/ntCount);
	      if (itset >=ntsetmx) { itset = ntsetmx-1;}
	      if (itsetold != itset) {

		fileOut->cd();
		TDatime T0x(tmptimerate);
		for (int ix=0; ix<nlayer; ix++) {
		  sprintf(name,"totalentry_l%i_itset%i", ix, itset);
		sprintf(title, "totalentry_l%i[%s]", ix, T0x.AsString());
		totalentry_set[ix][itset]=new TH2F(name, title, nstrip/2, -0.5, nstrip-0.5, nstrip/2, -0.5, nstrip-0.5);

	     	sprintf(name, "inefficiency_corx_l%i_itset%i", ix, itset);
		sprintf(title, "inefficiency_corx_l%i[%s]", ix, T0x.AsString());
		inefficiency_corx_set[ix][itset]=new TH2F(name, title, nstrip/2, -.5, nstrip -0.5, nstrip/2, -0.5, nstrip -0.5);

		sprintf(name, "triggereffi_x_l%i_itset%i", ix, itset);
		sprintf(title, "triggereffi_x_l%i[%s]", ix, T0x.AsString());
		triggereffi_x_set[ix][itset]=new TH2F(name, title, nstrip/2, -0.5, nstrip-0.5, nstrip/2, -0.5, nstrip-0.5);

		sprintf(name, "triggereffi_y_l%i_itset%i", ix, itset);
		sprintf(title, "triggereffi_y_l%i[%s]", ix, T0x.AsString());
		triggereffi_y_set[ix][itset]=new TH2F(name, title, nstrip/2, -0.5, nstrip-0.5, nstrip/2, -0.5, nstrip-0.5);
		}
		itsetold = itset;
	       	fileIn->cd();
	      }


	    }



	    //  cout<<iev<<"  "<<"PASSSSSSSSSSSS 1"<<endl;
	    if (iev%20000==1) {
#ifdef OLDFORMAT
	      cout <<"iev "<< iev<<" "<<event->ENum<<" "<<" "<<event->t.GetDate()<<" "<<event->t.GetTime()<<" "<<event->t.GetSec()<<" "<<event->t.GetSec()<<" "<<event->t.GetTime(0, nsec)<<" "<<nsec<<" "<<event->t.GetSec()-nsec<<endl;
#else
	      /*
	      cout <<"iev "<< iev<<" "<<event->ENum<<" "<<" "<<event->Evetime->GetDate()<<" "<<event->Evetime->GetTime()<<" "<<event->Evetime->GetSec()<<" "<<event->Evetime->GetSec()<<" "<<event->Evetime->GetTime(0, nsec)<<" "<<nsec<<" "<<event->Evetime->GetSec()-nsec<<endl;
	      */
#endif
	    }

	    //if(iev == 0) { init_time = event->Evetime->GetSec(); }
	      //if(iev == nentry-1) { end_time = event->Evetime->GetSec();
			  //total_time  =  end_time -init_time ;
			  //cout<<"total_time"<<"			"<<total_time<<endl; }
	    //	    nsec=event->Evetime->GetSec();
	    //	    double layerdx[nlayer], layerdy[nlayer];
	    //	    for (int ij=0; ij<nlayer; ij++) { layerdx[ij] = layerdy[ij]=-10;}

	    xslope= xinters= yslope= yinters= timexslope= timeyslope= timex2slope= timey2slope=-100.;
	    xchi2= ychi2= xt0chi2= yt0chi2=-100.;
	    nxstrip= nystrip=Nx=Ny=nxtime=nytime=ntxyla=0;

	    for (int ix=0; ix<nlayer; ix++) {
	      for (int iy=0; iy<2*nmxiter; iy++) {
		istime_xyreso[ix][iy]=false;
		is_xyreso[ix][iy]=false;
		xposEffPass[ix][iy] = false;
		yposEffPass[ix][iy] = false;
	      }
	    }

	    vector<int> xpts[nlayer];
	    vector<int> ypts[nlayer];

	    vector<int> xptsall[nlayer]; //For trigger criteria
	    vector<int> yptsall[nlayer];

	    fileOut->cd();

#ifdef OLDFORMAT
	    BID=0; EvtIndx = 0;
	    for(int jk=0;jk<nlayer;jk++) {
	      for (int kl=0; kl<nTDCpLayer; kl++) {
		xallhits[jk][kl] = yallhits[jk][kl] =-1;
	      }
	      xpts[jk].clear(); ypts[jk].clear();  xptsall[jk].clear(); yptsall[jk].clear();

	      BID = (( event->Evedata[EvtIndx]>>10)&0x0F);
	      if(event->Evedata[EvtIndx]&0x8000){// ychannel
		//cout<<"Wrong X side"<<endl;
	      } else { // x channel
                for(int kl=15; kl>=0; kl--){
		  if(event->Evedata[EvtIndx+1] >> kl & 1){
		    xptsall[BID].push_back(kl+16);
		  }
		}
                for(int kl=15; kl>=0; kl--){
		  if(event->Evedata[EvtIndx+2] >> kl & 1){
		    xptsall[BID].push_back(kl);
		  }
		}
	      }
	      EvtIndx+=3;
	    }

	    EvtIndx = 0;
            for(int jk=0;jk<nlayer;jk++) {
	      BID=(( event->Evedata[EvtIndx+48]>>10)&0x0F);
	      if(event->Evedata[EvtIndx+48]&0x8000) {// y channel
                for(int kl=15; kl>=0; kl--) {
		  if(event->Evedata[EvtIndx+49] >> kl & 1){
		    yptsall[BID].push_back(kl+16);
		  }
		}
                for(int kl=15; kl>=0; kl--){
		  if(event->Evedata[EvtIndx+50] >> kl & 1){
		    yptsall[BID].push_back(kl);
		  }
		}
	      } else { // y channel
		//        cout<<"Wrong Y side"<<endl;
	      }
	      EvtIndx+=3;
	    }
	    // Store total number of hits in X/Y layer irrespective of noise etc.
	    for (int iz=0; iz<nlayer; iz++) {
	      xallhits[iz][0] = (xptsall[iz].size()==1) ? xptsall[iz][0] : -10-xptsall[iz].size(); //avoid confusion with no hits and only hits at strip# 1
	      yallhits[iz][0] = (yptsall[iz].size()==1) ? yptsall[iz][0] : -10-yptsall[iz].size();
	    }

#else
	    vector<int> xptsalltdc[nlayer][nTDCpLayer]; //signal for each TDC
	    vector<int> yptsalltdc[nlayer][nTDCpLayer];
	    for(int jk=0;jk<nlayer-2;jk++) {
	      for (int kl=0; kl<nTDCpLayer; kl++) {
		xallhits[jk][kl] = yallhits[jk][kl] =-1;
		xptsalltdc[jk][kl].clear();  yptsalltdc[jk][kl].clear();
	      }
	      xpts[jk].clear(); ypts[jk].clear();  xptsall[jk].clear(); yptsall[jk].clear();
	      //printf("layer\t");
	      for(int kl=0; kl<nstrip; kl++) {
#ifdef BARC_EVTBLR
		if(event->xstriphits[jk]->TestBitNumber(kl)) {
#else
		if(event->xLayer[jk]->TestBitNumber(kl)) {
#endif
		    xptsall[jk].push_back(kl);
#if defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
		  xptsalltdc[jk][kl%8].push_back(kl);
#else
		  xptsalltdc[jk][int((nTDCpLayer*kl)/nstrip)].push_back(kl);
#endif
		}
               // cout<<endl;
#ifdef BARC_EVTBLR
		if(event->ystriphits[jk]->TestBitNumber(kl)) {
#else
		if(event->yLayer[jk]->TestBitNumber(kl)) {
#endif
		  yptsall[jk].push_back(kl);
#if defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
		   yptsalltdc[jk][kl%8].push_back(kl);
#else
		  yptsalltdc[jk][int((nTDCpLayer*kl)/nstrip)].push_back(kl);
#endif
		}
	      }
	      //cout<<endl;
	      //	      cout<<"iev "<<iev<<" "<< jk<<" "<<xptsall[jk].size()<<" "<<yptsall[jk].size()<<endl;
	    }





	    // Store total number of hits in X/Y layer irrespective of noise etc.
	    for (int iz=0; iz<nlayer; iz++) {
	      for (int itdc=0; itdc<nTDCpLayer; itdc++) {
		xallhits[iz][itdc] = (xptsalltdc[iz][itdc].size()==1) ? xptsalltdc[iz][itdc][0] : -10-xptsalltdc[iz][itdc].size(); //avoid confusion with no hits and only hits at strip# 1
		yallhits[iz][itdc] = (yptsalltdc[iz][itdc].size()==1) ? yptsalltdc[iz][itdc][0] : -10-yptsalltdc[iz][itdc].size();

		//	      cout <<"xallpts "<<iz<<" "<<xallhits[iz]<<" "<<xptsall[iz].size()<<" "<<yallhits[iz]<<" "<<yptsall[iz].size()<<endl;
	      }
	    }
#endif
	    ////////////////////////////////////////////
	    //
	    //  First clean up noisey layers and then
	    //  Strip efficiency and resolutions etc etc
	    //
	    ////////////////////////////////////////////
	    // cout<<iev<<"  "<<"PASSSSSSSSSSSS 2"<<endl;
	    xslope = 0;
	    xinters = 0;
	    yslope = 0;
	    yinters = 0;
	    ntrigX = ntrigY = 0;

	    for (int iz=0; iz<nlayer; iz++) {
	      xhits[iz] = yhits[iz] =0;
	      //Don't use it here	      if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && xptsall[iz].size()>0) {ntrigX++;} //18 01 2012
	      if (filloccu) {
		xlayer_allmult[iz]->Fill(xptsall[iz].size());
#ifndef MONTECARLO
		strp_xmult_set[iz][iset]->Fill(xptsall[iz].size());
#endif
	      }
	      for (int ix=0; ix<xptsall[iz].size(); ix++) {
		if (filloccu) {
		  xlayer_alloccu->Fill(nstrip*iz+xptsall[iz][ix]);
		  xlayer_occu[iz]->Fill(xptsall[iz][ix]);
#ifndef MONTECARLO
		  strp_count_set[iset]->Fill(iz, xptsall[iz][ix]);
#endif
		}
		bool failed=false;
		if (!posinuse[0][iz][xptsall[iz][ix]]) {failed=true;}
		//	if(iz==3 && xptsall[iz][ix]==57) {failed=true;}
		//	if(iz==7 && (xptsall[iz][ix]==34 || xptsall[iz][ix]==35 || xptsall[iz][ix]==36)) {failed=true;}
		//if(iz==8 && (xptsall[iz][ix]==42 || xptsall[iz][ix]==43)) {failed=true;}
		//	if(iz==9 && xptsall[iz][ix]==24) {failed=true;}


		// if(iz==1 && xptsall[iz][ix]>=25 && xptsall[iz][ix]<=31) {failed=true;}
				// if(iz==2 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=7)  {failed=true;}
				// if(iz==4 && xptsall[iz][ix]>=24 && xptsall[iz][ix]<=31)  {failed=true;}
				// if(iz==8 && (xptsall[iz][ix]==17 || xptsall[iz][ix]==30))  {failed=true;}
				// if(iz==10 && xptsall[iz][ix]==10)  {failed=true;}
				// if(iz==11 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=7)  {failed=true;}


		//	if(iz==3 && xptsall[iz][ix]>=49 && xptsall[iz][ix]<=50) {failed=true;}
		//	if(iz==4 && xptsall[iz][ix]>=11 && xptsall[iz][ix]<=13) {failed=true;}
		//	if(iz==4 && xptsall[iz][ix]>=52 && xptsall[iz][ix]<=54) {failed=true;}
		//	if(iz==5 && xptsall[iz][ix]>=8 && xptsall[iz][ix]<=9) {failed=true;}
		//	if(iz==8 && xptsall[iz][ix]>=51 && xptsall[iz][ix]<=53) {failed=true;}


		//	if(iz==4 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=nstrip) {failed=true;}
		//	if(iz==5 && xptsall[iz][ix]>=28 && xptsall[iz][ix]<=nstrip) {failed=true;}
		//if(iz==6 && xptsall[iz][ix]>=28 && xptsall[iz][ix]<=nstrip) {failed=true;}
		//if(iz==7 && xptsall[iz][ix]>=28 && xptsall[iz][ix]<=nstrip) {failed=true;}
		//if(iz==8 && xptsall[iz][ix]>=28 && xptsall[iz][ix]<=nstrip) {failed=true;}


	        //	if(xptsall[iz][ix]>=nstrip/2) { failed = true;}
		/*
		if(iz==0 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=1) {failed=true;}
		if(iz==1 && xptsall[iz][ix]>=1 && xptsall[iz][ix]<=3)  {failed=true;}
		if(iz==2 && xptsall[iz][ix]>=3 && xptsall[iz][ix]<=4)  {failed=true;}
		if(iz==3 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=2)  {failed=true;}
		if(iz==4 && (xptsall[iz][ix]==0 || xptsall[iz][ix]==37))  {failed=true;}
		if(iz==5 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=2)  {failed=true;}
		if(iz==5 && xptsall[iz][ix]==20 )  {failed=true;}
		if(iz==6 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=2)  {failed=true;}
		if(iz==7 && (xptsall[iz][ix]==0 || xptsall[iz][ix]==30))  {failed=true;}
		if(iz==8 && xptsall[iz][ix]==0 )  {failed=true;}
		if(iz==9 && (xptsall[iz][ix]>=0 && xptsall[iz][ix]<=1))  {failed=true;}
		if(iz==9 && xptsall[iz][ix]==22)  {failed=true;}
		if(iz==10 && (xptsall[iz][ix]==0 ))  {failed=true;}
		if(iz==11 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=3)  {failed=true;}
		if(iz==11 && xptsall[iz][ix]>=28 && xptsall[iz][ix]<=29)  {failed=true;}
		if(iz==11 && xptsall[iz][ix]>=56 && xptsall[iz][ix]<=59)  {failed=true;}
		*/

		//shift position
#if defined(SHIFT0)&&defined(NEWRPCDAQ1)
		if (iz==3) {
		  for (int xix=33; xix<40; xix++) {
		    if (xptsall[iz][ix]==xix) { xptsall[iz][ix]--;}
		  }
		}

		// if(iz==8 && xptsall[iz][ix]==53) {    //rejection of ystrip = 21 when ystrip 45 is having hit
		//   for(int ixx=0;ixx<xptsall[iz].size();ixx++){
		//     if(xptsall[iz][ixx]==51) {
		//       failed=true;
		//       break;
		//     }
		//   }
		// }
#endif
		if (!failed){
		  xpts[iz].push_back(xptsall[iz][ix]);
		  if (filloccu) {xlayer_alloccusel->Fill(nstrip*iz+xptsall[iz][ix]);}
		}

		//		if (failed) {xpts[iz].erase(xpts[iz].begin()+ix); ix--;} else {
		//		  if (filloccu) {xlayer_alloccusel->Fill(nstrip*iz+xpts[iz][ix]);}
		//		}
	      }
	      xhits[iz] = xpts[iz].size();
	      //	      cout <<"xhits[iz] "<< iz<<" "<<xhits[iz]<<endl;
	     // if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && yptsall[iz].size()>0) { ntrigY++;} //18 01 2012
	      if (filloccu) {
		ylayer_allmult[iz]->Fill(yptsall[iz].size());
#ifndef MONTECARLO
		strp_ymult_set[iz][iset]->Fill(yptsall[iz].size());
#endif
	      }
	      for (int iy=0; iy<yptsall[iz].size(); iy++) {
		// cout<<"yptsall"<<"	"<<iz<<"	"<<iy<<"	"<<yptsall[iz][iy]<<endl;
		if (filloccu) {
		  ylayer_occu[iz]->Fill(yptsall[iz][iy]);
		  ylayer_alloccu->Fill(nstrip*iz+yptsall[iz][iy]);
#ifndef MONTECARLO
		  strp_count_set[iset]->Fill(nlayer+iz, yptsall[iz][iy]);
#endif
		}

		bool failed=false;
		if (!posinuse[1][iz][yptsall[iz][iy]]) {failed=true;}
		//	if (iz==4 && (yptsall[iz][iy]==37 || yptsall[iz][iy]==38)) {failed=true;}
		//	if (iz==6 && yptsall[iz][iy]==21) {failed=true;}
		// if(iz==0 && yptsall[iz][iy]==21) {    //rejection of ystrip = 21 when ystrip 45 is having hit
		//   for(int iyy=iy+1;iyy<yptsall[iz].size();iyy++){
		//     if(yptsall[iz][iyy]==45) {
		//       failed=true;
		//       break;
		//     }
		//   }
		// }

		// if(iz==6 && yptsall[iz][iy]==62) {    //rejection of ystrip = 62 when ystrip 56 is having hit
		//   for(int iyy=0;iyy<yptsall[iz].size();iyy++){
		//     if(yptsall[iz][iyy]==56) {
		//       failed=true;
		//       break;
		//     }
		//   }
		// }
		//	if(iz==4 && yptsall[iz][iy]>=0 && yptsall[iz][iy]<=nstrip) {failed=true;}
		//	if(iz==5 && yptsall[iz][iy]>=31 && yptsall[iz][iy]<=nstrip) {failed=true;}
		//	if(iz==6 && yptsall[iz][iy]>=31 && yptsall[iz][iy]<=nstrip) {failed=true;}
		//	if(iz==7 && yptsall[iz][iy]>=31 && yptsall[iz][iy]<=nstrip) {failed=true;}
		//	if(iz==8 && yptsall[iz][iy]>=31 && yptsall[iz][iy]<=nstrip) {failed=true;}

		//	if(yptsall[iz][iy]>=nstrip/2) { failed = true;}
		/*
		if(iz==1 && yptsall[iz][iy]==4) {failed=true;}
		if(iz==2  && yptsall[iz][iy]==45) {failed=true;}
		if(iz==3 && yptsall[iz][iy]>=3 && yptsall[iz][iy]<=5)  {failed=true;}
		if(iz==4  && yptsall[iz][iy]==0) {failed=true;}
		if(iz==4  && yptsall[iz][iy]==28) {failed=true;}
		if(iz==6  && (yptsall[iz][iy]==23 ||yptsall[iz][iy]==57 )) {failed=true;}
		if(iz==7  && (yptsall[iz][iy]==0 ||yptsall[iz][iy]==20 ||yptsall[iz][iy]==23 )) {failed=true;}
		if(iz==8 && yptsall[iz][iy]>=53 && yptsall[iz][iy]<=55)  {failed=true;}
		if(iz==8  && yptsall[iz][iy]==58) {failed=true;}
		if(iz==11 && yptsall[iz][iy]>=10 && yptsall[iz][iy]<=11)  {failed=true;}
		if(iz==11 && yptsall[iz][iy]>=48 && yptsall[iz][iy]<=51)  {failed=true;}
		*/

		//cout<<"After"<<
		if (!failed) {
		  ypts[iz].push_back(yptsall[iz][iy]);
		  if (filloccu) {ylayer_alloccusel->Fill(nstrip*iz+yptsall[iz][iy]);}
		}
		//		if (failed) {ypts[iz].erase(ypts[iz].begin()+iy); iy--;} else {
		//		  if (filloccu) {ylayer_alloccusel->Fill(nstrip*iz+ypts[iz][iy]);}
		//		}
	      }
	      yhits[iz] = ypts[iz].size();
	    }

	    //Should not be only on MC code   if(ntrigX<4) continue;

	    //	    if (filloccu) { // (isfill) {}
	    if (isfill) {
	      int nrawoccu[2][nlayer]={0};
	      for (int ix=0; ix<nlayer; ix++) {
		for (int jkx=0; jkx<xlayer_occu[ix]->GetNbinsX(); jkx++) {
		  if (xlayer_occu[ix]->GetBinContent(jkx+1)>1.e-3) nrawoccu[0][ix]++;
		}
		for (int jkx=0; jkx<ylayer_occu[ix]->GetNbinsX(); jkx++) {
		  if (ylayer_occu[ix]->GetBinContent(jkx+1)>1.e-3) nrawoccu[1][ix]++;
		}
		nrawoccu[0][ix] -=3; //Maximum-3 hits is assumed to be large multiplicity layer
		nrawoccu[1][ix] -=3;
	      }

	      h_corr_layer_mult->Fill(-1.0, -1.0);
	      for (int ix=0; ix<nlayer; ix++) {
		if (xptsall[ix].size() >=nrawoccu[0][ix]) {
		  for (int iy=0; iy<nlayer; iy++) {
		    if (yptsall[iy].size() >=nrawoccu[1][ix]) {
		      h_corr_layer_mult->Fill(ix, iy+nlayer);
		    }
		  }

		  for (int ix1=ix+1; ix1<nlayer; ix1++) {
		    if (xptsall[ix1].size() >=nrawoccu[0][ix]) {
		      h_corr_layer_mult->Fill(ix, ix1);
		    }
		  }
		}
		if (yptsall[ix].size() >=nrawoccu[1][ix]) {
		  for (int ix1=ix+1; ix1<nlayer; ix1++) {
		    if (yptsall[ix1].size() >=nrawoccu[1][ix]) {
		      h_corr_layer_mult->Fill(ix+nlayer, ix1+nlayer);
		    }
		  }
		}
	      }
	    }



	    if (filloccu) {

              for (int ij=0; ij<nlayer; ij++) {
		xlayer_mult[ij]->Fill(xpts[ij].size());
		ylayer_mult[ij]->Fill(ypts[ij].size());

		rawhits_corr_xymul[ij]->Fill(xptsall[ij].size(), yptsall[ij].size());

		for (int ix=0; ix<xptsall[ij].size(); ix++) {
                  for (int iy=0; iy<yptsall[ij].size(); iy++) {
		    raw_occu[ij]->Fill(xptsall[ij][ix], yptsall[ij][iy]);
		    raw_occu_set[ij][iset]->Fill(xptsall[ij][ix], yptsall[ij][iy]);
		  }
		}
		for (int jk=ij+1; jk<nlayer; jk++) {
		  rawhits_xlay_corr_mul[ij][jk]->Fill(xptsall[ij].size(), xptsall[jk].size());
		  rawhits_ylay_corr_mul[ij][jk]->Fill(yptsall[ij].size(), yptsall[jk].size());
		}
	      }
	    }

#endif  //end of else of MONTECARLO

	    // bool TrgCheckL0  = (xptsall[1].size()>0 && xptsall[2].size()>0 && xptsall[9].size()>0 && xptsall[10].size()>0 )   ?  true : false;
	    ////   bool TrgCheckL11 =  (xpts[10].size()>0 && ypts[10].size()>0)   ?  true : false;

	    // //  if((!TrgCheckL0) || (!TrgCheckL11)) continue;
	    // if(!TrgCheckL0) continue;
	    //	    if((xptsall[0].size()==0) && (yptsall[0].size()==0)) continue;

            for (int ij=0; ij<nlayer; ij++) { Xdev[ij] = 100; Xpos[ij]=0.0;}
	    for (int ij=0; ij<nlayer; ij++) { Xdev1[ij] =Xdev2[ij] = 100; Xpos1[ij]= Xpos2[ij]=0.0;}

            for (int ij=0;ij<nlayer;ij++) {

	      if (xhits[ij]<=0 || xhits[ij]>nmxhits) {
		Xpos[ij]= -100;
	      } else {
		for (int ix=0; ix<xhits[ij]; ix++) {
		  if (ix<xhits[ij]-1 && abs(xpts[ij][ix]-xpts[ij][ix+1])>1) { Xpos[ij]=-100; break;}
		  Xpos[ij] +=xpts[ij][ix];
		}
		if (Xpos[ij]>=0.0) {
		  Xpos[ij]  = Xpos[ij]/xhits[ij] + 0.5 - xoff[ij];
#if defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
		  double tmptimexx[nmxtimehit-1] = {0.};
		  double tmptotxx[nmxtimehit-1] = {0.};
		  for(int ix=0;ix<min(nmxhits,xhits[ij]);ix++) {
		    if(event->vxtdc_l[ij][xpts[ij][ix]%8]->size()) {
		    tmptimexx[ix] = 1250.0 + 0.1*(int(event->vxtdc_l[ij][xpts[ij][ix]%8]->at(0)-event->tdc_ref_l[ij]));
		    tmptimexx[ix] -= xtoffset[ij][xpts[ij][ix]];
		    if(event->vxtdc_t[ij][xpts[ij][ix]%8]->size()){
		      tmptotxx[ix] = 0.1*int(event->vxtdc_t[ij][xpts[ij][ix]%8]->at(0))-0.1*int(event->vxtdc_l[ij][xpts[ij][ix]%8]->at(0));
		    }
		    }
		  }
		  /*
		    TOT and time correction for X-position
		  if(xhits[ij]>1 && xhits[ij]<=4) {
		    double lrange;
		     if(xhits[ij]==2)  { lrange = 5.;}
		    if(xhits[ij]==3)  { lrange = 2.5;}
		    else if(xhits[ij]==4)  { lrange = 1.5;}
		    if(abs(tmptimexx[xhits[ij]-1]-tmptimexx[0])<lrange) {
		      //  if(tmptimexx[xhits[ij]-1]>100 && tmptimexx[0]>100) {  Xpos[ij] += (xpos_xtime_corr[xhits[ij]-2][ij][0]+xpos_xtime_corr[xhits[ij]-2][ij][1]*(tmptimexx[xhits[ij]-1]-tmptimexx[0]));}
		    }
		    if((tmptotxx[xhits[ij]-1]/tmptotxx[0])>0.8 && (tmptotxx[xhits[ij]-1]/tmptotxx[0])<1.2 && tmptotxx[xhits[ij]-1]>3 && tmptotxx[0]>3) {Xpos[ij] += (xpos_xtot_corr[xhits[ij]-2][ij][0]+xpos_xtot_corr[xhits[ij]-2][ij][1]*(tmptotxx[xhits[ij]-1]/tmptotxx[0])); }
		  }
*/
		  Xpos[ij] -= cal_slope2(Xpos[ij], &align_xstr_xdev[ij][0]);
#endif
		  xxerr[ij] = xposerrsq[xhits[ij]-1][ij];
		}
	      }
	      if (isfill && Xpos[ij] >=0) xstrip_mult->Fill(6*ij+xhits[ij]);

	    //cout<<"Xpos/Ypos"<<"		"<<Xpos[ij]<<"		"<<Ypos[ij]<<endl;
	    }

	    for (int ij=0; ij<nlayer; ij++) { Ydev[ij] = 100; Ypos[ij]=0.0;}
	    for (int ij=0; ij<nlayer; ij++) { Ydev1[ij] = Ydev2[ij] = 100; Ypos1[ij]= Ypos2[ij]=0.0;}
            for (int ij=0;ij<nlayer;ij++) {
	      if (yhits[ij]<=0 || yhits[ij]>nmxhits) {
		Ypos[ij]=-100;
	      } else {
		for (int iy=0; iy<yhits[ij]; iy++) {
		  if (iy<yhits[ij]-1 && abs(ypts[ij][iy]-ypts[ij][iy+1])>1) { Ypos[ij]=-100; break;}
		  Ypos[ij] +=ypts[ij][iy];
		}
		if (Ypos[ij]>=0.0) {
		  Ypos[ij]  = Ypos[ij]/yhits[ij] + 0.5 - yoff[ij];
#if defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
		  double tmptimeyy[nmxtimehit-1] = {0.};
		  double tmptotyy[nmxtimehit-1] = {0.};
		  for(int ix=0;ix<min(nmxhits,yhits[ij]);ix++) {
		    if(event->vytdc_l[ij][ypts[ij][ix]%8]->size()){
		    tmptimeyy[ix] = 1250.0 + 0.1*(int(event->vytdc_l[ij][ypts[ij][ix]%8]->at(0)-event->tdc_ref_l[ij]));
		    tmptimeyy[ix] -= ytoffset[ij][ypts[ij][ix]];
		    if(event->vytdc_t[ij][ypts[ij][ix]%8]->size()){
		      tmptotyy[ix] = 0.1*int(event->vytdc_t[ij][ypts[ij][ix]%8]->at(0))-0.1*int(event->vytdc_l[ij][ypts[ij][ix]%8]->at(0));
		    }
		    }
		  }
		  /*
		    TOT and time correction for Y-position
		  if(yhits[ij]>1 && yhits[ij]<=4) {
		    double lrange;
		    if(yhits[ij]==2)  { lrange = 5.;}
		    else if(yhits[ij]==3)  { lrange = 2.5;}
		    else if(yhits[ij]==4)  { lrange = 1.5;}
		    if(abs(tmptimeyy[yhits[ij]-1]-tmptimeyy[0])<lrange) {
		      //	    if(tmptimeyy[yhits[ij]-1] > 100 && tmptimeyy[0]>100) {Ypos[ij] += (ypos_ytime_corr[yhits[ij]-2][ij][0]+ypos_ytime_corr[yhits[ij]-2][ij][1]*(tmptimeyy[yhits[ij]-1]-tmptimeyy[0]));}
		  }
		    if((tmptotyy[yhits[ij]-1]/tmptotyy[0])>0.8 && (tmptotyy[yhits[ij]-1]/tmptotyy[0])<1.2 && tmptotyy[yhits[ij]-1]>3 && tmptotyy[0]>3) {Ypos[ij] += (ypos_ytot_corr[yhits[ij]-2][ij][0]+ypos_ytot_corr[yhits[ij]-2][ij][1]*(tmptotyy[yhits[ij]-1]/tmptotyy[0])); }
		  }
*/
		  Ypos[ij] -= cal_slope2(Ypos[ij], &align_ystr_ydev[ij][0]);
#endif
		  yyerr[ij] = yposerrsq[yhits[ij]-1][ij];//+(20.*yposerrsq[yhits[ij]-1][ij])/100.);
		}
	      }
	      if (isfill && Ypos[ij] >=0) ystrip_mult->Fill(6*ij+yhits[ij]);
	    }

//#endif  //end of else of MONTECARLO
	    ntrigX=0; ntrigY=0;
	    for (int iz=0; iz<nlayer; iz++) {
	      if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && xptsall[iz].size()>0) {ntrigX++;} //18 01 2012
	      if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && yptsall[iz].size()>0) {ntrigY++;} //18 01 2012

	    }

	    if(isfill) {
	      // filecorOut->cd();
	      EvNum = iev;
	      XYhitcor = 0;
	      XYhitcor2 = 0;
	      int strp_cor[64]={0,1,   //0      //2
				21,45,
				1,2,3, //1     //7
				4,
				3,4,   //2     //10
				45,
				0,1,2,32, //3   //15
				3,4,5,
				0,37,  //4     //20
				0,28,
				0,1, 20,  //5   //25

				0,1,2, 58, //6   //29
				28,56,57, 62,
				0, 30,       //7  //35
				0,20,23,
				0,           //8   //39
				53,54,55,58,
				0,1,22,      //9  //46

				0,59,        //10   //48

				0,1,2,3,28,29,56,57,58,59, //11
				10,11,48,49,50,51};

	      int strp_cor2[48]={5,7,27,30,31,41,20,46,28,35,20,30,20,45,25,35,5,30,20,40,30,35,35,45,10,40,15,45,15,50,10,30,25,35,20,40,30,45,10,50,25,31,20,40,15,35,20,30};
	      ULong64_t tmphitcor=0;
	      ULong64_t tmphitcor2=0;
	      for(int ij=0;ij<nlayer;ij++) {
		XYhitall[ij] = 0;
		XYhitall[ij] += xptsall[ij].size();
		XYhitall[ij]<<=8;
		XYhitall[ij] += yptsall[ij].size();
		int istx,isty,iendx,iendy;
		if(ij==0) { istx = 0;iendx=1; isty=2;iendy=3;}
		else if(ij==1) { istx = 4;iendx=6; isty=7;iendy=7;}
		else if(ij==2) { istx = 8;iendx=9; isty=10;iendy=10;}
		else if(ij==3) { istx = 11;iendx=14; isty=15;iendy=17;}
		else if(ij==4) { istx = 18;iendx=19; isty=20;iendy=21;}
		else if(ij==5) { istx = 22;iendx=24; isty=0;iendy=-1;}
		else if(ij==6) { istx = 25;iendx=28; isty=29;iendy=32;}
		else if(ij==7) { istx = 33;iendx=34; isty=35;iendy=37;}
		else if(ij==8) { istx = 38;iendx=38; isty=39;iendy=42;}
		else if(ij==9) { istx = 43;iendx=45; isty=0;iendy=-1;}
		else if(ij==10) { istx = 46;iendx=47; isty=0;iendy=-1;}
	       else if(ij==11) { istx = 48;iendx=57; isty=58;iendy=63;}

		ULong64_t tmpxy=0;
		for(int strx=istx;strx<=iendx;strx++) {
		  ULong64_t tmpx =0;
		  for(int ii=0;ii<xptsall[ij].size();ii++) {
		    if(strp_cor[strx]==xptsall[ij][ii]){
		      tmpx = 1;
		      tmpx<<=strx;
		    }
		  }
		  tmpxy+=tmpx;
		}


		for(int stry=isty;stry<=iendy;stry++) {
		  ULong64_t tmpy =0;
		  for(int jj=0;jj<yptsall[ij].size();jj++) {
		    if(strp_cor[stry]==yptsall[ij][jj]){
		      tmpy  = 1;
		      tmpy<<=stry;
		    }
		  }
		  tmpxy+=tmpy;
		}
		tmphitcor += tmpxy;

		ULong64_t tmpxy2=0;
		for(int strx=ij*4;strx<ij*4+2;strx++) {
		  ULong64_t tmpx2 =0;
		  for(int ii=0;ii<xptsall[ij].size();ii++) {
		    if(strp_cor2[strx]==xptsall[ij][ii]){
		      tmpx2 = 1;
		      tmpx2<<=strx;
		    }
		  }
		  tmpxy2+=tmpx2;
		}


		for(int stry=ij*4+2;stry<ij*4+4;stry++) {
		  ULong64_t tmpy2 =0;
		  for(int jj=0;jj<yptsall[ij].size();jj++) {
		    if(strp_cor2[stry]==yptsall[ij][jj]){
		      tmpy2  = 1;
		      tmpy2<<=stry;
		    }
		  }
		  tmpxy2+=tmpy2;
		}

		tmphitcor2 += tmpxy2;


		if (xptsall[ij].size() >0) {
		  for (int ix=0; ix<xptsall[ij].size(); ix++) {
		    h_raw_xcorstrips[ij]->Fill(xptsall[ij][ix], -1.);
		    h_raw_xstrpnhits[ij]->Fill(xptsall[ij][ix], -1.);
		    h_raw_xystrpnhits[ij]->Fill(xptsall[ij][ix], -1.);
		    h_raw_xstrpnhits[ij]->Fill(xptsall[ij][ix],xptsall[ij].size());
		    h_raw_xystrpnhits[ij]->Fill(xptsall[ij][ix],yptsall[ij].size());
		    for(int ixx=0;ixx<xptsall[ij].size();ixx++) {

		      if(abs(xptsall[ij][ix]-xptsall[ij][ixx])>0) {
			h_raw_xcorstrips[ij]->Fill(xptsall[ij][ix],xptsall[ij][ixx]);
		      }
		    }
		  }
		}

		if (yptsall[ij].size() >0) {
		  for (int iy=0; iy<yptsall[ij].size(); iy++) {
		    h_raw_ystrpnhits[ij]->Fill(yptsall[ij][iy], -1.);
		    h_raw_ycorstrips[ij]->Fill(yptsall[ij][iy], -1.);
		    h_raw_yxstrpnhits[ij]->Fill(yptsall[ij][iy], -1.);
		    h_raw_ystrpnhits[ij]->Fill(yptsall[ij][iy],yptsall[ij].size());
		    h_raw_yxstrpnhits[ij]->Fill(yptsall[ij][iy],xptsall[ij].size());
		    for(int iyy=0;iyy<yptsall[ij].size();iyy++) { //choose all of them and then use maximum value before plot

		      if(abs(yptsall[ij][iy]-yptsall[ij][iyy])>0) {
			h_raw_ycorstrips[ij]->Fill(yptsall[ij][iy],yptsall[ij][iyy]);
		    }
		    }
		  }
		}
	      }
	      XYhitcor=tmphitcor;
	      XYhitcor2=tmphitcor2;

	      //	      T3->Fill();
	    }

	    //  fileOut->cd();
	    // fileIn->cd();

	    if(ntrigX<4 && ntrigY<4) {/*cout<<iev<<"  without trigger hits"<<endl;*/ nxytrig++; }

	    if(xptsall[trigly1].size()==0) { nx4++;}
	    if(xptsall[trigly2].size()==0) { nx5++;}
	    if(xptsall[trigly3].size()==0) { nx6++;}
	    if(xptsall[trigly4].size()==0) { nx7++;}
	    if(yptsall[trigly1].size()==0) { ny4++;}
	    if(yptsall[trigly2].size()==0) { ny5++;}
	    if(yptsall[trigly3].size()==0) { ny6++;}
	    if(yptsall[trigly4].size()==0) { ny7++;}
	    if(xptsall[trigly1].size()==0 && yptsall[trigly1].size()==0) { nxy4++;}
	    if(xptsall[trigly2].size()==0 && yptsall[trigly2].size()==0) { nxy5++;}
	    if(xptsall[trigly3].size()==0 && yptsall[trigly3].size()==0) { nxy6++;}
	    if(xptsall[trigly4].size()==0 && yptsall[trigly4].size()==0) { nxy7++;}
	    if(xptsall[trigly1].size()==0 || yptsall[trigly1].size()==0) { nxory4++;}
	    if(xptsall[trigly2].size()==0 || yptsall[trigly2].size()==0) { nxory5++;}
	    if(xptsall[trigly3].size()==0 || yptsall[trigly3].size()==0) { nxory6++;}
	    if(xptsall[trigly4].size()==0 || yptsall[trigly4].size()==0) { nxory7++;}

#ifdef NEWRPCDAQ1
	    if(filloccu){



	      for(int ij=0;ij<nlayer;ij++) {
		bool tdccountx =false; bool tdccounty=false;
		for(int jk=0;jk<8;jk++) {
		  if(event->vxtdc_l[ij][jk]->size()>0) { tdccountx=true; break;}
		  // if(event->vytdc_l[ij][jk]->size()>0) { tdccounty=true; break;}
		}
		for(int jk=0;jk<8;jk++) {
		  //if(event->vxtdc_l[ij][jk]->size()>0) { tdccountx=true; break;}
		   if(event->vytdc_l[ij][jk]->size()>0) { tdccounty=true; break;}
		}
		if(xptsall[ij].size()) {
		  ntxext[ij][0]++;
		  if(event->tdc_ref_l[ij]>0) {
		    ntxref[ij][0]++;
		    if(tdccountx) {
		      nttxext[ij][0]++;
		    }
		  }
		}

		if(yptsall[ij].size()) {
		  ntyext[ij][0]++;
		  if(event->tdc_ref_l[ij]>0) {
		    ntyref[ij][0]++;
		    if(tdccounty) {
		      nttyext[ij][0]++;
		    }
		  }
		}


		if(xptsall[ij].size()) {
		  ntxext[ij][1]++;
		  if(tdccountx) {
		    nttxext[ij][1]++;
		    if(event->tdc_ref_l[ij]>0) {
		      ntxref[ij][1]++;
		    }
		  }
		}

		if(yptsall[ij].size()) {
		  ntyext[ij][1]++;
		  if(tdccounty) {
		    nttyext[ij][1]++;
		    if(event->tdc_ref_l[ij]>0) {
		      ntyref[ij][1]++;
		    }
		  }
		}

		if(event->tdc_ref_l[ij]>0) {
		  ntxref[ij][2]++;
		  if(xptsall[ij].size()) {
		    ntxext[ij][2]++;
		    if(tdccountx) {
		      nttxext[ij][2]++;
		    }
		  }
		}

	        if(event->tdc_ref_l[ij]>0) {
		  ntyref[ij][2]++;
		  if(yptsall[ij].size()) {
		    ntyext[ij][2]++;
		    if(tdccounty) {
		      nttyext[ij][2]++;
		    }
		  }
		}

		if(event->tdc_ref_l[ij]>0) {
		  ntxref[ij][3]++;
		  if(tdccountx) {
		    nttxext[ij][3]++;
		    if(xptsall[ij].size()) {
		      ntxext[ij][3]++;
		    }
		  }
		}


		if(event->tdc_ref_l[ij]>0) {
		  ntyref[ij][3]++;
		  if(tdccounty) {
		    nttyext[ij][3]++;
		    if(yptsall[ij].size()) {
		      ntyext[ij][3]++;
		    }
		  }
		}


		if(tdccountx) {
		  nttxext[ij][4]++;
		  if(xptsall[ij].size()) {
		    ntxext[ij][4]++;
		    if(event->tdc_ref_l[ij]>0) {
		      ntxref[ij][4]++;
		    }
		  }
		}

		if(tdccounty) {
		  nttyext[ij][4]++;
		  if(yptsall[ij].size()) {
		    ntyext[ij][4]++;
		    if(event->tdc_ref_l[ij]>0) {
		      ntyref[ij][4]++;
		    }
		  }
		}
		if(tdccountx) {
		  nttxext[ij][5]++;
		  if(event->tdc_ref_l[ij]>0) {
		    ntxref[ij][5]++;
		    if(xptsall[ij].size()) {
		      ntxext[ij][5]++;
		    }
		  }
		}

		if(tdccounty) {
		  nttyext[ij][5]++;
		  if(event->tdc_ref_l[ij]>0) {
		    ntyref[ij][5]++;
		    if(yptsall[ij].size()) {
		      ntyext[ij][5]++;
		    }
		  }
		}
	      }
	      int ntdcref =0;
	      for(int ij=0;ij<nlayer;ij++) {
		if(event->tdc_ref_l[ij]<=0) {
		  ntdcref++;
		  h_corr_tdc_ref->Fill(ij,-1.);
		  for(int jk=0;jk<nlayer;jk++) {
		    if(abs(ij-jk)>0 && event->tdc_ref_l[jk]<=0) {
		      h_corr_tdc_ref->Fill(ij,jk);
		    }
		  }
		}
	      }

	      h_tdc_ref->Fill(ntdcref);

	    }
#endif



	    //Sort out hits, which can be used for fit
	    for (int ij=0;ij<nlayer;ij++) {
	      Xusedpos[ij] = (Xpos[ij]>-100 && xhits[ij]<=nmxusedhits) ? true : false; //Xpos[ij] : -101;
	    }


	    int tmpnx2=0;
	    for (int ix=0; ix<nlayer; ix++) { if (Xpos[ix]>-90) {tmpnx2++;} }
	    h_tmp2xndf->Fill(tmpnx2);

            Nx=0;
	    int nxfail = 0;
	    xchi2 = 0;
	    double xresol = 0;
	    double zval[nlayer], xext[nlayer], xextloc[nlayer], xexter[nlayer], xposinstr[nlayer];
	    double yext[nlayer], yextloc[nlayer], yexter[nlayer], yposinstr[nlayer];;
	    for (int ix=0; ix<nlayer; ix++) { yext[ix]= yextloc[ix] = yexter[ix] =yposinstr[ix] =  100;}
	    for (int ix=0; ix<nlayer; ix++) { zval[ix]=layerzpos[ix];}
	    for (int ix=0; ix<nlayer; ix++) { xext[ix]= xextloc[ix] = xexter[ix] =xposinstr[ix] =  100;}
	    // for (int ix=0; ix<nlayer; ix++) {


	      // }
	    StraightLineFit xposfit(1, zval, Xpos,  xxerr, Xusedpos,occulyr, occulyr, layfirst, laylast, xyPosDev);
	    xposfit.GetParameters(nxfail, xinters, xslope);
	    xposfit.GetFitValues(xext, Xdev, xexter);
#ifndef NEWRPCDAQ1
	    xposfit.GetChisqure(Nx,xchi2);
	    h_tmp3xndf->Fill(Nx);
	    GetXposInStrip(0, xext, yext, xoff, xposinstr, xextloc);
#endif
	    //Sort out hits, which can be used for fit
	    for (int ij=0;ij<nlayer;ij++) {
	      Yusedpos[ij] = (Ypos[ij]>-99 && yhits[ij]<=nmxusedhits) ? true : false; //Ypos[ij] :  -101;
	    }
            Ny=0;
	    int nyfail = 0;
	    ychi2 = 0;
	    double yresol = 0;
	    //  double yext[nlayer], yextloc[nlayer], yexter[nlayer], yposinstr[nlayer];;
	    // for (int ix=0; ix<nlayer; ix++) { yext[ix]= yextloc[ix] = yexter[ix] =yposinstr[ix] =  100;}

	    StraightLineFit yposfit(1, zval, Ypos,  yyerr, Yusedpos,occulyr, occulyr, layfirst, laylast, xyPosDev);
	    yposfit.GetParameters(nyfail, yinters, yslope);
	    yposfit.GetFitValues(yext, Ydev, yexter);
#ifndef NEWRPCDAQ1
	    //	    yposfit.GetError(errcst, errlin, errcov);
	    yposfit.GetChisqure(Ny, ychi2);
	    GetXposInStrip(1, yext, xext, yoff, yposinstr, yextloc);
#endif

#if defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
	    for (int ix=0; ix<nlayer; ix++) {
	      //	      cout <<"ix "<< ix<<" "<<Xpos[ix]<<" "<<Ypos[ix]<<" "<<align_ystr_xdev[ix][0]<<" "<<align_ystr_xdev[ix][1]<<" "<<yext[ix]<<endl;
	      Xpos[ix] -=cal_slope2(yext[ix], &align_ystr_xdev[ix][0]);
	      Ypos[ix] -=cal_slope2(xext[ix], &align_xstr_ydev[ix][0]);
	      //	      cout <<"ix "<< ix<<" "<<Xpos[ix]<<" "<<Ypos[ix]<<" "<<align_xstr_ydev[ix][0]<<" "<<align_xstr_ydev[ix][1]<<" "<<xext[ix]<<endl;
	    }

	    xposfit = StraightLineFit(1, zval, Xpos,  xxerr, Xusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
	    xposfit.GetParameters(nxfail, xinters, xslope);
	    xposfit.GetFitValues(xext, Xdev, xexter);
	    xposfit.GetChisqure(Nx,xchi2);
	    h_tmp3xndf->Fill(Nx);
	    GetXposInStrip(0, xext, yext, xoff, xposinstr, xextloc);

	    yposfit =StraightLineFit(1, zval, Ypos,  yyerr, Yusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
	    yposfit.GetParameters(nyfail, yinters, yslope);
	    yposfit.GetFitValues(yext, Ydev, yexter);
	    yposfit.GetChisqure(Ny, ychi2);
	    GetXposInStrip(1, yext, xext, yoff, yposinstr, yextloc);
#endif
	    xbina = 0;
	    isfid_l9_eff = false; //jim jim
	    if(ntcor==1 && occulyr==obslay) {
	      xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[obslay]),int(yextloc[obslay])));
	      // isfid_l9_eff = ((xextloc[occulyr]>0 && xextloc[occulyr]<lastXstrip && yextloc[occulyr]>0 && yextloc[occulyr]<lastYstrip) && (int(xextloc[occulyr]-4)%6 <5 && int(xextloc[occulyr]-4)%6 >1 && int(yextloc[occulyr]-3)%7 <6 && int(yextloc[occulyr]-3)%7 >1)) ? true:false;
	      bool cond1 = (int(xextloc[occulyr]) == 4 && int(yextloc[occulyr])==3) ? true:false;
	      bool cond2 = (int(xextloc[occulyr])>4 && int(yextloc[occulyr])>3 && int(xextloc[occulyr]-4)%6 < 5 && int(xextloc[occulyr]-4)%6 > 1 && int(yextloc[occulyr]-3)%7 < 6 && int(yextloc[occulyr]-3)%7 > 1) ? true:false;
	      isfid_l9_eff = ((xextloc[occulyr]>0 && xextloc[occulyr]<lastXstrip && yextloc[occulyr]>0 && yextloc[occulyr]<lastYstrip)/* && (cond1||cond2)*/) ? true:false;
	    } else {
		xbina = 0;
		isfid_l9_eff = false;
	      }
	      if(xbina>nsplit-1) xbina = nsplit-1;

	      nhits = 100*Nx + Ny;
	      // for(int ij=0;ij<nlayer;ij++) {
	      if(ntcor==1 && iiter==nmxiter-1) {
	      tmpxmult = xpts[occulyr].size();
	      tmpymult = ypts[occulyr].size();
	      for(int jk=0;jk<tmpxmult;jk++) { xtpulsewidth_mul[jk] =  rawtmpxtime[jk] =  tmpxpts[jk] = 10000;}
	      for(int jk=0;jk<tmpymult;jk++) { ytpulsewidth_mul[jk] =  rawtmpytime[jk] =  tmpypts[jk] = 10000;}
	      tmpxext = tmpyext = tmpxdev = tmpydev =  tmpxtdev = tmpytdev = xtpulsewidth = ytpulsewidth = tmpxtexter =  tmpytexter = tmpxexter = tmpyexter = 1000;
	      tmpxtime = tmpytime = tmpxtext = tmpyext = 10000;
	      // tmpxmult = tmpymult = -10;
	      layoccu = occulyr;
	      tmpxext =  xextloc[occulyr];
	      tmpyext =  yextloc[occulyr];
	      tmpxdev = Xdev[occulyr];
	      tmpydev = Ydev[occulyr];

	      tmpxexter = xexter[occulyr];
	      tmpyexter = yexter[occulyr];
	      tmpNx = Nx;
	      tmpNy = Ny;
	      tmpxchis2 = xchi2;
	      tmpychis2 = ychi2;
	      for(int jk=0;jk<tmpxmult;jk++) {
		//cout<<iev<<"   "<<jk<<"   "<<"    "<<endl;
		if(event->vxtdc_l[occulyr][xpts[occulyr][jk]%8]->size()) {
		  if(event->vxtdc_t[occulyr][xpts[occulyr][jk]%8]->size()) {

		    tmpxpts[jk] = xpts[occulyr][jk];
		    xtpulsewidth_mul[jk] = 0.1*int(event->vxtdc_t[occulyr][xpts[occulyr][jk]%8]->at(0)) -  0.1*int(event->vxtdc_l[occulyr][xpts[occulyr][jk]%8]->at(0));
		    rawtmpxtime[jk] = 1250. + 0.1*(int(event->vxtdc_l[occulyr][xpts[occulyr][jk]%8]->at(0))-int(event->tdc_ref_l[occulyr]));
		  }
		}
	      }
	      for(int jk=0;jk<tmpymult;jk++) {
		if(event->vytdc_l[occulyr][ypts[occulyr][jk]%8]->size()) {
		  if(event->vytdc_t[occulyr][ypts[occulyr][jk]%8]->size()) {
		    tmpypts[jk] = ypts[occulyr][jk];
		    ytpulsewidth_mul[jk] = 0.1*int(event->vytdc_t[occulyr][ypts[occulyr][jk]%8]->at(0)) -  0.1*int(event->vytdc_l[occulyr][ypts[occulyr][jk]%8]->at(0));
		    rawtmpytime[jk] = 1250. + 0.1*(int(event->vytdc_l[occulyr][ypts[occulyr][jk]%8]->at(0))-int(event->tdc_ref_l[occulyr]));
		  }
		}
	      }

	    }

	    for(int ij = 0;ij<nlayer;ij++) { Xpos1[ij] = Xpos2[ij]=Xpos[ij]; Ypos1[ij] = Ypos2[ij] = Ypos[ij];}
	    for (int ij=0;ij<nlayer;ij++) {
	      if(ij<=5) { Xusedpos1[ij] = (Xpos1[ij]>-100 && xhits[ij]<=nmxusedhits) ? true : false; } else { Xusedpos1[ij]=false;}
	    }

            Nx1=0;
	    int nxfail1 = 0;
	    xchi21 = 0;

	    double  xext1[nlayer], xextloc1[nlayer],xexter1[nlayer], xposinstr1[nlayer];

	    for (int ix=0; ix<nlayer; ix++) { xext1[ix]= xextloc1[ix] = xexter1[ix] = 100;}
	    StraightLineFit xposfit1(1, zval, Xpos1,  xxerr, Xusedpos1, occulyr, occulyr,0, 5, xyPosDev);
	    xposfit1.GetParameters(nxfail1, xinters1, xslope1);
	    //	    xposfit.GetError(errcst, errlin, errcov);
	    xposfit1.GetChisqure(Nx1, xchi21);

	    xposfit1.GetFitValues(xext1, Xdev1, xexter1);

	    for (int ij=0;ij<nlayer;ij++) {
	      if(ij<=5) { Yusedpos1[ij] = (Ypos1[ij]>-99 && yhits[ij]<=nmxusedhits) ? true : false; } else { Yusedpos1[ij]=false;}
	    }

            Ny1=0;
	    int nyfail1 = 0;
	    ychi21 = 0;

	    double yext1[nlayer], yextloc1[nlayer],yexter1[nlayer], yposinstr1[nlayer];;
            for (int ix=0; ix<nlayer; ix++) { yext1[ix]= yextloc1[ix] = yexter1[ix] =  100;}

	    StraightLineFit yposfit1(1, zval, Ypos1,  yyerr, Yusedpos1, occulyr, occulyr, 0, 5, xyPosDev);
	    yposfit1.GetParameters(nyfail1, yinters1, yslope1);
	    //	    yposfit.GetError(errcst, errlin, errcov);
	    yposfit1.GetChisqure(Ny1, ychi21);
	    yposfit1.GetFitValues(yext1, Ydev1, yexter1);

	    for (int ij=0;ij<nlayer;ij++) {
	      if(ij>5) { Xusedpos2[ij] = (Xpos2[ij]>-100 && xhits[ij]<=nmxusedhits) ? true : false; } else { Xusedpos2[ij]=false;}
	    }

	    Nx2=0;
	    int nxfail2 = 0;
	    xchi22 = 0;

	    double  xext2[nlayer], xextloc2[nlayer], xexter2[nlayer], xposinstr2[nlayer];

	    for (int ix=0; ix<nlayer; ix++) { xext2[ix]= xextloc2[ix] = xexter2[ix] = 100;}
	    StraightLineFit xposfit2(1, zval, Xpos2,  xxerr, Xusedpos2, occulyr, occulyr, 6, 11, xyPosDev);
	    xposfit2.GetParameters(nxfail2, xinters2, xslope2);
	    //	    xposfit.GetError(errcst, errlin, errcov);
	    xposfit2.GetChisqure(Nx2, xchi22);

	    xposfit2.GetFitValues(xext2, Xdev2, xexter2);

	    for (int ij=0;ij<nlayer;ij++) {
	      if(ij>5) { Yusedpos2[ij] = (Ypos2[ij]>-99 && yhits[ij]<=nmxusedhits) ? true : false; } else { Yusedpos2[ij]=false;}
	    }

            Ny2=0;
	    int nyfail2 = 0;
	    ychi22 = 0;

	    double yext2[nlayer],  yextloc2[nlayer], yexter2[nlayer], yposinstr2[nlayer];;
            for (int ix=0; ix<nlayer; ix++) { yext2[ix]=  yextloc2[ix] = yexter2[ix] =  100;}

	    StraightLineFit yposfit2(1, zval, Ypos2,  yyerr, Yusedpos2, occulyr, occulyr, 6, 11, xyPosDev);
	    yposfit2.GetParameters(nyfail2, yinters2, yslope2);
	    //	    yposfit.GetError(errcst, errlin, errcov);
	    yposfit2.GetChisqure(Ny2, ychi22);
	    yposfit2.GetFitValues(yext2, Ydev2, yexter2);

	    if (isfill) { //Correlated hits
		for(int ij=0;ij<nlayer;ij++) {
		  double expp = xextloc[ij];//-100;
		  double expdiff = 10000;
		  // for (int ix=0; ix<xpts[ij].size(); ix++) {
		  //   if (abs(xextloc[ij] - xpts[ij][ix])<1.0) {
		  //     double tmpexpdiff = abs(xextloc[ij] - xpts[ij][ix]);
		  //     if(tmpexpdiff < expdiff) {
		  //       expdiff =tmpexpdiff;	expp = xpts[ij][ix];// break;
		  //     }
		  //   }
		  // }
		  if (expp >-50) {
		    h_xmucorstrips[ij]->Fill(expp, -1.);
		    h_xymucorstrips[ij]->Fill(expp, -1.);
		    h_xmucornhits[ij]->Fill(expp,-1.);
		    h_xmucornhits[ij]->Fill(expp,xpts[ij].size());
		    h_xymucornhits[ij]->Fill(expp,-1.);
		    h_xymucornhits[ij]->Fill(expp,ypts[ij].size());
		    for (int ix=0; ix<xpts[ij].size(); ix++) {
		      if (abs(expp - xpts[ij][ix])>0) { //choose all of them and then use maximum value before plot
			h_xmucorstrips[ij]->Fill(expp, xpts[ij][ix]);
			//	xstrip_xdev[occulyr][iiterrs]->Fill(xepp[occulyr],Xdev[occulyr]);
		      }
		    }

		    for (int ix=0; ix<ypts[ij].size(); ix++) {
		      h_xymucorstrips[ij]->Fill(expp, ypts[ij][ix]);
		    }
		  }
		}

		for(int ij=0;ij<nlayer;ij++) {
		  double expp = yextloc[ij];//-100;
		  double expdiff = 10000;
		  // for (int ix=0; ix<ypts[ij].size(); ix++) {
		  //   if (abs(yextloc[ij] - ypts[ij][ix])<1.0) {
		  //     double tmpexpdiff = abs(yextloc[ij] - ypts[ij][ix]);
		  //     if(tmpexpdiff < expdiff) {
		  //       expdiff =tmpexpdiff;	expp = ypts[ij][ix];// break;
		  //     }
		  //   }
		  // }
		  if (expp >-50) {
		    h_ymucorstrips[ij]->Fill(expp, -1.);
		    h_yxmucorstrips[ij]->Fill(expp, -1.);
		    h_ymucornhits[ij]->Fill(expp,-1.);
		    h_ymucornhits[ij]->Fill(expp,ypts[ij].size());
		    h_yxmucornhits[ij]->Fill(expp,-1.);
		    h_yxmucornhits[ij]->Fill(expp,ypts[ij].size());

		    for (int ix=0; ix<ypts[ij].size(); ix++) {
		      if (abs(expp - ypts[ij][ix])>0) { //choose all of them and then use maximum value before plot
			h_ymucorstrips[ij]->Fill(expp, ypts[ij][ix]);
			//	xstrip_xdev[occulyr][iiterrs]->Fill(xepp[occulyr],Xdev[occulyr]);
		      }
		    }

		    for (int ix=0; ix<xpts[ij].size(); ix++) {
		      h_yxmucorstrips[ij]->Fill(expp, xpts[ij][ix]);
		    }
		  }
		}





	    } //if (isfill)


            if (Nx>= nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {//4Nov,2011
	      if (filloccu) {
		if (iev>0 && tmptimerate >=tmpoldmuposxrate+timebin) {
		  file_out <<"posxrate "<<datafile<<" "<<iev<<" "<<nallmuposxrate<<" "<<nmuposxrate<<endl;
		  muposxrate->Fill(nmuposxrate*(1.0/timebin));
		  muposxratex->Fill(nallmuposxrate, nmuposxrate*(1.0/timebin));
		  nmuposxrate = 0;
		  nallmuposxrate++;
		  tmpoldmuposxrate = tmptimerate;
		}
		nmuposxrate++;
	      }




	      pos_xslope[occulyr][iiter]->Fill(stripwidth*xslope);
	      if (occulyr>=nlayer) {
                for (int ij=0; ij<nlayer; ij++) {
		  if (abs(Xdev[ij])<1.5 && Xusedpos[ij]) {
		    for (int jk=0; jk<nlayer; jk++) {
		      if (abs(Xdev[jk])<1.5 && Xusedpos[ij]) {
			deltaposxcov[ij][jk] += Xdev[ij]*Xdev[jk];
			deltaposxCount[ij][jk] +=1 ;
		      }
		    }
		  }
#ifndef MONTECARLO
		  if (filloccu && Xusedpos[ij]) {xlayer_reso_set[ij][iset]->Fill(Xdev[ij]);}
#endif
		  // xstrip_xdev[ij][iiterrs]->Fill(xextloc[ij],Xdev[ij]);
		  if (abs(Xdev[ij])<6.0 && Xusedpos[ij] && xextloc[ij]>firstXstrip && xextloc[ij]<lastXstrip && (posinuse[0][ij][int(xextloc[ij]+0.5)]==true)) {
		    // cout<<"xxxxxxxxxxxxx"<<endl;
		    xlayer_reso[ij][iiterrs]->Fill(Xdev[ij]);
		    xlayer_reso_nTermSplit[ij][int(xextloc[ij]/20)][iiterrs]->Fill(Xdev[ij]);
		    //		     xstrip_xdev[ij]->Fill(xextloc[ij],Xdev[ij]);
		    // int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[ij]+0.5),int(yextloc[ij]+0.5)));//(nsplit/8)*int(xextloc[ij]/7.25)+int(yextloc[ij]/7.625);
		    //if(xbina>nsplit-1) xbina = nsplit-1;
		    if(isfid_l9_eff) { xlayer_reso_split[ij][xbina][iiterrs]->Fill(Xdev[ij]); }
		    is_xyreso[ij][iiterrs] = true;
		    xlayer_exterr[ij][iiterrs]->Fill(xexter[ij]);
		    xlayer_exterr_nTermSplit[ij][int(xextloc[ij]/20)][iiterrs]->Fill(xexter[ij]);
		    int ixx=min(xhits[ij], nmxhits);
		    if (ixx>0 && Xpos[ij]>1.0 && Xpos[ij]<nstrip-1) {
		      xlayer_reso_mul[ij][iiterrs][ixx-1]->Fill(Xdev[ij]);
		      if(isfid_l9_eff) { xlayer_reso_mul_split[ij][xbina][iiterrs][ixx-1]->Fill(Xdev[ij]);}
		    }
		  }
		}
		nTotallx++;
		for (int ij=0; ij<nlayer; ij++) {
		  if (abs(Xdev[ij])<1.5 && Xusedpos[ij]) {
		    for (int jk=0; jk<nlayer; jk++) {
		      if (abs(Xdev[jk])<1.5 && Xusedpos[jk]) {
			h_posxcoreff->Fill(ij,jk) ;
			posxcoreffCount +=1 ;
		      }
		    }
		  }
		}
	      } else {
		//	xstrip_xdev[occulyr][iiterrs]->Fill(xextloc[occulyr],Xdev[occulyr]);
		//	//int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[occulyr]+0.5),int(yextloc[occulyr]+0.5)));//(nsplit/8)*int(xextloc[occulyr]/7.25)+int(yextloc[occulyr]/7.625);
		//if(xbina>nsplit-1) xbina = nsplit-1;
		if (abs(Xdev[occulyr])<6.0 && Xusedpos[occulyr] && xextloc[occulyr]>firstXstrip && xextloc[occulyr]<lastXstrip && (posinuse[0][occulyr][int(xextloc[occulyr]+0.5)]==true)) {
		  //		  cout<<"yyyyyyyyyyyyyyyy"<<endl;
		  is_xyreso[occulyr][iiterrs] = true;
		  xlayer_reso[occulyr][iiterrs]->Fill(Xdev[occulyr]);
		  xlayer_reso_nTermSplit[occulyr][int(xextloc[occulyr]/20)][iiterrs]->Fill(Xdev[occulyr]);
		  if(isfid_l9_eff) { xlayer_reso_split[occulyr][xbina][iiterrs]->Fill(Xdev[occulyr]);}
		  xlayer_exterr[occulyr][iiterrs]->Fill(xexter[occulyr]);
		  xlayer_exterr_nTermSplit[occulyr][int(xextloc[occulyr]/20)][iiterrs]->Fill(xexter[occulyr]);

		  int ixx=min(xhits[occulyr], nmxhits);
		  if (ixx>0 && Xpos[occulyr]>1&&Xpos[occulyr]<nstrip-1 ) {
		    xlayer_reso_mul[occulyr][iiterrs][ixx-1]->Fill(Xdev[occulyr]);
		    if(isfid_l9_eff) {xlayer_reso_mul_split[occulyr][xbina][iiterrs][ixx-1]->Fill(Xdev[occulyr]);}
		  }
		}




		if (iiter==nmxiter-1) { //Correlated hits
		  double expp = -100;
		  double expdiff = 10000;
		  for (int ix=0; ix<xpts[occulyr].size(); ix++) {
		    if (abs(xextloc[occulyr] - xpts[occulyr][ix])<1.0) {
		      double tmpexpdiff = abs(xextloc[occulyr] - xpts[occulyr][ix]);
		      if(tmpexpdiff < expdiff) {
		        expdiff =tmpexpdiff;	expp = xpts[occulyr][ix];// break;
		      }
		    }
		  }




		  if (expp >-50) {
		    h_xcorstrips[occulyr]->Fill(expp, -1.);
		    h_xycorstrips[occulyr]->Fill(expp, -1.);
		    for (int ix=0; ix<xpts[occulyr].size(); ix++) {
		      if (abs(expp - xpts[occulyr][ix])>0) { //choose all of them and then use maximum value before plot
			h_xcorstrips[occulyr]->Fill(expp, xpts[occulyr][ix]);
			//	xstrip_xdev[occulyr][iiterrs]->Fill(xepp[occulyr],Xdev[occulyr]);
		      }
		    }

		    for (int ix=0; ix<ypts[occulyr].size(); ix++) {
		      h_xycorstrips[occulyr]->Fill(expp, ypts[occulyr][ix]);
		    }
		  }
		} //if (iiter==nmxiter-1)
	      }
	    }
            if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {//4Nov
	      if (filloccu) {
		if (iev>0 && tmptimerate >=tmpoldmuposyrate+timebin) {
		  file_out <<"posyrate "<<datafile<<" "<<iev<<" "<<nallmuposyrate<<" "<<nmuposyrate<<endl;
		  muposyrate->Fill(nmuposyrate*(1.0/timebin));
		  muposyratex->Fill(nallmuposyrate, nmuposyrate*(1.0/timebin));
		  nmuposyrate = 0;
		  nallmuposyrate++;
		  tmpoldmuposyrate = tmptimerate;
		}
		nmuposyrate++;
	      }

	      pos_yslope[occulyr][iiter]->Fill(stripwidth*yslope);
              if (occulyr>=nlayer) {
                for (int ij=0; ij<nlayer; ij++) {
		  if (abs(Ydev[ij])<1.5 && Yusedpos[ij]) {
		    for (int jk=0; jk<nlayer; jk++) {
		      if (abs(Ydev[jk])<1.5 && Yusedpos[ij]) {
			deltaposycov[ij][jk] += Ydev[ij]*Ydev[jk];
			deltaposyCount[ij][jk] +=1 ;
		      }
		    }
		  }
#ifndef MONTECARLO
		  if (filloccu && Yusedpos[ij]) {ylayer_reso_set[ij][iset]->Fill(Ydev[ij]);}
#endif
		  //  //int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[ij]+0.5),int(yextloc[ij]+0.5)));//(nsplit/8)*int(xextloc[ij]/7.25)+int(yextloc[ij]/7.625);
		  //   ystrip_ydev[ij][iiterrs]->Fill(yextloc[ij],Ydev[ij]);
		  //if(xbina>nsplit-1) xbina = nsplit-1;
                  if (abs(Ydev[ij])<6.0  && Yusedpos[ij] && yextloc[ij]>firstYstrip && yextloc[ij]<lastYstrip && posinuse[1][ij][int(yextloc[ij]+0.5)]==true) {
                    ylayer_reso[ij][iiterrs]->Fill(Ydev[ij]);
		    ylayer_reso_nTermSplit[ij][int(yextloc[ij]/20)][iiterrs]->Fill(Ydev[ij]);
		    if(isfid_l9_eff) {ylayer_reso_split[ij][xbina][iiterrs]->Fill(Ydev[ij]);}
		    if (is_xyreso[ij][iiterrs]) {xylayer_reso[ij][iiterrs]->Fill(Xdev[ij], Ydev[ij]);}
                    ylayer_exterr[ij][iiterrs]->Fill(yexter[ij]);
		    ylayer_exterr_nTermSplit[ij][int(yextloc[ij]/20)][iiterrs]->Fill(yexter[ij]);
		    int iyy = min(yhits[ij],nmxhits);
		    if (iyy>0 && Ypos[ij]>1.0 && Ypos[ij]<nstrip-1 ) {
		      ylayer_reso_mul[ij][iiterrs][iyy-1]->Fill(Ydev[ij]);
		      if(isfid_l9_eff) { ylayer_reso_mul_split[ij][xbina][iiterrs][iyy-1]->Fill(Ydev[ij]);}
		    }


		  }
                }
                nTotally++;
                 for (int ij=0; ij<nlayer; ij++) {
		  if (abs(Ydev[ij])<1.5 && Yusedpos[ij]) {
		    for (int jk=0; jk<nlayer; jk++) {
		      if (abs(Ydev[jk])<1.5 && Yusedpos[jk]) {
			h_posycoreff ->Fill(ij,jk) ;
			posycoreffCount +=1 ;
		      }
		    }
		  }
	     }
	      } else {
		//	ystrip_ydev[occulyr][iiterrs]->Fill(yextloc[occulyr],Ydev[occulyr]);
		//	//int xbina = (nsplit/8)*int(xextloc[occulyr]/7.25)+int(yextloc[occulyr]/7.625);
		//int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[occulyr]+0.5),int(yextloc[occulyr]+0.5)));//(nsplit/8)*int(xextloc[occulyr]/7.25)+int(yextloc[occulyr]/7.625);
		//if(xbina>nsplit-1) xbina = nsplit-1;
		if (abs(Ydev[occulyr])<6.0 && Yusedpos[occulyr] && yextloc[occulyr]>firstYstrip && yextloc[occulyr]<lastYstrip && posinuse[1][occulyr][int(yextloc[occulyr]+0.5)]==true) {
                  ylayer_reso[occulyr][iiterrs]->Fill(Ydev[occulyr]);
		  ylayer_reso_nTermSplit[occulyr][int(yextloc[occulyr]/20)][iiterrs]->Fill(Ydev[occulyr]);
		  if(xextloc[occulyr]>0 && xextloc[occulyr]<lastXstrip && yextloc[occulyr]>0 && yextloc[occulyr]<lastYstrip ) { ylayer_reso_split[occulyr][xbina][iiterrs]->Fill(Ydev[occulyr]);}
		  if (is_xyreso[occulyr][iiterrs]) {xylayer_reso[occulyr][iiterrs]->Fill(Xdev[occulyr], Ydev[occulyr]);}
                  ylayer_exterr[occulyr][iiterrs]->Fill(yexter[occulyr]);
		  ylayer_exterr_nTermSplit[occulyr][int(yextloc[occulyr]/20)][iiterrs]->Fill(yexter[occulyr]);
		  int iyy = min(yhits[occulyr],nmxhits);
		  if (iyy>0 && Ypos[occulyr]>1.0 && Ypos[occulyr]<nstrip-1 ) {
		    ylayer_reso_mul[occulyr][iiterrs][iyy-1]->Fill(Ydev[occulyr]);
		    if(isfid_l9_eff ) { ylayer_reso_mul_split[occulyr][xbina][iiterrs][iyy-1]->Fill(Ydev[occulyr]);}
		  }
		}

		if (iiter==nmxiter-1) { //Correlated hits
		  double expp = -100;
		  double expdiff = 10000;
		  for (int ix=0; ix<ypts[occulyr].size(); ix++) {
		    if (abs(yextloc[occulyr] - ypts[occulyr][ix])<1.0) {
		      double tmpexpdiff = abs(yextloc[occulyr] - ypts[occulyr][ix]);
		      if(tmpexpdiff < expdiff) {
			 expdiff =tmpexpdiff;  expp = ypts[occulyr][ix];// break;
		      }
		    }
		  }
		  if (expp >-50) {
		    h_ycorstrips[occulyr]->Fill(expp, -1.0);
		    h_yxcorstrips[occulyr]->Fill(expp, -1.0);
		    for (int ix=0; ix<ypts[occulyr].size(); ix++) {
		      if (abs(expp - ypts[occulyr][ix])>0) { //choose all of them and then use maximum value before plot
			h_ycorstrips[occulyr]->Fill(expp, ypts[occulyr][ix]);
		      }
		    }

		    for (int ix=0; ix<xpts[occulyr].size(); ix++) {
		      h_yxcorstrips[occulyr]->Fill(expp, xpts[occulyr][ix]);
		    }
		  }
		} //if (iiter==nmxiter-1)
              }
            }

            if (nxfail==0 && isfill) {
              h_chisqx->Fill(xchi2);
              if (Nx-2>0) {
		h_reduchisqx->Fill(xchi2/(Nx-2));
		double probx = TMath::Prob(xchi2, Nx-2);
		h_xprob->Fill(probx);
		int ibin = getbinid(Nx, nprob, probs);
		if (ibin>=0 && ibin<nprob) {h_xnprob[ibin]->Fill(probx);}
	      }
              h_xndf->Fill(Nx);
            }
            if (nyfail==0 && isfill) {
	      h_chisqy->Fill(ychi2);
	      if (Ny-2>0) {
		h_reduchisqy->Fill(ychi2/(Ny-2));
		double probx =TMath::Prob(ychi2, Ny-2);
		h_yprob->Fill(probx);
		int ibin = getbinid(Ny, nprob, probs);
		if (ibin>=0 && ibin<nprob) {h_ynprob[ibin]->Fill(probx);}
	      }
	      h_yndf->Fill(Ny);
	    }

	    //	    cout <<"nxy "<<iev<<" "<< Nx<<" "<<Ny<<" "<<xchi2<<" "<<ychi2<<" "<<nxfail<<" "<<nyfail<<endl;


	    double theta11 = acos(sqrt(1./(1+pow(stripwidth*xslope1,2.)+pow(stripwidth*yslope1,2.))));
	    double theta1 = (180./pival)*theta11;
	    double theta22 = acos(sqrt(1./(1+pow(stripwidth*xslope2,2.)+pow(stripwidth*yslope2,2.))));
	    double theta2 = (180./pival)*theta22;
	    double theta_diff_val1 = 0.;//fabs(theta1-theta2);

	    int  xext1l5 = int(stripwidth*xextloc1[5]); int yext1l5 = int(stripwidth*yext1[5]);
	    int xstp = int(xextloc1[5]); int ystp = int(yext1[5]);
	    double zval1[nlayer];
	    for(int ij=0;ij<nlayer;ij++) { zval1[ij] = zval[ij]/stripwidth;}
	    double R0 = sqrt(pow(xextloc1[5]-xextloc1[0],2.)+pow(yextloc1[5]-yextloc1[0],2.)+pow(zval1[5]-zval1[0],2.));
	    double A0 = (xextloc1[5]-xextloc1[0])/R0;
	    double B0 = (yextloc1[5]-yextloc1[0])/R0;
	    double C0 = (zval1[5]-zval1[0])/R0;
	    double R1 = sqrt(pow(xextloc2[11]-xextloc2[6],2.)+pow(yextloc2[11]-yextloc2[6],2.)+pow(zval1[11]-zval1[6],2.));
	    double A1 = (xextloc2[11]-xextloc2[6])/R1;
	    double B1 = (yextloc2[11]-yextloc2[6])/R1;
	    double C1 = (zval1[11]-zval1[6])/R1;
	    double theta_diff_val = 0.;//fabs(acos(A0*A1+B0*B1+C0*C1))*180./TMath::Pi();

	    if(Nx1>3 && Ny1>3 && Nx2>3 && Ny2>3 && ychi2/(Ny-2)<15 && nyfail==0 && xchi2/(Nx-2)<15 && nxfail==0 && Nx>=6 && Ny>=6 ) {

	      if (isfill) {

	       	theta_diff_val = abs(acos(A0*A1+B0*B1+C0*C1))*180./TMath::Pi();
		theta_diff_val1 = abs(theta1-theta2);
		//	cout<<"theta_diff"<<"    "<<theta_diff_val<<endl;
		//	if(theta_diff_val1<10.0) {pixel_scatang1[xstp][ystp]->Fill(theta_diff_val1); }
		//	cout<<"xxxxx"<<"    "<<xstp<<"   "<<ystp<<"    "<<theta_diff_val<<endl;
		if(xstp>=0. && xstp <28.0 && ystp>=0. && ystp <31.0){pixel_scatang[xstp][ystp]->Fill(theta_diff_val,1.0); }
		theta_diff_1->Fill(theta_diff_val1);
		theta_diff_2->Fill(theta_diff_val);
		if(theta_diff_val1>0.0 && theta_diff_val1<=.3) pixel_diff_theta[0]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>0.3 && theta_diff_val1<=.6) pixel_diff_theta[1]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>0.6 && theta_diff_val1<=.9) pixel_diff_theta[2]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>0.9 && theta_diff_val1<=1.2) pixel_diff_theta[3]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>1.2 && theta_diff_val1<=1.5) pixel_diff_theta[4]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>1.5 && theta_diff_val1<=1.8) pixel_diff_theta[5]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>1.8 && theta_diff_val1<=2.1) pixel_diff_theta[6]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>2.1 && theta_diff_val1<=2.4) pixel_diff_theta[7]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>2.4 && theta_diff_val1<=2.7) pixel_diff_theta[8]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>2.7 && theta_diff_val1<=3.0) pixel_diff_theta[9]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>3.0 && theta_diff_val1<=3.3) pixel_diff_theta[10]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>3.3 && theta_diff_val1<=3.6) pixel_diff_theta[11]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>3.6 && theta_diff_val1<=3.9) pixel_diff_theta[12]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>3.9 && theta_diff_val1<=4.2) pixel_diff_theta[13]->Fill(xext1l5,yext1l5);
		if(theta_diff_val1>4.2 && theta_diff_val1<=4.5) pixel_diff_theta[14]->Fill(xext1l5,yext1l5);
		if(theta2<10.0){ if(theta_diff_val>=0.5 && theta_diff_val<=2.5) pixel_diff_theta[15]->Fill(xext1l5,yext1l5); }
	      }
	    }

	    double fitthe1 = acos(sqrt(1./(1+pow(stripwidth*xslope,2.)+pow(stripwidth*yslope,2.)))); //10Nov//S= sqrt(dx^2+dy^2+dz^2)) theta= acos(height / S) height = dz;
	    double fitthe = (180./pival)*fitthe1; //acos(sqrt(1./(1+pow(stripwidth*xslope,2.)+pow(stripwidth*yslope,2.))));

#ifdef C217STRIP
	    double fitphi = atan2(-xslope, yslope);
#else
	    double fitphi = atan2(yslope, xslope);  // What is the direction
#endif

// 	    if (isfill) {
// 	      if (Ny>=nmnhits && ychi2/(Ny-2)<mxchisq && nyfail==0
// 		  && Nx>=nmnhits && xchi2/(Nx-2)<mxchisq && nxfail==0) {
// 		T3->Fill();
// 	      } else {
// 		EvNum *=-1;
// 		T3->Fill();
// 	      }
// 	    }

	    if (ntcor==1 && iiter==nmxiter-1) {
	      bool xexit=false;
	      bool yexit=false;
	      int nxc=0;
	      int nyc=0;
	      for (int ix=0; ix<xptsall[occulyr].size(); ix++) {
		if (abs(xextloc[occulyr] - xptsall[occulyr][ix]) <1.1) {
		  xexit = true; nxc++;
		}
	      }
	      nxc = xptsall[occulyr].size() - nxc;

	      for (int iy=0; iy<yptsall[occulyr].size(); iy++) {
		if (abs (yextloc[occulyr] - yptsall[occulyr][iy]) <1.1) {
		  yexit = true; nyc++;
		}
	      }
	      nyc = yptsall[occulyr].size() - nyc;

	      bool muselec[nMuselec];
	      muselec[0] = (Nx>=nmnhits  && ychi2/(Nx-2)<mxchisq && nxfail==0) ? true : false;
	      muselec[1] = (Ny>=nmnhits  && ychi2/(Ny-2)<mxchisq && nyfail==0) ? true : false;
	      muselec[2] = muselec[0] && muselec[1];
	      muselec[3] = muselec[0] && (!muselec[1]);
	      muselec[4] = (!muselec[0]) && (muselec[1]);
	      muselec[5] = (!muselec[0]) && (!muselec[1]);
	      muselec[7] = (!muselec[0]);
	      muselec[8] = (!muselec[1]);

	      bool hitselec[nEnv];
	      hitselec[0] = (xexit && yexit);
	      hitselec[1] = (xexit && (!yexit));
	      hitselec[2] = ((!xexit) && yexit);
	      hitselec[3] = ((!xexit) && (!yexit));

	      for (int ix=0; ix<nMuselec; ix++) {
		for (int iy=0; iy<nEnv; iy++) {
		  if (muselec[ix] && hitselec[iy]) {
		    mult2d_selecxy[occulyr][ix][iy]->Fill(nxc, nyc);
		  }
		}
	      }
	       }

            if (Ny>=nmnhits /*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {//4Nov
	      if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {
		if (filloccu) {
		  for (int iz=0; iz<nlayer; iz++) {
		    xlayer_allmumult[iz]->Fill(xptsall[iz].size());
		    ylayer_allmumult[iz]->Fill(yptsall[iz].size());
		  }
		}



		//       C217
		//       x1 x2 x3 ...............................x30   x31
		//       ^  ^  ^  ^  ^  ^  ^  ^  W  ^  ^  ^  ^  ^  ^  ^  ^
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y1
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y2
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//  S<-- |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.--> N
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y30
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y31
                //                               |
		//                               E

		pos_theta[occulyr][iiter]->Fill(fitthe);
		// sel_theta_mom->Fill(fitthe,momin);
		//  sel_phi_mom->Fill(fitphi*180./pival,momin);
		pos_phi[occulyr][iiter]->Fill(fitphi*180./pival);
		if (isfill) {
#ifdef MONTECARLO
#ifdef FULLG4SIM
		  sel_theta_mom->Fill(fitthe,abs(momin[0]));
		  sel_phi_mom->Fill(fitphi*180./pival,abs(momin[0]));
		  sel_theta_rec_gen->Fill(60.*(thegen[0]),fitthe);
#endif
#endif
		  costhe[0]->Fill(fitthe, 1.0);
		  phiang[0]->Fill(fitphi, 1.0);

		  if ((xpts[trigly1].size()>0 &&
		      xpts[trigly2].size()>0 &&
		      xpts[trigly3].size()>0 &&
		      xpts[trigly4].size()>0)||
		      (ypts[trigly1].size()>0 &&
		      ypts[trigly2].size()>0 &&
		      ypts[trigly3].size()>0 &&
		       ypts[trigly4].size()>0)){
		    costhe[1]->Fill(fitthe, 1.0);
		    phiang[1]->Fill(fitphi, 1.0);

		    if (Xpos[trigly1] >-99 && Xpos[trigly2] >-99 &&
			Xpos[trigly3] >-99 && Xpos[trigly4] >-99) {
		      costhe[2]->Fill(fitthe, 1.0);
		      phiang[2]->Fill(fitphi, 1.0);

		      if (Ypos[trigly1] >-99 && Ypos[trigly2] >-99 &&
			  Ypos[trigly3] >-99 && Ypos[trigly4] >-99) {
			costhe[3]->Fill(fitthe, 1.0);
			phiang[3]->Fill(fitphi, 1.0);

			if (abs(Xdev[trigly1]) < 1.5 && Xusedpos[trigly1] && abs(Ydev[trigly1]) < 1.5 &&
			    abs(Xdev[trigly2]) < 1.5 && Xusedpos[trigly2] && abs(Ydev[trigly2]) < 1.5 &&
			    abs(Xdev[trigly3]) < 1.5 && Xusedpos[trigly3] && abs(Ydev[trigly3]) < 1.5 &&
			    abs(Xdev[trigly4]) < 1.5 && Xusedpos[trigly4] && abs(Ydev[trigly4]) < 1.5) {
			  costhe[4]->Fill(fitthe, 1.0);
			  phiang[4]->Fill(fitphi, 1.0);
			  costhe[5]->Fill(cos(fitthe1),1.0);//10Nov
			  zen = fitthe;
			}
		      }
		    }
		  }
		} // if (isfill)
	      } // if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0)
            } // if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0)

	    //9th July 2016 : redundant costheta distributions for "n" check
	    if (isfill && nxfail==0 && nyfail==0) {
	      bool passx[ntrkselcrit]={0};
	      int nxyval[7] = {3,4,5,6,7,8,9};//{6,7,8,9,10,11,12};
	      double xychisval[20]= {2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0};
	      for(int nxych =0;nxych<20;nxych++) {
		for(int nxy =0;nxy<7;nxy++) {
		  if (Nx>=nxyval[nxy] && Ny>=nxyval[nxy] && xchi2/(Nx-2)<xychisval[nxych] && ychi2/(Ny-2)<xychisval[nxych]) {
		    passx[nxy+nxych*7]=1;}
		}
	      }

	      //  if(fitphi<0.) { fitphi= fitphi+2*pival;}

	      for (int jk=0; jk<ntrkselcrit; jk++) {
		if (passx[jk]) {
		  //zenithang[0][jk]->Fill(fitthe, 1.0);
		  //  sel_theta_phi_reco[0][jk]->Fill(fitphi*180./pival,fitthe);
		  //zenithang[0][jk]->Fill(fitthe, 1.0);
		  // azimuthang[0][jk]->Fill(fitphi*180./pival, 1.0);

		  if ((xpts[trigly1].size()>0 &&
		       xpts[trigly2].size()>0 &&
		       xpts[trigly3].size()>0 &&
		       xpts[trigly4].size()>0)||
		      (ypts[trigly1].size()>0 &&
		       ypts[trigly2].size()>0 &&
		       ypts[trigly3].size()>0 &&
		       ypts[trigly4].size()>0)) {
		    if ((Xpos[trigly1] >-99  &&
			 Xpos[trigly2] >-99 &&
			 Xpos[trigly3] >-99  &&
			 Xpos[trigly4] >-99 ) ||
			(Ypos[trigly1] >-99  &&
			 Ypos[trigly2] >-99 &&
			 Ypos[trigly3] >-99  &&
			 Ypos[trigly4] >-99 ))	{



		      sel_theta_phi_reco[0][jk]->Fill(fitphi*180./pival,fitthe);
		      zenithang[0][jk]->Fill(fitthe, 1.0);
		      azimuthang[0][jk]->Fill(fitphi*180./pival, 1.0);

		      //  if(abs(momin[0])>0.25) {
		      sel_theta_phi_reco[1][jk]->Fill(fitphi*180./pival,fitthe);
		      zenithang[1][jk]->Fill(fitthe, 1.0);
		      azimuthang[1][jk]->Fill(fitphi*180./pival, 1.0);
		      //  }
		      //  if(abs(momin[0])>0.45) {
		      sel_theta_phi_reco[2][jk]->Fill(fitphi*180./pival,fitthe);
		      zenithang[2][jk]->Fill(fitthe, 1.0);
		      azimuthang[2][jk]->Fill(fitphi*180./pival, 1.0);
		      //  }
		      //  if(abs(momin[0])>0.65) {
		      sel_theta_phi_reco[3][jk]->Fill(fitphi*180./pival,fitthe);
		      zenithang[3][jk]->Fill(fitthe, 1.0);
		      azimuthang[3][jk]->Fill(fitphi*180./pival, 1.0);
		      //  }
		      //  if(abs(momin[0])>0.85) {
		      sel_theta_phi_reco[4][jk]->Fill(fitphi*180./pival,fitthe);
		      zenithang[4][jk]->Fill(fitthe, 1.0);
		      azimuthang[4][jk]->Fill(fitphi*180./pival, 1.0);
		      //  }
		      // if(abs(momin[0])>1.05) {
		      sel_theta_phi_reco[5][jk]->Fill(fitphi*180./pival,fitthe);
		      zenithang[5][jk]->Fill(fitthe, 1.0);
		      azimuthang[5][jk]->Fill(fitphi*180./pival, 1.0);
		      // }


		    //   if (Xpos[trigly1] >-99 &&  Xpos[trigly2] >-99 &&
		    //	Xpos[trigly3] >-99 &&  Xpos[trigly3] >-99  ) {

		      //  sel_theta_phi_reco[2][jk]->Fill(fitphi*180./pival,fitthe);
		      // zenithang[2][jk]->Fill(fitthe, 1.0);
		      //azimuthang[2][jk]->Fill(fitphi*180./pival, 1.0);

		         if (Ypos[trigly1] >-99 && Ypos[trigly2] >-99 &&
		      	  Ypos[trigly3] >-99 && Ypos[trigly4] >-99) {

			//	if (Xpos[1] >-99 && Xpos[2] >-99 &&
			//	Xpos[9] >-99 && Xpos[10] >-99) {
			//  if (Ypos[1] >-99 && Ypos[2] >-99 &&
			//	  Ypos[9] >-99 && Ypos[9] >-99) {
			//	cout<<"       Phi value            "<<double(fitphi)*180./pival<<endl;

			   //  sel_theta_phi_reco[3][jk]->Fill(fitphi*180./pival,fitthe);
			   //zenithang[3][jk]->Fill(fitthe, 1.0);
			   //azimuthang[3][jk]->Fill(fitphi*180./pival, 1.0);


		      // //  if (abs(momin[0]) > 0.25) {
			// zenithang[1][jk]->Fill(fitthe, 1.0);
			// azimuthang[1][jk]->Fill(fitphi*180./pival, 1.0);
			// //  if (abs(momin[0]) > 0.35) {
			// zenithang[2][jk]->Fill(fitthe, 1.0);
			// azimuthang[2][jk]->Fill(fitphi*180./pival, 1.0);
			// //   if (abs(momin[0]) > 0.45) {
			// zenithang[3][jk]->Fill(fitthe, 1.0);
			// azimuthang[3][jk]->Fill(fitphi*180./pival, 1.0);
			// //    if (abs(momin[0]) > 0.55) {
			// zenithang[4][jk]->Fill(fitthe, 1.0);
			// azimuthang[4][jk]->Fill(fitphi*180./pival, 1.0);
			// //   if (abs(momin[0]) > 0.65) {
			// zenithang[5][jk]->Fill(fitthe, 1.0);
			// azimuthang[5][jk]->Fill(fitphi*180./pival, 1.0);
			// //}
			// //}
			// //}
			// //}
			// //}

			//	if(fitphi<0.) { fitphi= fitphi+2*pival;}
			/*
			if(fitphi*180./pival >=0. && fitphi*180./pival <45.0) { zenithang_azimuth_8[0][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=45.0 &&  fitphi*180./pival <90.0) {zenithang_azimuth_8[1][jk]->Fill(fitthe, 1.0); }
			else if(fitphi*180./pival >=90.0 &&  fitphi*180./pival <135.0) {zenithang_azimuth_8[2][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=135.0 &&  fitphi*180./pival <180.0){ zenithang_azimuth_8[3][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=180.0 &&  fitphi*180./pival <225.0){ zenithang_azimuth_8[4][jk]->Fill(fitthe, 1.0); }
			else if(fitphi*180./pival >=225.0 &&  fitphi*180./pival <270.0){ zenithang_azimuth_8[5][jk]->Fill(fitthe, 1.0); }
			else if(fitphi*180./pival >=270.0 &&  fitphi*180./pival <315.0){ zenithang_azimuth_8[6][jk]->Fill(fitthe, 1.0); }
			else if(fitphi*180./pival >=315.0 &&  fitphi*180./pival <360.0) {zenithang_azimuth_8[7][jk]->Fill(fitthe, 1.0);}

		    // if(fitphi*180./pival >=0. && fitphi*180./pival <45.0


			if(fitphi*180./pival >=0. && fitphi*180./pival <22.5){ zenithang_azimuth[0][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=22.5 && fitphi*180./pival <45.0) {zenithang_azimuth[1][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=45.0 && fitphi*180./pival <67.5) {zenithang_azimuth[2][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=67.5 && fitphi*180./pival <90.0) {zenithang_azimuth[3][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=90.0 && fitphi*180./pival <112.5) {zenithang_azimuth[4][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=112.5 && fitphi*180./pival <135.0) {zenithang_azimuth[5][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=135.0 && fitphi*180./pival <157.5) {zenithang_azimuth[6][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=157.5 && fitphi*180./pival <180.0) {zenithang_azimuth[7][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=180.0 && fitphi*180./pival <202.5) {zenithang_azimuth[8][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=202.5 && fitphi*180./pival <225.0) {zenithang_azimuth[9][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=225.0 && fitphi*180./pival <247.5) {zenithang_azimuth[10][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=247.5 && fitphi*180./pival <270.0) {zenithang_azimuth[11][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=270.0 && fitphi*180./pival <292.5) {zenithang_azimuth[12][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=292.5 && fitphi*180./pival <315.0){ zenithang_azimuth[13][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=315.0 && fitphi*180./pival <337.5) {zenithang_azimuth[14][jk]->Fill(fitthe, 1.0);}
			else if(fitphi*180./pival >=337.5 && fitphi*180./pival <360.0) {zenithang_azimuth[15][jk]->Fill(fitthe, 1.0);}
			*/
			// if(fitphi*180./pival >=348.75 || fitphi*180./pival <11.25){ zenithang_azimuth[0][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=11.25 && fitphi*180./pival <33.75) {zenithang_azimuth[1][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=33.75 && fitphi*180./pival <56.25) {zenithang_azimuth[2][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=56.25 && fitphi*180./pival <78.75) {zenithang_azimuth[3][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=78.75 && fitphi*180./pival <101.25) {zenithang_azimuth[4][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=101.25 && fitphi*180./pival <123.75) {zenithang_azimuth[5][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=123.75 && fitphi*180./pival <146.25) {zenithang_azimuth[6][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=146.25 && fitphi*180./pival <168.75) {zenithang_azimuth[7][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=168.75 && fitphi*180./pival <191.25) {zenithang_azimuth[8][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=191.25 && fitphi*180./pival <213.75) {zenithang_azimuth[9][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=213.75 && fitphi*180./pival <236.25) {zenithang_azimuth[10][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=236.25 && fitphi*180./pival <258.75) {zenithang_azimuth[11][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=258.75 && fitphi*180./pival <281.25) {zenithang_azimuth[12][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=281.25 && fitphi*180./pival <303.75){ zenithang_azimuth[13][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >= 303.75&& fitphi*180./pival <326.25) {zenithang_azimuth[14][jk]->Fill(fitthe, 1.0);}
			// else if(fitphi*180./pival >=326.25 && fitphi*180./pival <348.75) {zenithang_azimuth[15][jk]->Fill(fitthe, 1.0);}


			//  zenithang[2][jk]->Fill(fitthe, 1.0);
			//  if (Ypos[trigly1] >-99 && Ypos[trigly2] >-99 &&
			// Ypos[trigly3] >-99 && Ypos[trigly4] >-99) {
			//zenithang[3][jk]->Fill(fitthe, 1.0);
			if (abs(Xdev[trigly1]) < 1.5 && Xusedpos[trigly1]  &&
			    abs(Xdev[trigly2]) < 1.5 && Xusedpos[trigly2]  &&
			    abs(Xdev[trigly3]) < 1.5 && Xusedpos[trigly3]  &&
			    abs(Xdev[trigly4]) < 1.5 && Xusedpos[trigly4] ) {
			  // zenithang[3][jk]->Fill(fitthe, 1.0);
			  //  sel_theta_phi_reco[4][jk]->Fill(fitphi*180./pival,fitthe);
			  //zenithang[4][jk]->Fill(fitthe, 1.0);
			  //azimuthang[4][jk]->Fill(fitphi*180./pival, 1.0);
			  //zenithang[4][jk]->Fill(fitthe, 1.0);
			  if(  abs(Ydev[trigly1]) < 1.5 && Yusedpos[trigly1]  &&
			       abs(Ydev[trigly2]) < 1.5 && Yusedpos[trigly2]  &&
			       abs(Ydev[trigly3]) < 1.5 && Yusedpos[trigly3]  &&
			       abs(Ydev[trigly4]) < 1.5 && Yusedpos[trigly4]) {
			    // zenithang[5][jk]->Fill(fitthe, 1.0);

			    //sel_theta_phi_reco[5][jk]->Fill(fitphi*180./pival,fitthe);
			    //zenithang[5][jk]->Fill(fitthe, 1.0);
			    //azimuthang[5][jk]->Fill(fitphi*180./pival, 1.0);

			    //  }
			    //	}
			  }
			}
		         }
		    }
		  }
		}
	      } // for (int jk=0; jk<ntrkselcrit; jk++)
	    } // if (isfill && nxfail==0 && nyfail==0)

	    if (isfill) {
	      nTotalp++;
	      for (int ix=0; ix<nlayer; ix++) {
		if (Xpos[ix]>=-1 && Xpos[ix]<=nstrip) {
		  for (int iy=ix+1; iy<nlayer; iy++) {
		    if (Xpos[iy]>=-1 && Xpos[iy]<=nstrip) {
		      h_xcorhits->Fill(ix, iy);
		    }
		  }
		  for (int iy=0; iy<nlayer; iy++) {
		    if (Ypos[iy]>=-1 && Ypos[iy]<=nstrip) {
		      h_xycorhits->Fill(ix, iy);
		    }
		  }
		}
	      }

	      for (int ix=0; ix<nlayer-1; ix++) {
		if (Ypos[ix]>=-1 && Ypos[ix]<=nstrip) {
		  for (int iy=ix+1; iy<nlayer; iy++) {
		    if (Ypos[iy]>=-1 && Ypos[iy]<=nstrip) {
		      h_ycorhits->Fill(ix, iy);
		    }
		  }
		}
	      }
	    }

	    //15th Oct
	    //	    if (isfiducial) {
	    if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {
	      if (occulyr>=nlayer) {
		for (int ij=0; ij<nlayer; ij++) {
		  if ( xexter[ij]<0.5) {
		    strp_xmul[ij][iiterrs]->Fill(xposinstr[ij], xhits[ij]);
		  }
		}
	      } else {
		strp_xmul[occulyr][iiterrs]->Fill(xposinstr[occulyr], xhits[occulyr]);
	      }
	    }
	    if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {
	      if (occulyr>=nlayer) {
		for (int ij=0; ij<nlayer; ij++) {
		  if ( yexter[ij]<0.5) {
		    strp_ymul[ij][iiterrs]->Fill(yposinstr[ij], yhits[ij]);
		  }
		}
	      } else {
		strp_ymul[occulyr][iiterrs]->Fill(yposinstr[occulyr], yhits[occulyr]);
	      }
	    }
	    //	    }


            if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0 &&
                Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0 ) {//4Nov
	      //position resolution as a function of position and multiplicity
	      int istrt = (ntcor==0) ? 0 : occulyr;
	      int iend = (ntcor==0) ? nlayer : occulyr+1;
	      bool whichFill = (ntcor==0) ? true : false;
	      for (int ij=istrt; ij<iend; ij++) {
		if ( yexter[ij]<0.4 &&  xexter[ij]<0.4) {
		  int ixx = xpts[ij].size()-1; if (ixx>nstr_posmx-1) ixx=nstr_posmx-1;
		  int iyy = ypts[ij].size()-1; if (iyy>nstr_posmx-1) iyy=nstr_posmx-1;
		  if (iiter==nmxiter-1) {
		    if (whichFill) {
		      if (ixx>=0) {xpos_xdev[ij][ixx]->Fill(xposinstr[ij], Xdev[ij]);}
		      if (iyy>=0) {ypos_ydev[ij][iyy]->Fill(yposinstr[ij], Ydev[ij]);}
		    } else {
		      if (ixx>=0) {ypos_xdev[ij][ixx]->Fill(yposinstr[ij], Xdev[ij]);}
		      if (iyy>=0) {xpos_ydev[ij][iyy]->Fill(xposinstr[ij], Ydev[ij]);}
		    }
		  }

		  if (ixx>=0 && ixx<nstr_posmx) {
		    if (whichFill) {
		      xstr_xdev[ij][iiter][ixx]->Fill(xextloc[ij], Xdev[ij]);
		    } else {
		      //  cout<<ij<<"   "<<iiter<<"   "<<ixx<<"   "<<ystr_xdev[ij][iiter][ixx]->GetName()<<"   "<<yextloc[ij]<<"    "<<Xdev[ij]<<endl;
		      xstr_xdev_occu[ij][iiter][ixx]->Fill(xextloc[ij], Xdev[ij]);
		      ystr_xdev[ij][iiter][ixx]->Fill(yextloc[ij], Xdev[ij]);
		    }
		    if (ixx<nstr_posmx-1) {
		      bool spass=true;
		      int ixz=int(xextloc[ij]+0.5);
		      for (int iz=0; iz<deadchannel[0][ij].size(); iz++) {
			if (ixz==deadchannel[0][ij][iz]) { spass=false; break;}
		      }
		      if (spass) {
			if (whichFill) {
			  prof_xstr_xdev[ij][iiter]->Fill(xextloc[ij], Xdev[ij]);
			} else {
			  prof_ystr_xdev[ij][iiter]->Fill(yextloc[ij], Xdev[ij]);
			}
		      }
		    }
		  }

		  if (iyy>=0 && iyy<nstr_posmx) {
		    if (whichFill) {
		      ystr_ydev[ij][iiter][iyy]->Fill(yextloc[ij], Ydev[ij]);
		    } else {
		      ystr_ydev_occu[ij][iiter][iyy]->Fill(yextloc[ij], Ydev[ij]);
		      xstr_ydev[ij][iiter][iyy]->Fill(xextloc[ij], Ydev[ij]);
		    }
		    if (iyy<nstr_posmx-1) {
		      bool spass=true;
		      int ixz=int(yextloc[ij]+0.5);
		      for (int iz=0; iz<deadchannel[1][ij].size(); iz++) {
			if (ixz==deadchannel[1][ij][iz]) { spass=false; break;}
		      }
		      if(spass) {
			if (whichFill) {
			  prof_ystr_ydev[ij][iiter]->Fill(yextloc[ij], Ydev[ij]);
			} else {
			  prof_xstr_ydev[ij][iiter]->Fill(xextloc[ij], Ydev[ij]);
			}
		      }
		    }
		  }
		}
	      }

	      if (occulyr>=nlayer) {
                for (int ij=0; ij<nlayer; ij++) {
		  //int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[ij]+0.5),int(yextloc[ij]+0.5)));//(nsplit/8)*int(xextloc[ij]/7.25)+int(yextloc[ij]/7.625);
		  //if(xbina>nsplit-1) xbina = nsplit-1;
		  if ( yexter[ij]<0.4 &&  xexter[ij]<0.4) {

		    bool  isfiducial = (int(xextloc[ij])>1 && int(xextloc[ij])<58 && int(yextloc[ij])>1 && int(yextloc[ij])<61) ? 1 : 0;
		    double gap=accprng+effirng; //find_max_gap(ij);
                    if (abs(xposinstr[ij])+xexter[ij]<accprng && abs(yposinstr[ij])+yexter[ij]<accprng) {
#ifdef ISEFFICIENCY
		      totalentry[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]); //jamessantosh
		      fine_totalentry[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]); //jamessantosh
#ifndef MONTECARLO
		      if(filloccu) { totalentry_set[ij][itset]->Fill(xextloc[ij],yextloc[ij]); }
#endif
#endif
		      if (isfiducial) {

			//			strp_xmul[ij][iiterrs]->Fill(xposinstr[ij], xhits[ij]);
			//			strp_ymul[ij][iiterrs]->Fill(yposinstr[ij], yhits[ij]);
			if (Xusedpos[ij]) {
			  strp_xmulsim[ij][iiterrs]->Fill(xposinstr[ij], xhits[ij]);
			  if(isfid_l9_eff) { strp_xmulsim_split[ij][xbina][iiterrs]->Fill(xposinstr[ij], xhits[ij]);}

			} else {
			  strp_xmulsim[ij][iiterrs]->Fill(xposinstr[ij], 0.0);
			  if(isfid_l9_eff ) { strp_xmulsim_split[ij][xbina][iiterrs]->Fill(xposinstr[ij], 0.0);}

			}
			if (Yusedpos[ij]) {
			  strp_ymulsim[ij][iiterrs]->Fill(yposinstr[ij], yhits[ij]);
			  if(isfid_l9_eff) { strp_ymulsim_split[ij][xbina][iiterrs]->Fill(yposinstr[ij], yhits[ij]);}
		} else {
			  strp_ymulsim[ij][iiterrs]->Fill(yposinstr[ij], 0.0);
			  if(isfid_l9_eff ) { strp_ymulsim_split[ij][xbina][iiterrs]->Fill(yposinstr[ij], 0.);}
		}

#ifdef ISEFFICIENCY
			efficiency_xallpixel[ij][iiterrs]->Fill(xposinstr[ij]);
			efficiency_yallpixel[ij][iiterrs]->Fill(yposinstr[ij]);
			if (abs(Xdev[ij])<accprng && Xusedpos[ij]) {
			      if (isfiducial) { efficiency_xpixel[ij][iiterrs]->Fill(xposinstr[ij]);}
			}
			 if (abs(Ydev[ij])<accprng && Yusedpos[ij]) {
			      if (isfiducial) { efficiency_ypixel[ij][iiterrs]->Fill(yposinstr[ij]);}
			    }




#endif
		      }
#ifdef ISEFFICIENCY
		      if (xptsall[ij].size()>0) { triggereffi_xevt[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);}//jamessantosh
		      if (yptsall[ij].size()>0) { triggereffi_yevt[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);}//jamessantosh

		      for (int ix=0; ix<xptsall[ij].size(); ix++) {
			if (abs(xextloc[ij]- xptsall[ij][ix])<gap) {
			  triggereffi_x[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);//jamessantosh
			  fine_triggereffi_x[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);//jamessantosh
#ifndef MONTECARLO
			  if(filloccu) { triggereffi_x_set[ij][itset]->Fill(xextloc[ij],yextloc[ij]);  }
#endif
			  break;
			}
		      }

		      for (int iy=0; iy<yptsall[ij].size(); iy++) {
			if (abs(yextloc[ij]- yptsall[ij][iy])<gap) {
			  triggereffi_y[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);//jamessantosh
			  fine_triggereffi_y[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);//jamessantosh
#ifndef MONTECARLO
			  if(filloccu) { triggereffi_y_set[ij][itset]->Fill(xextloc[ij],yextloc[ij]);  }
#endif
			  break;
			}
		      }
		      //  if(ij==0) { cout<<Xdev[ij]<<"   "<<Ydev[ij]<<endl; }
		      bool tmpxptsall = false;
		      bool tmpyptsall = false;
		      for(int ix=0;ix<xptsall[ij].size();ix++) {
			tmpxptsall =  abs(xextloc[ij]- xptsall[ij][ix])<gap ? true:false;
			if(tmpxptsall) break;
		      }
		      for(int ix=0;ix<yptsall[ij].size();ix++) {
			tmpyptsall =  abs(yextloc[ij]- yptsall[ij][ix])<gap ? true:false;
			if(tmpyptsall) break;
		      }


		      if(!tmpxptsall) {
			if(!tmpyptsall) {
			  inefficiencytrue_corx[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);
			}
		      }
		      if(!tmpyptsall) {
			if(!tmpxptsall) {
			 inefficiencytrue_cory[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);
			}
		      }
		      if (abs(Xdev[ij])>accprng || (!Xusedpos[ij])) {
			if (abs(Ydev[ij])>accprng || (!Yusedpos[ij])) {
			  // if(ij==0) cout<<"PassY"<<"   "<<abs(Ydev[ij])<<endl;
			  inefficiency_corx[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);//jamessantosh
#ifndef MONTECARLO
			  if(filloccu) { inefficiency_corx_set[ij][itset]->Fill(xextloc[ij],yextloc[ij]); }

#endif
			} else {
			  inefficiency_uncx[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);//jamessantosh
			}
			xposEffPass[ij][iiterrs]=true;
		      }

		      if (abs(Ydev[ij])>accprng || (!Yusedpos[ij])) {
			if (abs(Xdev[ij])>accprng || (!Xusedpos[ij])) {
			  inefficiency_cory[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);//jamessantosh
			} else {
			  inefficiency_uncy[ij][iiterrs]->Fill(xextloc[ij],yextloc[ij]);//jamessantosh
			}
			yposEffPass[ij][iiterrs]=true;
		      }

#endif
		    } // if (abs(xposinstr[ij])+xexter[ij]<accprng && abs(yposinstr[ij])+yexter[ij]<accprng)
		  } // if ( yexter[ij]<0.4 &&  xexter[ij]<0.4)
		} // for (int ij=0; ij<nlayer; ij++)
	      } else { // if (occulyr>=nlayer)
	      //int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[occulyr]+0.5),int(yextloc[occulyr]+0.5)));//(nsplit/8)*int(xextloc[occulyr]/7.25)+int(yextloc[occulyr]/7.625);
	      if ( yexter[occulyr]<0.4 &&  xexter[occulyr]<0.4) {

		 //	 if(occulyr<=3||occulyr>=9) {
		 bool isfiducial = (int(xextloc[occulyr])>1 && int(xextloc[occulyr])<lastXstrip-1 &&
				     int(yextloc[occulyr])>1 && int(yextloc[occulyr])<lastYstrip-1) ? 1 : 0;
		  // } else {
		  //	  isfiducial = (int(xextloc[occulyr])>1 && int(xextloc[occulyr])<27 &&
		  //		     int(yextloc[occulyr])>1 && int(yextloc[occulyr])<30) ? 1 : 0;
		  //	 }
		  double gap=accprng+effirng; // find_max_gap(occulyr);

		  if (abs(xposinstr[occulyr])+xexter[occulyr]<accprng && abs(yposinstr[occulyr])+yexter[occulyr]<accprng) {
#ifdef ISEFFICIENCY
		    totalentry[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		    fine_totalentry[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		    if (xptsall[occulyr].size()>0) { triggereffi_xevt[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);}
		    if (yptsall[occulyr].size()>0) { triggereffi_yevt[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);}

		    for (int ix=0; ix<xptsall[occulyr].size(); ix++) {
		      if (abs(xextloc[occulyr]- xptsall[occulyr][ix])<gap) {
			triggereffi_x[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
			fine_triggereffi_x[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
			break;
		      }
		    }

		    for (int iy=0; iy<yptsall[occulyr].size(); iy++) {
		      if (abs(yextloc[occulyr]- yptsall[occulyr][iy])<gap) {
			triggereffi_y[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
			fine_triggereffi_y[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
			break;
		      }
		    }
#endif
		    if (isfiducial) {
		      //		      strp_xmul[occulyr][iiterrs]->Fill(xposinstr[occulyr], xhits[occulyr]);
		      //		      strp_ymul[occulyr][iiterrs]->Fill(yposinstr[occulyr], yhits[occulyr]);

		      if (Xusedpos[occulyr]) {
			strp_xmulsim[occulyr][iiterrs]->Fill(xposinstr[occulyr], xhits[occulyr]);
			if(isfid_l9_eff ) { strp_xmulsim_split[occulyr][xbina][iiterrs]->Fill(xposinstr[occulyr], xhits[occulyr]);}

		      } else {
			strp_xmulsim[occulyr][iiterrs]->Fill(xposinstr[occulyr], 0.0);
			if(isfid_l9_eff ) { strp_xmulsim_split[occulyr][xbina][iiterrs]->Fill(xposinstr[occulyr],0.0);}

		      }
		      if (Yusedpos[occulyr]) {
			strp_ymulsim[occulyr][iiterrs]->Fill(yposinstr[occulyr], yhits[occulyr]);
			if(isfid_l9_eff ) { strp_ymulsim_split[occulyr][xbina][iiterrs]->Fill(yposinstr[occulyr], yhits[occulyr]);}


		      } else {
			strp_ymulsim[occulyr][iiterrs]->Fill(yposinstr[occulyr], 0.0);
			if(isfid_l9_eff) { strp_ymulsim_split[occulyr][xbina][iiterrs]->Fill(yposinstr[occulyr], 0.);}

		      }

#ifdef ISEFFICIENCY
		      efficiency_xallpixel[occulyr][iiterrs]->Fill(xposinstr[occulyr]);
		      efficiency_yallpixel[occulyr][iiterrs]->Fill(yposinstr[occulyr]);
#endif
		    }
#ifdef ISEFFICIENCY


		    bool tmpxptsall = false;
		    bool tmpyptsall = false;
		    for(int ix=0;ix<xptsall[occulyr].size();ix++) {
		      tmpxptsall =  abs(xextloc[occulyr]- xptsall[occulyr][ix])<gap ? true:false;
		      if(tmpxptsall) break;
		    }
		    for(int ix=0;ix<yptsall[occulyr].size();ix++) {
		      tmpyptsall =  abs(yextloc[occulyr]- yptsall[occulyr][ix])<gap ? true:false;
		      if(tmpyptsall)  break;
		    }


		    if(!tmpxptsall) {
		      if(!tmpyptsall) {
			inefficiencytrue_corx[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		      }
		    }
		    if(!tmpyptsall) {
		      if(!tmpxptsall) {
			inefficiencytrue_cory[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		      }
		    }



		    if (abs(Xdev[occulyr])>gap || (!Xusedpos[occulyr])) {
		     // inefficiency_uncx[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		      if (abs(Ydev[occulyr])>gap || (!Yusedpos[occulyr])) {

			inefficiency_corx[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		      } else {
			inefficiency_uncx[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		      }
		      xposEffPass[occulyr][iiterrs]=true;
		    }
		    if (abs(Xdev[occulyr])<gap && Xusedpos[occulyr]) {
		      if (isfiducial) { efficiency_xpixel[occulyr][iiterrs]->Fill(xposinstr[occulyr]);}
		    }
		    if (abs(Ydev[occulyr])>gap || (!Yusedpos[occulyr])) {
		     // inefficiency_uncy[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		      if (abs(Xdev[occulyr])>gap || (!Xusedpos[occulyr])) {
			inefficiency_cory[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		      } else {
			inefficiency_uncy[occulyr][iiterrs]->Fill(xextloc[occulyr],yextloc[occulyr]);
		      }
		      yposEffPass[occulyr][iiterrs]=true;
		    }
		    if (abs(Ydev[occulyr])<gap && Yusedpos[occulyr]) {
		      if (isfiducial) { efficiency_ypixel[occulyr][iiterrs]->Fill(yposinstr[occulyr]);}
		    }
#endif
		  } //if (abs(dx)+xexter[occulyr]<accprng && abs(dy)+yexter[occulyr]<accprng)
		} //if ( yexter[occulyr]<0.4 &&  xexter[occulyr]<0.4)
	      } //else of  if (occulyr>=nlayer)
	    } //if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0 ....)


	    // cout<<iev<<"  "<<"PASSSSSSSSSSSS 3"<<endl;

            // in the root structure; otherwise it will fill after time fit
	    //GMA            if ( !isTiming && isfill &&  nxfail==0 &&  nyfail==0 && Nx>nmnhits && Ny>nmnhits && xchi2/(Ny-2) >mxchisq && ychi2/(Nx-2) >mxchisq) T2->Fill();
	    // if ((!isTiming) && isfill &&  nxfail==0 &&  nyfail==0) T2->Fill(); //07/09/2011
	    //	    T2->Fill();


            // Now filling T2 contains only proper track fit data
	    nxtime=0;
	    nytime=0;
	    ntxyla=0;
	    //	    cout <<"nX "<<Nx<<" "<<Ny<<" "<<xchi2<<" "<<ychi2<<" "<<nmnhits<<" "<<mxchisq<<" "<<nxfail<<" "<<nyfail<<endl;
	    //	    if (isTiming&& Nx>=nmnhits/*-ntcor*/ && Ny>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && ychi2/(Ny-2)<mxchisq && nxfail==0 && nyfail==0) {
	    if (isTiming) {
              //////////////////////////////////////////////
              //                                          //
              //  Timing informations and directionality  //
              //                                          //
              //////////////////////////////////////////////

	      //Store time informations
#ifndef MONTECARLO

	      int istartx[nlayer] = {0}; int iendx[nlayer] = {0};
	      int istarty[nlayer] = {0}; int iendy[nlayer] = {0};
	      istartx[0] = 190; iendx[0] = 210;istarty[0] = 190; iendy[0] = 210 ;
	      istartx[1] = 195; iendx[1] = 215;istarty[1] = 200; iendy[1] = 220;
	      istartx[2] = 195; iendx[2] = 215;istarty[2] = 195; iendy[2] = 215;
	      istartx[3] = 190; iendx[3] = 210;istarty[3] = 191; iendy[3] = 211;
	      istartx[4] = 200; iendx[4] = 220;istarty[4] = 190; iendy[4] = 210 ;
	      istartx[5] = 190; iendx[5] = 210;istarty[5] = 190; iendy[5] = 210;
	      istartx[6] = 190; iendx[6] = 210;istarty[6] = 185; iendy[6] = 205;
	      istartx[7] = 183; iendx[7] = 203;istarty[7] = 180; iendy[7] = 200;
	      istartx[8] = 183; iendx[8] = 203;istarty[8] = 180; iendy[8] = 200;
	      istartx[9] = 180; iendx[9] = 200;istarty[9] = 180; iendy[9] = 200;




		 for(int jk=0;jk<nlayer;jk++) {
		   int nxtdchitmult = 0;
		   int nytdchitmult = 0;
		   for (int itdc=0; itdc<nTDCpLayer; itdc++) {

		  int jkrd = (itdc==0) ? jk : jk + 32;
		  int jkfl = (itdc==0) ? jk : jk + 12;
		  //x-side tdc data
		  //		cout <<"time "<<iev<<" "<<jk<<" "<<event->tdcdata[jk] <<" "<<event->tdcdata[jk+16]<<endl;

		  int tmpxxtime1 = -100;
#if defined(ONETIMEFORMAT) && defined(OLDFORMAT)
		  int tmpxxtime = event->Tdc[jk];
#elif defined(NEWRPCDAQ1)
		  int tmpxxtime = -100;
		  int xvectsize = event->vxtdc_l[jk][itdc]->size();
		  nxtdchitmult += xvectsize;

		  // for(int jkl = 0;jkl<xvectsize;jkl++) {

		  if(event->vxtdc_l[jk][itdc]->size()){
		    //  if(jk==3 || jk==5) {6000.+int(event->vxtdc_l[jk][itdc]->at(0)-event->tdc_ref_l[jk])-440;}
		    //   else {
		    tmpxxtime = 6000.+int(event->vxtdc_l[jk][itdc]->at(0)-event->tdc_ref_l[jk]);
		    //}
		    tmpxxtime1 = 12500.+int(event->vxtdc_l[jk][itdc]->at(0)-event->tdc_ref_l[jk]);
		    h_time_layerprx[jk]->Fill(tmpxxtime);
		    h_time_layerprx1[jk]->Fill(tmpxxtime1);
		    // }
		    //   cout<<"X --- xtdc"<<"   "<<iev<<"  "<<xvectsize<<"  "<<jk<<"   "<<itdc<<"   "<<tmpxxtime<<"   "<<int(event->vxtdc_l[jk][itdc]->at(0)-event->tdc_ref_l[jk])<<endl;
		    if (filloccu) {
		      ntxtime2[jk][itdc]++;
		      int ncxpt=0;;
		      if(xpts[jk].size()>0){
			for(int kk=0;kk<xpts[jk].size();kk++) {
			  for(int ii=0;ii<8;ii++) {
			    if(xpts[jk][kk]==itdc+8*ii) { ncxpt++;  }

			  }
			}
			if(ncxpt>0) { ntotxext2[jk][itdc]++;}
		      }
		    }
		  }
#elif defined(BARC_EVTBLR)
		  int tmpxxtime = -100;
		  if(event->xtdc_l[jk][itdc]->GetEntries()) {
		    TDCInfo *xtdc = (TDCInfo*)event->xtdc_l[jk][itdc]->At(0);
		    tmpxxtime = 6000.+int(xtdc->GetValue()-event->tdc_ref_l[jk]);
		    tmpxxtime1 = 12500. +int(xtdc->GetValue()-event->tdc_ref_l[jk]);
		  }

#else
		  int tmpxxtime = event->tdcdata[jkrd];
#endif

#if defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
		  initxtime[jk][itdc]=(tmpxxtime<50.) ? -100.0 : 0.1*tmpxxtime;
#else
		  initxtime[jk][itdc]=(tmpxxtime<50.) ? -100.0 : 0.1*tmpxxtime;
#endif
		  if(filloccu) {
		    if (jk==0 && itdc==0) {
		      time_layerstrip->Fill(-1., -1.);
		      time_layerstripcpy->Fill(-1., -1.);
		    } //For total count
		    //printf("inside tdcdate read\n");
		    time_layer->Fill(jkfl, tmpxxtime);

		    //  if(itdc==0 && jk==3) rawtimex->Fill(tmpxxtime);

		    if (xallhits[jk][itdc]>=0 && xallhits[jk][itdc]<nstrip) {
		      time_xraw[jk][xallhits[jk][itdc]]->Fill(tmpxxtime);
		      time_layerstrip->Fill(nstrip*jk+xallhits[jk][itdc],tmpxxtime);
		      time_layerstripcpy->Fill(nstrip*jk+xallhits[jk][itdc],tmpxxtime1);

#ifdef TIMESLOPE
		      if (initxtime[jk][itdc]>0 && abs(xslope)<1.0&&abs(yslope)<1.0) {
			time_X_L7_corr[jk]->Fill(xextloc[7], initxtime[jk][itdc]);
			time_X_corr[jk][itdc]->Fill(xextloc[jk], initxtime[jk][itdc]);
			time_xslope[jk][itdc]->Fill(yextloc[7], initxtime[jk][itdc]);
		      }
#endif
		    }
		  }

		  jkrd +=16;
		  jkfl +=25; // one gap after nTDCpLayer*nlayer
		  int tmpyytime1 = -100;
		  //y-side tdc data
#if defined(ONETIMEFORMAT) && defined(OLDFORMAT)
		  int tmpyytime = event->Tdc[jk+16];
#elif defined(NEWRPCDAQ1)
		  int yvectsize = event->vytdc_l[jk][itdc]->size();

		  // for(int jkl = 0;jkl<xvectsize;jkl++) {
		  int tmpyytime= -100;

		  if(event->vytdc_l[jk][itdc]->size()) {
		  //if(jk==3 || jk==50) {tmpyytime = 6000.+int(event->vytdc_l[jk][itdc]->at(0)-event->tdc_ref_l[jk])-440;}
		  // else {
		    tmpyytime = 6000.+int(event->vytdc_l[jk][itdc]->at(0)-event->tdc_ref_l[jk]);
		    //}
		    // tmpyytime = 6000.+int(event->vytdc_l[jk][itdc]->at(0)-event->tdc_ref_l[jk]);
		    tmpyytime1 = 12500.+int(event->vytdc_l[jk][itdc]->at(0)-event->tdc_ref_l[jk]);

		    h_time_layerpry[jk]->Fill(tmpyytime);
		    h_time_layerpry1[jk]->Fill(tmpyytime1);
		    //    cout<<"Y --- ytdc"<<"   "<<iev<<"  "<<yvectsize<<"  "<<jk<<"   "<<itdc<<"   "<<tmpyytime<<"   "<<int(event->vytdc_l[jk][itdc]->at(0)-event->tdc_ref_l[jk])<<endl;
		    if (filloccu) {
		      ntytime2[jk][itdc]++;
		      int ncypt=0;;
		      if(ypts[jk].size()>0){
			for(int kk=0;kk<ypts[jk].size();kk++) {
			  for(int ii=0;ii<8;ii++) {
			    if(ypts[jk][kk]==itdc+8*ii) { ncypt++;  }

			  }
			}
			if(ncypt>0) { ntotyext2[jk][itdc]++;}
		      }
		    }
		  }
		  //	  cout<<iev<<"   "<<jk<<"  "<<itdc<<"   "<<xvectsize<<"   "<<yvectsize<<"  "<<tmpxxtime<<"  "<<tmpyytime<<endl;

#elif defined(BARC_EVTBLR)
		  int tmpyytime = -100;

		  if(event->ytdc_l[jk][itdc]->GetEntries()) {
		    TDCInfo *ytdc = (TDCInfo*)event->ytdc_l[jk][itdc]->At(0);
		    tmpyytime = 6000.+int(ytdc->GetValue()-event->tdc_ref_l[jk]);
		    tmpyytime1 = 12500.+int(ytdc->GetValue()-event->tdc_ref_l[jk]);
		  }

#else
		  int tmpyytime = event->tdcdata[jkrd];
#endif
#if defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
		  initytime[jk][itdc]=(tmpyytime<50.0) ? -100.0 : 0.1*tmpyytime;
#else
		  initytime[jk][itdc]=(tmpyytime<50.0) ? -100.0 : 0.1*tmpyytime;
#endif
		  //		if (jk==4) timesy[jk] -=25.0;
		  if(filloccu) {
		    time_layer->Fill(jkfl, tmpyytime);
		    if (yallhits[jk][itdc]>=0 && yallhits[jk][itdc]<nstrip) {
		      time_yraw[jk][yallhits[jk][itdc]]->Fill(tmpyytime);
		      time_layerstrip->Fill(nstrip*(jk+nlayer)+yallhits[jk][itdc], tmpyytime);
		      time_layerstripcpy->Fill(nstrip*(jk+nlayer)+yallhits[jk][itdc], tmpyytime1);
#ifdef TIMESLOPE
		      if (initytime[jk][itdc]>0 && abs(xslope)<1.0&&abs(yslope)<1.0) {
			time_Y_corr[jk][itdc]->Fill(yextloc[jk], initytime[jk][itdc]);
			time_yslope[jk][itdc]->Fill(xextloc[jk], initytime[jk][itdc]);
		      }
#endif
		    }
		  }
		   } // for (int itdc=0; itdc<nTDCpLayer; itdc++)
#ifdef NEWRPCDAQ1
		   if (filloccu) {
		     /*
		     h_xtstrpmult_xtdchitmult[jk]->Fill(xptsall[jk].size(),nxtdchitmult);
		     h_ytstrpmult_ytdchitmult[jk]->Fill(yptsall[jk].size(),nytdchitmult);
		     h_xtstrpmult_ytdchitmult[jk]->Fill(xptsall[jk].size(),nytdchitmult);
		     h_ytstrpmult_xtdchitmult[jk]->Fill(yptsall[jk].size(),nxtdchitmult);
		     */
		     if( xptsall[jk].size()==1 &&  yptsall[jk].size()==1) {
		       if(initxtime[jk][xptsall[jk][0]%8]>istartx[jk] && initxtime[jk][xptsall[jk][0]%8]<iendx[jk] ) {
			 xlayer_timeoccu[jk]->Fill(xptsall[jk][0],yptsall[jk][0]);
		       }
		        if(initytime[jk][yptsall[jk][0]%8]>istarty[jk] && initytime[jk][yptsall[jk][0]%8]<iendy[jk]) {
			 ylayer_timeoccu[jk]->Fill(xptsall[jk][0],yptsall[jk][0]);
		       }
		       if(initxtime[jk][xptsall[jk][0]%8]>istartx[jk] && initxtime[jk][xptsall[jk][0]%8]<iendx[jk] &&
			  initytime[jk][yptsall[jk][0]%8]>istarty[jk] && initytime[jk][yptsall[jk][0]%8]<iendy[jk]) {
			 xylayer_timeoccu[jk]->Fill(xptsall[jk][0],yptsall[jk][0]);
		       }

		     }
		     if(event->tdc_ref_l[jk]<=0 && jk ==trigly1) {ntxytref[0]++;}
		     if(event->tdc_ref_l[jk]<=0 && jk ==trigly2) {ntxytref[1]++;}
		     if(event->tdc_ref_l[jk]<=0 && jk ==trigly3) {ntxytref[2]++;}
		     if(event->tdc_ref_l[jk]<=0 && jk ==trigly4) {ntxytref[3]++;}
		     int xtdccount =0; int ytdccount=0;

		     for(int tdc = 0;tdc<8;tdc++){
		       if(event->vxtdc_l[jk][tdc]->size()>0) { xtdccount++;}
		       if(event->vytdc_l[jk][tdc]->size()>0) { ytdccount++;}
		     }
		     if(xtdccount==0 && ytdccount==0 && jk==trigly1){ ntxytime2[0]++;}
		     if(xtdccount==0 && ytdccount==0 && jk==trigly2){ ntxytime2[1]++;}
		     if(xtdccount==0 && ytdccount==0 && jk==trigly3){ ntxytime2[2]++;}
		     if(xtdccount==0 && ytdccount==0 && jk==trigly4){ ntxytime2[3]++;}
		     if(event->tdc_ref_l[jk]>0) {
		       ntxyttref[jk]++;
		       if(xtdccount>0) { ntxreftdc[jk]++;}
		       if(ytdccount>0) { ntyreftdc[jk]++;}

		     }

		   }
#endif

		 } //  for(int jk=0;jk<nlayer;jk++)

#endif // ifndef MONTECARLO




		 //	      cout<<"timing "<<iev<<" "<<ntiming++<<endl;
		 //	      ntiming++;

		 for (int ij=0; ij<nlayer; ij++) {
		dist[ij]= -100.; istrxtime[ij] = istrytime[ij] = -1;
		passxtime[ij] = passytime[ij] = false;
		passxtmy[ij] = passytmx[ij] = true;
#ifdef MONTECARLO
		rawxtime[ij] = rawxtime1[ij] = rawxtime2[ij] = rawxtime3[ij] = timesx[ij] = xtime[ij] ;
		rawytime[ij] = rawytime1[ij] = rawytime2[ij] = rawytime3[ij] = timesy[ij] = ytime[ij] ;
		xusedtime[ij] = yusedtime[ij] = true;
#else
		rawxtime[ij] = rawxtime0[ij] = rawxtime1[ij] = rawxtime2[ij] = rawxtime3[ij] = timesx[ij] = xtime[ij] =-100;
		rawytime[ij] = rawytime0[ij] = rawytime1[ij] = rawytime2[ij] = rawytime3[ij] = timesy[ij] = ytime[ij] =-100;
		xusedtime[ij] = yusedtime[ij] = false;
#endif
		xtdev[ij] = 100;
		ytdev[ij] = 100;
	      }


#ifdef CAUDATA

	      while(event->Evetime->GetSec() >= CauData[cauentry].CauTime){
		if(cauentry < cau_tree->GetEntries()-2)cauentry = cauentry+2 ;
	      }

	      	for(int kl =0;kl<nlayer;kl++) {
		  if(cauentry==0) {
		      if(kl==1 || kl==4 || kl==5 || kl==6 || kl==7 ) {cau_calib1[kl] = cau_calib2[kl]= CauData[cauentry+1].CauVal[kl];}
		      if(kl==8 ||kl==9) {cau_calib1[kl] = cau_calib2[kl]= CauData[cauentry].CauVal[kl];}
		  }else {
		if(kl==1 || kl==4 || kl==5 || kl==6 || kl==7 ) {
		  cau_calib2[kl] = CauData[cauentry+1].CauVal[kl];
		  cau_calib1[kl] = CauData[cauentry-1].CauVal[kl];
		}
		if(kl==8 ||kl==9) {
		  cau_calib2[kl]  = CauData[cauentry].CauVal[kl];
		  cau_calib1[kl]  = CauData[cauentry-2].CauVal[kl];
		}
		  }
		  //  cout<<"CAU Values"<<"    "<<kl<<"    "<<CauData[cauentry].CENum<<"  "<<cau_calib1[kl]<<"   "<<cau_calib2[kl]<<endl;

		}

	      /*	for(int kl =0;kl<nlayer;kl++) {
		  if(cauentry==1) {
		      if(kl==1 || kl==4 || kl==5 || kl==6 || kl==7 ) {cau_calib1[kl] = cau_calib2[kl]= CauData[cauentry-1].CauVal[kl];}
		      if(kl==8 ||kl==9) {cau_calib1[kl] = cau_calib2[kl]= CauData[cauentry].CauVal[kl];}
		  }else {
		if(kl==1 || kl==4 || kl==5 || kl==6 || kl==7 ) {
		  cau_calib1[kl] = CauData[cauentry-3].CauVal[kl];
		  cau_calib2[kl] = CauData[cauentry-1].CauVal[kl];
		}
		if(kl==8 ||kl==9) {
		  cau_calib1[kl]  = CauData[cauentry-2].CauVal[kl];
		  cau_calib2[kl]  = CauData[cauentry].CauVal[kl];
		}
		  }
		  cout<<"CAU Values"<<"    "<<kl<<"    "<<CauData[cauentry].CENum<<"  "<<cau_calib1[kl]<<"   "<<cau_calib2[kl]<<endl;

		}
	      */

#elif defined(CAUDATA1)
		int tmpncalib=0;
		while(event->Evetime->GetSec() >= CauData[cauentry].CauTime){
		  if(cauentry < cau_tree->GetEntries()-2)cauentry = cauentry+2 ;
		  tmpncalib++;
		  if (tmpncalib>10000) {
		    cout<<"Wrong CALUNIT time "<<iev<<" "<<event->Evetime->GetSec()<<" "<<CauData[cauentry].CauTime<<endl;
		    break;

		  }
		}

	      	for(int kl =0;kl<nlayer;kl++) {
		  if(cauentry==1) {
		    if(kl==0 || kl==1 || kl==2 || kl==3 || kl==4 || kl==5 || kl==6|| kl==7) {cau_calib1[kl] = cau_calib2[kl]= CauData[cauentry-1].CauVal[kl];}
		    if(kl==8 ||kl==9 || kl==10 || kl==11) {cau_calib1[kl] = cau_calib2[kl]= CauData[cauentry].CauVal[kl];}
		  }else {
		    if(kl==1 || kl==2 || kl==3 || kl==4 || kl==5 || kl==6 || kl==7 ) {
		      cau_calib1[kl] = CauData[cauentry-3].CauVal[kl];
		      cau_calib2[kl] = CauData[cauentry-1].CauVal[kl];
		    }
		    if(kl==8 ||kl==9 || kl==10 || kl==11) {
		      cau_calib2[kl]  = CauData[cauentry].CauVal[kl];
		      cau_calib1[kl]  = CauData[cauentry-2].CauVal[kl];
		    }
		  }
		  //  cout<<"CAU Values"<<"    "<<kl<<"    "<<CauData[cauentry].CENum<<"  "<<cau_calib1[kl]<<"   "<<cau_calib2[kl]<<endl;

		}
		//	cout<<iev<<"  "<<"PASSSSSSSSSSS4"<<endl;
#endif
	      xval=-100; yval=-100;
              int init=-1;
	      double initzpos=0;
	      bool xbandcond[nlayer][nband];
	      bool ybandcond[nlayer][nband];
	      bool xbandzone[nlayer][nband+1];
	      bool ybandzone[nlayer][nband+1];
	      bool isButton = false;// button jim
	      if(ntcor==1) {
		for(int ix=0;ix<nxybutton;ix++) {
		  for(int iy=0;iy<nxybutton;iy++) {
		    if(abs(xextloc[occulyr]-button_xpos[occulyr][ix])<1. && abs(yextloc[occulyr]-button_ypos[occulyr][iy])<1. ) { isButton = true; break;}
		  }
		  if(isButton) { break;}
		}
	    }
	      for( int ij=0;ij<nlayer-2;ij++) {
		if(isButton) { break; } //button jim
		//		if(abs(Xdev[ij]) < 2.0 && abs(Ydev[ij]) < 2.0) {
		//		//		if(abs(Xdev[ij]) < xyPosDev && abs(Ydev[ij]) < xyPosDev) {
		//		  if (xpts[ij].size()==0 || xpts[ij].size() >nmxhits || ypts[ij].size()==0 || ypts[ij].size() >nmxhits) continue;
		for(int kk=0;kk<nband;kk++) { xbandcond[ij][kk] = ybandcond[ij][kk] = false;}
		for(int kk=0;kk<nband+1;kk++) { xbandzone[ij][kk] = ybandzone[ij][kk] = false;}
		if (init<0) {
		  xval = xextloc[ij];
		  yval = yextloc[ij];
		  dist[ij] = 0.0;
		  init = ij;
		  initzpos = layerzpos[ij];
		} else {
		  dist[ij] = sqrt( pow((xextloc[ij] - xval)*stripwidth, 2.) +
				   pow((yextloc[ij] - yval)*stripwidth, 2.) +
				   pow(layerzpos[ij] - initzpos, 2.));
		  //                                     pow((ij - init)*layergap, 2.));
		}

		bool xfidul = (xextloc[ij]+0.5 >-1 && xextloc[ij]+0.5<nstrip+1) ? 1 : 0; //GMA
		bool yfidul = (yextloc[ij]+0.5 >-1 && yextloc[ij]+0.5<nstrip+1) ? 1 : 0; //GMA
		bool isButtonall = false; //button jim
		for(int ix=0;ix<nxybutton;ix++) {
		  for(int iy=0;iy<nxybutton;iy++) {
		    if(abs(xextloc[ij]-button_xpos[ij][ix])<1. && abs(yextloc[ij]-button_ypos[ij][iy])<1. ) { isButtonall = true; break;}
		  }
		  if(isButtonall) { break; }
		}

		//if(ij==occulyr){isButtonall=false;}  // jim jim added to remove the layer under study to include also button

		if(!isButtonall && abs(Xdev[ij]) < 2.0 && abs(Ydev[ij]) < 2.0 &&
		   xpts[ij].size()>0 && xpts[ij].size() <=nmxhits && ypts[ij].size()>0 && ypts[ij].size() <=nmxhits) { //button jim here !isButtonall

		  //calcualte strips where signal come earlier and the value of offset
		  double tshft=1000.0;
		  bool tpass=true;
		  passxtime[ij] = true;
		  passxtmy[ij] = true; // (yextloc[ij] > firstYstrip && yextloc[ij] < lastYstrip) ? true : false;

		  int istr = istrxtime[ij] = xpts[ij][0];
		  //		  cout <<"istr = "<<istrxtime[ij]<<" "<<xpts[ij][0]<<endl;
#ifndef MONTECARLO
		  for (int ix=0; ix<xpts[ij].size(); ix++) {
		    if ((!timeinuse[0][ij][xpts[ij][ix]]) ||
			xpts[ij][ix]< firstXstrip ||
			xpts[ij][ix]> lastXstrip) {
		      tpass= passxtime[ij] = false;
		    }

		    //		    if (xtoffset[ij][xpts[ij][ix]]<tshft) {
		    //		      tshft =xtoffset[ij][xpts[ij][ix]]; istr = istrxtime[ij] = xpts[ij][ix];
		    //		    }
		  }
		  //  cout<<"1 Passsssssssssssss"<<endl;
#ifdef ONETIMEFORMAT
		  timesx[ij] = initxtime[ij][0];
#elif defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
		  //cau_tree->GetEntries();

		  double stripdelay[8] = {0.4,0.54,0.82,1.09,1.23,1.5,1.77,2.04};
		  int stripno[8] = {0,8,16,24,32,40,48,56};

		  // timesx[ij] = -100.;
		  // if(abs(Xdev[ij])<1.5) {

		    int ixxx[nmxtimehit];
		    for (int xix=0; xix<min(nmxtimehit, int(xpts[ij].size())); xix++) {
		      ixxx[xix] = xpts[ij][xix]%8;
#ifdef SHIFT0
		      if (ij==3 && xpts[ij][xix]>=32 && xpts[ij][xix]<40) {
			ixxx[xix]++;
			if (ixxx[xix]<0 || ixxx[xix]>=8) { ixxx[xix]=0;}
		      }
#endif
		    }
		    if (filloccu && xpts[ij].size()==1) {
		      ntotxext[ij][ixxx[0]]++;
#ifdef NEWRPCDAQ1
		      if(event->vxtdc_l[ij][ixxx[0]]->size()) {
#else
			if(event->xtdc_l[ij][ixxx[0]]->GetEntries()) {
#endif
			  ntxtime[ij][ixxx[0]]++;
			}
		      }

		      double inittime = 11000;
		      double initwidth = 0;
		      double xpulsewidth= 0.;
		      for (int xix=0; xix<min(nmxtimehit, int(xpts[ij].size())); xix++) {
#ifdef NEWRPCDAQ1
			if (event->vxtdc_l[ij][ixxx[xix]]->size()) {
			  double tmptime;
			  //if(ij==3 || ij==5) { //jim jim
			  //tmptime = 1250.0 + 0.1*(int(event->vxtdc_l[ij][ixxx[xix]]->at(0)-event->tdc_ref_l[ij]))-44.; }
			  //else {
			    tmptime = 1250.0 + 0.1*(int(event->vxtdc_l[ij][ixxx[xix]]->at(0)-event->tdc_ref_l[ij]));
			    //}
//- (5.2*stripdelay[int(xpts[ij][xix]/8)]); //+0.1*(cau_calib1[ij]+cau_calib2[ij])/4. ;

			  double tmpwidth=0;
			  if(event->vxtdc_t[ij][ixxx[xix]]->size()) {
			    tmpwidth = 0.1*int(event->vxtdc_t[ij][ixxx[xix]]->at(0)) -  0.1*int(event->vxtdc_l[ij][ixxx[xix]]->at(0));
			  }
#else
			//	cout<<"XXXXXXXXXXXX"<<endl;
			if (event->xtdc_l[ij][ixxx[xix]]->GetEntries()) {
			  //  cout<<ij<<"      "<<"XXXXXXXXXXXX"<<endl;
			  TDCInfo *xxtdc = (TDCInfo*)event->xtdc_l[ij][ixxx[xix]]->At(0);
			  double tmptime = 400.0 + 0.1*(int(xxtdc->GetValue()-event->tdc_ref_l[ij]));
#endif


			  if(tmptime < inittime) {
			    initwidth = xpulsewidth = tmpwidth; inittime = tmptime; tshft =xtoffset[ij][xpts[ij][xix]]; istr = istrxtime[ij] = xpts[ij][xix];
			}
		      }
		    }

		    if(inittime <10000) {
#ifdef CAUDATA1
		      timesx[ij] = inittime /*- 5.2*stripdelay[int(istr/8)] */+0.1*(cau_calib1[ij]+cau_calib2[ij])/4.;
#else
		      timesx[ij] = inittime;// - 5.2*stripdelay[int(istr/8)];
		      if(ntcor==1 && iiter==nmxiter-1 && ij==occulyr) { xtpulsewidth = xpulsewidth;}
#endif
		    }
		    // }

#else
		  timesx[ij] = (xextloc[ij] < nstrip/2-0.5) ? initxtime[ij][0] : initxtime[ij][1];
#endif
		  //  cout<<"2 Passsssssssssssss"<<endl;
		  // cout<<"timesX"<<"    "<<timesx[ij]<<endl;
		  //cout<<iev<<" X  "<<ij<<"  "<<timesx[ij]<<endl;
		  rawxtime[ij] = timesx[ij];

		     // cout<<"timesY"<<"    "<<timesy[ij]<<endl;
		  timesx[ij] -=timeoffsetx[ij];

		  rawxtime0[ij] = timesx[ij];
		  if(ntcor==1 && iiter == nmxiter-1 && ij==occulyr) { xtimecpy[ij][0]->Fill(rawxtime0[ij]);}
		  timesx[ij] -=slope_path*yext[ij]; //(5./32.)*yext[ij];

		  rawxtime1[ij] = timesx[ij];
		  // if(ij!=occulyr) { timesx[ij] -=tshft; }
		   timesx[ij] -=tshft;
		  rawxtime2[ij] = timesx[ij];
		  if(ntcor==1 && iiter == nmxiter-1 && ij==occulyr) { xtimecpy[ij][1]->Fill(rawxtime2[ij]);}
		  istr = int(yextloc[ij]+0.5);
		  if (istr<0) istr=0;
		  if (istr>=nstrip) istr = nstrip-1;
		  double dx = yextloc[ij]-istr;

		  // Linear extrapolation using only two points
		  // if(ij!=occulyr) {
		  if ((istr==0 && dx<=0.0) || (istr==nstrip-1 && dx>=0.0)) {
 		    timesx[ij] -=xt_slope_cor[ij][istrxtime[ij]][istr];
 		  } else if (dx>0) {
 		    timesx[ij] -=(1-dx)*xt_slope_cor[ij][istrxtime[ij]][istr]+dx*xt_slope_cor[ij][istrxtime[ij]][istr+1];
 		  } else {
 		    timesx[ij] -=abs(dx)*xt_slope_cor[ij][istrxtime[ij]][istr-1]+(1-abs(dx))*xt_slope_cor[ij][istrxtime[ij]][istr];
 		  }
		  //}
		  //  rawxtime[ij] =  rawxtime0[ij] =  rawxtime1[ij] = rawxtime2[ij] = rawxtime3[ij] = timesx[ij];





		  if (xfidul && (xpts[ij].size() >0 && xpts[ij].size()<=nmxtimehit) && (usealltime || (ij==occulyr) || (tpass && passxtmy[ij]))) { //GMA
		    //  if(ij!=occulyr) {

		    int isiz = max(0, min(nmxtimehit,int(xpts[ij].size()))-1);
		    xtime[ij] = timesx[ij];// - parabolic(xposinstr[ij], strpos_vs_time[isiz][ij]); //jim jim enable parabolic
		    //  }else {
		    // xtime[ij] = timesx[ij];
		      //  }
		    //if(ij!=occulyr) {
		       //if(xpulsewidth>23) { xpulsewidth=23;}
		       ///*This one enable*/if(xextloc[ij]>0 && xextloc[ij]<nstrip && (ystr_vs_xtot_corr[ij]->GetBinContent(ystr_vs_xtot_corr[ij]->FindBin(yextloc[ij],xpulsewidth))) > 1.) { xtime[ij] -=  (ystr_vs_xtot_corr[ij]->GetBinContent(ystr_vs_xtot_corr[ij]->FindBin(yextloc[ij],xpulsewidth))-15.);} //jim jim enable for tot pixelwise combined
		      //if(yextloc[ij]>=0 && xpulsewidth>=0 && yextloc[ij]<60  && xextloc[ij]>=0 && xextloc[ij]<60) {
		      //if(xextloc[ij]>0 && xextloc[ij]<nstrip && (ystr_vs_xtot_corr[ij][int(xextloc[ij]/20)]->GetBinContent(ystr_vs_xtot_corr[ij][int(xextloc[ij]/20)]->FindBin(yextloc[ij],xpulsewidth))) > 1.) { xtime[ij] -=  (ystr_vs_xtot_corr[ij][int(xextloc[ij]/20)]->GetBinContent(ystr_vs_xtot_corr[ij][int(xextloc[ij]/20)]->FindBin(yextloc[ij],xpulsewidth))-15.);} //jim jim enable for tot pixelwise 0-20, 20-40, 40-60 tot
		      //}
			//if(/*xpulsewidth<23 && */xextloc[ij]>0 && xextloc[ij]<nstrip && (ystr_vs_xtot_corr[ij][int(xextloc[ij]/20)]->GetBinContent(ystr_vs_xtot_corr[ij][int(xextloc[ij]/20)]->FindBin(xpulsewidth))) > 1.) { xtime[ij] -=  (ystr_vs_xtot_corr[ij][int(xextloc[ij]/20)]->GetBinContent(ystr_vs_xtot_corr[ij][int(xextloc[ij]/20)]->FindBin(xpulsewidth))-15.);} //jim jim enable for tot
		    //}
		    rawxtime3[ij] = xtime[ij] + slope_path*yext[ij]; //undo the pathlenght correction
		    xusedtime[ij] = (xpts[ij].size()<=nmxusedtimehit) ? true : false; //  xtime[ij] : -101;
		     if(ntcor==1 && iiter == nmxiter-1 && ij==occulyr) { xtimecpy[ij][2]->Fill(xtime[ij]);}
		  } // else {rawxtime[ij] = rawxtime0[ij] = rawxtime1[ij] = rawxtime2[ij] = rawxtime3[ij] =-90;}

#endif

		  tshft=1000.0;
		  tpass=true;
		  passytime[ij] = true;
		  passytmx[ij] = true; //(xextloc[ij] > firstXstrip && xextloc[ij] < lastXstrip) ? true : false;
		  istr = istrytime[ij] = ypts[ij][0];
#ifndef MONTECARLO
		  for (int iy=0; iy<ypts[ij].size(); iy++) {
		    if ((!timeinuse[1][ij][ypts[ij][iy]]) ||
			ypts[ij][iy]< firstYstrip ||
			ypts[ij][iy]> lastYstrip) {
		      tpass = passytime[ij] = false;
		    }
		  }

#ifdef ONETIMEFORMAT
		  timesy[ij] = initytime[ij][0];
#elif defined(NEWRPCDAQ1) || defined(BARC_EVTBLR)
		  // timesy[ij] = -100.;
		  // if(abs(Ydev[ij])<1.5) {

		    int iyyy[nmxtimehit];


		    inittime = 11000;
		   double ypulsewidth = 0.;
                   for (int yiy=0; yiy<min(nmxtimehit, int(ypts[ij].size())); yiy++) {
                     iyyy[yiy] = ypts[ij][yiy]%8;
#ifdef NEWRPCDAQ1
                     if (event->vytdc_l[ij][iyyy[yiy]]->size()) {

		       double tmptime;
		       //if(ij==3 || ij==5) { //jim jim
		       //tmptime = 1250.0 + 0.1*(int(event->vytdc_l[ij][iyyy[yiy]]->at(0)-event->tdc_ref_l[ij]))-44.;}
		       //else  {
			 tmptime = 1250.0 + 0.1*(int(event->vytdc_l[ij][iyyy[yiy]]->at(0)-event->tdc_ref_l[ij]));
			 //}//- (5.2*stripdelay[int(ypts[ij][yiy]/8)]); //+0.1*(cau_calib1[ij]+cau_calib2[ij])/4. ;
			double tmpwidth;
			if(event->vytdc_t[ij][iyyy[yiy]]->size()) {
			    tmpwidth = 0.1*int(event->vytdc_t[ij][iyyy[yiy]]->at(0)) -  0.1*int(event->vytdc_l[ij][iyyy[yiy]]->at(0));

			  }
#else
		       if (event->ytdc_l[ij][iyyy[yiy]]->GetEntries()) {
			  TDCInfo *yytdc = (TDCInfo*)event->ytdc_l[ij][iyyy[yiy]]->At(0);
			  double tmptime = 300.0 + 0.1*(int(yytdc->GetValue()-event->tdc_ref_l[ij]));
#endif

                       if (tmptime < inittime) {
                         initwidth = ypulsewidth = tmpwidth; inittime = tmptime; tshft =ytoffset[ij][ypts[ij][yiy]]; istr = istrytime[ij] = ypts[ij][yiy];

		       }
		       }
		     }
		   if (inittime <10000) {
#ifdef CAUDATA1
		     timesy[ij] = inittime /*- (5.2*stripdelay[int(istr/8)])*/ +0.1*(cau_calib1[ij]+cau_calib2[ij])/4.;
#else
		     timesy[ij] = inittime ;//- (5.2*stripdelay[int(istr/8)]);
		     if(ntcor==1 && iiter==nmxiter-1 && ij==occulyr) { ytpulsewidth = ypulsewidth;}
#endif
		   }


		    if (filloccu && ypts[ij].size()==1) {
		      ntotyext[ij][iyyy[0]]++;
#ifdef NEWRPCDAQ1
		      if(event->vytdc_l[ij][iyyy[0]]->size()) {
#else
			if(event->ytdc_l[ij][iyyy[0]]->GetEntries()) {
#endif
			  ntytime[ij][iyyy[0]]++;
		      }
		    }
		      // }
#else
		  timesy[ij] = (yextloc[ij] < nstrip/2-0.5) ? initytime[ij][0] : initytime[ij][1];
#endif

		  //  cout<<"timesY"<<"    "<<timesy[ij]<<endl;
		  rawytime[ij] = timesy[ij];

		     //  cout<<iev<<" Y  "<<ij<<"  "<<timesy[ij]<<endl;
		  timesy[ij] -=timeoffsety[ij];

		  rawytime0[ij] = timesy[ij];
		  if(ntcor==1 && iiter == nmxiter-1 && ij==occulyr) {
		    ytimecpy[ij][0]->Fill(rawytime0[ij]);
		    xytimecpy[ij][0]->Fill(rawxtime0[ij]-rawytime0[ij]);
		    for(int ii=0;ii<nTermResSplit;ii++) {
		      if(int(xextloc[ij])/20 ==ii && int(yextloc[ij])/20 ==ii) {
		      xytimecpy_nTermSplit[ij][ii][0]->Fill(rawxtime0[ij]-rawytime0[ij]); }
		    }
		  }
		  timesy[ij] -=ytimeshift; // shift y-time to match with x-time
#ifdef C217STRIP
		  timesy[ij] -= 5. - slope_path*xext[ij];
#else
		  timesy[ij] -= slope_path*xext[ij]; // Here both X/Y will follow same convention 5. - slope_path*xext[ij]; //(5./32.)*xext[ij];
#endif
		  rawytime1[ij] = timesy[ij];

		  // if(ij!=occulyr) { timesy[ij] -=tshft; }
		  timesy[ij] -=tshft;
		  rawytime2[ij] = timesy[ij];
		  if(ntcor==1 && iiter == nmxiter-1 && ij==occulyr) {
		    ytimecpy[ij][1]->Fill(rawytime2[ij]);
		    xytimecpy[ij][1]->Fill(rawxtime2[ij]-rawytime2[ij]);
		    for(int ii=0;ii<nTermResSplit;ii++) {

		    if(int(xextloc[ij])/20 ==ii && int(yextloc[ij])/20 ==ii) {
		      xytimecpy_nTermSplit[ij][ii][1]->Fill(rawxtime2[ij]-rawytime2[ij]); }
		    }
		  }
		  istr = int(xextloc[ij]+0.5);
		  if (istr<0) istr=0;
		  if (istr>=nstrip) istr = nstrip-1;

		  double dy = xextloc[ij]-istr;

		  // Linear extrapolation using only two points
		  // if(ij!=occulyr) {
		  if ((istr==0 && dy<=0.0) || (istr==nstrip-1 && dy>=0.0)) {
 		    timesy[ij] -=yt_slope_cor[ij][istrytime[ij]][istr];
 		  } else if (dy>0) {
 		    timesy[ij] -=(1-dy)*yt_slope_cor[ij][istrytime[ij]][istr]+dy*yt_slope_cor[ij][istrytime[ij]][istr+1];
 		  } else {
 		    timesy[ij] -=abs(dy)*yt_slope_cor[ij][istrytime[ij]][istr-1]+(1-abs(dy))*yt_slope_cor[ij][istrytime[ij]][istr];
 		  }
		  //	  }


		  if (yfidul && (ypts[ij].size() >0 && ypts[ij].size()<=nmxtimehit) && (usealltime || (ij==occulyr) || (tpass && passytmx[ij]))) { //GMA
		    //  if(ij!=occulyr){
		    int isiz = max(0, min(nmxtimehit,int(ypts[ij].size()))-1)+nmxtimehit;
		    ytime[ij] = timesy[ij];// - parabolic(yposinstr[ij], strpos_vs_time[isiz][ij]); //jim jim enable parabolic
		    //		    cout<<"yyyy"<<"   "<<isiz<<"   "<<ij<<"   "<<ytime[ij]<<"   "<<timesy[ij]<<"   "<<parabolic(yposinstr[ij], strpos_vs_time[isiz][ij])<<"   "<<strpos_vs_time[isiz][ij][0]<<"   "<<strpos_vs_time[isiz][ij][1]<<"   "<<strpos_vs_time[isiz][ij][2]<<endl;
		    // } else {
		    //  ytime[ij] = timesy[ij];
		      //}
		    //if(ij!=occulyr) {
			//if(ypulsewidth>23) {ypulsewidth=23;}
		     ///*Jim for tot correction*/ if(yextloc[ij]>0 && yextloc[ij]<nstrip && (xstr_vs_ytot_corr[ij]->GetBinContent(xstr_vs_ytot_corr[ij]->FindBin(xextloc[ij],ypulsewidth)))>1.) { ytime[ij] -=  (xstr_vs_ytot_corr[ij]->GetBinContent(xstr_vs_ytot_corr[ij]->FindBin(xextloc[ij],ypulsewidth))-15.);} //jim jim enable this for tot pixelwise combined
		    //if(xextloc[ij]>=0 && ypulsewidth>=0 && xextloc[ij] <60 && yextloc[ij]>=0 && yextloc[ij]<60) {
		      // if(yextloc[ij]>0 && yextloc[ij]<nstrip && (xstr_vs_ytot_corr[ij][int(yextloc[ij]/20)]->GetBinContent(xstr_vs_ytot_corr[ij][int(yextloc[ij]/20)]->FindBin(xextloc[ij],ypulsewidth)))>1.) { ytime[ij] -=  (xstr_vs_ytot_corr[ij][int(yextloc[ij]/20)]->GetBinContent(xstr_vs_ytot_corr[ij][int(yextloc[ij]/20)]->FindBin(xextloc[ij],ypulsewidth))-15.);} //jim jim enable this for tot pixelwise 0-20, 20-40, 40-60
		      //}

			//if(/*ypulsewidth<23 &&*/ yextloc[ij]>0 && yextloc[ij]<nstrip && (xstr_vs_ytot_corr[ij][int(yextloc[ij]/20)]->GetBinContent(xstr_vs_ytot_corr[ij][int(yextloc[ij]/20)]->FindBin(ypulsewidth)))>1.) { ytime[ij] -=  (xstr_vs_ytot_corr[ij][int(yextloc[ij]/20)]->GetBinContent(xstr_vs_ytot_corr[ij][int(yextloc[ij]/20)]->FindBin(ypulsewidth))-15.); } //jim jim enable this for tot
		      //}
#ifdef C217STRIP
		    rawytime3[ij] = ytime[ij] + 5. - slope_path*xext[ij];
#else
		    rawytime3[ij] = ytime[ij] + slope_path*xext[ij];
#endif
		    yusedtime[ij] = (ypts[ij].size()<=nmxusedtimehit) ?  true : false; //ytime[ij] : -101;
		    // if (ypts[ij].size()>nmxusedtimehit && yusedtime[ij]) cout <<"ypts "<<ij<<" "<< ypts[ij].size()<<" "<<nmxusedtimehit<<" "<<int(yusedtime[ij])<<endl;
		    if(ntcor==1 && iiter == nmxiter-1 && ij==occulyr) {
		    ytimecpy[ij][2]->Fill(ytime[ij]);
		    xytimecpy[ij][2]->Fill(xtime[ij]-ytime[ij]);
		    for(int ii=0;ii<nTermResSplit;ii++) {
		      if(int(xextloc[ij])/20 ==ii && int(yextloc[ij])/20 ==ii) {
			xytimecpy_nTermSplit[ij][ii][2]->Fill(xtime[ij]-ytime[ij]); }
		    }
		    }
		  } // else { rawytime[ij] = rawytime0[ij] = rawytime1[ij] = rawytime2[ij] = rawytime3[ij] =-90;}

#endif


		  if (filloccu) {
		    h_xrawcorhits->Fill(-1.0, -1.0);

		    for (int ix=0; ix<nlayer; ix++) {
		      if (xptsall[ix].size()==0) continue;
		      for (int iy=ix+1; iy<nlayer; iy++) {
			if (yptsall[iy].size()==0) continue;
			h_xrawcorhits->Fill(ix, iy);
			if (rawxtime[ix]>-50 && rawxtime[iy]>-50) {
			  h_xtrawcorhits->Fill(ix, iy);
			}
		      }
		      for (int iy=0; iy<nlayer; iy++) {
			if (yptsall[iy].size()==0) continue;
			h_xyrawcorhits->Fill(ix, iy);
			if (rawxtime[ix]>-50 && rawytime[iy]>-50) {
			  h_xytrawcorhits->Fill(ix, iy);
			}
		      }
		    }

		    for (int ix=0; ix<nlayer-1; ix++) {
		      if (yptsall[ix].size()==0) continue;
		      for (int iy=ix+1; iy<nlayer; iy++) {
			if (yptsall[iy].size()==0) continue;
			h_yrawcorhits->Fill(ix, iy);
			if (rawytime[ix]>-50 && rawytime[iy]>-50) {
			  h_ytrawcorhits->Fill(ix, iy);
			}
		      }
		    }
		  }

//  cout<<"eve no "<<"   "<<iev<<"   "<<ij<<"  "<<xext[ij]<<"   "<<Xpos[ij]<<"   "<<rawxtime[ij]<<"   "<<xtime[ij] <<"  "<<yext[ij]<<"    "<<Ypos[ij]<<"   "<<rawytime[ij]<<"  "<<ytime[ij]<< endl;
		  int ifirst=0;//4;//1;//0; //was used one, while Layer-0 was not active
		  double dtime = sqrt (pow((xextloc[ij] - xextloc[ifirst])*stripwidth, 2.) +
				       pow((yextloc[ij] - yextloc[ifirst])*stripwidth, 2.) +
				       pow(layerzpos[ij] - layerzpos[ifirst], 2.0))/cval;
		  //				       pow(ij*layergap, 2.))/cval;

		  double dtime2 = sqrt (pow((xextloc[ij] - xextloc[nlayer-2])*stripwidth, 2.) +
					pow((yextloc[ij] - yextloc[nlayer-2])*stripwidth, 2.) +
					pow(layerzpos[ij] - layerzpos[nlayer-2], 2.0))/cval;
		  //		pow((nlayer-1-ij)*layergap, 2.))/cval;
		  //  cout<<"XXXXXXXXXXX"<<"   "<<ij<<"   "<<iiter<<"   "<<dtime<<"   "<<xusedtime[ij]<<"  "<<xusedtime[ifirst]<<"  "<<rawxtime[ij] - rawxtime[ifirst]<<endl;
		  if (firstiter==0 && ntcor==0) {
		    if (init==ifirst && xpts[ifirst].size()>0 && xpts[ifirst].size()<=nmxhits && xtime[ifirst]>-50 && xpts[ij].size()>0 && xpts[ij].size()<=nmxhits && xtime[ij]>-50) {

		      if ((usealltime || timeinuse[0][ij][istrxtime[ij]]) &&  xusedtime[ij] && xusedtime[ifirst]) {
			timex_shift[ij][iiter+1]->Fill(xtime[ij] - xtime[ifirst] + dtime);
			timex_2dshift[ij][iiter+1]->Fill(xextloc[ij], xtime[ij] - xtime[ifirst] + dtime);
			if (iiter==0) {
			  timex_shift[ij][0]->Fill(rawxtime[ij] - rawxtime[ifirst] + dtime);
			  // cout<<"timediff"<<"    "<<rawxtime[ij] - rawxtime[ifirst] + dtime<<endl;
			  timex_2dshift[ij][0]->Fill(xextloc[ij], rawxtime[ij] - rawxtime[ifirst] + dtime);
			}
		      }
		    }

		    if (init==ifirst && ypts[ifirst].size()>0 && ypts[ifirst].size()<=nmxhits && ytime[ifirst]>-50 && ypts[ij].size()>0 && ypts[ij].size()<=nmxhits && ytime[ij]>-50) {
		      if ((usealltime || timeinuse[1][ij][istrytime[ij]]) &&  yusedtime[ij] && yusedtime[ifirst]) {
			timey_shift[ij][iiter+1]->Fill(ytime[ij] - ytime[ifirst] + dtime);
			timey_2dshift[ij][iiter+1]->Fill(yextloc[ij], ytime[ij] - ytime[ifirst] + dtime);
			if (iiter==0) {
			  timey_shift[ij][0]->Fill(rawytime[ij] - rawytime[ifirst] + dtime);
			  timey_2dshift[ij][0]->Fill(yextloc[ij], rawytime[ij] - rawytime[ifirst] + dtime);
			}
		      }
		    }
		  }

		  if (isfill) {
		    if (xpts[ij].size()==1 && timesx[ij]>0 && passxtime[ij] && ypts[ij].size()>0) { //GMA
		      if (xpts[0].size() ==1 && xpts[0][0] >=13 && xpts[0][0] <=15 &&
			  ypts[0].size() ==1 && ypts[0][0] >=13 && ypts[0][0] <=15 ) {

			timex_fy[ij]->Fill(xpts[ij][0]*nstrip + ypts[ij][0], timesx[ij] - timesx[0] + dtime);
			timex_pfy[ij]->Fill(xpts[ij][0]*nstrip + ypts[ij][0], timesx[ij] - timesx[0] + dtime);
			timex_pfy[ij]->Fill(nstrip*nstrip + ypts[ij][0], timesx[ij] - timesx[0] + dtime);

		      }

		      if (xpts[nlayer-1].size() ==1 && xpts[nlayer-1][0] >=13 && xpts[nlayer-1][0] <=15 &&
			  ypts[nlayer-1].size() ==1 && ypts[nlayer-1][0] >=13 && ypts[nlayer-1][0] <=15 ) {

			timex_fy2[ij]->Fill(xpts[ij][0]*nstrip + ypts[ij][0], timesx[ij] - timesx[nlayer-1] - dtime2);
			timex_pfy2[ij]->Fill(xpts[ij][0]*nstrip + ypts[ij][0], timesx[ij] - timesx[nlayer-1] - dtime2);
			timex_pfy2[ij]->Fill(nstrip*nstrip + ypts[ij][0], timesx[ij] - timesx[nlayer-1] - dtime2);
		      }
		    } // if (xpts[ij].size()==1 && timesx[ij]>0 && passxtime[ij])

		    if (ypts[ij].size()==1 && timesy[ij]>0.0 && passytime[ij]  && xpts[ij].size()>0) { //GMA
		      if (xpts[0].size() ==1 && xpts[0][0] >=13 && xpts[0][0] <=15 &&
			  ypts[0].size() ==1 && ypts[0][0] >=13 && ypts[0][0] <=15 ) {

			timey_fx[ij]->Fill(ypts[ij][0]*nstrip + xpts[ij][0], timesy[ij] - timesy[0] + dtime);
			timey_pfx[ij]->Fill(ypts[ij][0]*nstrip + xpts[ij][0], timesy[ij] - timesy[0] + dtime);
			timey_pfx[ij]->Fill(nstrip*nstrip + xpts[ij][0], timesy[ij] - timesy[0] + dtime);
		      }

		      if (xpts[nlayer-1].size() ==1 && xpts[nlayer-1][0] >=13 && xpts[nlayer-1][0] <=15 &&
			  ypts[nlayer-1].size() ==1 && ypts[nlayer-1][0] >=13 && ypts[nlayer-1][0] <=15 ) {

			timey_fx2[ij]->Fill(ypts[ij][0]*nstrip + xpts[ij][0], timesy[ij] - timesy[nlayer-1] - dtime2);
			timey_pfx2[ij]->Fill(ypts[ij][0]*nstrip + xpts[ij][0], timesy[ij] - timesy[nlayer-1] - dtime2);
			timey_pfx2[ij]->Fill(nstrip*nstrip + xpts[ij][0], timesy[ij] - timesy[nlayer-1] - dtime2);
		      }
		    } // if (ypts[ij].size()==1 && timesy[ij]>0.0 && passytime[ij] )
		  } // if (isfill)
		} else { //if(abs(Xdev[ij]) < xyPosDev && abs(Ydev[ij]) < xyPosDev)
                  xusedtime[ij] = yusedtime[ij] = false;
		}
	      } //for(int ij=0;ij<nlayer;ij++)
	      //  cout<<"XXXXXXXXXXXXXXXXXXXXXXX"<<endl;
	      //  cout<<"\n";

	      int tmpxtent[nlayer];
	      int tmpytent[nlayer];

	      timexslope = 0;
              xc0inters = 0;
	      xt0chi2 = 1.e+6;
              nxtfail = 0;

	      timeyslope = 0;
              yc0inters = 0;
	      yt0chi2 = 1.e+6;
              nytfail = 0;

	      for (int ix=0; ix<nlayer; ix++) {xtext[ix]= xtexter[ix] = 1000;}

	      // Xtime fit
	      int iTimeSlopeConst = 1; //((isTimeCorOrReso) && ntcor==1) ? 0 : -1; //-1 : 0; //jim jim
	      //was used to check correlation	      int iTimeSlopeConst = (iiter<nmxiter-1 && ntcor==1) ? 0 : -1;
	      StraightLineFit xtimefit(iTimeSlopeConst, dist, xtime,  timeserrx2, xusedtime, occulyr, occulyr, xtcorstr, xtcorend, float(1000.0));
	      xtimefit.GetParameters(nxtfail, xc0inters, timexslope);
	      //	    xtimefit.GetError(errcst, errlin, errcov);
	      xtimefit.GetChisqure(nxtime, xt0chi2);
	      xtimefit.GetFitValues(xtext, xtdev, xtexter);
	      //	      cout<<"nxfail "<<iev<<" "<<iiterrs<<" "<< nxtfail<<" "<<xc0inters<<" "<<timexslope<<" "<<nxtime<<" "<< xt0chi2<<endl;
	      //	      for (int ixx=0; ixx<nlayer; ixx++) {
	      //		cout<<"ixx "<<setw(6)<<iev<<setw(6)<<" "<<ixx<<setw(6)<<" "<<dist[ixx]<<setw(6)<<" "<<xtime[ixx]<<setw(6)<<" "<<xtext[ixx]<<setw(6)<<" "<<xtdev[ixx]<<setw(6)<<" "<<timeserrx2[ixx]<<endl;
	      //	      }


	      //	tmpytdev[ij] = ytdev[ij];
	      if(ntcor==1 && iiter == nmxiter-1) {
	      tmpxtdev = xtdev[occulyr];
	      tmpnxtime = nxtime;
	      tmpxtime = xtime[occulyr];
	      tmpxtext = xtext[occulyr];
	      }

	      if (nxtfail==0 && isfill) {
		h_tchisqx->Fill(xt0chi2);
		if (nxtime-2>0) {
		  h_treduchisqx->Fill(xt0chi2/(nxtime-2));
		  double probx =TMath::Prob(xt0chi2, nxtime-2);
		  h_xtprob->Fill(probx);
		  xtprob_vs_tchisqx->Fill(probx, xt0chi2);
		  xtprob_vs_ndf->Fill(probx, nxtime);
		  tchisqx_vs_ndf->Fill(xt0chi2, nxtime);
		  xtprob_vs_treduchisqx->Fill(probx, xt0chi2/(nxtime-2));
		  int ibin = getbinid(nxtime, nprob, probs);
		  if (ibin>=0 && ibin<nprob) {h_xtnprob[ibin]->Fill(probx);}
		}
		h_txndf->Fill(nxtime);
	      }

	      if ((isfill || occulyr==nlayer) && nxtfail==0 && nxtime >nmnhits && xt0chi2/(nxtime-2)<mxtimechisq) {
		for(int ij=0;ij<nlayer;ij++){

		  if (dist[ij]<0 || xtime[ij] <-50) continue;
		  double dt1 = xtdev[ij]-biasinxtime[ij];
		  if (abs(dt1)>3*sigmarpc || (!xusedtime[ij])) continue;

		  if (isfill) {
		    for(int jk=0;jk<nlayer;jk++){

		      if (dist[jk]<0 || xtime[jk] <-50) continue;
		      //                                      if(jk != ij )continue;
		      double dt2 = xtdev[jk]-biasinxtime[jk];
		      if (abs(dt2) >3*sigmarpc || (!xusedtime[jk])) continue;
		      deltatcov[ij][jk] += dt1*dt2;
		      deltatCount[ij][jk] +=1 ;
		    }
		  }
		  //140923
		  if (occulyr==nlayer) {
		    if (passxtime[ij]) {
		      ntotxtime++; totxtime +=xtime[ij];
		      for(int jk=ij+1;jk<nlayer;jk++){
			if (dist[jk]<0 || xtime[jk] <-50 || abs(xtdev[jk]-biasinxtime[jk])>4*sigmarpc || (!xusedtime[jk]) || (!passxtime[jk])) continue;
			timex_correl[ij][jk][iiter+1]->Fill(xtime[jk]-xtime[ij] - (dist[ij]-dist[jk])/cval);
			if (iiter==0) { timex_correl[ij][jk][0]->Fill(rawxtime[jk]-rawxtime[ij] - (dist[ij]-dist[jk])/cval);}
		      } // for(jk=ij+1;jk<nlayer;jk++)
		    }
		  }
		}// for(int ij=0;ij<nlayer;ij++)
	      } //if ((isfill || occulyr==nlayer) && nxtfail==0 && nxtime >nmnhits && xt0chi2/(nxtime-2)<mxtimechisq)

	      //	      cout <<"nxtime "<<nxtime<<" "<<nxtfail<<" "<<occulyr<<endl;
              if (nxtime>=nmnhits/*-ntcor*/ && nxtfail==0/* && xt0chi2/(nxtime-2)<2*/) {//4Nov For resolution 8layers are considered //jim jim tighter criteria of xhi2/ndf for better alignment

		//Was it twice 150101
		//Or has memory overflow !!!!!!!!!
		if (xt0chi2/(nxtime-2)<mxtimechisq) {
		  if (filloccu) {
		    if (iev>0 && tmptimerate >=tmpoldmutimexrate+timebin) {
		      file_out <<"timexrate "<<datafile<<" "<<iev<<" "<<nallmutimexrate<<" "<<nmutimexrate<<endl;
		      //		      cout <<"timexrate "<<datafile<<" "<<iev<<" "<<nallmutimexrate<<" "<<nmutimexrate<<endl;
		      mutimexrate->Fill(nmutimexrate*(1.0/timebin));
		      mutimexratex->Fill(nallmutimexrate, nmutimexrate*(1.0/timebin));
		      nmutimexrate = 0;
		      nallmutimexrate++;
		      tmpoldmutimexrate = tmptimerate;
		    }
		    nmutimexrate++;
		  }

		  dir_cx[occulyr][iiter]->Fill(-cval*timexslope);
		  int ibin = getbinid(nxtime, nprob, probs);
		  if (ibin>=0 && ibin<nprob) {dir_cxlay[occulyr][iiter][ibin]->Fill(-cval*timexslope);
		    //if(ibin>4 && (-cval*timexslope >2.0)){
		    //cout<<"dir_bias"<<"	"<<iev<<"	"<<ibin<<"	"<<-cval*timexslope<<endl;}
}

		  dir_cx2[occulyr][iiter]->Fill(xtimefit.GetSlope2());
		  dir_c0x[occulyr][iiter]->Fill(xc0inters);
		}
		if (isfill) dir_cxchi->Fill(xt0chi2, -1./timexslope/cval);

		if (occulyr >=nlayer) {

		  for (int ij=0; ij<nlayer-1; ij++) { //review1
		    for(int jk=ij+1; jk<nlayer; jk++) {
		      //cval is in cm/ns (29.979)
		      //timeshift_xreso[ij]->Fill((xtime[ij]-dist[ij]*timexslope)-(xtime[ij+1]-dist[ij+1]*timexslope));   //version1
		      //timeshift_xreso[ij]->Fill((xtime[ij])-(xtime[ij+1]));                                             //version2
		      if (dist[ij]<0 || xtime[ij] <-50 || dist[jk]<0 || xtime[jk]<-50) continue;
		      if (passxtime[ij] && passxtmy[ij] && passxtime[jk] && passxtmy[jk]) {
			if (xusedtime[ij] && xusedtime[jk]) {
			  timeshift_xreso[ij][jk]->Fill((xtime[ij]+dist[ij]*cval)-(xtime[jk]+dist[jk]*cval));   //version3
			}
		      }                                                                                                    //These three condirions i put for version 4
		    }
		  }

		  for (int ij=0; ij<nlayer; ij++) {
		    //int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[ij]+0.5),int(yextloc[ij]+0.5)));//(nsplit/8)*int(xextloc[ij]/7.25)+int(yextloc[ij]/7.625);
		    //if(xbina>nsplit-1) xbina = nsplit-1;
		    if (dist[ij]<0 || xtime[ij] <-50) continue;
		    if (xusedtime[ij]) {
		      if (istrxtime[ij]>=0 && istrxtime[ij] <nstrip && passxtmy[ij]) {
			time_xstrreso[ij][istrxtime[ij]][iiterrs]->Fill(xtdev[ij]);
			tmptime_xreso[ij][iiterrs]->Fill(xtdev[ij]);
		      }
		      xtime_exterr[ij][iiterrs]->Fill(xtexter[ij]);
		      xtime_exterr_vs_ndf[ij][iiterrs]->Fill(xtexter[ij], nxtime);
		      xtime_exterr_nTermSplit[ij][int(xextloc[ij]/20)][iiterrs]->Fill(xtexter[ij]);
		    }
                    if (passxtime[ij] && passxtmy[ij]) {
		      if (xusedtime[ij]) {
#ifndef MONTECARLO
			if (filloccu) { time_xreso_set[ij][iset]->Fill(xtdev[ij]);}
#endif
			time_xreso[ij][iiterrs]->Fill(xtdev[ij]); //jim jim
			time_xreso_vs_ndf[ij][iiterrs]->Fill(xtdev[ij],nxtime); //jim jim
			time_xreso_nTermSplit[ij][int(xextloc[ij]/20)][iiterrs]->Fill(xtdev[ij]);
			if(isfid_l9_eff) { time_xreso_split[ij][xbina][iiterrs]->Fill(xtdev[ij]);}
			istime_xyreso[ij][iiterrs] = true;
		      }  // Filling and saving. Not fitted.
		      //		      int ixx = max(0,min(xhits[ij],nmxtimehit)-1);
		      int ixx = min(xhits[ij],nmxtimehit)-1;
		      //		      if (ixx>0) cout<<"ixx "<<ixx<<" "<<xtdev[ij]<<endl;
		      if (ixx>=0) {
			time_mulxreso[ij][iiterrs][ixx]->Fill(xtdev[ij]);
			if(isfid_l9_eff ) {time_mulxreso_split[ij][xbina][iiterrs][ixx]->Fill(xtdev[ij]);}
			/*
			int ybinext = int(yext[ij]/5)+1;

			  if(yext[ij]>1 && yext[ij]<20) {
			    if(xbandcond[0]) {
			    time_mulxreso_band[ij][iiterrs][ixx][0]->Fill(xtdev[ij]);
			    }else if(xcondband[1] || xcondband[2]) {
			     time_mulxreso_band[ij][iiterrs][ixx][ybinext]->Fill(xtdev[ij]);
			    }else if(xbandcond[3]) {
			    time_mulxreso_band[ij][iiterrs][ixx][5]->Fill(xtdev[ij]);
			    }
			  }else if(yext[ij]>=20 && yext[ij]<40) {
			    if(xbandcond[0]) {
			      time_mulxreso_band[ij][iiterrs][ixx][6]->Fill(xtdev[ij]);
			    }else if(xbandcond[1]) {
			       time_mulxreso_band[ij][iiterrs][ixx][7]->Fill(xtdev[ij]);
			    }else if(xbandcond[2]) {
			       time_mulxreso_band[ij][iiterrs][ixx][8]->Fill(xtdev[ij]);
			    }else if(xbandcond[3]) {
			       time_mulxreso_band[ij][iiterrs][ixx][9]->Fill(xtdev[ij]);
			    }
			  }else if(yext[ij]>=40 && yext[ij]<=62) {
			    if(xcondband[0] || xcondband[1]) {
			      time_mulxreso_band[ij][iiterrs][ixx][int(yext[ij]/5)+2]->Fill(xtdev[ij]);
			    }else if(xcondband[2] || xcondband[3]) {
			     time_mulxreso_band[ij][iiterrs][ixx][int(yext[ij]/5)+2+5]->Fill(xtdev[ij]);
			    }
			  }
			  */
		      }

		    }

#ifdef ISEFFICIENCY
		    if (xposEffPass[ij][iiterrs]) {
		      total_xt[ij][iiterrs]->Fill(istrxtime[ij], istrytime[ij]); //11/12/2011
		      if((!xusedtime[ij]) || abs(xtdev[ij])>15*time_xrms[ij]) {
			inefficiency_xt[ij][iiterrs]->Fill(istrxtime[ij], istrytime[ij]);
		      }
		    }
#endif
		  }
		} else {

		  //Muon Multiplicity Filling  - X-side
		  muon_xmul_eq[occulyr][iiterrs][nxtime][0]->Fill(xpts[occulyr].size());  //Just like that no criteria applied
		  int muonTrue=0;  //Position Criteria
		  int mulWithTime=0;
		  int mulWithBothPosTime=0;
		  for (int ix=0; ix<xpts[occulyr].size(); ix++) {
		    if(ix<xpts[occulyr].size()-1 &&  abs(xpts[occulyr][ix]-xpts[occulyr][ix+1])>1) {muonTrue= -1 ; break;}
		    if (abs(xextloc[occulyr]- xpts[occulyr][ix])<1) {muonTrue=1;}
		  }
		  for (int ix=0; ix<xpts[occulyr].size(); ix++) {
		    int isNoise=0;
		    for(int iy=0; iy<event->vxtdc_l[occulyr][xpts[occulyr][ix]%8]->size();iy++) {
		      double tmptime;
		      tmptime = 1250.0 + 0.1*(int(event->vxtdc_l[occulyr][xpts[occulyr][ix]%8]->at(iy)-event->tdc_ref_l[occulyr]));
		      if(abs(tmptime-xtext[occulyr])>20) {isNoise=1;}
		    }
		    if(isNoise==0) {mulWithTime++;}
		    if(isNoise==0 && muonTrue==1) {mulWithBothPosTime++;}
		  }
		  if(muonTrue<=0) {muon_xmul_eq[occulyr][iiterrs][nxtime][1]->Fill(0);}
		  if(muonTrue==1) {muon_xmul_eq[occulyr][iiterrs][nxtime][1]->Fill(xpts[occulyr].size());}

		  muon_xmul_eq[occulyr][iiterrs][nxtime][2]->Fill(mulWithTime);
		  muon_xmul_eq[occulyr][iiterrs][nxtime][3]->Fill(mulWithBothPosTime);
		  if(nxtime>=6) {
		    muon_xmul_gr[occulyr][iiterrs][6][0]->Fill(xpts[occulyr].size());
		    if(muonTrue==1) {muon_xmul_gr[occulyr][iiterrs][6][1]->Fill(xpts[occulyr].size());}
		    muon_xmul_gr[occulyr][iiterrs][6][2]->Fill(mulWithTime);
		    muon_xmul_gr[occulyr][iiterrs][6][3]->Fill(mulWithBothPosTime);
		  }
		  if(nxtime>=7) {
		    muon_xmul_gr[occulyr][iiterrs][7][0]->Fill(xpts[occulyr].size());
		    if(muonTrue==1) {muon_xmul_gr[occulyr][iiterrs][7][1]->Fill(xpts[occulyr].size());}
		    muon_xmul_gr[occulyr][iiterrs][7][2]->Fill(mulWithTime);
		    muon_xmul_gr[occulyr][iiterrs][7][3]->Fill(mulWithBothPosTime);
		  }
		  if(nxtime>=8) {
		    muon_xmul_gr[occulyr][iiterrs][8][0]->Fill(xpts[occulyr].size());
		    if(muonTrue==1) {muon_xmul_gr[occulyr][iiterrs][8][1]->Fill(xpts[occulyr].size());}
		    muon_xmul_gr[occulyr][iiterrs][8][2]->Fill(mulWithTime);
		    muon_xmul_gr[occulyr][iiterrs][8][3]->Fill(mulWithBothPosTime);
		  }




		  if (dist[occulyr] >=0 && xtime[occulyr] >-50.0) {
		    if (passxtime[occulyr]) {

		      if (abs(xtdev[occulyr]) >8.0) {
			//			file_out<<"xocculyrx "<<occulyr<<" "<<xtdev[occulyr]<<" "<<xext[occulyr]<<" "<<yext[occulyr]<<" "<<Nx<<" "<<Ny<<" "<<nxtime<<" "<<xtime[occulyr]<<" "<<xtext[occulyr]<<" "<<xt0chi2<<" "<<xhits[occulyr]<<" "<<yhits[occulyr]<<endl;
			if (abs(xtdev[occulyr]) <30.0) {
			  if (xtdev[occulyr]>0) {
			    nxypos_xytdev[0]->Fill(xextloc[occulyr], yextloc[occulyr]);
			  } else {
			    nxypos_xytdev[1]->Fill(xextloc[occulyr], yextloc[occulyr]);
			  }
			}
		      }

		      if (xusedtime[occulyr]) {
			xstr_xtdev[occulyr][iiter]->Fill(xextloc[occulyr], xtdev[occulyr]);
			widthx_timex[occulyr][iiter]->Fill(xtdev[occulyr], widthx[occulyr]); //160315

			ystr_xtdev[occulyr][iiter]->Fill(yextloc[occulyr], xtdev[occulyr]);

			if (abs(xtdev[occulyr])<20) {
			  nxystr_xtdev[occulyr][iiter]->Fill(xextloc[occulyr], yextloc[occulyr], 1.);
			  xystr_xtdev[occulyr][iiter]->Fill(xextloc[occulyr], yextloc[occulyr], xtdev[occulyr]+20.0);
			}
		      }
		      if (iiter==nmxiter-1) {
			double diff = -30;
			for (int jkl=0; jkl<6; jkl++) {
			  switch(jkl) {
			  case 0 : diff = rawxtime0[occulyr]-xtext[occulyr]; break;
			  case 1 : diff = rawxtime1[occulyr]-xtext[occulyr]; break;
			  case 2 : diff = rawxtime2[occulyr]-xtext[occulyr]; break;
			  case 3 : diff = timesx[occulyr]-xtext[occulyr]; break;
			  case 4 : diff = rawxtime3[occulyr]-xtext[occulyr]; break;
			  default : diff = xtime[occulyr]-xtext[occulyr]; break;
			  }
			  if (abs(diff)<20 && xusedtime[occulyr]) {
			    nxystr_xtdev[occulyr][nmxiter+jkl]->Fill(xextloc[occulyr], yextloc[occulyr], 1.);
			    xystr_xtdev[occulyr][nmxiter+jkl]->Fill(xextloc[occulyr], yextloc[occulyr], diff+20.0);
			  }
			}

			int nsize=xpts[occulyr].size();
			if (nsize>=1 && nsize<=nmxtimehit) {
			  //			  xpos_xdev[occulyr][nsize-1]->Fill(xposinstr[occulyr], Xdev[occulyr]);
			  //			  ypos_xdev[occulyr][nsize-1]->Fill(yposinstr[occulyr], Xdev[occulyr]);
			  xpos_xtdev[occulyr][nsize-1]->Fill(xposinstr[occulyr], xtdev[occulyr]);
			  ypos_xtdev[occulyr][nsize-1]->Fill(yposinstr[occulyr], xtdev[occulyr]);
			  xpos_xtdev_str[occulyr][nsize-1]->Fill(xposinstr[occulyr], diff);
			  for(int ij=0;ij<nband;ij++) {
			    if(xbandcond[occulyr][ij]) { xpos_xtdev_band_str[occulyr][ij][nsize-1]->Fill(xposinstr[occulyr], diff);}

			  }
			  xpos_xtdev_glb[occulyr][nsize-1]->Fill(xextloc[occulyr], diff);
			  ypos_xtdev_str[occulyr][nsize-1]->Fill(yposinstr[occulyr], diff);
			  ypos_xtdev_glb[occulyr][nsize-1]->Fill(yextloc[occulyr], diff);

			}
		      }
		    }

		    //		    if (xpts[occulyr].size()==1 && istrxtime[occulyr]>=0 && istrxtime[occulyr] <nstrip) {
		    if (xusedtime[occulyr]) {
		      if (istrxtime[occulyr]>=0 && istrxtime[occulyr] <nstrip) {
#ifdef TIMESLOPE
			time_xslope_pr[occulyr][istrxtime[occulyr]][iiter]->Fill( yextloc[occulyr], xtdev[occulyr]);
#endif
			if (passxtmy[occulyr]) {
			  time_xstrreso[occulyr][istrxtime[occulyr]][iiterrs]->Fill(xtdev[occulyr]);
			  tmptime_xreso[occulyr][iiterrs]->Fill(xtdev[occulyr]);
			}
		      }
		      xtime_exterr[occulyr][iiterrs]->Fill(xtexter[occulyr]);
		      xtime_exterr_vs_ndf[occulyr][iiterrs]->Fill(xtexter[occulyr], nxtime);
		      xtime_exterr_nTermSplit[occulyr][int(xextloc[occulyr]/20)][iiterrs]->Fill(xtexter[occulyr]);
		    }
		   //int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[occulyr]+0.5),int(yextloc[occulyr]+0.5)));//(nsplit/8)*int(xextloc[occulyr]/7.25)+int(yextloc[occulyr]/7.625);
		   //if(xbina>nsplit-1) xbina = nsplit-1;



		    if (passxtime[occulyr] && passxtmy[occulyr]) {
		      if (xusedtime[occulyr]) {
			time_xreso[occulyr][iiterrs]->Fill(xtdev[occulyr]);

			time_xreso_chi2_ndf[occulyr][iiterrs][0]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][0]->Fill(xtexter[occulyr]);
			if(xt0chi2/(nxtime-2)<20){time_xreso_chi2_ndf[occulyr][iiterrs][1]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][1]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<10){time_xreso_chi2_ndf[occulyr][iiterrs][2]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][2]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<9){time_xreso_chi2_ndf[occulyr][iiterrs][3]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][3]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<8){time_xreso_chi2_ndf[occulyr][iiterrs][4]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][4]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<7){time_xreso_chi2_ndf[occulyr][iiterrs][5]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][5]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<6){time_xreso_chi2_ndf[occulyr][iiterrs][6]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][6]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<5){time_xreso_chi2_ndf[occulyr][iiterrs][7]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][7]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<4){time_xreso_chi2_ndf[occulyr][iiterrs][8]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][8]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<3){time_xreso_chi2_ndf[occulyr][iiterrs][9]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][9]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<2){time_xreso_chi2_ndf[occulyr][iiterrs][10]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][10]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<1.5){time_xreso_chi2_ndf[occulyr][iiterrs][11]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][11]->Fill(xtexter[occulyr]);}
			if(xt0chi2/(nxtime-2)<1){time_xreso_chi2_ndf[occulyr][iiterrs][12]->Fill(xtdev[occulyr]); xtime_exterr_chi2_ndf[occulyr][iiterrs][12]->Fill(xtexter[occulyr]);}


			time_xreso_vs_ndf[occulyr][iiterrs]->Fill(xtdev[occulyr],nxtime);
			time_xreso_nTermSplit[occulyr][int(xextloc[occulyr]/20)][iiterrs]->Fill(xtdev[occulyr]);
			if(isfid_l9_eff ) { time_xreso_split[occulyr][xbina][iiterrs]->Fill(xtdev[occulyr]);}
			istime_xyreso[occulyr][iiterrs] = true;
		      }
		      int ixx = min(xhits[occulyr],nmxtimehit)-1;
		      if (ixx>=0) {time_mulxreso[occulyr][iiterrs][ixx]->Fill(xtdev[occulyr]);
			if(isfid_l9_eff ) {time_mulxreso_split[occulyr][xbina][iiterrs][ixx]->Fill(xtdev[occulyr]);}
			int ybinext = int(yext[occulyr]/5)+1;

			if(xextloc[occulyr]>0 && xextloc[occulyr]<nstrip-4 && yextloc[occulyr]>0 && yextloc[occulyr]<nstrip-4) {
			  // cout<<"ssssss"<<"   "<<iev<<"   "<<endl;
			  for(int ij=0;ij<nband;ij++) {
			    if(xbandcond[occulyr][ij]) {
			      //  cout<<iev<<"   "<<occulyr<<"   "<<ij<<endl;
			      time_mulxreso_band[occulyr][iiterrs][ixx][ij]->Fill(xtdev[occulyr]);
			      timex_mulxyreso_band[occulyr][iiterrs][int(xextloc[occulyr]+0.5)/2][int(yextloc[occulyr]+0.5)/2][ij]->Fill(xtdev[occulyr]);
			    }
			  }
			}

			/*
			if(yext[occulyr]>1 && yext[occulyr]<20) {
			    if(xbandcond[0]) {
			    time_mulxreso_band[occulyr][iiterrs][ixx][0]->Fill(xtdev[occulyr]);
			    }else if(xbandcond[1] || xbandcond[2]) {
			     time_mulxreso_band[occulyr][iiterrs][ixx][ybinext]->Fill(xtdev[occulyr]);
			    }else if(xbandcond[3]) {
			    time_mulxreso_band[occulyr][iiterrs][ixx][5]->Fill(xtdev[occulyr]);
			    }
			  }else if(yext[occulyr]>=20 && yext[occulyr]<40) {
			    if(xbandcond[0]) {
			      time_mulxreso_band[occulyr][iiterrs][ixx][6]->Fill(xtdev[occulyr]);
			    }else if(xbandcond[1]) {
			       time_mulxreso_band[occulyr][iiterrs][ixx][7]->Fill(xtdev[occulyr]);
			    }else if(xbandcond[2]) {
			       time_mulxreso_band[occulyr][iiterrs][ixx][8]->Fill(xtdev[occulyr]);
			    }else if(xbandcond[3]) {
			       time_mulxreso_band[occulyr][iiterrs][ixx][9]->Fill(xtdev[occulyr]);
			    }
			  }else if(yext[occulyr]>=40 && yext[occulyr]<=62) {
			    if(xbandcond[0] || xbandcond[1]) {
			      time_mulxreso_band[occulyr][iiterrs][ixx][int(yext[occulyr]/5)+2]->Fill(xtdev[occulyr]);
			    }else if(xbandcond[2] || xbandcond[3]) {
			     time_mulxreso_band[occulyr][iiterrs][ixx][int(yext[occulyr]/5)+7]->Fill(xtdev[occulyr]);
			    }
			  }
			*/
		      }

		    }
#ifdef ISEFFICIENCY
		    if (xposEffPass[occulyr][iiterrs]) {
		      total_xt[occulyr][iiterrs]->Fill(istrxtime[occulyr], istrytime[occulyr]); //11/12/2011
		      if((!xusedtime[occulyr]) || abs(xtdev[occulyr])>15*time_xrms[occulyr]) {
			inefficiency_xt[occulyr][iiterrs]->Fill(istrxtime[occulyr], istrytime[occulyr]); //150112
		      }
		    }
#endif
		  } //if (dist[occulyr] >=0 && xtime[occulyr] >-50.0)

		  if (iiter==nmxiter-1) { //Correlated hits
		    if (abs(timesx[occulyr] - xtext[occulyr])<10.0) {
		      double expp = -100;
		      double expdiff = 10000;
		      for (int ix=0; ix<xpts[occulyr].size(); ix++) {
			if (abs(xextloc[occulyr] - xpts[occulyr][ix])<1.0) {
			  double tmpexpdiff = abs(xextloc[occulyr] - xpts[occulyr][ix]);
			  if(tmpexpdiff < expdiff) {
			    expdiff =tmpexpdiff;  expp = xpts[occulyr][ix];// break;
			  }
			}
		      }

		      if (expp >-50) {
			h_xtcorstrips[occulyr]->Fill(expp, -1.0);
			h_xytcorstrips[occulyr]->Fill(expp, -1.0);
			for (int ix=0; ix<xpts[occulyr].size(); ix++) {
			  if (abs(expp - xpts[occulyr][ix])>0) { //choose all of them and then use maximum value before plot
			    h_xtcorstrips[occulyr]->Fill(expp, xpts[occulyr][ix]);
			  }
			}

			for (int ix=0; ix<ypts[occulyr].size(); ix++) {
			  h_xytcorstrips[occulyr]->Fill(expp, ypts[occulyr][ix]);
			}
		      }
		    }
		  } //if (iiter==nmxiter-1)
		} //else of if (occulyr >=nlayer)
	      } // if (nxtime>=nmnhits/*-ntcor*/ && nxtfail==0)

	      for (int ix=0; ix<nlayer; ix++) { ytext[ix]= ytexter[ix] = 100;}

	      int tmpntxy = 0;
#if defined(MADURAIINTM2) || defined(MADURAIINTM3) //excluding 4 NINO layer
	      StraightLineFit ytimefit(iTimeSlopeConst, dist, ytime, timeserry2, yusedtime, occulyr, 4, ytcorstr, ytcorend, float(1000.0));
#else
	      StraightLineFit ytimefit(iTimeSlopeConst, dist, ytime, timeserry2, yusedtime,occulyr,occulyr, ytcorstr, ytcorend, float(1000.0));
#endif
	      ytimefit.GetParameters(nytfail, yc0inters, timeyslope);
	      //	    ytimefit.GetError(errcst, errlin, errcov);
	      ytimefit.GetChisqure(nytime, yt0chi2);
	      ytimefit.GetFitValues(ytext, ytdev, ytexter);

	      //150101		  tmpytent[ij]= 0;

	      // for(int ij=0;ij<nlayer;ij++) {
	      if(ntcor==1 && iiter == nmxiter-1) {
	      tmpytdev = ytdev[occulyr];
	      tmpnytime = nytime;
	      tmpytime = ytime[occulyr];
	      tmpytext = ytext[occulyr];
	      }
	      //
	      if (nytfail==0 && isfill) {
		h_tchisqy->Fill(yt0chi2);
		if (nytime-2>0) {
		  h_treduchisqy->Fill(yt0chi2/(nytime-2));
		  double probx =TMath::Prob(yt0chi2, nytime-2);
		  h_ytprob->Fill(probx);
		  ytprob_vs_tchisqy->Fill(probx, yt0chi2);
		  ytprob_vs_ndf->Fill(probx, nytime);
		  tchisqy_vs_ndf->Fill(yt0chi2, nytime);
		  ytprob_vs_treduchisqy->Fill(probx, yt0chi2/(nytime-2));
		  int ibin = getbinid(nytime, nprob, probs);
		  if (ibin>=0 && ibin<nprob) {h_ytnprob[ibin]->Fill(probx);}
		}
		h_tyndf->Fill(nytime);
	      }

	      if ((isfill || occulyr==nlayer) && nytfail==0 && nytime >nmnhits && yt0chi2/(nytime-2)<mxtimechisq) {
		for(int ij=0;ij<nlayer;ij++) {
		  if (dist[ij]<0 || ytime[ij]<-50) continue;
		  double dt1 = ytdev[ij]-biasinytime[ij];
		  if (abs(dt1)>4*sigmarpc || (!yusedtime[ij])) continue;

		  //140923
		  //		  cout <<"occulyr "<< occulyr<<endl;
		  if (occulyr==nlayer) {
		    if (passytime[ij]) {
		      ntotytime++; totytime +=ytime[ij];
		      for(int jk=ij+1;jk<nlayer;jk++){
			if (dist[jk]<0 || ytime[jk]<-50 || abs(ytdev[jk]-biasinytime[jk])>4*sigmarpc || (!yusedtime[jk])  || (!passytime[jk])) continue;
			timey_correl[ij][jk][iiter+1]->Fill(ytime[jk]-ytime[ij] - (dist[ij]-dist[jk])/cval);
			if (iiter==0) { timey_correl[ij][jk][0]->Fill(rawytime[jk]-rawytime[ij] - (dist[ij]-dist[jk])/cval);}
		      } // for(jk=ij+1;jk<nlayer;jk++)
		      //		      cout<<"passxtime[ij] "<<passxtime[ij]<<endl;
		      if (passxtime[ij]) {
			//X-Y correlation

			if (nxtfail==0 && nxtime >nmnhits && xt0chi2/(nxtime-2)<mxtimechisq && xtime[ij] >-50 &&
			    abs(xtdev[ij]-biasinxtime[ij])<4*sigmarpc && xusedtime[ij]) {
			  int ipix = 4*int((4*istrxtime[ij])/nstrip) + int((4*istrytime[ij])/nstrip);
			  if (ipix<0 || ipix>15) cout <<" ipix "<<ipix<<" "<<istrxtime[ij]<<" "<<istrytime[ij]<<endl;

			  timexy_correl[ij][ipix][iiter+1]->Fill(ytime[ij]-xtime[ij]);
			  timexy_correl[ij][npixel][iiter+1]->Fill(ytime[ij]-xtime[ij]);
			  if (iiter==0) {
			    timexy_correl[ij][ipix][0]->Fill(rawytime[ij]-rawxtime[ij]);
			    timexy_correl[ij][npixel][0]->Fill(rawytime[ij]-rawxtime[ij]);
			  }
			  //			  cout <<"iiter "<<itter<<endl;
			  if (iiter==0) { //compare timing of all combination of X/Y strips
			    int ixx = istrxtime[ij];
			    int iyy = istrytime[ij];
			    if (ixx>=0 && ixx <nstrip && iyy>=0 && iyy <nstrip) {
			      //			      cout <<"ijxy "<<ij<<" "<<ixx<<" "<<iyy<<" "<<rawytime0[ij]-rawxtime0[ij]<<" "<<rawytime3[ij]-rawxtime3[ij]<<endl;
			      indtimexy_correl[ij][ixx][iyy][0]->Fill(rawytime0[ij]-rawxtime0[ij]);
			      indtimexy_correl[ij][ixx][iyy][1]->Fill(rawytime3[ij]-rawxtime3[ij]);
#ifdef C217STRIP
			      indtimexy_prof[ij][0]->Fill(ixx+iyy, rawytime0[ij]-rawxtime0[ij]);
			      indtimexy_prof[ij][1]->Fill(ixx+iyy, rawytime3[ij]-rawxtime3[ij]);
#else
			      indtimexy_prof[ij][0]->Fill(iyy-ixx+nstrip, rawytime0[ij]-rawxtime0[ij]);
			      indtimexy_prof[ij][1]->Fill(iyy-ixx+nstrip, rawytime3[ij]-rawxtime3[ij]);
#endif
			    } else {
			      cout <<"wrong istr "<< ij<<" "<<ixx<<" "<<iyy<<endl;
			    }
			  }
			}
		      } // if (passxtime[ij])
		    } // if (passytime[ij])
		  } //if (occulyr==nlayer)
		} // ij
	      } // if ((isfill || occulyr==nlayer) && nytfail==0 && nytime >nmnhits && yt0chi2/(nytime-2)<mxtimechisq)

	      ntxyla = (xtimefit.GetLayerIds()<<nlayer) +  ytimefit.GetLayerIds();

	      if (isfill) {
		nTotalt++;
		for (int ix=0; ix<nlayer; ix++) {
		  if (Xpos[ix]>=-1 && Xpos[ix]<=nstrip && abs(xtdev[ix])<5 && xusedtime[ix]) {
		    for (int iy=ix+1; iy<nlayer; iy++) {
		      if (Xpos[iy]>=-1 && Xpos[iy]<=nstrip && abs(xtdev[iy])<5 && xusedtime[iy]) {
			h_xtcorhits->Fill(ix, iy);
		      }
		    }
		    for (int iy=0; iy<nlayer; iy++) {
		      if (Ypos[iy]>=-1 && Ypos[iy]<=nstrip && abs(ytdev[iy])<5 && yusedtime[iy]) {
			h_xytcorhits->Fill(ix, iy);
		      }
		    }
		  }
		}
		for (int ix=0; ix<nlayer-1; ix++) {
		  if (Ypos[ix]>=-1 && Ypos[ix]<=nstrip && abs(ytdev[ix])<5 && yusedtime[ix]) {
		    for (int iy=ix+1; iy<nlayer; iy++) {
		      if (Ypos[iy]>=-1 && Ypos[iy]<=nstrip && abs(ytdev[iy])<5  && yusedtime[iy]) {
			h_ytcorhits->Fill(ix, iy);
		      }
		    }
		  }
		}
	      }

              if (nytime>=nmnhits/*-ntcor*/ && nytfail==0/* && yt0chi2/(nytime-2)<2*/) {//4Nov        //jim jim tightening criteria of chi2/ndf for better alignment correction
		if (yt0chi2/(nytime-2)<mxtimechisq) {

		  if (filloccu) {
		    if (iev>0 && tmptimerate >=tmpoldmutimeyrate+timebin) {
		      file_out <<"timeyrate "<<datafile<<" "<<iev<<" "<<nallmutimeyrate<<" "<<nmutimeyrate<<endl;
		      mutimeyrate->Fill(nmutimeyrate*(1.0/timebin));
		      mutimeyratex->Fill(nallmutimeyrate, nmutimeyrate*(1.0/timebin));
		      nmutimeyrate = 0;
		      nallmutimeyrate++;
		      tmpoldmutimeyrate = tmptimerate;
		    }
		    nmutimeyrate++;
		  }

		  dir_cy[occulyr][iiter]->Fill(-cval*timeyslope);
		  int ibin = getbinid(nytime, nprob, probs);
		  if (ibin>=0 && ibin<nprob) {dir_cylay[occulyr][iiter][ibin]->Fill(-cval*timeyslope);}
		  dir_cy2[occulyr][iiter]->Fill(ytimefit.GetSlope2());
		  dir_c0y[occulyr][iiter]->Fill(yc0inters);
		}

		if (isfill) dir_cychi->Fill(yt0chi2, -1./timeyslope/cval);
		//		cout<<"nxtime "<<nxtime<<" "<< nmnhits<<" "<<nxtfail<<" "<<isfill<<endl;
		if (nxtime>=nmnhits/*-ntcor*/ && nxtfail==0) {
		  if (filloccu) {
		    for (int iz=0; iz<nlayer; iz++) {
		      xlayer_allmutimemult[iz]->Fill(xptsall[iz].size());
		      ylayer_allmutimemult[iz]->Fill(yptsall[iz].size());
		    }
		  }

		  double xt = (timexslope <0) ? log10(1+abs(1./timexslope/cval)) : -log10(1+abs(1./timexslope/cval));
		  double yt = (timeyslope <0) ? log10(1+abs(1./timeyslope/cval)) : -log10(1+abs(1./timeyslope/cval));

		  dir_cxy->Fill(xt, yt);
		}

		if (occulyr >=nlayer) {

		  for (int ij=0; ij<nlayer-1; ij++) { //review1
		    for(int jk=ij+1; jk<nlayer; jk++) {
		      //cval is in cm/ns (29.979)
		      if (dist[ij]<0 || ytime[ij] <-50 || dist[jk]<0 || ytime[jk] <-50) continue;
		      if (passytime[ij] && passytmx[ij] && passytime[jk] && passytmx[jk]) {
			if (yusedtime[ij] && yusedtime[jk]) {
			  timeshift_yreso[ij][jk]->Fill((ytime[ij]-dist[ij]*timeyslope)-(ytime[jk]-dist[jk]*timeyslope));
			  timeshift_xyreso[ij][jk]->Fill(((xtime[ij]-dist[ij]*timexslope)+(ytime[ij]-dist[ij]*timeyslope))/2.0-((xtime[jk]-dist[jk]*timexslope)+(ytime[jk]-dist[jk]*timeyslope))/2.0);
			}
		      }
		    }
		  }

                  for (int ij=0; ij<nlayer; ij++) {
		    if (dist[ij]<0 || ytime[ij] <-50) continue;
		    if (yusedtime[ij]) {
		      if (istrytime[ij]>=0 && istrytime[ij] <nstrip && passytmx[ij]) {
			time_ystrreso[ij][istrytime[ij]][iiterrs]->Fill(ytdev[ij]);
			tmptime_yreso[ij][iiterrs]->Fill(ytdev[ij]);
		      }
		      ytime_exterr[ij][iiterrs]->Fill(ytexter[ij]);
		      ytime_exterr_vs_ndf[ij][iiterrs]->Fill(ytexter[ij], nytime);
		      ytime_exterr_nTermSplit[ij][int(yextloc[ij]/20)][iiterrs]->Fill(ytexter[ij]);
		    }
		    //int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[ij]+0.5),int(yextloc[ij]+0.5)));//(nsplit/8)*int(xextloc[ij]/7.25)+int(yextloc[ij]/7.625);
		    if (passytime[ij] && passytmx[ij]) {
		      if (yusedtime[ij]) {
#ifndef MONTECARLO
			if (filloccu) {time_yreso_set[ij][iset]->Fill(ytdev[ij]);}
#endif
			time_yreso[ij][iiterrs]->Fill(ytdev[ij]);
			time_yreso_vs_ndf[ij][iiterrs]->Fill(ytdev[ij],nytime);
			time_yreso_nTermSplit[ij][int(yextloc[ij]/20)][iiterrs]->Fill(ytdev[ij]);
			if(isfid_l9_eff ) { time_yreso_split[ij][xbina][iiterrs]->Fill(ytdev[ij]);}
			if (istime_xyreso[ij][iiterrs]) {
			  time_xyreso[ij][iiterrs]->Fill(xtdev[ij], ytdev[ij]);
			}
		      }
		      int iyy = min(yhits[ij],nmxtimehit)-1;
		      if (iyy>=0) {time_mulyreso[ij][iiterrs][iyy]->Fill(ytdev[ij]);
			if(isfid_l9_eff) {time_mulyreso_split[ij][xbina][iiterrs][iyy]->Fill(ytdev[ij]);}
			/*
			int xbinext = int(xext[ij]/5)+1;

			  if(xext[ij]>1 && xext[ij]<20) {
			    if(ybandcond[0]) {
			    time_mulyreso_band[ij][iiterrs][iyy][0]->Fill(ytdev[ij]);
			    }else if(ycondband[1] || ycondband[2]) {
			     time_mulyreso_band[ij][iiterrs][iyy][xbinext]->Fill(ytdev[ij]);
			    }else if(ybandcond[3]) {
			    time_mulyreso_band[ij][iiterrs][iyy][5]->Fill(ytdev[ij]);
			    }
			  }else if(xext[ij]>=20 && xext[ij]<40) {
			    if(ybandcond[0]) {
			      time_mulyreso_band[ij][iiterrs][iyy][6]->Fill(ytdev[ij]);
			    }else if(xbandcond[1]) {
			       time_mulyreso_band[ij][iiterrs][iyy][7]->Fill(ytdev[ij]);
			    }else if(ybandcond[2]) {
			       time_mulyreso_band[ij][iiterrs][iyy][8]->Fill(ytdev[ij]);
			    }else if(ybandcond[3]) {
			       time_mulyreso_band[ij][iiterrs][iyy][9]->Fill(ytdev[ij]);
			    }
			  }else if(xext[ij]>=40 && xext[ij]<=62) {
			    if(ycondband[0] || ycondband[1]) {
			      time_mulyreso_band[ij][iiterrs][iyy][int(xext[ij]/5)+2]->Fill(ytdev[ij]);
			    }else if(ycondband[2] || ycondband[3]) {
			     time_mulyreso_band[ij][iiterrs][iyy][int(xext[ij]/5)+2+5]->Fill(ytdev[ij]);
			    }
			  }
			*/
		      }

		    }
#ifdef ISEFFICIENCY
		    if (yposEffPass[ij][iiterrs]) {
		      total_yt[ij][iiterrs]->Fill(istrxtime[ij], istrytime[ij]); //11/12/2011
		      if((!yusedtime[ij]) || abs(ytdev[ij])>15*time_yrms[ij]) {
			inefficiency_yt[ij][iiterrs]->Fill(istrxtime[ij], istrytime[ij]);
		      }
		    }
#endif
		  }
		} else {

		  //Muon Multiplicity Filling - YSide
		  muon_ymul_eq[occulyr][iiterrs][nytime][0]->Fill(ypts[occulyr].size());  //Just like that no criteria applied
		  int muonTrue=0;  //Position Criteria
		  int mulWithTime=0;
		  int mulWithBothPosTime=0;
		  for (int ix=0; ix<ypts[occulyr].size(); ix++) {
		    if (abs(yextloc[occulyr]- ypts[occulyr][ix])<1) {muonTrue=1;}
		  }
		  for (int ix=0; ix<ypts[occulyr].size(); ix++) {
		    int isNoise=0;
		    for(int iy=0; iy<event->vytdc_l[occulyr][ypts[occulyr][ix]%8]->size();iy++) {
		      double tmptime;
		      tmptime = 1250.0 + 0.1*(int(event->vytdc_l[occulyr][ypts[occulyr][ix]%8]->at(iy)-event->tdc_ref_l[occulyr]));
		      if(abs(tmptime-ytext[occulyr])>20) {isNoise=1;}
		    }
		    if(isNoise==0) {mulWithTime++;}
		    if(isNoise==0 && muonTrue==1) {mulWithBothPosTime++;}
		  }
		  if(muonTrue==1) {muon_ymul_eq[occulyr][iiterrs][nytime][1]->Fill(ypts[occulyr].size());}
		  muon_ymul_eq[occulyr][iiterrs][nytime][2]->Fill(mulWithTime);
		  muon_ymul_eq[occulyr][iiterrs][nytime][3]->Fill(mulWithBothPosTime);
		  if(nytime>=6) {
		    muon_ymul_gr[occulyr][iiterrs][6][0]->Fill(ypts[occulyr].size());
		    if(muonTrue==1) {muon_ymul_gr[occulyr][iiterrs][6][1]->Fill(ypts[occulyr].size());}
		    muon_ymul_gr[occulyr][iiterrs][6][2]->Fill(mulWithTime);
		    muon_ymul_gr[occulyr][iiterrs][6][3]->Fill(mulWithBothPosTime);
		  }
		  if(nytime>=7) {
		    muon_ymul_gr[occulyr][iiterrs][7][0]->Fill(ypts[occulyr].size());
		    if(muonTrue==1) {muon_ymul_gr[occulyr][iiterrs][7][1]->Fill(ypts[occulyr].size());}
		    muon_ymul_gr[occulyr][iiterrs][7][2]->Fill(mulWithTime);
		    muon_ymul_gr[occulyr][iiterrs][7][3]->Fill(mulWithBothPosTime);
		  }
		  if(nytime>=8) {
		    muon_ymul_gr[occulyr][iiterrs][8][0]->Fill(ypts[occulyr].size());
		    if(muonTrue==1) {muon_ymul_gr[occulyr][iiterrs][8][1]->Fill(ypts[occulyr].size());}
		    muon_ymul_gr[occulyr][iiterrs][8][2]->Fill(mulWithTime);
		    muon_ymul_gr[occulyr][iiterrs][8][3]->Fill(mulWithBothPosTime);
		  }



		  if (dist[occulyr]>=0 && ytime[occulyr] >-50) {
		    if (passytime[occulyr]) {
		      if (xtime[occulyr] >-50 && passxtime[occulyr]) {
			if ((abs(xtdev[occulyr])>5.0 || abs(ytdev[occulyr])>5.0) &&
			    abs(xtdev[occulyr])<30.0 && Ny>=9 && abs(ytdev[occulyr])<30.0 && Nx>=9) {
			  //			  double xx = xext[occulyr]*3.0+3.5;
			  //			  double yy = yext[occulyr]*3.0+3.5;

			  if (xextloc[occulyr]<1 || xextloc[occulyr]>60) {
			    xtdev_ytdev[0]->Fill(xtdev[occulyr], ytdev[occulyr]);
			  } else if (yextloc[occulyr]<1 || yextloc[occulyr]>61) {
			    xtdev_ytdev[1]->Fill(xtdev[occulyr], ytdev[occulyr]);
			  } else if ((abs(xextloc[occulyr]-5.0)<1.0 || abs(xextloc[occulyr]-11.0)<1.0 || abs(xextloc[occulyr]-16.0)<1.0 || abs(xextloc[occulyr]-22.0)<1.0 || abs(xextloc[occulyr]-28.0)<1.0) && (abs(yextloc[occulyr]-5.0)<1.0 || abs(yextloc[occulyr]-11.0)<1.0 || abs(yextloc[occulyr]-17.0)<1.0 || abs(yextloc[occulyr]-23.0)<1.0 || abs(yextloc[occulyr]-29)<1.0)) {

			    xtdev_ytdev[2]->Fill(xtdev[occulyr], ytdev[occulyr]);
			  } else {
			    xtdev_ytdev[3]->Fill(xtdev[occulyr], ytdev[occulyr]);
			  }
			}
		      }

		      if (abs(ytdev[occulyr])>8.0) {
			//			file_out<<"yocculyry "<<occulyr<<" "<<ytdev[occulyr]<<" "<<xext[occulyr]<<" "<<yext[occulyr]<<" "<<Nx<<" "<<Ny<<" "<<nytime<<" "<<ytime[occulyr]<<" "<<ytext[occulyr]<<" "<<yt0chi2<<" "<<xhits[occulyr]<<" "<<yhits[occulyr]<<" "<<xtdev[occulyr]<<endl;
			if (abs(ytdev[occulyr]) <30.0) {
			  if (ytdev[occulyr]>0) {
			    nxypos_xytdev[2]->Fill(xextloc[occulyr], yextloc[occulyr]);
			  } else {
			    nxypos_xytdev[3]->Fill(xextloc[occulyr], yextloc[occulyr]);
			  }

			  if (xtime[occulyr] >-50 && passxtime[occulyr]) {


			    if (abs(xtdev[occulyr])>8.0 &&abs(xtdev[occulyr])<30.0) {
			      if (xtdev[occulyr]>0) {
				if (ytdev[occulyr]>0) {
				  nxypos_xytdev[4]->Fill(xextloc[occulyr], yextloc[occulyr]);
				} else {
				  nxypos_xytdev[5]->Fill(xextloc[occulyr], yextloc[occulyr]);
				}
			      } else {
				if (ytdev[occulyr]>0) {
				  nxypos_xytdev[6]->Fill(xextloc[occulyr], yextloc[occulyr]);
				} else {
				  nxypos_xytdev[7]->Fill(xextloc[occulyr], yextloc[occulyr]);
				}
			      }
			    }
			  }
			}
		      }
		      if (yusedtime[occulyr]) {
			widthx_timex[occulyr][iiter]->Fill(ytdev[occulyr], widthy[occulyr]);
			xstr_ytdev[occulyr][iiter]->Fill(xextloc[occulyr], ytdev[occulyr]);
			ystr_ytdev[occulyr][iiter]->Fill(yextloc[occulyr], ytdev[occulyr]);

			if (abs(ytdev[occulyr])<20) {
			  nxystr_ytdev[occulyr][iiter]->Fill(xextloc[occulyr], yextloc[occulyr], 1.);
			  xystr_ytdev[occulyr][iiter]->Fill(xextloc[occulyr], yextloc[occulyr], ytdev[occulyr]+20.0);
			}
		      }

		      if (iiter==nmxiter-1) {
			double diff = -30;
			for (int jkl=0; jkl<6; jkl++) {
			  switch(jkl) {
			  case 0 : diff = rawytime0[occulyr]-ytext[occulyr]; break;
			  case 1 : diff = rawytime1[occulyr]-ytext[occulyr]; break;
			  case 2 : diff = rawytime2[occulyr]-ytext[occulyr]; break;
			  case 3 : diff = timesy[occulyr]-ytext[occulyr]; break;
			  case 4 : diff = rawytime3[occulyr]-ytext[occulyr]; break;
			  default : diff = ytime[occulyr]-ytext[occulyr]; break;
			  }
			  if (abs(diff)<20 && yusedtime[occulyr]) {
			    nxystr_ytdev[occulyr][nmxiter+jkl]->Fill(xextloc[occulyr], yextloc[occulyr], 1.);
			    xystr_ytdev[occulyr][nmxiter+jkl]->Fill(xextloc[occulyr], yextloc[occulyr], diff+20.0);
			  }
			}

			int nsize=ypts[occulyr].size();
			if (nsize>=1 && nsize<=nmxtimehit) {
			  xpos_ytdev[occulyr][nsize-1]->Fill(xposinstr[occulyr], ytdev[occulyr]);
			  ypos_ytdev[occulyr][nsize-1]->Fill(yposinstr[occulyr], ytdev[occulyr]);

			  ypos_ytdev_str[occulyr][nsize-1]->Fill(yposinstr[occulyr], diff);
			  for(int ij=0;ij<nband;ij++) {
			  if(ybandcond[occulyr][ij]) { ypos_ytdev_band_str[occulyr][ij][nsize-1]->Fill(yposinstr[occulyr], diff);}
			  }
			  ypos_ytdev_glb[occulyr][nsize-1]->Fill(yextloc[occulyr], diff);
			  xpos_ytdev_str[occulyr][nsize-1]->Fill(xposinstr[occulyr], diff);
			  xpos_ytdev_glb[occulyr][nsize-1]->Fill(xextloc[occulyr], diff);
			}
		      }
		    }

		    //		    if (ypts[occulyr].size()==1 && istrytime[occulyr]>=0 && istrytime[occulyr] <nstrip) {
		    if (yusedtime[occulyr]) {
		      if (istrytime[occulyr]>=0 && istrytime[occulyr] <nstrip) {
#ifdef TIMESLOPE
			time_yslope_pr[occulyr][istrytime[occulyr]][iiter]->Fill(xextloc[occulyr], ytdev[occulyr]);
#endif
			//if (ypts[occulyr].size()>nmxusedtimehit) cout <<"occulyr "<<occulyr<<" "<<ypts[occulyr].size()<<" "<<nmxusedtimehit<<" "<<ytdev[occulyr]<<endl;
			if (passytmx[occulyr]) {
			  time_ystrreso[occulyr][istrytime[occulyr]][iiterrs]->Fill(ytdev[occulyr]);
			  tmptime_yreso[occulyr][iiterrs]->Fill(ytdev[occulyr]);
			}
		      }
		      ytime_exterr[occulyr][iiterrs]->Fill(ytexter[occulyr]);
		      ytime_exterr_vs_ndf[occulyr][iiterrs]->Fill(ytexter[occulyr], nytime);
		      ytime_exterr_nTermSplit[occulyr][int(yextloc[occulyr]/20)][iiterrs]->Fill(ytexter[occulyr]);
		    }
		    //int xbina = int(10.*efficiency_x_l9->GetBinContent(int(xextloc[occulyr]+0.5),int(yextloc[occulyr]+0.5)));//(nsplit/8)*int(xextloc[occulyr]/7.25)+int(yextloc[occulyr]/7.625);
		    //if(xbina>nsplit-1) xbina = nsplit-1;



		    if (passytime[occulyr] && passytmx[occulyr]) {
		      if (yusedtime[occulyr]) {
			time_yreso[occulyr][iiterrs]->Fill(ytdev[occulyr]);

			time_yreso_chi2_ndf[occulyr][iiterrs][0]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][0]->Fill(ytexter[occulyr]);
			if(yt0chi2/(nytime-2)<20){time_yreso_chi2_ndf[occulyr][iiterrs][1]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][1]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<10){time_yreso_chi2_ndf[occulyr][iiterrs][2]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][2]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<9){time_yreso_chi2_ndf[occulyr][iiterrs][3]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][3]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<8){time_yreso_chi2_ndf[occulyr][iiterrs][4]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][4]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<7){time_yreso_chi2_ndf[occulyr][iiterrs][5]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][5]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<6){time_yreso_chi2_ndf[occulyr][iiterrs][6]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][6]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<5){time_yreso_chi2_ndf[occulyr][iiterrs][7]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][7]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<4){time_yreso_chi2_ndf[occulyr][iiterrs][8]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][8]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<3){time_yreso_chi2_ndf[occulyr][iiterrs][9]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][9]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<2){time_yreso_chi2_ndf[occulyr][iiterrs][10]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][10]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<1.5){time_yreso_chi2_ndf[occulyr][iiterrs][11]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][11]->Fill(ytexter[occulyr]);}
			if(yt0chi2/(nytime-2)<1){time_yreso_chi2_ndf[occulyr][iiterrs][12]->Fill(ytdev[occulyr]); ytime_exterr_chi2_ndf[occulyr][iiterrs][12]->Fill(ytexter[occulyr]);}

			if(passxtime[occulyr] && passxtmy[occulyr] && xusedtime[occulyr]) {time_xtdev_ytdev[occulyr][iiterrs]->Fill(xtdev[occulyr],ytdev[occulyr]);}
			time_yreso_vs_ndf[occulyr][iiterrs]->Fill(ytdev[occulyr],nytime);
			time_yreso_nTermSplit[occulyr][int(yextloc[occulyr]/20)][iiterrs]->Fill(ytdev[occulyr]);
			if(isfid_l9_eff ) { time_yreso_split[occulyr][xbina][iiterrs]->Fill(ytdev[occulyr]);}
			if (istime_xyreso[occulyr][iiterrs]) {
			  time_xyreso[occulyr][iiterrs]->Fill(xtdev[occulyr], ytdev[occulyr]);
			}
		      }
		      int iyy = min(yhits[occulyr],nmxtimehit)-1;
		      if (iyy>=0) {
			time_mulyreso[occulyr][iiterrs][iyy]->Fill(ytdev[occulyr]);
			if(isfid_l9_eff ) {time_mulyreso_split[occulyr][xbina][iiterrs][iyy]->Fill(ytdev[occulyr]);}
			int xbinext = int(xext[occulyr]/5)+1;

			if(xextloc[occulyr]>0 && xextloc[occulyr]<nstrip-4 && yextloc[occulyr]>0 && yextloc[occulyr]<nstrip-4) {
			  for(int ij=0;ij<nband;ij++) {
			    if(ybandcond[occulyr][ij]) {
			      time_mulyreso_band[occulyr][iiterrs][iyy][ij]->Fill(ytdev[occulyr]);
			      timey_mulxyreso_band[occulyr][iiterrs][int(xextloc[occulyr]+0.5)/2][int(yextloc[occulyr]+0.5)/2][ij]->Fill(ytdev[occulyr]);
			    }
			  }
			}
			/*
			  if(xext[occulyr]>1 && xext[occulyr]<20) {
			    if(ybandcond[0]) {
			    time_mulyreso_band[occulyr][iiterrs][iyy][0]->Fill(ytdev[occulyr]);
			    }else if(ybandcond[1] || ybandcond[2]) {
			     time_mulyreso_band[occulyr][iiterrs][iyy][xbinext]->Fill(ytdev[occulyr]);
			    }else if(ybandcond[3]) {
			    time_mulyreso_band[occulyr][iiterrs][iyy][5]->Fill(ytdev[occulyr]);
			    }
			  }else if(xext[occulyr]>=20 && xext[occulyr]<40) {
			    if(ybandcond[0]) {
			      time_mulyreso_band[occulyr][iiterrs][iyy][6]->Fill(ytdev[occulyr]);
			    }else if(ybandcond[1]) {
			       time_mulyreso_band[occulyr][iiterrs][iyy][7]->Fill(ytdev[occulyr]);
			    }else if(ybandcond[2]) {
			       time_mulyreso_band[occulyr][iiterrs][iyy][8]->Fill(ytdev[occulyr]);
			    }else if(ybandcond[3]) {
			       time_mulyreso_band[occulyr][iiterrs][iyy][9]->Fill(ytdev[occulyr]);
			    }
			  }else if(xext[occulyr]>=40 && xext[occulyr]<=62) {
			    if((ybandcond[0]) || (ybandcond[1])) {
			      time_mulyreso_band[occulyr][iiterrs][iyy][int(xext[occulyr]/5)+2]->Fill(ytdev[occulyr]);
			    }else if(ybandcond[2] || ybandcond[3]) {
			     time_mulyreso_band[occulyr][iiterrs][iyy][int(xext[occulyr]/5)+2+5]->Fill(ytdev[occulyr]);
			    }
			  }
			*/
		      }

		    }
#ifdef ISEFFICIENCY
		    if (yposEffPass[occulyr][iiterrs]) {
		      total_yt[occulyr][iiterrs]->Fill(istrxtime[occulyr], istrytime[occulyr]); //11/12/2011
		      if((!yusedtime[occulyr]) || abs(ytdev[occulyr])<15*time_yrms[occulyr]) {
			inefficiency_yt[occulyr][iiterrs]->Fill(istrxtime[occulyr], istrytime[occulyr]);
		      }
		    }
#endif
		  } //if (dist[occulyr] >=0 && xtime[occulyr] >-50.0)
		  if (iiter==nmxiter-1) { //Correlated hits
		    if (abs(timesy[occulyr] - ytext[occulyr])<10.0) {
		      double expp = -100;
		      double expdiff = 10000;
		      for (int ix=0; ix<ypts[occulyr].size(); ix++) {
			if (abs(yextloc[occulyr] - ypts[occulyr][ix])<1.0) {
			  double tmpexpdiff = abs(yextloc[occulyr] - ypts[occulyr][ix]);
			  if(tmpexpdiff < expdiff) {
			    expdiff =tmpexpdiff;  expp = ypts[occulyr][ix]; //break;
			  }
			}
		      }

		      if (expp >-50) {
			h_ytcorstrips[occulyr]->Fill(expp, -1.0);
			h_yxtcorstrips[occulyr]->Fill(expp, -1.0);
			for (int ix=0; ix<ypts[occulyr].size(); ix++) {
			  if (abs(expp - ypts[occulyr][ix])>0) {//choose all of them and then use maximum value before plot
			    h_ytcorstrips[occulyr]->Fill(expp, ypts[occulyr][ix]);
			  }
			}

			for (int ix=0; ix<xpts[occulyr].size(); ix++) {
			  h_yxtcorstrips[occulyr]->Fill(expp, xpts[occulyr][ix]);
			}
		      }
		    }
		  } //if (iiter==nmxiter-1)
		} // else of if (occulyr >=nlayer)
              } //if (nytime>=nmnhits/*-ntcor*/ && nytfail==0)
	    } // if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0 && nxfail==0 && nyfail==0); Position fit cut
	    if (isTiming && isfill &&  nxfail==0 &&  nyfail==0 && nxtime>2 && nytime>2) {

	      if (isfill && nytfail==0 && nytime >nmnhits && yt0chi2/(nytime-2)<mxtimechisq) {
		if (nxtfail==0 && nxtime >nmnhits && xt0chi2/(nxtime-2)<mxtimechisq ) {
		  for(int ij=0;ij<nlayer;ij++){

		    if (dist[ij]<0 || xtime[ij] <-50 || (!xusedtime[ij])) continue;
		    double dt1 = xtdev[ij]-biasinxtime[ij];

		    if (abs(dt1)>3*sigmarpc) continue;
		    for(int jk=0;jk<nlayer;jk++){

		      if (dist[jk]<0 || xtime[jk] <-50 || (!xusedtime[jk])) continue;
		      double dt2 = xtdev[jk]-biasinxtime[jk];
		      if (abs(dt2) >3*sigmarpc) continue;
		      deltatcov2[ij][jk] += dt1*dt2;
		      deltatCount2[ij][jk] +=1 ;
		  }
		  }// ij
		}

		for(int ij=0;ij<nlayer;ij++) {
		  if (dist[ij]<0 || ytime[ij] <-50 || (!yusedtime[ij])) continue;
		  double dt1 = ytdev[ij]-biasinytime[ij];
		  if (abs(dt1)>3*sigmarpc || (!yusedtime[ij])) continue;
		  for(int jk=0;jk<nlayer;jk++){

		    if (dist[jk]<0 || ytime[jk] <-50 || (!yusedtime[jk]) ) continue;
		    //                                      if(jk != ij )continue;
		    double dt2 = ytdev[jk]-biasinytime[jk];
		    if (abs(dt2) >3*sigmarpc) continue;
		    deltatcovy[ij][jk] += dt1*dt2;
		    deltatCounty[ij][jk] +=1 ;
		  }
		}// ij
	      } //if (isfill && nytfail==0 && nytime >nmnhits && yt0chi2/(nytime-2)<mxtimechisq)

	      txslop = -1./timexslope/cval;
	      tyslop = -1./timeyslope/cval;

	      //  if (nxtime>3 && nxtfail==0 && nytime>3 && nytfail==0) T2->Fill();

	      //	      timexinters= xc0inters;
	      //	      timeyinters= yc0inters;
	      //	      errtimexinters = sqrt(errcst_tx);
	      errtimexslope =sqrt(errlin_tx);

	      //	      errtimeyinters = sqrt(errcst_ty);
	      //	      errtimeyslope =sqrt(errlin_ty);
	      // Filling rootupule
	      //	      cout<<"nxtime "<< nytime<<endl;
	      //if (nxtime>11 || nytime>11 || gRandom->Uniform()<1.e-5) T2->Fill(); //07/09/2011
            }

	    if (ntcor==0) {
	      passed_strip[iiter]->Fill(Nx, 1.);
	      passed_strip[iiter]->Fill(nlayer+2+Ny, 1.);
	      passed_strip[iiter]->Fill(2*(nlayer+2)+nxtime, 1.);
	      passed_strip[iiter]->Fill(3*(nlayer+2)+nytime, 1.);
	    }
	   // cout<<iev<<"		"<<"XXXXXXXXXXXXXXXXXXXXXXX"<<endl;
	    //--------------------- Don't change anything here --------------
	    if((iev%100000)==1)
	      { cout<<"Processed " <<ntotal<<" ("<<iev<<") events so far for iteration# "<< iiter<<" occu "<<occulyr<<endl; }

	    //#if ndefined(MONTECARLO) && defined(NOISE_SIM)
#if defined(NOISE_SIM)
	    //#ifndef MONTECARLO
	    //Generate root files for noise and subsequently use for simulation
	    if (isfill && nxfail==0 && Nx >7 && xchi2<50 && nyfail==0 && Ny >7 && ychi2<50) {
	      NoisefileOut->cd();

	      //Copy isolated events
	      bool isxnoise = false;
	      bool isynoise = false;
	      vector<int> xx2ptsall[nlayer];
	      vector<int> yy2ptsall[nlayer];
	      nseltot++;

	      isxnoise = true;
	      for (int ij=0; ij<nlayer; ij++) {
		for (int ix=0; ix<xptsall[ij].size(); ix++) {
		  xlayer_seloccu[ij]->Fill(xptsall[ij][ix]);
		  if (xptsall[ij].size()<4) {xlayer_sel2occu[ij]->Fill(xptsall[ij][ix]);}

		  if (abs(xextloc[ij] - xptsall[ij][ix])>2.5) {
		    xx2ptsall[ij].push_back(xptsall[ij][ix]);
		    xlayer_noiseoccu[ij]->Fill(xptsall[ij][ix]);
		  }
		}
	      }

	      isynoise = true;
	      for (int ij=0; ij<nlayer; ij++) {
		for (int iy=0; iy<yptsall[ij].size(); iy++) {
		  ylayer_seloccu[ij]->Fill(yptsall[ij][iy]);
		  if (yptsall[ij].size()<4) {ylayer_sel2occu[ij]->Fill(yptsall[ij][iy]);}

		  if (abs(yextloc[ij] - yptsall[ij][iy])>2.5) {
		    yy2ptsall[ij].push_back(yptsall[ij][iy]);
		    ylayer_noiseoccu[ij]->Fill(yptsall[ij][iy]);
		  }
		}
	      }

	      for (int ij=0; ij<nlayer; ij++) {
		for (int ix=0; ix<xptsall[ij].size(); ix++) {
		  for (int iy=0; iy<yptsall[ij].size(); iy++) {
		    raw_seloccu[ij]->Fill(xptsall[ij][ix], yptsall[ij][iy]);
		  }
		}

		for (int ix=0; ix<xx2ptsall[ij].size(); ix++) {
		  for (int iy=0; iy<yy2ptsall[ij].size(); iy++) {
		    raw_noiseoccu[ij]->Fill(xx2ptsall[ij][ix], yy2ptsall[ij][iy]);
		  }
		}
	      }

	      //Store hits after removing muon hits, simple noise
	      if (isxnoise && isynoise) {

		int nxmulall = 0; int nymulall = 0;
		//	      int xhitsall[nlayer]={0};
		for (int ij=0; ij<nlayer; ij++) {
		  nxmulall += xptsall[ij].size();
		  nymulall += yptsall[ij].size();
		}
		n_multall[0]->Fill(nxmulall);
		n_multall[1]->Fill(nymulall);

		for (int ij=0; ij<nlayer; ij++) {
		  //for (int jk=0; jk<=xptsall[ij].size(); jk++) {
		  total_vs_indmulall[0][ij]->Fill(nxmulall, xptsall[ij].size()); //check this
		  //}
		  //for (int jk=0; jk<=yptsall[ij].size(); jk++) {
		  total_vs_indmulall[1][ij]->Fill(nymulall,yptsall[ij].size()); //check this
		  //}
		}

		for (int ij=0; ij<nlayer; ij++) {
		  layer_multall[0][ij]->Fill(xptsall[ij].size());
		  layer_multall[1][ij]->Fill(yptsall[ij].size());

		  for (int jk=0; jk<xptsall[ij].size(); jk++) { // X
		    layer_hitsall[0][ij]->Fill(xptsall[ij][jk]);
		    hist_correlall[0][ij]->Fill(xptsall[ij][jk], xptsall[ij][jk]);
		    for (int kl=jk+1; kl<xptsall[ij].size(); kl++) {
		      hist_correlall[0][ij]->Fill(xptsall[ij][jk], xptsall[ij][kl]);
		    }
		  }

		  for (int jk=0; jk<yptsall[ij].size(); jk++) { // Y
		    layer_hitsall[1][ij]->Fill(yptsall[ij][jk]);
		    hist_correlall[1][ij]->Fill(yptsall[ij][jk], yptsall[ij][jk]);
		    for (int kl=jk+1; kl<yptsall[ij].size(); kl++) {
		      hist_correlall[1][ij]->Fill(yptsall[ij][jk], yptsall[ij][kl]);
		    }
		  }

		  for (int jk=ij+1; jk<nlayer; jk++) {
		    mult_correlall[0][jk]->Fill(xptsall[ij].size(), xptsall[jk].size());
		    mult_correlall[1][jk]->Fill(yptsall[ij].size(), yptsall[jk].size());
		  }

		  if (xptsall[ij].size()>0) {layer_timeall[0][ij]->Fill(timesx[ij]);}
		  if (yptsall[ij].size()>0) {layer_timeall[1][ij]->Fill(timesy[ij]);}
		}

		//Hits excluding muon trajectory

		nxmulall = nymulall = 0;
		ntot_uncsim++;
		//	      for (int ij=0; ij<nlayer; ij++) { xhitsall[ij]=0;}
		for (int ij=0; ij<nlayer; ij++) {
		  //		xhits[ij] = xptsall[ij].size();
		  //		yhits[ij] = yptsall[ij].size();
		  //		nmulall += xhitsall[ij] = xptsall[ij].size()  + yptsall[ij].size();
		  nxmulall += xx2ptsall[ij].size();
		  nymulall += yy2ptsall[ij].size();
		}

		n_mult[0]->Fill(nxmulall);
		n_mult[1]->Fill(nymulall);

		for (int ij=0; ij<nlayer; ij++) {
		  //for (int jk=0; jk<=xx2ptsall[ij].size(); jk++) {
		  total_vs_indmul[0][ij]->Fill(nxmulall, xx2ptsall[ij].size()); //check this
		  //}
		  //for (int jk=0; jk<=yy2ptsall[ij].size(); jk++) {
		  total_vs_indmul[1][ij]->Fill(nymulall, yy2ptsall[ij].size()); //check this
		  //}
		}

		for (int ij=0; ij<nlayer; ij++) {
		  layer_mult[0][ij]->Fill(xx2ptsall[ij].size());
		  layer_mult[1][ij]->Fill(yy2ptsall[ij].size());

		  for (int jk=0; jk<xx2ptsall[ij].size(); jk++) { // X&&Y
		    layer_hits[0][ij]->Fill(xx2ptsall[ij][jk]);
		    hist_correl[0][ij]->Fill(xx2ptsall[ij][jk], xx2ptsall[ij][jk]);
		    for (int kl=jk+1; kl<xx2ptsall[ij].size(); kl++) {
		      hist_correl[0][ij]->Fill(xx2ptsall[ij][jk], xx2ptsall[ij][kl]);
		    }
		  }
		  for (int jk=0; jk<yy2ptsall[ij].size(); jk++) { // X&&Y
		    layer_hits[1][ij]->Fill(yy2ptsall[ij][jk]);
		    hist_correl[1][ij]->Fill(yy2ptsall[ij][jk], yy2ptsall[ij][jk]);
		    for (int kl=jk+1; kl<yy2ptsall[ij].size(); kl++) {
		      hist_correl[1][ij]->Fill(yy2ptsall[ij][jk], yy2ptsall[ij][kl]);
		    }
		  }

		  for (int jk=ij+1; jk<nlayer; jk++) {
		    mult_correl[0][jk]->Fill(xx2ptsall[ij].size(), xx2ptsall[jk].size());
		    mult_correl[1][jk]->Fill(yy2ptsall[ij].size(), yy2ptsall[jk].size());
		  }

		  if (xx2ptsall[ij].size()>0) {layer_time[0][ij]->Fill(timesx[ij]);}
		  if (yy2ptsall[ij].size()>0) {layer_time[1][ij]->Fill(timesy[ij]);}
		}
	      }
	    }
	    //#endif  // #if ndefined(MONTECARLO) && defined(NOISE_SIM)
#endif

	  //if(iev == 0) { init_time = event->Evetime->GetSec(); }
	      //if(iev == nentry-1) { end_time = event->Evetime->GetSec();
			  //total_time  =  end_time -init_time ;
			  //cout<<"total_time"<<"			"<<total_time<<endl; }
	    //if(ntcor==1 && iiter == nmxiter-1) {T4->Fill();}// jim jim enable line
	}  // EVENT loop

	  if (filloccu) {
	    cout<<"missing hits from X trigger layer"<<"   "<<nentry<<"    "<<nx4<<"  "<<nx5<<"   "<<nx6<<"  "<<nx7<<endl;
	    cout<<"missing hits from Y trigger layer"<<"   "<<nentry<<"    "<<ny4<<"  "<<ny5<<"   "<<ny6<<"  "<<ny7<<endl;
	    cout<<"missing hits from X and Y trigger layer"<<"   "<<nentry<<"    "<<nxy4<<"  "<<nxy5<<"   "<<nxy6<<"  "<<nxy7<<endl;
	    cout<<"missing hits from X or Y trigger layer"<<"   "<<nentry<<"    "<<nxory4<<"  "<<nxory5<<"   "<<nxory6<<"  "<<nxory7<<endl;
	    cout<<"missing hits from X and Y final trigger"<<"   "<<nentry<<"    "<<nxytrig<<endl;
	    cout<<"missing TDC from X and Y trigger layer"<<"    "<<nentry<<"    "<<ntxytime2[0]<<"   "<<ntxytime2[1]<<"   "<<ntxytime2[2]<<"   "<<ntxytime2[3]<<endl;
	    cout<<"missing TDC ref trig layer"<<"    "<<nentry<<"    "<<ntxytref[0]<<"   "<<ntxytref[1]<<"   "<<ntxytref[2]<<"   "<<ntxytref[3]<<endl;

	    cout<<"6 combination of counter"<<endl;
	    cout<<"Latch"<<"    "<<"Ref"<<"   "<<"TDC"<<endl;
	    for(int jj=0;jj<nlayer;jj++) {
	      cout<<jj<<"   "<<ntxext[jj][0]<<"  "<<ntxref[jj][0]<<"   "<<nttxext[jj][0]<<"   "<<ntyext[jj][0]<<"  "<<ntyref[jj][0]<<"   "<<nttyext[jj][0]<<endl;
	    }
	    cout<<endl;
	    cout<<"Latch"<<"    "<<"TDC"<<"   "<<"Ref"<<endl;
	     for(int jj=0;jj<nlayer;jj++) {
	      cout<<jj<<"   "<<ntxext[jj][1]<<"  "<<nttxext[jj][1]<<"   "<<ntxref[jj][1]<<"   "<<ntyext[jj][1]<<"  "<<nttyext[jj][1]<<"   "<<ntyref[jj][1]<<endl;
	    }
	     cout<<endl;
	     cout<<"Ref"<<"    "<<"Latch"<<"   "<<"TDC"<<endl;
	     for(int jj=0;jj<nlayer;jj++) {
	      cout<<jj<<"   "<<ntxref[jj][2]<<"  "<<ntxext[jj][2]<<"   "<<nttxext[jj][2]<<"   "<<ntyref[jj][2]<<"  "<<ntyext[jj][2]<<"   "<<nttyext[jj][2]<<endl;
	    }

	      cout<<endl;
	      cout<<"Ref"<<"    "<<"TDC"<<"   "<<"Latch"<<endl;
	       for(int jj=0;jj<nlayer;jj++) {
	      cout<<jj<<"   "<<ntxref[jj][3]<<"  "<<nttxext[jj][3]<<"   "<<ntxext[jj][3]<<"   "<<ntyref[jj][3]<<"  "<<nttyext[jj][3]<<"   "<<ntyext[jj][3]<<endl;
	    }
	       cout<<endl;
	       cout<<"TDC"<<"    "<<"Latch"<<"   "<<"Ref"<<endl;
	       for(int jj=0;jj<nlayer;jj++) {
	      cout<<jj<<"   "<<nttxext[jj][4]<<"  "<<ntxext[jj][4]<<"   "<<ntxref[jj][4]<<"   "<<nttyext[jj][4]<<"  "<<ntyext[jj][4]<<"   "<<ntyref[jj][4]<<endl;
	    }
	       cout<<endl;
	        cout<<"TDC"<<"    "<<"Ref"<<"   "<<"Latch"<<endl;

		 for(int jj=0;jj<nlayer;jj++) {
		   cout<<jj<<"   "<<nttxext[jj][5]<<"  "<<ntxref[jj][5]<<"   "<<ntxext[jj][5]<<"   "<<nttyext[jj][5]<<"  "<<ntyref[jj][5]<<"   "<<ntyext[jj][5]<<endl;
	    }
		 cout<<endl;

	    cout<<"TDC reference for all layers"<<endl;
	    for(int jj=0;jj<nlayer;jj++) {
	       for(int kk=0;kk<8;kk++) {
		 cout<<jj<<"    "<<kk<<"    "<<ntotxext[jj][kk]<<"    "<<ntxxref[jj][kk]<<"   "<<ntotyext[jj][kk]<<"    "<<ntyyref[jj][kk]<<endl;
	       }
	       }

	    cout<<"TDC reference vs TDC hit for all layers"<<endl;
	    for(int jj=0;jj<nlayer;jj++) {
	      cout<<jj<<"    "<<ntxyttref[jj]<<"   "<<ntxreftdc[jj]<<"    "<<ntyreftdc[jj]<<endl;
	    }
	    cout<<"latched vs tdc counts"<<endl;
	    for(int jj=0;jj<nlayer;jj++) {
	      for(int kk=0;kk<8;kk++) {
		cout<<jj<<"   "<<kk<<"   "<<ntotxext[jj][kk]<<"    "<<ntxtime[jj][kk]<<"   "<<ntotyext[jj][kk]<<"    "<<ntytime[jj][kk]<<endl;
	      }
	    }

	    cout<<"tdc vs latched"<<endl;
	    for(int jj=0;jj<nlayer;jj++) {
	      for(int kk=0;kk<8;kk++) {
		cout<<jj<<"   "<<kk<<"   "<<ntxtime2[jj][kk]<<"    "<<ntotxext2[jj][kk]<<"   "<<ntytime2[jj][kk]<<"   "<<ntotyext2[jj][kk]<<"    "<<endl;
	      }
	    }
	  }



	  fileIn->cd();
#ifdef MONTECARLO
	  delete T1;
#else

	  // #elif DLEVEL
	  delete event_tree;
#if CAUDATA1
	  delete cau_tree;
#endif
#endif
	  delete fileIn;
	} // while(!(file_db.eof()))

	if (ntotal<1) ntotal=1;
	if (nseltot<1) nseltot = 1;
	file_db.close();

	if(isfill) {
	  for(int ij=0;ij<nstrip/2;ij++) {
	    for(int jk=0;jk<nstrip/2;jk++) {
	      double scatmean1 = pixel_scatang1[ij][jk]->GetMean();
	      // double scatmed1 =  Median(pixel_scatang1[ij][jk]);
	       pixel_scatmean1->SetBinContent(ij,jk,scatmean1);

	      double scatmean = pixel_scatang[ij][jk]->GetMean();
	      // double scatmed =  Median(pixel_scatang[ij][jk]);
	       pixel_scatmean->SetBinContent(ij,jk,scatmean);
	    }
	  }
	}



	//#if ndefined(MONTECARLO) && defined(NOISE_SIM)
#if defined(NOISE_SIM)
	//#ifndef MONTECARLO
	if (isfill) {
	  NoisefileOut->cd();
	  // Normalise layer multiplity for each individual value of total multiplicity
	  for (int ixy=0; ixy<2; ixy++) {
	    for (int ij=0; ij<nlayer; ij++) {
	      for (int jk=0; jk<total_vs_indmulall[ixy][ij]->GetNbinsX(); jk++) {
		double ntt = 0.0;
		for (int kl=0; kl<total_vs_indmulall[ixy][ij]->GetNbinsY(); kl++) {
		  ntt += total_vs_indmulall[ixy][ij]->GetBinContent(jk+1,kl+1);
		}
		//		double ntt = total_vs_indmulall[ixy][ij]->GetBinContent(jk+1,1);
		if (ntt<1.0) ntt=1.0;
		for (int kl=0; kl<total_vs_indmulall[ixy][ij]->GetNbinsY(); kl++) {
		  total_vs_indmulall[ixy][ij]->SetBinContent(jk+1, kl+1, total_vs_indmulall[ixy][ij]->GetBinContent(jk+1,kl+1)/ntt); // (max(1, ntot_uncsim))/*ntt*/);
		}
	      }

	      for (int jk=0; jk<total_vs_indmul[ixy][ij]->GetNbinsX(); jk++) {
		double ntt = 0;
		for (int kl=0; kl<total_vs_indmul[ixy][ij]->GetNbinsY(); kl++) {
		  ntt += total_vs_indmul[ixy][ij]->GetBinContent(jk+1,kl+1);
		}
		if (ntt<1.0) ntt=1.0;
		for (int kl=0; kl<total_vs_indmul[ixy][ij]->GetNbinsY(); kl++) {
		  total_vs_indmul[ixy][ij]->SetBinContent(jk+1, kl+1, total_vs_indmul[ixy][ij]->GetBinContent(jk+1,kl+1)/ntt); //(max(1, ntot_uncsim))/*ntt*/);
		}
	      }
	    }
	  }

	  double scal = 1./(max(1, ntot_uncsim));
	  double scalunc = 1./(max(1, ntot_uncsim));
	  cout <<"scal "<<scal<<endl;

	  for (int ixy=0; ixy<2; ixy++) {
	    n_multall[ixy]->Scale(scal);
	    n_mult[ixy]->Scale(scalunc);

	    for (int ij=0; ij<nlayer; ij++) {
	      layer_hitsall[ixy][ij]->Scale(scal);
	      layer_multall[ixy][ij]->Scale(scal);
	      layer_timeall[ixy][ij]->Scale(scal);

	      mult_correlall[ixy][ij]->Scale(scal);
	      hist_correlall[ixy][ij]->Scale(scal);

	      layer_hits[ixy][ij]->Scale(scalunc);
	      layer_mult[ixy][ij]->Scale(scalunc);
	      layer_time[ixy][ij]->Scale(scalunc);

	      mult_correl[ixy][ij]->Scale(scalunc);
	      hist_correl[ixy][ij]->Scale(scalunc);
	    }
	  }
	  NoisefileOut->Write();
	  NoisefileOut->Close();
	}
	//#endif
#endif
	ps.NewPage();
	fileOut->cd();
	if ((filloccu) || (isfill)) {
	  //	  filloccu = false;

	  TF1* fitxy[nlayer]={0};
	  TH1F* xy_occu[nlayer]={0};

	  double xycount[2][nlayer]={0}; //Sum of X and Y strip for all layer

	  gStyle->SetOptStat(0);
	  gStyle->SetStatW(.36);
	  gStyle->SetStatH(.20);
	  gStyle->SetStatY(.99);
	  gStyle->SetStatX(.99);
	  gStyle->SetOptFit(0);
	  gStyle->SetOptLogy(1);
	  gStyle->SetPadBottomMargin(0.11);
	  gStyle->SetPadTopMargin(0.07);
	  gStyle->SetPadLeftMargin(0.10);
	  gStyle->SetPadRightMargin(0.11);
	  gStyle->SetPaintTextFormat("g");
	  gStyle->SetOptTitle(1);

	  for (int ixy=0; ixy<8; ixy++) {
	    if ((!filloccu) && ixy<=1) continue;
	    if ((!isfill) && ixy>1) continue;

	    ps.NewPage();
	    gStyle->SetPadRightMargin(0.01);
	    gStyle->SetOptLogy(0);

	    TCanvas *c4=new TCanvas ("c4","Strip occupancy",500,700);
	    c4->Divide(3,4);

	    int nxx = (ixy<2) ? max(1, ntotal) : max(1, nseltot);
	    double ascl=100./nxx;

            for (int il=0; il<nlayer; il++) {
	      switch (ixy) {
	      case 0 : xy_occu[il] = (TH1F*)xlayer_occu[il]->Clone(); break;
	      case 1 : xy_occu[il] = (TH1F*)ylayer_occu[il]->Clone(); break;
	      case 2 : xy_occu[il] = (TH1F*)xlayer_seloccu[il]->Clone(); break;
	      case 3 : xy_occu[il] = (TH1F*)ylayer_seloccu[il]->Clone(); break;
	      case 4 : xy_occu[il] = (TH1F*)xlayer_sel2occu[il]->Clone(); break;
	      case 5 : xy_occu[il] = (TH1F*)ylayer_sel2occu[il]->Clone(); break;
	      case 6 : xy_occu[il] = (TH1F*)xlayer_noiseoccu[il]->Clone(); break;
	      case 7 : xy_occu[il] = (TH1F*)ylayer_noiseoccu[il]->Clone(); break;
	      default : xy_occu[il] = (TH1F*)xlayer_occu[il]->Clone(); break;
	      }

	      xy_occu[il]->Scale(ascl);
	      // store strips with dead/close to channels
	      if (ixy<=1) {
		for (int jkx=0; jkx<xy_occu[il]->GetNbinsX(); jkx++) {
		  if (xy_occu[il]->GetBinContent(jkx+1)<1.e-3) {deadchannel[ixy][il].push_back(jkx);}
		}
	      }




	      c4->cd(il+1);

	      sprintf(name, "fitxy_%i_%i", il, ixy);
#ifdef ONEMBYONEM
	      fitxy[il] = new TF1(name, fitspec, 0, nusedstrip, 4);
#else
#if defined(MADURAIINTM1)
	      if (ixy==0) {
		fitxy[il] = new TF1(name, fitspec, 0, nusedstrip/2, 4);
	      } else {
		if (il>=5 && il<=8) {
		  fitxy[il] = new TF1(name, fitspec, 0, nusedstrip/2, 4);
		} else {
		  fitxy[il] = new TF1(name, fitspec, 0, (ixy==0) ? lastXstrip : lastYstrip, 4);
		}
	      }
#endif
#ifdef MADURAIINTM2
	      if (ixy==0) {
		fitxy[il] = new TF1(name, fitspec, 0, nusedstrip/2, 4);
	      } else {
		if (il>=5 && il<=8) {
		  fitxy[il] = new TF1(name, fitspec, 0, nusedstrip/2, 4);
		} else {
		  fitxy[il] = new TF1(name, fitspec, 0, (ixy==0) ? lastXstrip : lastYstrip, 4);
		}
	      }
#endif
#ifdef MADURAIINTM3
	      if (il>=5 && il<=8) {
		fitxy[il] = new TF1(name, fitspec, 1, nusedstrip/2-2, 4);
	      } else {
		fitxy[il] = new TF1(name, fitspec, 1, (ixy==0) ? lastXstrip : lastYstrip, 4);
	      }
#endif
#ifdef MADURAIINTM4
	      if(ixy<2) {fitxy[il] = new TF1(name, fitspec, 1, (ixy==0 ) ? lastXstrip : lastYstrip, 4); }
#endif

#endif
	      double hgh ;//= 0.5*xy_occu[il]->GetMaximum();
	      double parx[4] ={hgh, 0., hgh, 1.};
	      if(ixy<2) {
		hgh = 0.5*xy_occu[il]->GetMaximum();
		fitxy[il]->SetParLimits(0, 0.0, 2.0*hgh);
		fitxy[il]->SetParLimits(1, -1.0, 1.0);
		fitxy[il]->SetParLimits(2, 0.5*hgh, 4.0*hgh);

		fitxy[il]->FixParameter(3, 1.0);
		fitxy[il]->SetParameters(parx);
	      }

	      xy_occu[il]->SetMarkerStyle(24);
	      xy_occu[il]->SetMarkerSize(0.8);
	      xy_occu[il]->GetXaxis()->SetTitle("Strip No");
	      xy_occu[il]->GetXaxis()->SetTitleOffset(.8);
	      xy_occu[il]->GetXaxis()->SetTitleSize(.06);
	      xy_occu[il]->GetXaxis()->CenterTitle();
	      xy_occu[il]->GetYaxis()->SetTitle("Occupancy (%)");
	      xy_occu[il]->GetYaxis()->SetTitleOffset(.8);
	      xy_occu[il]->GetYaxis()->SetTitleSize(.06);
	      xy_occu[il]->GetYaxis()->CenterTitle();

	      xy_occu[il]->GetXaxis()->SetLabelSize(0.06);
	      xy_occu[il]->GetYaxis()->SetLabelSize(0.06);
	      xy_occu[il]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
	  //     if(ixy<2) {
    //           xy_occu[il]->Fit(name, "WW:BRQ");
	  //     fitxy[il]->GetParameters(parx);
    //
	  //     double xx[2]={0,0};
    //           for (int istr=0; istr<nusedstrip; istr++) {
		// xx[0] = xy_occu[il]->GetBinCenter(istr+1);
		// double yy = xy_occu[il]->GetBinContent(istr+1);
		// // 		cout <<"ixy "<< ixy<<" "<<il<<" "<<ntotal<<" "<<istr<<" "<<xx[0]<<" "<<yy<<" "<<xy_occu[il]->GetBinContent(istr+1)<<endl;
		// xycount[ixy][il] +=yy;
    //
	  //     }
	  //     } else {
		// xy_occu[il]->Draw("");
	  //     }

      xy_occu[il]->Draw("HIST");

	    }

	    c4->Update();
	    if (c4) { delete c4; c4=0;}

            for (int il=0; il<nlayer; il++) {
	      if (fitxy[il]) {delete fitxy[il]; fitxy[il]=0;}
	      if (xy_occu[il]) {delete xy_occu[il]; xy_occu[il]=0;}
	    }
	  } // for (int ixy=0; ixy<2; ixy++)
	}

	if (isfill) {
	  ps.NewPage();
	  gStyle->SetOptTitle(1);
	  gStyle->SetPadLeftMargin(0.05);
	  gStyle->SetPadRightMargin(0.12);
	  gStyle->SetPadTopMargin(0.10);
	  gStyle->SetPadBottomMargin(0.12);

	  gStyle->SetOptLogy(0);
	  gStyle->SetOptLogz(0);
	  //	  gStyle->SetOptStat(1110);
	  //	  gStyle->SetOptFit(100);
	  gStyle->SetStatW(.36); //40);
	  gStyle->SetStatH(.28); //30);
	  gStyle->SetTitleFontSize(0.07);
	  gStyle->SetOptStat(0);
	  gStyle->SetOptFit(0);

	  TCanvas *c4z=new TCanvas ("c4z","Strip occupancy",500,700);
	  c4z->Divide(3,4);

	   for (int ij=0; ij<nlayer; ij++) {
	     c4z->cd(ij+1);
	     gPad->SetLeftMargin(0.12);
	     gPad->SetRightMargin(0.15);
	     gPad->SetBottomMargin(0.12);
	     gPad->SetTopMargin(0.10);
	     raw_seloccu[ij]->Scale(100./nseltot);
	     raw_seloccu[ij]->GetXaxis()->SetLabelSize(.07);
	     raw_seloccu[ij]->GetXaxis()->SetTitle("X-strip");
	     raw_seloccu[ij]->GetXaxis()->CenterTitle();
	     raw_seloccu[ij]->GetXaxis()->SetTitleSize(0.07);
	     raw_seloccu[ij]->GetXaxis()->SetTitleOffset(.8);

	     raw_seloccu[ij]->GetYaxis()->SetTitle("Y-strip");
	     raw_seloccu[ij]->GetYaxis()->CenterTitle();
	     raw_seloccu[ij]->GetYaxis()->SetTitleSize(0.06);
	     raw_seloccu[ij]->GetYaxis()->SetTitleOffset(.78);
	     raw_seloccu[ij]->GetYaxis()->SetLabelSize(.06);
	     raw_seloccu[ij]->GetXaxis()->SetRangeUser(-1., nusedstrip);
	     raw_seloccu[ij]->GetYaxis()->SetRangeUser(-1., nusedstrip);
	     raw_seloccu[ij]->GetZaxis()->SetLabelSize(.06);
	    raw_seloccu[ij]->Draw("colz");
	   }
	   c4z->Update();

	   ps.NewPage();

	   for (int ij=0; ij<nlayer; ij++) {
	     c4z->cd(ij+1);
	     gPad->SetLeftMargin(0.12);
	     gPad->SetRightMargin(0.15);
	     gPad->SetBottomMargin(0.12);
	     gPad->SetTopMargin(0.10);
	     raw_noiseoccu[ij]->Scale(100./nseltot);
	     raw_noiseoccu[ij]->GetXaxis()->SetLabelSize(.07);
	     raw_noiseoccu[ij]->GetXaxis()->SetTitle("X-strip");
	     raw_noiseoccu[ij]->GetXaxis()->CenterTitle();
	     raw_noiseoccu[ij]->GetXaxis()->SetTitleSize(0.07);
	     raw_noiseoccu[ij]->GetXaxis()->SetTitleOffset(.8);

	     raw_noiseoccu[ij]->GetYaxis()->SetTitle("Y-strip");
	     raw_noiseoccu[ij]->GetYaxis()->CenterTitle();
	     raw_noiseoccu[ij]->GetYaxis()->SetTitleSize(0.06);
	     raw_noiseoccu[ij]->GetYaxis()->SetTitleOffset(.78);
	     raw_noiseoccu[ij]->GetYaxis()->SetLabelSize(.06);
	     raw_noiseoccu[ij]->GetXaxis()->SetRangeUser(-1., nusedstrip);
	     raw_noiseoccu[ij]->GetYaxis()->SetRangeUser(-1., nusedstrip);
	     raw_noiseoccu[ij]->GetZaxis()->SetLabelSize(.06);
	     raw_noiseoccu[ij]->Draw("colz");
	   }
	   c4z->Update();
	   if (c4z) { delete c4z; c4z=0;}
	}

	if (filloccu) {
	  filloccu = false;
	  ps.NewPage();
	  gStyle->SetOptStat(1110);
	  gStyle->SetOptTitle(1);
	  gStyle->SetOptLogy(1);
	  gStyle->SetStatW(.40); //40);
	  gStyle->SetStatH(.24); //30);
	  gStyle->SetTitleFontSize(0.07);
	  gStyle->SetPadLeftMargin(0.09);
	  gStyle->SetPadBottomMargin(0.06);
	  gStyle->SetPadTopMargin(0.09); //(0.03);
	  gStyle->SetPadRightMargin(0.01);
	  TCanvas* c4c = new TCanvas("c4c", "c4c", 700, 900);
	  c4c->Divide(3,4);

	  //  ps.NewPage();
	  gStyle->SetOptStat(1100);
	  gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    xlayer_allmult[ij]->Scale(1./ntotal);
	    xlayer_allmult[ij]->GetXaxis()->SetLabelSize(0.054);
	    xlayer_allmult[ij]->GetYaxis()->SetLabelSize(0.054);
	    xlayer_allmult[ij]->GetYaxis()->SetRangeUser(1.e-5, 0.8);
	    xlayer_allmult[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
	    xlayer_allmult[ij]->SetLineColor(1); xlayer_allmult[ij]->Draw();
	  }

	  c4c->Update();
	  gStyle->SetStatY(.87); gStyle->SetStatTextColor(2);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    ylayer_allmult[ij]->Scale(1./ntotal);
	    ylayer_allmult[ij]->SetLineColor(2); ylayer_allmult[ij]->Draw("sames");
	  }
	  c4c->Update();
	  gStyle->SetStatY(.75); gStyle->SetStatTextColor(3);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    xlayer_mult[ij]->Scale(1./ntotal);
	    xlayer_mult[ij]->SetLineColor(3); xlayer_mult[ij]->Draw("sames");
	  }
	  c4c->Update();
	  gStyle->SetStatY(.64); gStyle->SetStatTextColor(4);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    ylayer_mult[ij]->Scale(1./ntotal);
	    ylayer_mult[ij]->SetLineColor(4); ylayer_mult[ij]->Draw("sames");
	  }

	  c4c->Update();

	  ps.NewPage();
	  gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    xlayer_allmumult[ij]->Scale(1./ntotal);
	    xlayer_allmumult[ij]->GetXaxis()->SetLabelSize(0.054);
	    xlayer_allmumult[ij]->GetYaxis()->SetLabelSize(0.054);
	    xlayer_allmumult[ij]->GetYaxis()->SetRangeUser(1.e-5, 0.8);
	    xlayer_allmumult[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
	    xlayer_allmumult[ij]->SetLineColor(1); xlayer_allmumult[ij]->Draw();
	  }
	  c4c->Update();
	  gStyle->SetStatY(.87); gStyle->SetStatTextColor(2);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    ylayer_allmumult[ij]->Scale(1./ntotal);
	    ylayer_allmumult[ij]->SetLineColor(2); ylayer_allmumult[ij]->Draw("sames");
	  }
	  c4c->Update();
	  gStyle->SetStatY(.75); gStyle->SetStatTextColor(3);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    xlayer_allmutimemult[ij]->Scale(1./ntotal);
	    xlayer_allmutimemult[ij]->SetLineColor(3); xlayer_allmutimemult[ij]->Draw("sames");
	  }
	  c4c->Update();
	  gStyle->SetStatY(.64); gStyle->SetStatTextColor(4);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    ylayer_allmutimemult[ij]->Scale(1./ntotal);
	    ylayer_allmutimemult[ij]->SetLineColor(4); ylayer_allmutimemult[ij]->Draw("sames");
	  }

	  c4c->Update();

	  ps.NewPage();
	  gStyle->SetStatTextColor(1); gStyle->SetStatY(.99);

	  gStyle->SetOptStat(1110);
	  gStyle->SetStatW(.36);
	  gStyle->SetStatH(.20);
	  gStyle->SetStatY(.99);
	  gStyle->SetStatX(.99);
	  gStyle->SetOptFit(0);
	  gStyle->SetOptLogy(1);
	  gStyle->SetPadBottomMargin(0.12);
	  gStyle->SetPadTopMargin(0.09);
	  gStyle->SetPadLeftMargin(0.12);
	  gStyle->SetPadRightMargin(0.01);
	  gStyle->SetPaintTextFormat("g");
	  gStyle->SetOptTitle(1);

#ifndef MONTECARLO
	  if (plot_level>90) {
	    TCanvas *c44=new TCanvas ("c44","Strip occupancy",500,700);
	    c44->Divide(3,4);
	    for (int ix=0; ix<=iset; ix++) {
	      if (ix>0) ps.NewPage();
	      gStyle->SetStatY(0.99); gStyle->SetStatTextColor(1);

	      for (int il=0; il<nlayer; il++) {
		c44->cd(il+1);
		strp_xmult_set[il][ix]->SetLineColor(1);
		strp_xmult_set[il][ix]->Scale(1./nCount);
		strp_xmult_set[il][ix]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
		strp_xmult_set[il][ix]->GetXaxis()->SetTitle("Multiplicity");
		strp_xmult_set[il][ix]->GetYaxis()->SetTitle("Normalised Entry");
		strp_xmult_set[il][ix]->GetXaxis()->SetTitleOffset(.75);
		strp_xmult_set[il][ix]->GetXaxis()->SetTitleSize(.06);
		strp_xmult_set[il][ix]->GetXaxis()->CenterTitle();
		strp_xmult_set[il][ix]->GetXaxis()->SetLabelSize(.06);
		strp_xmult_set[il][ix]->GetXaxis()->SetLabelOffset(-0.001);
		strp_xmult_set[il][ix]->GetYaxis()->SetTitleOffset(0.9);
		strp_xmult_set[il][ix]->GetYaxis()->SetTitleSize(.06);
		strp_xmult_set[il][ix]->GetYaxis()->CenterTitle();
		strp_xmult_set[il][ix]->GetYaxis()->SetLabelSize(.05);
		strp_xmult_set[il][ix]->GetYaxis()->SetLabelOffset(-0.001);
		strp_xmult_set[il][ix]->Draw();
	      }
	      c44->Update();
	      gStyle->SetStatY(0.80); gStyle->SetStatTextColor(2);
	      for (int il=0; il<nlayer; il++) {
		c44->cd(il+1);
		strp_ymult_set[il][ix]->SetLineColor(2);
		strp_ymult_set[il][ix]->Scale(1./nCount);
		strp_ymult_set[il][ix]->Draw("sames");
	      }
	      c44->Update();
	    }

	    //	  ps.NewPage();
	    gStyle->SetOptStat(0);  gStyle->SetTitleFontSize(0.05);
	    gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
	    for (int ix=0; ix<=iset; ix++) {
	      ps.NewPage();
	      c44->cd();
	      gPad->SetLogy(0);
	      gPad->SetLeftMargin(0.05);
	      gPad->SetRightMargin(0.07);
	      gStyle->SetPadTopMargin(0.05);

	      strp_count_set[ix]->Scale(100./nCount);
	      strp_count_set[ix]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
	      strp_count_set[ix]->GetXaxis()->SetLabelSize(.040);
	      strp_count_set[ix]->GetXaxis()->SetLabelOffset(0.003);
	      strp_count_set[ix]->GetYaxis()->SetLabelSize(.030);
	      strp_count_set[ix]->GetYaxis()->SetLabelOffset(0.003);
	      strp_count_set[ix]->SetMarkerSize(0.75);
	      strp_count_set[ix]->GetZaxis()->SetLabelSize(.030);
	      for (int ij=0; ij<nusedstrip; ij++) {
		strp_count_set[ix]->GetYaxis()->SetBinLabel(ij+1, labels[ij]);
	      }
	      strp_count_set[ix]->GetYaxis()->LabelsOption("h");

	      for (int ij=0; ij<2*nlayer; ij++) {
		strp_count_set[ix]->GetXaxis()->SetBinLabel(ij+1, xlabels[ij]);
	      }
	      strp_count_set[ix]->GetXaxis()->LabelsOption("v");
	      strp_count_set[ix]->Draw("colz:text25");
	      c44->Update();
	    }
	    c44->Clear();
	  }
#endif
	  ps.NewPage();
	  gStyle->SetOptTitle(1);
	  gPad->SetLeftMargin(0.05);
	  gStyle->SetPadRightMargin(0.12);
	  gStyle->SetPadTopMargin(0.10);
	  gStyle->SetPadBottomMargin(0.12);

	  gStyle->SetOptLogy(0);
	  gStyle->SetOptLogz(0);
	  //	  gStyle->SetOptStat(1110);
	  //	  gStyle->SetOptFit(100);
	  gStyle->SetStatW(.36); //40);
	  gStyle->SetStatH(.28); //30);
	  gStyle->SetTitleFontSize(0.07);

	  TCanvas* c4a = new TCanvas("c4a", "c4a", 700, 900);
	  c4a->Divide(3,4);
#ifndef MONTECARLO
	  if (plot_level>90) {
	    for (int iyy=0; iyy<4; iyy++) {
	      if (!isTiming && iyy>=2) continue;
	      for (int ix=0; ix<=iset; ix++) {
		TH1F* histxx[nlayer]={0};
		ps.NewPage();
		for (int ij=0; ij<nlayer; ij++) {
		  c4a->cd(ij+1);
		  switch(iyy) {
		  case 0 : histxx[ij] = (TH1F*)xlayer_reso_set[ij][ix]->Clone(); break;
		  case 1 : histxx[ij] = (TH1F*)ylayer_reso_set[ij][ix]->Clone(); break;
		  case 2 : histxx[ij] = (TH1F*)time_xreso_set[ij][ix]->Clone(); break;
		  case 3 : histxx[ij] = (TH1F*)time_yreso_set[ij][ix]->Clone(); break;
		  default : histxx[ij] = (TH1F*)xlayer_reso_set[ij][ix]->Clone(); break;
		  }
		  histxx[ij]->GetXaxis()->SetLabelSize(0.07);
		  histxx[ij]->GetXaxis()->SetTitle(histxx[ij]->GetName());
		  histxx[ij]->GetXaxis()->CenterTitle();
		  histxx[ij]->GetXaxis()->SetTitleSize(0.07);
		  histxx[ij]->GetXaxis()->SetTitleOffset(.9);

		  histxx[ij]->GetYaxis()->SetLabelSize(0.07);
		  TFitResultPtr ptr = histxx[ij]->Fit("gaus", "SQ");
		  Int_t fitStatus = ptr;
		  if (fitStatus==0) {
		    latex.DrawLatex(0.60, 0.84,Form("#scale[0.5]{%g/%i}", int(100*ptr->Chi2())/100., ptr->Ndf()));
		    latex.DrawLatex(0.72, 0.76,Form("#scale[0.5]{%g}", int(1000*ptr->Parameter(1))/1000.));
		    latex.DrawLatex(0.72, 0.68,Form("#scale[0.5]{%g}", int(1000*ptr->Parameter(2))/1000.));
		  }
		  latex.DrawLatex(0.72, 0.60, Form("#scale[0.5]{%g}", histxx[ij]->Integral()));
		  latex.DrawLatex(0.72, 0.52,Form("#scale[0.5]{%g}", int(1000*histxx[ij]->GetMean())/1000.));
		  latex.DrawLatex(0.72, 0.44,Form("#scale[0.5]{%g}", int(1000*histxx[ij]->GetRMS())/1000.));


		}
		c4a->Update();
		for (int ij=0; ij<nlayer; ij++) {
		  if (histxx[ij]) { delete histxx[ij]; histxx[ij]=0;}
		}
	      }
	    }
	  }
#endif

	  gStyle->SetOptStat(0);
	  gStyle->SetOptFit(0);

	  for (int ij=0; ij<nlayer; ij++) {
	    c4a->cd(ij+1);
	    raw_occu[ij]->Scale(100./ntotal);
	    raw_occu[ij]->GetXaxis()->SetLabelSize(.07);
	    raw_occu[ij]->GetXaxis()->SetTitle("X-strip");
	    raw_occu[ij]->GetXaxis()->CenterTitle();
	    raw_occu[ij]->GetXaxis()->SetTitleSize(0.07);
	    raw_occu[ij]->GetXaxis()->SetTitleOffset(.8);

	    raw_occu[ij]->GetYaxis()->SetTitle("Y-strip");
	    raw_occu[ij]->GetYaxis()->CenterTitle();
	    raw_occu[ij]->GetYaxis()->SetTitleSize(0.07);
	    raw_occu[ij]->GetYaxis()->SetTitleOffset(.78);
	    raw_occu[ij]->GetYaxis()->SetLabelSize(.07);
	    raw_occu[ij]->GetXaxis()->SetRangeUser(-1., nusedstrip);
	    raw_occu[ij]->GetYaxis()->SetRangeUser(-1., nusedstrip);
	    //	    raw_occu[ij]->GetZaxis()->SetLabelSize(.07);
	    raw_occu[ij]->Draw("colz");
	  }
	  c4a->Update();



	  if (plot_level>90) {
	    for (int ix=0; ix<iset; ix++) {
	      ps.NewPage();
	      for (int ij=0; ij<nlayer; ij++) {
		c4a->cd(ij+1);
		raw_occu_set[ij][ix]->Scale(100./nCount);
		raw_occu_set[ij][ix]->GetXaxis()->SetLabelSize(.07);
		raw_occu_set[ij][ix]->GetXaxis()->SetTitle("X-strip");
		raw_occu_set[ij][ix]->GetXaxis()->CenterTitle();
		raw_occu_set[ij][ix]->GetXaxis()->SetTitleSize(0.07);
		raw_occu_set[ij][ix]->GetXaxis()->SetTitleOffset(.78);

		raw_occu_set[ij][ix]->GetYaxis()->SetTitle("Y-strip");
		raw_occu_set[ij][ix]->GetYaxis()->CenterTitle();
		raw_occu_set[ij][ix]->GetYaxis()->SetTitleSize(0.07);
		raw_occu_set[ij][ix]->GetYaxis()->SetTitleOffset(.7);
		raw_occu_set[ij][ix]->GetYaxis()->SetLabelSize(.07);

		raw_occu_set[ij][ix]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
		raw_occu_set[ij][ix]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
		raw_occu_set[ij][ix]->Draw("colz");
	      }
	      c4a->Update();
	    }
	  }
	  ps.NewPage();
	  for (int ij=0; ij<nlayer; ij++) {
	    c4a->cd(ij+1);
	    gPad->SetLogz(1); //GMA
	    rawhits_corr_xymul[ij]->Scale(100./ntotal);
	    rawhits_corr_xymul[ij]->GetXaxis()->SetLabelSize(.07);
	    rawhits_corr_xymul[ij]->GetYaxis()->SetLabelSize(.07);
	    rawhits_corr_xymul[ij]->GetXaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
	    rawhits_corr_xymul[ij]->GetYaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
	    rawhits_corr_xymul[ij]->Draw("colz");
	  }
	  c4a->Update();
	  if(plot_level>90) {
	  int icol=0;
	  for (int ij=0; ij<nlayer; ij++) {
	    for (int jk=ij+1; jk<nlayer; jk++) {
	      icol++;
	      if (icol==1) {ps.NewPage();}
	      c4a->cd(2*icol-1);
	      gPad->SetLogz(1); //GMA
	      rawhits_xlay_corr_mul[ij][jk]->Scale(100./ntotal);
	      rawhits_xlay_corr_mul[ij][jk]->GetXaxis()->SetLabelSize(.07);
	      rawhits_xlay_corr_mul[ij][jk]->GetYaxis()->SetLabelSize(.07);
	      rawhits_xlay_corr_mul[ij][jk]->GetXaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
	      rawhits_xlay_corr_mul[ij][jk]->GetYaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
	      rawhits_xlay_corr_mul[ij][jk]->Draw("colz");

	      c4a->cd(2*icol);
	      gPad->SetLogz(1); //GMA
	      rawhits_ylay_corr_mul[ij][jk]->Scale(100./ntotal);
	      rawhits_ylay_corr_mul[ij][jk]->GetXaxis()->SetLabelSize(.07);
	      rawhits_ylay_corr_mul[ij][jk]->GetYaxis()->SetLabelSize(.07);
	      rawhits_ylay_corr_mul[ij][jk]->GetXaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
	      rawhits_ylay_corr_mul[ij][jk]->GetYaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
	      rawhits_ylay_corr_mul[ij][jk]->Draw("colz");
	      if (icol==6) {c4a->Update(); icol=0;}
	    }
	  }
	  if (icol!=0) { c4a->Update();}
	  }






	  /*
	  // Remove noisy layer (X or Y)
	  for (int il=0; il<nlayer; il++) {
	    if (xycount[0][il] > 2*xycount[1][il]) {
	      for (int istr=0; istr<nstrip; istr++) {
		posinuse[0][il][istr]=false;
	      }
	    }

	    if (xycount[1][il] > 2*xycount[0][il]) {
	      for (int istr=0; istr<nstrip; istr++) {
		posinuse[1][il][istr]=false;
	      }
	    }
	  }
	  */
#ifdef TIMESLOPE
	  //#ifndef MONTECARLO
	  if (isTiming) {

	    ps.NewPage();
	    gStyle->SetOptTitle(1);
	    gStyle->SetTitleOffset(-0.06,"XYZ");
	    gStyle->SetTitleFillColor(10);
	    gStyle->SetTitleFontSize(0.06);
	    gStyle->SetPadTopMargin(0.08);
	    gStyle->SetPadRightMargin(0.15);
	    gStyle->SetPadLeftMargin(0.12);
	    gStyle->SetPadBottomMargin(0.11);

	    TCanvas* c3a = new TCanvas("c3a", "c3a", 500, 700);
	    c3a->Divide(3,4);
	    for (int itdc=0; itdc<nTDCpLayer; itdc++) {
	       if (itdc>0) {ps.NewPage();}
	      for (int il=0; il<nlayer; il++) {
		c3a->cd(il+1);
		time_xslope[il][itdc]->GetXaxis()->SetLabelSize(.07);
		time_xslope[il][itdc]->GetYaxis()->SetLabelSize(.07);
		time_xslope[il][itdc]->GetXaxis()->SetTitle("Strip No");
		time_xslope[il][itdc]->GetXaxis()->SetTitleOffset(.85);
		time_xslope[il][itdc]->GetXaxis()->SetTitleSize(.06);
		time_xslope[il][itdc]->GetXaxis()->CenterTitle();
		time_xslope[il][itdc]->GetYaxis()->SetTitle("Measured Time (ns)");
		time_xslope[il][itdc]->GetYaxis()->SetTitleOffset(.8);
		time_xslope[il][itdc]->GetYaxis()->SetTitleSize(.06);
		time_xslope[il][itdc]->GetYaxis()->CenterTitle();

		time_xslope[il][itdc]->Draw("colz");
		time_xslope[il][itdc]->ProfileX()->Draw("same");
	      }
	      c3a->Update();

	      ps.NewPage();
	      for (int il=0; il<nlayer; il++) {
		c3a->cd(il+1);
		time_yslope[il][itdc]->GetXaxis()->SetLabelSize(.07);
		time_yslope[il][itdc]->GetYaxis()->SetLabelSize(.07);
		time_yslope[il][itdc]->GetXaxis()->SetTitle("Strip No");
		time_yslope[il][itdc]->GetXaxis()->SetTitleOffset(.6);
		time_yslope[il][itdc]->GetXaxis()->SetTitleSize(.06);
		time_yslope[il][itdc]->GetXaxis()->CenterTitle();
		time_yslope[il][itdc]->GetYaxis()->SetTitle("Occupancy (%)");
		time_yslope[il][itdc]->GetYaxis()->SetTitleOffset(.6);
		time_yslope[il][itdc]->GetYaxis()->SetTitleSize(.06);
		time_yslope[il][itdc]->Draw("colz");
		time_yslope[il][itdc]->ProfileX()->Draw("same");
	      }
	      c3a->Update();
	    }
	    delete c3a;

	    gStyle->SetOptTitle(1);
	    gStyle->SetTitleOffset(-0.06,"XYZ");
	    gStyle->SetTitleFillColor(10);
	    gStyle->SetTitleFontSize(0.06);
	    gStyle->SetPadTopMargin(0.08);
	    gStyle->SetPadRightMargin(0.15);
	    gStyle->SetPadLeftMargin(0.12);
	    gStyle->SetPadBottomMargin(0.09);

	    ps.NewPage();
	    TCanvas* c7 = new TCanvas("c7", "c7", 700, 900);
	    c7->Divide(1,2);

	    TH2F* tmp2d = (TH2F*)time_layerstrip->Clone("hist2d");
	    c7->cd(1); gPad->SetLogz(1); tmp2d->Draw("colz");

	    TH1F* tmp1d = (TH1F*)tmp2d->ProjectionX("zeroTime", 0, 1);
	    TH1F* tmp1dall= (TH1F*)tmp2d->ProjectionX("allTime", 0, tmp2d->GetNbinsY());
	    tmp1d->Divide(tmp1dall);
	    //	    tmp1d->Scale(1./max(1., time_layerstrip->GetBinContent(0,0)));
	    tmp1d->GetYaxis()->SetTitle("Fraction of missing time");
	    c7->cd(2);gPad->SetLogy(1);tmp1d->GetYaxis()->SetRangeUser(0.0001,1.);tmp1d->Draw();

	    c7->Update();

	     gStyle->SetOptTitle(1);
	    gStyle->SetTitleOffset(-0.06,"XYZ");
	    gStyle->SetTitleFillColor(10);
	    gStyle->SetTitleFontSize(0.06);
	    gStyle->SetPadTopMargin(0.08);
	    gStyle->SetPadRightMargin(0.15);
	    gStyle->SetPadLeftMargin(0.12);
	    gStyle->SetPadBottomMargin(0.09);

	    ps.NewPage();
	    TCanvas* c7a = new TCanvas("c7a", "c7a", 700, 900);
	    c7a->Divide(1,2);
	    time_layerstrip->GetYaxis()->SetRangeUser(2500., 4500.0);
	    TH2F* tmpp2dx = (TH2F*)time_layerstrip->Clone("hist2dx");
	    TH2F* tmpp2dy = (TH2F*)time_layerstrip->Clone("hist2dy");
	    c7a->cd(1); gPad->SetLogz(1); tmpp2dx->GetXaxis()->SetRangeUser(-0.5,nstrip*nlayer-0.5);tmpp2dx->Draw("colz");
	    c7a->cd(2); gPad->SetLogz(1); tmpp2dy->GetXaxis()->SetRangeUser(nstrip*nlayer-0.5,2*nlayer*nstrip-0.5);tmpp2dy->Draw("colz");
	    c7a->Update();




// 	    // for (int ij=0; ij<nlayer; ij++) {
// // 	      for (int jk=0; jk<nstrip; jk++) {
// // 		if (time_xraw[ij][jk]->GetEntries() >50.0) {
// // 		  for (int kl=time_xraw[ij][jk]->GetNbinsX(); kl>10; kl--) {
// // 		    if (time_xraw[ij][jk]->GetBinContent(kl)>1) {
// // 		      time_xraw2d->Fill(ij, jk,
// // 					(time_xraw[ij][jk]->GetBinContent(kl)
// // 					 +time_xraw[ij][jk]->GetBinContent(kl-1))
// // 					/max(1., time_xraw[ij][jk]->GetBinContent(kl-9)
// // 					     +time_xraw[ij][jk]->GetBinContent(kl-10)));
// // 		      break;
// // 		    }
// // 		  }
// // 		} else {
// // 		  file_out<<"xstrip_tentry "<<ij<<" "<<jk<<" "<<time_xraw[ij][jk]->GetEntries()<<endl;
// // 		}

// // 		if (time_yraw[ij][jk]->GetEntries() >50.0) {
// // 		  for (int kl=time_yraw[ij][jk]->GetNbinsX(); kl>10; kl--) {
// // 		    if (time_yraw[ij][jk]->GetBinContent(kl)>1) {
// // 		      time_yraw2d->Fill(ij, jk,
// // 					(time_yraw[ij][jk]->GetBinContent(kl)
// // 					 +time_yraw[ij][jk]->GetBinContent(kl-1))
// // 					/max(1., time_yraw[ij][jk]->GetBinContent(kl-9)
// // 					     +time_yraw[ij][jk]->GetBinContent(kl-10)));
// // 		      break;
// // 		    }
// // 		  }
// // 		} else {
// // 		  file_out<<"ystrip_tentry "<<ij<<" "<<jk<<" "<<time_yraw[ij][jk]->GetEntries()<<endl;
// // 		}
// // 	      }
// // 	    }

// // 	    c7->cd(2);// c7->Draw();

// // 	    TPad* d1 = new TPad("d1", "d1", 0., 0., 1., 1.);
// // 	    d1->Draw();
// // 	    d1->Divide(2,1,1.e-10, 1.e-10,0.0);
// // 	    TPad* dd1 = (TPad*) d1->GetPad(1); dd1->SetPad(0.01, 0.01, 0.49, 0.99); dd1->cd();
// // 	    time_xraw2d->GetXaxis()->SetLabelOffset(0.001);
// // 	    time_xraw2d->SetMaximum(1.0); time_xraw2d->Draw("colz");
// // 	    TPad* dd2 = (TPad*) d1->GetPad(2); dd2->SetPad(0.51, 0.01, 0.99, 0.99); dd2->cd();
// // 	    time_yraw2d->GetXaxis()->SetLabelOffset(0.001);
// // 	    time_yraw2d->SetMaximum(1.0); time_yraw2d->Draw("colz");
// // 	    c7->Update();

	    for (int ij=0; ij<nlayer; ij++) {
	      for (int jk=0; jk<nstrip; jk++) {
		for (int kl=0; kl<nstrip; kl++) {
		  for (int lm=0; lm<ntimecor; lm++) {
		    //		    cout <<" ijjkkllm "<< ij<< " "<<jk<<" "<<kl<<" "<<lm<<endl;
		    double ent = indtimexy_correl[ij][jk][kl][lm]->GetEntries();
		    if (ent >10) {
		      double mean = indtimexy_correl[ij][jk][kl][lm]->GetMean();
		      double rms = indtimexy_correl[ij][jk][kl][lm]->GetRMS();

		      indtimexy_cormean[ij][lm]->Fill(jk, kl, mean+20);
		      indtimexy_corrms[ij][lm]->Fill(jk, kl, rms);

		      if (indtimexy_correl[ij][jk][kl][lm]->GetEntries()>3) {
			//			cout <<"0ijjkkllm "<< ij<< " "<<jk<<" "<<kl<<" "<<lm<<" "<<indtimexy_correl[ij][jk][kl][lm]->GetEntries()<<endl;
			TFitResultPtr ptr = indtimexy_correl[ij][jk][kl][lm]->Fit("gaus","SQ0"); //,"",istr, iend);
			Int_t fitStatus = ptr;
			if (fitStatus==0) {
			  mean = ptr->Parameter(1);
			  rms = ptr->Parameter(2);
			  //			  cout <<"status code "<<ptr<<" "<<indtimexy_correl[ij][jk][kl][lm]->GetEntries()<<" "<< ij<< " "<<jk<<" "<<kl<<" "<<lm<<" "<<mean<<" "<<rms<< endl;
			  //			cout <<"1ijjkkllm "<< ij<< " "<<jk<<" "<<kl<<" "<<lm<<" "<<indtimexy_correl[ij][jk][kl][lm]->GetEntries()<<" "<<mean<<" "<<rms<<endl;
			  indtimexy_fitmean[ij][lm]->Fill(jk, kl, mean+20);
			  indtimexy_fitrms[ij][lm]->Fill(jk, kl, rms);
			}
		      }
		    }
		  }
		}
	      }
	    }

	    //	    ps.NewPage();
	    gStyle->SetTitleFontSize(0.07);
	    for (int ixx=0; ixx<4; ixx++) {
	      for (int lm=0; lm<ntimecor; lm++) {
		TH2F* histxx[nlayer]={0};
		ps.NewPage();
		for (int ij=0; ij<nlayer; ij++) {
		  c4a->cd(ij+1);
		  switch(ixx) {
		  case 0 : histxx[ij] = (TH2F*)indtimexy_cormean[ij][lm]->Clone(); break;
		  case 1 : histxx[ij] = (TH2F*)indtimexy_corrms[ij][lm]->Clone(); break;
		  case 2 : histxx[ij] = (TH2F*)indtimexy_fitmean[ij][lm]->Clone(); break;
		  case 3 : histxx[ij] = (TH2F*)indtimexy_fitrms[ij][lm]->Clone(); break;
		  default : histxx[ij] = (TH2F*)indtimexy_cormean[ij][lm]->Clone(); break;
		  }

		  histxx[ij]->GetXaxis()->SetLabelSize(.07);
		  histxx[ij]->GetXaxis()->SetTitle("X-strip");
		  histxx[ij]->GetXaxis()->CenterTitle();
		  histxx[ij]->GetXaxis()->SetTitleSize(0.07);
		  histxx[ij]->GetXaxis()->SetTitleOffset(.9);

		  histxx[ij]->GetYaxis()->SetTitle("Y-strip");
		  histxx[ij]->GetYaxis()->CenterTitle();
		  histxx[ij]->GetYaxis()->SetTitleSize(0.07);
		  histxx[ij]->GetYaxis()->SetTitleOffset(.9);
		  histxx[ij]->GetYaxis()->SetLabelSize(.07);


		  double amn = histxx[ij]->GetMinimum();
		  double amx = histxx[ij]->GetMaximum();
		if (amx >amn+0.5) {
		  switch(ixx) {
		  case 0 : ;
		  case 2 :
		    histxx[ij]->SetMaximum(min(30.0, amx));
		    histxx[ij]->SetMinimum(max(10.0, amn)); break;
		  case 1 :
		    histxx[ij]->SetMaximum(min(3.0, amx));
		    histxx[ij]->SetMinimum(max(0.2, amn)); break;
		  case 3 :
		    histxx[ij]->SetMaximum(min(2.0, amx));
		    histxx[ij]->SetMinimum(max(0.2, amn)); break;
		  default : break;
		  }
		}
		histxx[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
		histxx[ij]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
		histxx[ij]->Draw("colz");
		}

		c4a->Update();
		for (int ij=0; ij<nlayer; ij++) {
		  if (histxx[ij]) { delete histxx[ij]; histxx[ij]=0;}
		}
	      }
	    }

	    ps.NewPage();
	    gStyle->SetOptFit(101);
	    gStyle->SetStatW(.30);
	    gStyle->SetStatH(.20);
	    for (int lm=0; lm<ntimecor; lm++) {
	      ps.NewPage();

	      for (int ij=0; ij<nlayer; ij++) {
		c4a->cd(ij+1);
		int istr=1;
		int iend =63;
		int nbin = indtimexy_prof[ij][lm]->GetNbinsX();
		for (int ix=nbin/2; ix>=0; ix--) {
		  if (indtimexy_prof[ij][lm]->GetBinError(ix+1) <1.e-4) {
		    istr = ix+1; break;
		  }
		}

		for (int ix=nbin/2; ix<nbin; ix++) {
		  if (indtimexy_prof[ij][lm]->GetBinError(ix+1) <1.e-4) {
		    iend = ix-1; break;
		  }
		}

		if (istr<nbin/2 && iend>nbin/2 ) {
		  indtimexy_prof[ij][lm]->GetXaxis()->SetLabelSize(.07);
		  indtimexy_prof[ij][lm]->GetYaxis()->SetLabelSize(.07);
		  /*TFitResultPtr ptr = */ indtimexy_prof[ij][lm]->Fit("pol2","SQ","",istr, iend);
		}
	      }
	      c4a->Update();
	    }
	  } // if (isTiming)
	  //#endif //MONTECARLO
#endif //ISTIMING
	} // if (filloccu)

	if (occulyr==nlayer) {
	  totxtime /=max(1,ntotxtime);
	  totytime /=max(1,ntotytime);
	  ytimeshift += totytime - totxtime;
    if (iiter == nmxiter-1) {
	    //file_out_h<<"double ytimeshift ="<< iiter <<" "<< ytimeshift<<" "<<totytime<<" "<<totxtime<<endl;
      file_out_h<<"double ytimeshift =" << ytimeshift<<endl;
    }
	}

	double ncx = 0;
	double ncy = 0;

	if (isTiming) {
	  if (firstiter==0 && ntcor==0) {

	    firstiter=1;

	    for (int ij=1; ij<nlayer; ij++) {
	      ps.NewPage();
	      gStyle->SetOptTitle(1);
	      gStyle->SetOptStat(111110);
	      gStyle->SetOptFit(101);
	      gStyle->SetOptLogy(1);
	      gStyle->SetTitleFontSize(0.06);
	      gStyle->SetPadBottomMargin(0.10);
	      gStyle->SetPadTopMargin(0.07);
	      gStyle->SetPadLeftMargin(0.10);
	      gStyle->SetPadRightMargin(0.02);
	      gStyle->SetStatW(.20);
	      gStyle->SetStatH(.14);

	      TCanvas *c3=new TCanvas ("c3","Time Residual",500,700);
	      c3->Divide(2,2);
	      TF1* fittx[2]={0};
	      TF1* fitty[2]={0};

	      for (int ix=0; ix<2; ix++) {

		c3->cd(2*ix+1);
		int iyy=(ix==0) ? 0 : iiter+1;
		sprintf(name, "fittx_%i", ix);
		fittx[ix] = new TF1(name, gausX, timex_shift[ij][iyy]->GetXaxis()->GetXmin(), timex_shift[ij][iyy]->GetXaxis()->GetXmax(), 3);
		double  parx[3]={timex_shift[ij][iyy]->GetMaximum(), timex_shift[ij][iyy]->GetMean(), timex_shift[ij][iyy]->GetRMS()};
		fittx[ix]->SetParameters(parx);

		timex_shift[ij][iyy]->GetXaxis()->SetTitle("#Deltat (ns)");
		timex_shift[ij][iyy]->GetXaxis()->SetTitleOffset(.7);
		timex_shift[ij][iyy]->GetXaxis()->SetTitleSize(.06);
		timex_shift[ij][iyy]->GetXaxis()->CenterTitle();
		timex_shift[ij][iyy]->GetXaxis()->SetLabelSize(.05);
		timex_shift[ij][iyy]->GetXaxis()->SetLabelOffset(-0.001);

		timex_shift[ij][iyy]->Fit(fittx[ix],"Q");

		if (isTimeCorOrReso && ix>0 && isalign >0 && iiter<nmxiter-isLast &&
		    timex_shift[ij][iyy]->Integral() >nmnentry &&
		    abs(timex_shift[ij][iyy]->GetMean() - fittx[ix]->GetParameter(1)) <
		    timex_shift[ij][iyy]->GetRMS()) {
		  timeoffsetx[ij] +=timex_shift[ij][iyy]->GetMean(); //fittx[ix]->GetParameter(1);

		  //		  file_out<<"timeoffsetx["<<ij<<"]= "<<ntcor<<" "<<timeoffsetx[ij]<<" "<<timex_shift[ij][iyy]->GetMean()<<endl;

		}

		c3->cd(2*ix+2);
		sprintf(name, "fitty_%i", ix);
		fitty[ix] = new TF1(name, gausX, timey_shift[ij][iyy]->GetXaxis()->GetXmin(), timey_shift[ij][iyy]->GetXaxis()->GetXmax(), 3);
		double pary[3]={timey_shift[ij][iyy]->GetMaximum(), timey_shift[ij][iyy]->GetMean(), timey_shift[ij][iyy]->GetRMS()};
		fitty[ix]->SetParameters(pary);
		timey_shift[ij][iyy]->GetXaxis()->SetTitle("#Deltat (ns)");
		timey_shift[ij][iyy]->GetXaxis()->CenterTitle();
		timey_shift[ij][iyy]->GetXaxis()->SetTitleOffset(.7);
		timey_shift[ij][iyy]->GetXaxis()->SetTitleSize(.06);
		timey_shift[ij][iyy]->GetXaxis()->SetLabelSize(.05);
		timey_shift[ij][iyy]->GetXaxis()->SetLabelOffset(-0.001);
		timey_shift[ij][iyy]->Fit(fitty[ix],"Q");

		if (ix>0) {
		  if (isTimeCorOrReso && isalign>0 && iiter<nmxiter-isLast &&
		      timey_shift[ij][iyy]->Integral() >nmnentry &&
		      abs(timey_shift[ij][iyy]->GetMean() - fitty[ix]->GetParameter(1)) <
		      timey_shift[ij][iyy]->GetRMS()) {
		    timeoffsety[ij] +=timey_shift[ij][iyy]->GetMean(); //fitty[ix]->GetParameter(1);
		    //		    file_out<<"timeoffsety["<<ij<<"]= "<<ntcor<<" "<<timeoffsety[ij]<<" "<<timex_shift[ij][iyy]->GetMean()<<endl;
		  }

		  file_out <<"global time "<<iiter<<" "<<ij<<" "<< parx[0]<<" "<<parx[1]<<" "<<timeoffsetx[ij]<<" "<<parx[2]<< " y "
			   << pary[0]<<" "<<pary[1]<<" "<<timeoffsety[ij]<<" "<<pary[2]<<endl;

		  file_out<<"fit "<<iiter<<" "<<ij<<" "<<fittx[ix]->GetParameter(0)<<" "<<fittx[ix]->GetParameter(1)<<" "<<fittx[ix]->GetParameter(2)<<" "<<fitty[ix]->GetParameter(0)<<" "<<fitty[ix]->GetParameter(1)<<" "<<fitty[ix]->GetParameter(2)<<endl;
		}

	      } //  for (int ix=0; ix<2; ix++)
	      c3->Update();
	      if (c3) { delete c3; c3=0;}
	      for (int ix=0; ix<2; ix++) {
		if (fittx[ix]) { delete fittx[ix]; fittx[ix]=0;}
		if (fitty[ix]) { delete fitty[ix]; fitty[ix]=0;}
	      }
	    } //for (int ij=1; ij<nlayer; ij++)


	    for (int ij=1; ij<nlayer; ij++) {
	      ps.NewPage();

	      gStyle->SetOptStat(0);
	      gStyle->SetOptFit(0);
	      gStyle->SetOptLogy(0);
	      gStyle->SetTitleFontSize(0.06);
	      gStyle->SetPadBottomMargin(0.10);
	      gStyle->SetPadTopMargin(0.07);
	      gStyle->SetPadLeftMargin(0.08);
	      gStyle->SetPadRightMargin(0.12);

	      TCanvas *c3=new TCanvas ("c3","Time Residual",500,700);
	      c3->Divide(2,2);

	      for (int ix=0; ix<2; ix++) {

		c3->cd(2*ix+1);
		int iyy=(ix==0) ? 0 : iiter+1;
		timex_2dshift[ij][iyy]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
		timex_2dshift[ij][iyy]->GetXaxis()->SetTitle("Strip No");
		timex_2dshift[ij][iyy]->GetXaxis()->SetTitleOffset(.7);
		timex_2dshift[ij][iyy]->GetXaxis()->SetTitleSize(.06);
		timex_2dshift[ij][iyy]->GetXaxis()->CenterTitle();
		timex_2dshift[ij][iyy]->GetXaxis()->SetLabelSize(.05);
		timex_2dshift[ij][iyy]->GetXaxis()->SetLabelOffset(-0.001);

		timex_2dshift[ij][iyy]->GetYaxis()->SetTitle("#Deltat (ns)");
		timex_2dshift[ij][iyy]->GetYaxis()->SetTitleOffset(.7);
		timex_2dshift[ij][iyy]->GetYaxis()->SetTitleSize(.06);
		timex_2dshift[ij][iyy]->GetYaxis()->CenterTitle();
		timex_2dshift[ij][iyy]->GetYaxis()->SetLabelSize(.05);
		timex_2dshift[ij][iyy]->GetZaxis()->SetLabelSize(.04);
		timex_2dshift[ij][iyy]->Draw("colz");

		c3->cd(2*ix+2);
		timey_2dshift[ij][iyy]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
		timey_2dshift[ij][iyy]->GetXaxis()->SetTitle("Strip No");
		timey_2dshift[ij][iyy]->GetXaxis()->CenterTitle();
		timey_2dshift[ij][iyy]->GetXaxis()->SetTitleOffset(.7);
		timey_2dshift[ij][iyy]->GetXaxis()->SetTitleSize(.06);
		timey_2dshift[ij][iyy]->GetXaxis()->SetLabelSize(.05);
		timey_2dshift[ij][iyy]->GetXaxis()->SetLabelOffset(-0.001);

		timey_2dshift[ij][iyy]->GetYaxis()->SetTitle("#Deltat (ns)");
		timey_2dshift[ij][iyy]->GetYaxis()->CenterTitle();
		timey_2dshift[ij][iyy]->GetYaxis()->SetTitleOffset(.7);
		timey_2dshift[ij][iyy]->GetYaxis()->SetTitleSize(.06);
		timey_2dshift[ij][iyy]->GetYaxis()->SetLabelSize(.05);
		timey_2dshift[ij][iyy]->GetZaxis()->SetLabelSize(.04);
		timey_2dshift[ij][iyy]->Draw("colz");

	      } //  for (int ix=0; ix<2; ix++)
	      c3->Update();
	      if (c3) { delete c3; c3=0;}
	    } //for (int ij=1; ij<nlayer; ij++)



	  } // if (firstiter==0 && ntcor==0)

	  //	  double ncx = 0;
	  //	  double ncy = 0;

	  for (int ij=0; ij<dir_cx[occulyr][iiter]->GetNbinsX(); ij++) {
	    if (dir_cx[occulyr][iiter]->GetBinCenter(ij)<0.0) {
	      ncx +=dir_cx[occulyr][iiter]->GetBinContent(ij);
	      ncy +=dir_cy[occulyr][iiter]->GetBinContent(ij);
	    } else { break;}
	  }
	  file_out <<" occulyr "<< occulyr<<" "<< iiter <<" xc "<<dir_cx[occulyr][iiter]->Integral()<<" "<<dir_cx[occulyr][iiter]->GetBinContent(0)<<" "<<dir_cx[occulyr][iiter]->GetEntries()<<" "<<dir_cx[occulyr][iiter]->GetMean()<<" "<<dir_cx[occulyr][iiter]->GetRMS()<<" yc "<<dir_cy[occulyr][iiter]->Integral()<<" "<<dir_cy[occulyr][iiter]->GetBinContent(0)<<" "<<dir_cy[occulyr][iiter]->GetEntries()<<" "<<dir_cy[occulyr][iiter]->GetMean()<<" "<<dir_cy[occulyr][iiter]->GetRMS()<<" ncx "<<ncx<<" ncy "<<ncy<<endl;

        } //if (isTiming)


	ps.NewPage();
	gStyle->SetOptStat(111110);
	gStyle->SetOptFit(101);
	gStyle->SetOptLogy(1);

	gStyle->SetPadBottomMargin(0.11);
	gStyle->SetPadTopMargin(0.09);
	gStyle->SetPadLeftMargin(0.07);
	gStyle->SetPadRightMargin(0.02);
	gStyle->SetOptTitle(1);
	gStyle->SetTitleFontSize(0.07);
	gStyle->SetStatW(.30);
	gStyle->SetStatH(.14);
	gStyle->SetStatY(.99);
	gStyle->SetStatX(.99);

	TCanvas *c2a=new TCanvas ("c2a","Slope",500,700);
	if (isTiming) {c2a->Divide(2,3);} else {c2a->Divide(2,2);}
	const int nxplot=6;
	//	latex.SetTextSize(0.08);
	TH1F* histax[nxplot]={0};
	for (int ix=0; ix<nxplot; ix++) {
	  if (!isTiming && ix>=4) continue;
	  switch(ix) {
	  case 0 : histax[ix] = (TH1F*)pos_xslope[occulyr][iiter]->Clone(); break;
	  case 1 : histax[ix] = (TH1F*)pos_yslope[occulyr][iiter]->Clone(); break;
	  case 2 : histax[ix] = (TH1F*)pos_theta[occulyr][iiter]->Clone(); break;
	  case 3 : histax[ix] = (TH1F*)pos_phi[occulyr][iiter]->Clone(); break;
	  case 4 : histax[ix] = (TH1F*)dir_cx[occulyr][iiter]->Clone(); break;
	  case 5 : histax[ix] = (TH1F*)dir_cy[occulyr][iiter]->Clone(); break;
	  default : histax[ix] = (TH1F*)pos_xslope[occulyr][iiter]->Clone(); break;
	  }
	  c2a->cd(ix+1);
	  histax[ix]->GetXaxis()->CenterTitle();
	  histax[ix]->GetXaxis()->SetTitleOffset(0.7);
	  histax[ix]->GetXaxis()->SetTitleSize(0.06);
	  histax[ix]->GetXaxis()->SetLabelOffset(-0.01);
	  histax[ix]->GetXaxis()->SetLabelSize(0.06);
	  histax[ix]->Draw();
	  //	  if (ix==4) {latex.DrawLatex(0.16, 0.66,Form("#font[12]#color[2]{%g}", int(100000*ncx/max(1.,dir_cx[occulyr][iiter]->GetEntries()))/1000.));}
	  //	  if (ix==5) {latex.DrawLatex(0.16, 0.66,Form("#scale[1.2]#color[2]{%g}", int(100000*ncy/max(1.,dir_cy[occulyr][iiter]->GetEntries()))/1000.));}
	  if (ix==4) {latex.DrawLatex(0.12, 0.66,Form("#scale[0.6]{%g\%}", int(1000000*ncx/max(1.,dir_cx[occulyr][iiter]->GetEntries()))/10000.));}
	  if (ix==5) {latex.DrawLatex(0.12, 0.66,Form("#scale[0.6]{%g\%}", int(1000000*ncy/max(1.,dir_cy[occulyr][iiter]->GetEntries()))/10000.));}

	}
	c2a->Update();
	if (c2a) { delete c2a; c2a=0;}
	for (int ix=0; ix<nxplot; ix++) {
	  if (histax[ix]) { delete histax[ix]; histax[ix]=0;}
	}

	if (isTiming) {
	  ps.NewPage();
	  TCanvas* c4c = new TCanvas("c4c", "c4c", 700, 900);
	  c4c->Divide(3,4);
	  TH1F* histbx[2*nprob]={0};
	  for (int ix=0; ix<2*nprob; ix++) {
	    if (ix<nprob) {
	      histbx[ix] = (TH1F*)dir_cxlay[occulyr][iiter][ix]->Clone();
	    } else {
	      histbx[ix] = (TH1F*)dir_cylay[occulyr][iiter][ix-nprob]->Clone();
	    }
	    c4c->cd(ix+1);
	    int ncx=0;
	    for (int ij=0; ij<histbx[ix]->GetNbinsX(); ij++) {
	      if (histbx[ix]->GetBinCenter(ij)<0.0) {
		ncx +=histbx[ix]->GetBinContent(ij);
	      } else { break;}
	    }
	    histbx[ix]->GetXaxis()->CenterTitle();
	    histbx[ix]->GetXaxis()->SetTitleOffset(0.8);
	    histbx[ix]->GetXaxis()->SetTitleSize(0.07);
	    histbx[ix]->GetXaxis()->SetLabelOffset(-0.01);
	    histbx[ix]->GetXaxis()->SetLabelSize(0.07);
	    histbx[ix]->Draw();
	    latex.DrawLatex(0.12, 0.66,Form("#scale[0.8]{%g\%}", int(1.e6*ncx/max(1.,histbx[ix]->GetEntries()))/1.e4));
	  }
	  c4c->Update();
	  if (c4c) {delete c4c; c4c=0;}
	  for (int ix=0; ix<2*nprob; ix++) {
	    if (histbx[ix]) { delete histbx[ix]; histbx[ix]=0;}
	  }
	} // if (isTiming)

	int istr = (occulyr >= nlayer) ? 0 : occulyr;
	int iend = (occulyr >= nlayer) ? nlayer : occulyr+1;

        for (int iocc = istr; iocc<iend; iocc++) {
	  double absWidth=3.0;//6.0;//3.0;//2.5;//min(2.5,max(1.9, 0.15*(nmxiter -iiter)));
	  double widthScale=1.5;//2.0;//1.0;//min(0.70,max(0.40, 0.05*(nmxiter -iiter))); //Initially one can allow very large asymmetric distribution and and gradually make it smaller

	  ps.NewPage();
	  gStyle->SetOptStat(1110);
	  gStyle->SetOptFit(101);
	  gStyle->SetOptLogy(1);

	  gStyle->SetPadBottomMargin(0.11);
	  gStyle->SetPadTopMargin(0.09);
	  gStyle->SetPadLeftMargin(0.09);
	  gStyle->SetPadRightMargin(0.02);
	  gStyle->SetOptTitle(1);
	  gStyle->SetTitleFontSize(0.07);
	  gStyle->SetStatW(.20);
	  gStyle->SetStatH(.12);
	  gStyle->SetStatY(.99);
	  gStyle->SetStatX(.99);

	  TCanvas *c2=new TCanvas ("c2","Residue",500,700);
	  if (isTiming) { c2->Divide(2,3);} else {c2->Divide(1,3);}

	  c2->cd(1);

	  double alw= (ntcor==0) ? -2.05 : -2.8; //xlayer_reso[iocc][iiterrs]->GetXaxis()->GetXmin();
	  double ahg= (ntcor==0) ?  2.05 :  2.8; //xlayer_reso[iocc][iiterrs]->GetXaxis()->GetXmax();
	  TF1* fitfx = new TF1("fitfx", gausX, alw, ahg, 3);

	  double  parx[3]={xlayer_reso[iocc][iiterrs]->GetMaximum(), xlayer_reso[iocc][iiterrs]->GetMean(), max(0.15,0.7*xlayer_reso[iocc][iiterrs]->GetRMS())};
	  fitfx->SetParameters(parx);
	  //	    fitfx->SetParLimits(0, 0.17*parx[0], 1.5*parx[0]);
	  fitfx->SetParLimits(1, parx[1]-1., parx[1]+1.);
	  fitfx->SetParLimits(2, 0.12, 1.5*parx[2]);

	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitle("X-residues (pitch)");
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->CenterTitle();
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);

	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
	  xlayer_reso[iocc][iiterrs]->GetYaxis()->SetLabelSize(0.06);
          xlayer_reso[iocc][iiterrs]->Fit(fitfx, "BRQ");

	  c2->cd(2);

	  TF1* fitfy = new TF1("fitfy", gausX, alw, ahg, 3);
	  double  pary[3]={ylayer_reso[iocc][iiterrs]->GetMaximum(), ylayer_reso[iocc][iiterrs]->GetMean(), max(0.15,0.7*ylayer_reso[iocc][iiterrs]->GetRMS())};
	  fitfy->SetParameters(pary);
	  //	    fitfy->SetParLimits(0, 0.15*pary[0], 1.5*pary[0]);
	  fitfy->SetParLimits(1, pary[1]-1., pary[1]+1.);
	  //	  fitfy->SetParLimits(2, 0.12, 1.5*pary[2]);

	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitle("Y-residues (pitch)");
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->CenterTitle();
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);

	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
	  ylayer_reso[iocc][iiterrs]->GetYaxis()->SetLabelSize(0.06);

          ylayer_reso[iocc][iiterrs]->Fit(fitfy, "BRQ");

	  //time resolution only for plot, but not for any calculation
	  alw= (ntcor==0) ? -7.0 :-10.0;
	  ahg= (ntcor==0) ?  7.0 : 10.0;
	  TF1* fittfx = new TF1("fittfx", gausX, alw, ahg, 3);
	  TF1* fittfy = new TF1("fittfy", gausX, alw, ahg, 3);
	  if (isTiming) {
	    c2->cd(3);
	    //	    TF1* fittfx = new TF1("fittfx", gausX, alw, ahg, 3);

	    double  partx[3]={time_xreso[iocc][iiterrs]->GetMaximum(), time_xreso[iocc][iiterrs]->GetMean(), max(0.15,0.7*time_xreso[iocc][iiterrs]->GetRMS())};
	    fittfx->SetParameters(partx);
	    //	    fittfx->SetParLimits(0, 0.17*partx[0], 1.5*partx[0]);
	    fittfx->SetParLimits(1, partx[1]-1., partx[1]+1.);
	    fittfx->SetParLimits(2, 0.6, 1.5*partx[2]);

	    time_xreso[iocc][iiterrs]->GetXaxis()->SetTitle("X - #Deltat (ns)");
	    time_xreso[iocc][iiterrs]->GetXaxis()->CenterTitle();

	    time_xreso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	    time_xreso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
	    time_xreso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	    time_xreso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
	    time_xreso[iocc][iiterrs]->GetYaxis()->SetLabelSize(0.06);
	    time_xreso[iocc][iiterrs]->Fit(fittfx, "BRQ");
	    time_xrms[iocc] = fittfx->GetParameter(2);
	    c2->cd(4);
	    //	    TF1* fittfy = new TF1("fittfy", gausX, alw, ahg, 3);

	    double  party[3]={time_yreso[iocc][iiterrs]->GetMaximum(), time_yreso[iocc][iiterrs]->GetMean(), max(0.15,0.7*time_yreso[iocc][iiterrs]->GetRMS())};
	    fittfy->SetParameters(party);
	    //	    fittfy->SetParLimits(0, 0.17*party[0], 1.5*party[0]);
	    fittfy->SetParLimits(1, party[1]-1., party[1]+1.);
	    fittfy->SetParLimits(2, 0.6, 1.5*party[2]);

	    time_yreso[iocc][iiterrs]->GetXaxis()->SetTitle("Y - #Deltat (ns)");
	    time_yreso[iocc][iiterrs]->GetXaxis()->CenterTitle();
	    time_yreso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	    time_yreso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
	    time_yreso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	    time_yreso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
	    time_yreso[iocc][iiterrs]->GetYaxis()->SetLabelSize(0.06);

	    time_yreso[iocc][iiterrs]->Fit(fittfy, "RQ");
	    time_yrms[iocc] = fittfy->GetParameter(2);
	    file_out <<"lay "<< iocc<<" it "<< iiterrs <<" X-sh "<< fitfx->GetParameter(1)<<" "<<fitfx->GetParameter(2)<<" Y-sh "<< fitfy->GetParameter(1)<<" "<<fitfy->GetParameter(2)<<" X-sh "<< fittfx->GetParameter(1)<<" "<<fittfx->GetParameter(2)<<" "<<time_xreso[iocc][iiterrs]->GetRMS()<<" Y-sh "<< fittfy->GetParameter(1)<<" "<<fittfy->GetParameter(2)<<" "<<time_yreso[iocc][iiterrs]->GetRMS()<<endl;

	    c2->cd(5);
	    tmptime_xreso[iocc][iiterrs]->GetXaxis()->SetTitle("tmpX - #Deltat (ns)");
	    tmptime_xreso[iocc][iiterrs]->GetXaxis()->CenterTitle();

	    tmptime_xreso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	    tmptime_xreso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
	    tmptime_xreso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	    tmptime_xreso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
	    tmptime_xreso[iocc][iiterrs]->GetYaxis()->SetLabelSize(0.06);
	    tmptime_xreso[iocc][iiterrs]->Fit("gaus", "Q");

	    c2->cd(6);
	    tmptime_yreso[iocc][iiterrs]->GetXaxis()->SetTitle("tmpY - #Deltat (ns)");
	    tmptime_yreso[iocc][iiterrs]->GetXaxis()->CenterTitle();
	    tmptime_yreso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	    tmptime_yreso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
	    tmptime_yreso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	    tmptime_yreso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
	    tmptime_yreso[iocc][iiterrs]->GetYaxis()->SetLabelSize(0.06);
	    tmptime_yreso[iocc][iiterrs]->Fit("gaus", "Q");

	  } //if (isTiming)

	  c2->Update();
	  if (c2) { delete c2; c2=0;}
	  // //	  if (isCombinedOffset) {
// 	  if (iiter<nmxiter-1 && ntcor==0) {
// 	    if (fabs(xlayer_reso[iocc][iiterrs]->GetMean()-
// 		     fitfx->GetParameter(1))
// 		< xlayer_reso[iocc][iiterrs]->GetRMS()) {
// 	      xoff[iocc] += fitfx->GetParameter(1);
// 	    }

// 	    if (fabs(ylayer_reso[iocc][iiterrs]->GetMean()-
// 		     fitfy->GetParameter(1))
// 		< ylayer_reso[iocc][iiterrs]->GetRMS()) {
// 	      yoff[iocc] += fitfy->GetParameter(1);
// 	    }
// 	  }
// 	    //	  }

	  xrms[iocc] = fitfx->GetParameter(2);
	  yrms[iocc] = fitfy->GetParameter(2);

	  delete fitfx; fitfx=0;
	  delete fitfy; fitfy=0;
	  delete fittfx; fittfx=0;
	  delete fittfy; fittfy=0;

	  if (ntcor==1 && iiter<nmxiter-isLast) {//Don't correct for last iteration //150112
	    bias_inpos_xreso2[iocc] = (xlayer_exterr[iocc][iiterrs]->GetMean())*(xlayer_exterr[iocc][iiterrs]->GetMean());
	    bias_inpos_yreso2[iocc] = (ylayer_exterr[iocc][iiterrs]->GetMean())*(ylayer_exterr[iocc][iiterrs]->GetMean());

	    for (int ix=0; ix<nmxhits; ix++) {
	      //	      pos_xrms[iocc][ix] = xlayer_reso_mul[iocc][iiterrs][ix]->GetRMS(); //Get value while this is used in fit
	      //	      pos_yrms[iocc][ix] = ylayer_reso_mul[iocc][iiterrs][ix]->GetRMS();

	      TFitResultPtr ptrx = xlayer_reso_mul[iocc][iiterrs][ix]->Fit("gaus","SQ0");
	      TFitResultPtr ptry = ylayer_reso_mul[iocc][iiterrs][ix]->Fit("gaus","SQ0");
	      Int_t fitStatus = ptrx;
	      pos_xrms[iocc][ix] = (fitStatus==0 && xlayer_reso_mul[iocc][iiterrs][ix]->GetEntries()>3) ? ptrx->Parameter(2) : 0.0;
	      fitStatus = ptry;
	      pos_yrms[iocc][ix] = (fitStatus==0 && ylayer_reso_mul[iocc][iiterrs][ix]->GetEntries()>3) ? ptry->Parameter(2) : 0.0;
	      //	      cout <<" pos_xrms[iocc][ix] "<<pos_xrms[iocc][ix]<<" "<<xlayer_reso_mul[iocc][iiterrs][ix]->GetRMS()<<" "
	      //		   <<pos_yrms[iocc][ix]<<" "<<ylayer_reso_mul[iocc][iiterrs][ix]->GetRMS()<<endl;
	      // 18th Feb 2015, for estimation of resolution
	      xposerrsq[ix][iocc] =  pos_xrms[iocc][ix]*pos_xrms[iocc][ix] - bias_inpos_xreso2[iocc];
	      yposerrsq[ix][iocc] =  pos_yrms[iocc][ix]*pos_yrms[iocc][ix] - bias_inpos_yreso2[iocc];

	      if (xposerrsq[ix][iocc]<0.02) xposerrsq[ix][iocc]=0.02; //resolution should be more than 0.15 strip width
	      if (yposerrsq[ix][iocc]<0.02) yposerrsq[ix][iocc]=0.02;
	    }
	  }
#ifdef NEWRPCDAQ1

	  //	  if ((!isCombinedOffset)) {
	  // if (ntcor == 1 && iiter<nmxiter-isLast) {//Don't correct for last iteration //150112
	    	  if (iiter<nmxiter-isLast) {//Don't correct for last iteration //150112
	      TF1* fitx = new TF1("fitex", cal_slope, 0.5, 59.5, nPosAlignPar);
	      double  parx[nPosAlignPar]={0,0,0};
	      if (ntcor==0)  {
	      	for (int ix=0; ix<nPosAlignPar; ix++) { parx[ix] = align_xstr_xdev[iocc][ix];}
	      	fitx->SetParameters(parx);
	      	prof_xstr_xdev[iocc][iiter]->GetXaxis()->SetRangeUser(2,57);
	      	TFitResultPtr ptrx = prof_xstr_xdev[iocc][iiter]->Fit(fitx,"SQ0");
	      	Int_t fitStatus = ptrx;
	      	if (fitStatus==0) {
	      	  align_xstr_xdev[iocc][0] +=ptrx->Parameter(0);
	      	  align_xstr_xdev[iocc][1] +=ptrx->Parameter(1);
	      	  align_xstr_xdev[iocc][2] +=ptrx->Parameter(2);
	      	}

	      	for (int ix=0; ix<nPosAlignPar; ix++) { parx[ix] = align_ystr_ydev[iocc][ix];}
	      	fitx->SetParameters(parx);
	      	prof_ystr_ydev[iocc][iiter]->GetXaxis()->SetRangeUser(2,59);
	      	ptrx = prof_ystr_ydev[iocc][iiter]->Fit(fitx,"SQ0");
	      	fitStatus = ptrx;
	      	if (fitStatus==0) {
	      	  align_ystr_ydev[iocc][0] +=ptrx->Parameter(0);
	      	  align_ystr_ydev[iocc][1] +=ptrx->Parameter(1);
	      	  align_ystr_ydev[iocc][2] +=ptrx->Parameter(2);
	      	}
	       } else {
		for (int ix=0; ix<nPosAlignPar; ix++) { parx[ix] = align_xstr_ydev[iocc][ix];}
		fitx->SetParameters(parx);
		prof_xstr_ydev[iocc][iiter]->GetXaxis()->SetRangeUser(2,57);
		 TFitResultPtr ptrx = prof_xstr_ydev[iocc][iiter]->Fit(fitx,"SQ0");
		 Int_t fitStatus = ptrx;
		if (fitStatus==0) {
		  align_xstr_ydev[iocc][0] +=ptrx->Parameter(0);
		  align_xstr_ydev[iocc][1] +=ptrx->Parameter(1);
		  align_xstr_ydev[iocc][2] +=ptrx->Parameter(2);
		}
		for (int ix=0; ix<nPosAlignPar; ix++) { parx[ix] = align_ystr_xdev[iocc][ix];}
		fitx->SetParameters(parx);
		prof_ystr_xdev[iocc][iiter]->GetXaxis()->SetRangeUser(2,59);
		ptrx = prof_ystr_xdev[iocc][iiter]->Fit(fitx,"SQ0");
		fitStatus = ptrx;
		if (fitStatus==0) {
		  align_ystr_xdev[iocc][0] +=ptrx->Parameter(0);
		  align_ystr_xdev[iocc][1] +=ptrx->Parameter(1);
		  align_ystr_xdev[iocc][2] +=ptrx->Parameter(2);
		}
		//	      }
		for (int ix=0; ix<nPosAlignPar; ix++) {
		  file_out<<"alignment "<<iiter<<" "<<iocc<<" "<<ix<<" "
			  <<align_xstr_xdev[iocc][ix]<<" "
			  <<align_xstr_ydev[iocc][ix]<<" "
			  <<align_ystr_ydev[iocc][ix]<<" "
			  <<align_ystr_xdev[iocc][ix]<<endl;
		}
	      }
	  } //if (iiter<nmxiter-1)
	    //	  } //if ((!isCombinedOffset))
#endif

	  if (isTiming) {
	    if (ntcor==1 && iiter<nmxiter-isLast) { //Don't correct for last iteration
	      // uncomment all Draw/Fit/c2x
	      // ps.NewPage();
// 	      gStyle->SetOptLogy(0);
// 	      gStyle->SetOptStat(0);
// 	      gStyle->SetOptFit(101);
// 	      gStyle->SetPadLeftMargin(0.08);
// 	      gStyle->SetPadRightMargin(0.12);
// 	      TCanvas *c2x=new TCanvas ("c2","Residue",500,700);
// 	      c2x->Divide(1,2);

	      TProfile* str_tdev[2]={0};
	      TH1F* str_tdev1d[2]={0};
	      for (int ixy=0; ixy<2; ixy++) {
		double meanx[nstrip]={0};
		double weightx[nstrip]={0};
		double sum=0;
		double sumwt=0;
		if (ixy==0) {
		  file_outstr <<"double xtoffystr[nlayer][nstrip] = {"<<endl;
		} else {
		  file_outstr <<"double ytoffxstr[nlayer][nstrip] = {"<<endl;
		}
		//		c2x->cd(ixy+1);
		if (ixy==0) {
		  //		  ystr_xtdev[iocc][iiter]->Draw("colz");
		  str_tdev[ixy] = (TProfile*)ystr_xtdev[iocc][iiter]->ProfileX("xxx");
		  str_tdev1d[ixy] = (TH1F*)ystr_xtdev[iocc][iiter]->ProjectionX("xxx1d");
		} else {
		  //		  xstr_ytdev[iocc][iiter]->Draw("colz");
		  str_tdev[ixy] = (TProfile*)xstr_ytdev[iocc][iiter]->ProfileX("yyy");
		  str_tdev1d[ixy] = (TH1F*)xstr_ytdev[iocc][iiter]->ProjectionX("yyy1d");
		}
		//		str_tdev[ixy]->SetMarkerSize(0.36);
		//		str_tdev[ixy]->SetMarkerStyle(24);
		//		str_tdev[ixy]->Fit("pol4","Q","same");
		int nbin = min(nstrip,str_tdev[ixy]->GetNbinsX());
		for (int jk=0; jk<nbin; jk++) {
		  meanx[jk] =str_tdev[ixy]->GetBinContent(jk+1);
		  weightx[jk] =str_tdev1d[ixy]->GetBinContent(jk+1);
		  sumwt +=weightx[jk];
		  sum +=meanx[jk]*weightx[jk];
		}
		sum /=max(1.0,sumwt);

		file_outstr<<"//toffstr "<<nbin<<" "<<ixy<<" "<<iocc<<" "<<iiter<<" "<<sum<<endl;
		if (isTimeCorOrReso) {
		  if (ixy==0) {
		    for (int jk=0; jk<nbin; jk++) {
		      xtoffystr[iocc][jk] +=meanx[jk]-sum;
		      file_outstr<<xtoffystr[iocc][jk]<<", ";
		      if ((jk+1)%8==0) file_outstr<<endl;
		    }
		  } else {
		    for (int jk=0; jk<nbin; jk++) {
		      ytoffxstr[iocc][jk] +=meanx[jk]-sum;
		      file_outstr<<ytoffxstr[iocc][jk]<<", ";
		      if ((jk+1)%8==0) file_outstr<<endl;
		    }
		  }
		} //if (isTimeCorOrReso)
	      } // for (int ixy=0; ixy<2; ixy++)
	      file_outstr <<"};"<<endl;
	      //	      c2x->Update();
	      //	      if (c2x) { delete c2x; c2x=0;}

	      for (int ixy=0; ixy<2; ixy++) {
		if (str_tdev[ixy]) { delete str_tdev[ixy]; str_tdev[ixy]=0;}
		if (str_tdev1d[ixy]) { delete str_tdev1d[ixy]; str_tdev1d[ixy]=0;}
	      }
	    }

	    //	    ps.NewPage();
	    // gStyle->SetOptTitle(0);
// 	    gStyle->SetOptStat(0);
// 	    gStyle->SetOptFit(0);
// 	    gStyle->SetOptLogy(1);
// 	    gStyle->SetPadTopMargin(.001);
// 	    gStyle->SetPadBottomMargin(0.001); //0.07
// 	    gStyle->SetPadLeftMargin(0.001);
// 	    gStyle->SetPadRightMargin(0.001);
//	    TCanvas* c1 = new TCanvas("c1", "c1", 700, 900);
//	    c1->Divide(8,8,1.e-6, 1.e-6);

	    TH1F* time_shift[nstrip][2]={0};
	    double fitmean[nstrip][2]={0};
	    double fitrms[nstrip][2]={0};
	    double fitchi[nstrip][2]={0};
	    double statmean[nstrip][2]={0};

	    TF1* fity[nstrip][2]={0};

	    const int nsgpr=3;
	    double fitres[nsgpr];
	    double parerr[nsgpr];
	    double fchisq;
	    ps.NewPage();
	    gStyle->SetOptTitle(0);
	    gStyle->SetOptStat(0);
	    gStyle->SetOptFit(0);
	    gStyle->SetOptLogy(1);
	    gStyle->SetPadTopMargin(.001);
	    gStyle->SetPadBottomMargin(0.001); //0.07
	    gStyle->SetPadLeftMargin(0.001);
	    gStyle->SetPadRightMargin(0.001);

#ifdef ONEMBYONEM
	    TCanvas* c1 = new TCanvas("c1", "c1", 700, 900);
	    c1->Divide(8,8,1.e-6, 1.e-6);
#endif
            for (int ij=0; ij<2; ij++) {
#ifndef ONEMBYONEM
	      TCanvas* c1 = new TCanvas("c1", "c1", 700, 900);
	      c1->Divide(8,8,1.e-6, 1.e-6);
#endif

              for (int jk=0; jk<nstrip; jk++) {
#ifdef ONEMBYONEM
		c1->cd(nstrip*ij+jk+1);
#else
		c1->cd(jk+1);
#endif
		if (ij==0) {
		  time_shift[jk][ij] = (TH1F*)time_xstrreso[iocc][jk][iiterrs]->Clone();
		} else {
		  time_shift[jk][ij] = (TH1F*)time_ystrreso[iocc][jk][iiterrs]->Clone();
		}

		if (time_shift[jk][ij]->Integral()>2) {
		  time_shift[jk][ij]->GetXaxis()->SetLabelSize(.07);

		  double  par[3]={time_shift[jk][ij]->GetMaximum(),time_shift[jk][ij]->GetMean(), time_shift[jk][ij]->GetRMS()};

		  int nbinx = time_shift[jk][ij]->GetNbinsX();
		  int nlow=0;
                  for (int kl=0; kl<nbinx; kl++) {
		    if (time_shift[jk][ij]->GetBinContent(kl+1) >0) {
		      nlow=kl; break;
		    }
		  }
		  int nhig = 0;
                  for (int kl=nbinx; kl>0; kl--) {
		    if (time_shift[jk][ij]->GetBinContent(kl+1) >0) {
		      nhig=kl; break;
		    }
		  }
		  float amean = 0.5*(nlow + nhig);
		  nlow = int(TMath::Max(0., nlow - 1.0*(amean - nlow)));
		  nhig = int(TMath::Min(nbinx-1., nhig + 1.0*(nhig - amean)));

		  nchannel = 0;
                  for (int kl=nlow; kl<=nhig; kl++) {
		    if (nchannel <nmxchn) {
		      m_data[nchannel] = time_shift[jk][ij]->GetBinContent(kl+1);
		      m_xpos[nchannel] = time_shift[jk][ij]->GetBinCenter(kl+1);
		      nchannel++;
		    }
		  }
		  double alw= time_shift[jk][ij]->GetBinCenter(nlow+1);
		  double ahg = time_shift[jk][ij]->GetBinCenter(nhig+1);

		  //		  time_shift[jk][ij]->GetXaxis()->SetRangeUser(alw, ahg);
		  //		  time_shift[jk][ij]->GetXaxis()->SetLabelSize(.1);
		  time_shift[jk][ij]->Draw();

		  TMinuit *gMinuit = new TMinuit(nsgpr);
		  gMinuit->SetPrintLevel(-1);

		  TString hname[nsgpr] = {"height", "mean", "rms"};

		  int nmx = time_shift[jk][ij]->GetMaximumBin();
		  double hgh = 0.35*(time_shift[jk][ij]->GetBinContent(nmx-1) +
				     time_shift[jk][ij]->GetBinContent(nmx) +
				     time_shift[jk][ij]->GetBinContent(nmx+1));
		  double strt[nsgpr] = {hgh,time_shift[jk][ij]->GetMean(), max(0.6, min(3.0, 0.9*time_shift[jk][ij]->GetRMS()))};

		  double alow[nsgpr] = {0.5*strt[0], strt[1]-2.0, max(0.4*strt[2],0.5)};
		  double ahig[nsgpr] = {2.0*strt[0], strt[1]+2.0, min(1.3*strt[2]+0.1,3.5)};
		  double step[nsgpr] = {0.5, 0.01, 0.01};

		  gMinuit->SetFCN(fcnsg);

		  double arglist[10];
		  int ierflg = 0;
		  arglist[0] =  1 ;
		  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

                  for (int kl=0; kl<nsgpr; kl++) {
		    gMinuit->mnparm(kl, hname[kl], strt[kl], step[kl], alow[kl], ahig[kl],ierflg);
		  }

		  arglist[0] = 0;
		  //	      gMinuit->mnexcm("MIGRAD", arglist, 0, ierflg);
		  gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);

		  arglist[0] = 0;
		  gMinuit->mnexcm("IMPROVE", arglist, 0, ierflg);

		  TString chnam;
		  double parv,err,xlo,xup, plerr, mierr, eparab, gcc;
		  int iuit;

                  for (int kl=0; kl<nsgpr; kl++) {
		    gMinuit->mnpout(kl, chnam, parv, err, xlo, xup, iuit);
		    gMinuit->mnerrs(kl, plerr, mierr, eparab, gcc);
		    fitres[kl] = parv;
		    parerr[kl] = err;
		  }
		  double  fedm, errdef;
		  int  nparx, istat, fitndof;
		  gMinuit->mnstat(fchisq, fedm, errdef, fitndof, nparx, istat);

		  //		  cout <<"chisq "<< fchisq<<" "<<fitndof<<" "<<nparx<<" "<<istat<<" "<<fitres[0]<<" "<<fitres[1]<<" "<<fitres[2]<<" "<<time_shift[jk][ij]->GetMean()<<" "<<time_shift[jk][ij]->GetRMS()<<" "<<time_shift[jk][ij]->GetEntries()<<endl;
		  if (istat==0) fchisq =100000000.0;
		  //		  time_shift[jk][ij]->Fit("gaus");

		  sprintf(name, "fity_%i_%i", ij, jk);
		  //	fity[jk][ij] = new TF1(name, gausX, alw, ahg,  3);
		  //	fity[jk][ij]->SetParameters(par);
		  //	time_shift[jk][ij]->Fit(name);

		  //	time_shift[jk][ij]->Draw();

		  fity[jk][ij] = new TF1(name, gausX, alw, ahg, 3);
		  fity[jk][ij]->SetParameters(fitres);
		  fity[jk][ij]->SetLineColor(2);
		  fity[jk][ij]->SetLineWidth(1);
		  fity[jk][ij]->Draw("same");

		  fitmean[jk][ij] = fitres[1]; // fity[jk][ij]->GetParameter(1);
		  fitrms[jk][ij] = fitres[2];
		  fitchi[jk][ij] = fchisq;

		  statmean[jk][ij] = time_shift[jk][ij]->GetMean();

		  latex.DrawLatex(0.32, 0.36,Form("%g", int(1000*fitres[1])/1000.));
		  latex.DrawLatex(0.32, 0.26,Form("%g", int(1000*statmean[jk][ij])/1000.));
		  latex.DrawLatex(0.32, 0.16,Form("%g", int(1000*fitres[2])/1000.));
		  latex.DrawLatex(0.32, 0.06,Form("%g", int(1000*time_shift[jk][ij]->GetRMS())/1000.));
		  delete gMinuit; gMinuit=0;
		} // if (time_shift[jk][ij]->GetEntries()>15)
	      } // for (int jk=0; jk<nstrip; jk++)
#ifndef ONEMBYONEM
	      c1->Update();
	      if (c1) { delete c1; c1=0;}
	      if (ij==0) {ps.NewPage();}
#endif
	    } //for (int ij=0; ij<2; ij++)
#ifdef ONEMBYONEM
	      c1->Update();
	      if (c1) { delete c1; c1=0;}
#endif

            for (int ij=0; ij<2; ij++) {
              for (int jk=0; jk<nstrip; jk++) {
		if (time_shift[jk][ij]) { delete time_shift[jk][ij]; time_shift[jk][ij]=0;}
		if (fity[jk][ij]) { delete fity[jk][ij]; fity[jk][ij]=0;}
	      }
	    }
	    file_out <<"lay3 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<" "<<laye<<" "<<occulyr<<" "<<iocc<<" "<<nlayer<<endl;


#ifdef 	TIMESLOPE
	    if (ntcor==1 && iiter<nmxiter-isLast) { //Don't correct for last iteration

	      //a	      ps.NewPage();
	      //a	      gStyle->SetOptLogy(0);
	      //a	      TCanvas* c1x = new TCanvas("c1x", "c1x", 700, 900);
	      //a	      c1x->Divide(8,8,1.e-6, 1.e-6);

	      TH2F* time_slope_pr[2][nstrip]={0};
	      TProfile* time_slope_proj[2][nstrip]={0};
	      TH1F* time_slope_1d[2][nstrip]={0};
	      TF1* time_slope_prfit[2][nstrip]={0};

	      for (int ixy=0; ixy<2; ixy++) {
		//	      for (int ixy=1; ixy>=0; ixy--) {

		file_outstr<<"//slope_cor "<<ixy<<" "<<iocc<<" "<<iiter<<endl;
		for (int str=0; str<nstrip; str++) {

		  double meanx[nstrip]={0};
		  double weightx[nstrip]={0};
		  double sum=0;
		  double sumwt=0;

		  if (ixy==0) {
		    time_slope_pr[ixy][str] = (TH2F*) time_xslope_pr[iocc][str][iiter]->Clone();
		  } else {
		    time_slope_pr[ixy][str] = (TH2F*) time_yslope_pr[iocc][str][iiter]->Clone();
		  }
		  time_slope_pr[ixy][str]->GetXaxis()->SetLabelSize(.08);
		  time_slope_pr[ixy][str]->GetYaxis()->SetLabelSize(.08);
		  //a		  c1x->cd(ixy*nstrip+str+1);
		  //a		  if (iiter==0) { time_slope_pr[ixy][str]->Draw("colz");} //Fit("pol1");
		  sprintf(name, "xxx_%i_%i", ixy, str);
		  time_slope_proj[ixy][str] =(TProfile*)time_slope_pr[ixy][str]->ProfileX(name);
		  time_slope_proj[ixy][str]->SetMarkerSize(0.4);
		  time_slope_proj[ixy][str]->SetMarkerStyle(24);

		  //a		  time_slope_proj[ixy][str]->Draw((iiter==0)?"same":"");
		  sprintf(name, "xxx1d_%i_%i", ixy, str);
		  time_slope_1d[ixy][str] =(TH1F*)time_slope_pr[ixy][str]->ProjectionX(name);

		  int nbin = min(nstrip,time_slope_1d[ixy][str]->GetNbinsX());
		  const int nmnhit=3;
		  for (int jk=0; jk<nbin; jk++) {
		    meanx[jk] =time_slope_proj[ixy][str]->GetBinContent(jk+1);
		    weightx[jk] =time_slope_1d[ixy][str]->GetBinContent(jk+1);
		    if (weightx[jk]>nmnhit) {
		      sumwt +=weightx[jk];
		      sum +=meanx[jk]*weightx[jk];
		    }
		  }
		  sum /=max(1.0,sumwt);
		  file_outstr<<"//slope_corstr "<<ixy<<" "<<nbin<<" "<<str<<" "<<sum<<" "<<statmean[str][ixy]<<" "<<time_slope_pr[ixy][str]->GetEntries()<<endl;
		  if (isTimeCorOrReso && iiter<nmxiter-isLast) {
		    if (ixy==0) {
		      for (int jk=0; jk<nbin; jk++) {
			if (weightx[jk]>nmnhit) {xt_slope_cor[iocc][str][jk] +=meanx[jk]-sum;} //statmean[str][ixy];} //sum;}
			file_outstr<<xt_slope_cor[iocc][str][jk]<<", ";
			if ((jk+1)%8==0) file_outstr<<endl;
		      }
		    } else {
		      for (int jk=0; jk<nbin; jk++) {
			if (weightx[jk]>nmnhit) {yt_slope_cor[iocc][str][jk] +=meanx[jk]-sum;} //statmean[str][ixy];} // sum;}
			file_outstr<<yt_slope_cor[iocc][str][jk]<<", ";
			if ((jk+1)%8==0) file_outstr<<endl;
		      }
		    }
		  }

 		  double par[nDelayPar]={0,0,0,0,0,0};
 		  if (time_slope_pr[ixy][str]->GetEntries()>120) {
 		    int istr=1;
 		    int iend =31;
 		    int nbin = time_slope_pr[ixy][str]->GetNbinsX();
 		    for (int ix=nbin/2; ix>=0; ix--) {
 		      if (time_slope_proj[ixy][str]->GetBinError(ix+1) <1.e-4) {
 			istr = ix+1; break;
 		      }
 		    }

 		    for (int ix=nbin/2; ix<nbin; ix++) {
 		      if (time_slope_proj[ixy][str]->GetBinError(ix+1) <1.e-4) {
 			iend = ix-1; break;
 		      }
 		    }

 		    if (istr<nbin/2 && iend>nbin/2 ) {
 		      TFitResultPtr ptr = time_slope_proj[ixy][str]->Fit("pol4","SQ0","",istr, iend);
		      Int_t fitStatus = ptr;
		      if (fitStatus==0) {
			par[0] = ptr->Parameter(0);
			par[1] = ptr->Parameter(1);
			par[2] = ptr->Parameter(2);
			par[3] = ptr->Parameter(3);
			par[4] = ptr->Parameter(4);

			file_outstr <<" ixy "<< ixy<<" "<<iocc<<" "<<str<<" "
				    <<par[0]<<" " <<par[1]<<" " <<par[2]<<" " <<par[3]<<" " <<par[4]<<" "<<
			  par[0]+par[1]*16+par[2]*256+par[3]*4096+par[4]*65536<<" "<<sum<<endl; //central strip
		      }
 		    } else {
 		      file_outstr <<" ixx "<< ixy<<" "<<iocc<<" "<<str<<" "<<par[0]<<" " <<par[1]<<" " <<par[2]<<" " <<par[3]<<" " <<par[4]<<" "<<par[5]<<" "<<sum<<endl;
 		    }
 		  } else {
 		    file_outstr <<" iyy "<< ixy<<" "<<iocc<<" "<<str<<" "<<par[0]<<" " <<par[1]<<" " <<par[2]<<" " <<par[3]<<" " <<par[4]<<" "<<par[5]<<" "<<sum<<endl;
 		  }

 		  time_slope_proj[ixy][str]->SetMarkerSize(0.6);
 		  time_slope_proj[ixy][str]->SetMarkerStyle(24);
		  //a 		  time_slope_proj[ixy][str]->Draw("same");
 		  sprintf(name, "pedfun_%i",str);
 		  time_slope_prfit[ixy][str] = new TF1(name, polfunc,0.,31.,nDelayPar-1);
 		  time_slope_prfit[ixy][str]->SetParameters(par);
 		  time_slope_prfit[ixy][str]->SetLineColor(2);
 		  time_slope_prfit[ixy][str]->SetLineWidth(1);
		  //a 		  time_slope_prfit[ixy][str]->Draw("same");

 		  //		  time_slope_pr[ixy][str]->ProfileX()->Draw("same");

		} // for (int str=0; str=nstrip; str++)
	      } // for (int ixy=0; ixy<2; ixy++)
	      //a	      c1x->Update();

	      for (int ixy=0; ixy<2; ixy++) {
		for (int str=0; str<nstrip; str++) {
		  if (time_slope_pr[ixy][str]) {
		    delete time_slope_pr[ixy][str]; time_slope_pr[ixy][str]=0;
		  }
		  if (time_slope_proj[ixy][str]) {
		    delete time_slope_proj[ixy][str]; time_slope_proj[ixy][str]=0;
		  }
		  if (time_slope_prfit[ixy][str]) {
		    delete time_slope_prfit[ixy][str]; time_slope_prfit[ixy][str]=0;
		  }
		}
	      } // for (int ixy=0; ixy<2; ixy++)
	      //a	      if (c1x) { delete c1x; c1x=0;}
	    } // if (ntcor==1)
#endif


	    //	    if (iiter<nmxiter && iocc<nlayer) {

	    if (ntcor==1 && iiter<nmxiter-isLast) { //Don't correct for last iteration //150112
	      bias_intime_xreso2[iocc] = (xtime_exterr[iocc][iiterrs]->GetMean())*(xtime_exterr[iocc][iiterrs]->GetMean())+(xtime_exterr[iocc][iiterrs]->GetRMS())*(xtime_exterr[iocc][iiterrs]->GetRMS());//jim
	      bias_intime_yreso2[iocc] = (ytime_exterr[iocc][iiterrs]->GetMean())*(ytime_exterr[iocc][iiterrs]->GetMean())+(ytime_exterr[iocc][iiterrs]->GetRMS())*(ytime_exterr[iocc][iiterrs]->GetRMS());//jim

	      TFitResultPtr ptrtimexSigma1 = time_xreso[iocc][iiterrs]->Fit("gaus","SQ");
	      TFitResultPtr ptrtimeySigma1 = time_yreso[iocc][iiterrs]->Fit("gaus","SQ");

	      timeserrx2[iocc] =  (ptrtimexSigma1->Parameter(2)*ptrtimexSigma1->Parameter(2)) - bias_intime_xreso2[iocc];
	      timeserry2[iocc] =  (ptrtimeySigma1->Parameter(2)*ptrtimeySigma1->Parameter(2)) - bias_intime_yreso2[iocc];

	      //	      timeserrx2[iocc] =  (time_xreso[iocc][iiterrs]->GetRMS())*time_xreso[iocc][iiterrs]->GetRMS() - bias_intime_xreso2[iocc];     //jim jim
	      //	      timeserry2[iocc] =  (time_yreso[iocc][iiterrs]->GetRMS())*time_yreso[iocc][iiterrs]->GetRMS() - bias_intime_yreso2[iocc];     //jim jim

	      if (timeserrx2[iocc]<0.25) timeserrx2[iocc]=0.25; //resolution should be more than 0.5 ns
	      if (timeserry2[iocc]<0.25) timeserry2[iocc]=0.25;
	    }

	    if (iiter<nmxiter && iocc<nlayer) {
	      int nbinxx = min(time_xstrreso[iocc][0][iiterrs]->GetNbinsX(),time_ystrreso[iocc][0][iiterrs]->GetNbinsX()) +1;
              for (int jk=0; jk<nstrip; jk++) {
		//	timeinuse[0][iocc][jk] = timeinuse[1][iocc][jk] = false;
		if (isalign >0 ) {
		  double ncont= time_xstrreso[iocc][jk][iiterrs]->Integral();
		  double tmpent= max(1.0, time_xstrreso[iocc][jk][iiterrs]->GetEntries());
		  double nunder=time_xstrreso[iocc][jk][iiterrs]->GetBinContent(0)/tmpent;
		  double nover=time_xstrreso[iocc][jk][iiterrs]->GetBinContent(nbinxx)/tmpent;
		  double statmean = time_xstrreso[iocc][jk][iiterrs]->GetMean();
		  double statrms = time_xstrreso[iocc][jk][iiterrs]->GetRMS();
		  double diff = abs(time_xstrreso[iocc][jk][iiterrs]->GetMean()-fitmean[jk][0]);
		  double width = sqrt(max(0.25, fitrms[jk][0]*fitrms[jk][0] -  bias_intime_xreso2[iocc] )); //time_xstrreso[iocc][jk][iiterrs]->GetRMS();

		  if (ncont >nmnentry && iocc>=xtcorstr && iocc<=xtcorend && ntcor==1) {
		    if (iiter<nmxiter-isLast) { //Do not update for last iteration
		      timeinuse[0][iocc][jk] = true;
		      if (width<absWidth && diff<widthScale*width) {
			if (isTimeCorOrReso) {
			  xtoffset[iocc][jk] +=fitmean[jk][0];//statmean; // time_xstrreso[iocc][jk][iiterrs]->GetMean(); // fitmean[jk][0];
			  //			} else if (diff<0.50*width) {
			  //			  xtoffset[iocc][jk] +=0.5*(statmean+fitmean[jk][0]);
			}
		      } else {
			timeinuse[0][iocc][jk] = false;
		      }
		    }
		    correction_xtime[iocc]->Fill(jk, iiter, fitmean[jk][0]); //(statmean);
		    fitted_rms_xtime[iocc]->Fill(jk, iiter, fitrms[jk][0]);

		    shift_time_mnft[iiter]->Fill(nstrip*iocc+jk, statmean-fitmean[jk][0]);
		    statmean_time[iiter]->Fill(nstrip*iocc+jk, statmean);
		    statrms_time[iiter]->Fill(nstrip*iocc+jk, statrms);
		    statskew_time[iiter]->Fill(nstrip*iocc+jk, time_xstrreso[iocc][jk][iiterrs]->GetSkewness());
		    statkurt_time[iiter]->Fill(nstrip*iocc+jk, time_xstrreso[iocc][jk][iiterrs]->GetKurtosis());

		    time_offset[iiter]->Fill(nstrip*iocc+jk, xtoffset[iocc][jk]);
		    time_entry[iiter]->Fill(nstrip*iocc+jk, ncont);
		    time_underflow[iiter]->Fill(nstrip*iocc+jk, nunder);
		    time_overflow[iiter]->Fill(nstrip*iocc+jk, nover);

		    rms_time[iiter]->Fill(nstrip*iocc+jk, fitrms[jk][0]);
		    if (!timeinuse[0][iocc][jk]) rms_timeused[iiter]->Fill(nstrip*iocc+jk, fitrms[jk][0]);

		    //////////////////////////////////
		    shift_time_mnftx[iiter]->Fill(statmean-fitmean[jk][0]);
		    statmean_timex[iiter]->Fill(statmean);
		    statrms_timex[iiter]->Fill(statrms);
		    statskew_timex[iiter]->Fill(time_xstrreso[iocc][jk][iiterrs]->GetSkewness());
		    statkurt_timex[iiter]->Fill(time_xstrreso[iocc][jk][iiterrs]->GetKurtosis());

		    if (abs(xtoffset[iocc][jk])>1.e-6) time_offsetx[iiter]->Fill(xtoffset[iocc][jk]);
		    rms_timex[iiter]->Fill(fitrms[jk][0]);
		    if (!timeinuse[0][iocc][jk]) rms_timeusedx[iiter]->Fill(fitrms[jk][0]);
		  }

		  ncont=time_ystrreso[iocc][jk][iiterrs]->Integral();
		  tmpent= max(1.0, time_ystrreso[iocc][jk][iiterrs]->GetEntries());
		  nunder=time_ystrreso[iocc][jk][iiterrs]->GetBinContent(0)/tmpent;
		  nover=time_ystrreso[iocc][jk][iiterrs]->GetBinContent(nbinxx)/tmpent;

		  statmean = time_ystrreso[iocc][jk][iiterrs]->GetMean();
		  statrms = time_ystrreso[iocc][jk][iiterrs]->GetRMS();
		  diff = abs(time_ystrreso[iocc][jk][iiterrs]->GetMean()-fitmean[jk][1]);
		  width = sqrt(max(0.25, fitrms[jk][1]*fitrms[jk][1] -  bias_intime_yreso2[iocc]));

		  if (ncont>nmnentry && iocc>=ytcorstr && iocc<=ytcorend && ntcor==1) {
		    if (iiter<nmxiter-isLast) { //Do not update for last iteration
		      timeinuse[1][iocc][jk] = true;
		      if (width<absWidth && diff<widthScale*width) {
			if (isTimeCorOrReso) {
			  ytoffset[iocc][jk] +=fitmean[jk][1];//statmean; // time_ystrreso[iocc][jk][iiterrs]->GetMean(); // fitmean[jk][1];
			//			} else if (diff<0.50*width) {
			//			  ytoffset[iocc][jk] +=0.5*(statmean+fitmean[jk][1]);
			}
		      } else {
			timeinuse[1][iocc][jk] = false;
			//			  ytoffset[iocc][jk] +=statmean;
		      }
		    }

		    correction_ytime[iocc]->Fill(jk, iiter, fitmean[jk][1]); //(statmean);
		    fitted_rms_ytime[iocc]->Fill(jk, iiter, fitrms[jk][1]);

		    shift_time_mnft[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, statmean-fitmean[jk][1]);
		    statmean_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, statmean);
		    statrms_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, statrms);
		    statskew_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, time_ystrreso[iocc][jk][iiterrs]->GetSkewness());
		    statkurt_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, time_ystrreso[iocc][jk][iiterrs]->GetKurtosis());

		    time_offset[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, ytoffset[iocc][jk]);
		    time_entry[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, ncont);
		    time_underflow[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, nunder);
		    time_overflow[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, nover);

		    rms_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, fitrms[jk][1]);
		    if (!timeinuse[1][iocc][jk]) rms_timeused[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, fitrms[jk][1]);

		    shift_time_mnfty[iiter]->Fill(statmean-fitmean[jk][1]);
		    statmean_timey[iiter]->Fill(statmean);
		    statrms_timey[iiter]->Fill(statrms);
		    statskew_timey[iiter]->Fill(time_ystrreso[iocc][jk][iiterrs]->GetSkewness());
		    statkurt_timey[iiter]->Fill(time_ystrreso[iocc][jk][iiterrs]->GetKurtosis());

		    if (abs(ytoffset[iocc][jk])>1.e-6) time_offsety[iiter]->Fill(ytoffset[iocc][jk]);
		    rms_timey[iiter]->Fill(fitrms[jk][1]);
		    if (!timeinuse[1][iocc][jk]) rms_timeusedy[iiter]->Fill(fitrms[jk][1]);
		  }
		} //if (isalign >0 )
		if (ntcor==1 && iiter<nmxiter-isLast) { //Don't update with last iteration
		//		  xtmean[iocc][jk] =fitmean[jk][0];
		//		  ytmean[iocc][jk] =fitmean[jk][1];

		  xtrms[iocc][jk] =fitrms[jk][0];
		  ytrms[iocc][jk] =fitrms[jk][1];
		}
	      } //for (int jk=0; jk<nstrip; jk++)

	      time_mean_reso->Fill(iocc, iiterrs, time_xreso[iocc][iiterrs]->GetMean());
	      time_rms_reso->Fill(iocc, iiterrs, time_xreso[iocc][iiterrs]->GetRMS());
	      time_mean_reso->Fill(nlayer+iocc, iiterrs, time_yreso[iocc][iiterrs]->GetMean());
	      time_rms_reso->Fill(nlayer+iocc, iiterrs, time_yreso[iocc][iiterrs]->GetRMS());

	      time_corrms_reso->Fill(iocc, iiterrs, sqrt(timeserrx2[iocc]));
	      time_corrms_reso->Fill(nlayer+iocc, iiterrs, sqrt(timeserry2[iocc]));

	      time_exterr_reso->Fill(iocc, iiterrs, sqrt(bias_intime_xreso2[iocc]));
	      time_exterr_reso->Fill(nlayer+iocc, iiterrs, sqrt(bias_intime_yreso2[iocc]));

	      file_out <<"time+pos " <<iiterrs<<" iocc "<<iocc<<" "
		       << setw(6)<<xoff[iocc]<<"+-"<<setw(6)<<xrms[iocc]<<" "
		       << setw(6)<<yoff[iocc]<<"+-"<<setw(6)<<yrms[iocc]<<" X "
		       << setw(6)<<time_xrms[iocc]<<" "<<xtime_exterr[iocc][iiterrs]->GetMean()<<" "
		       << setw(6)<<time_xreso[iocc][iiterrs]->GetMean()<<"+-"
		       << setw(6)<<time_xreso[iocc][iiterrs]->GetRMS()<<" "<<sqrt(timeserrx2[iocc])<<" Y "
		       << setw(6)<<time_yrms[iocc]<<" "<<ytime_exterr[iocc][iiterrs]->GetMean()<<" "
		       << setw(6)<<time_yreso[iocc][iiterrs]->GetMean()<<"+-"
		       << setw(6)<<time_yreso[iocc][iiterrs]->GetRMS()<<" "<<sqrt(timeserry2[iocc])<<endl;

	      file_out <<"// X-rms "<<iocc<<" "<< iiter<<endl;
              for (int jk=0; jk<nstrip; jk++) {
		file_out <<xtrms[iocc][jk]<<", ";
		if ((jk+1)%8==0) file_out<<endl;
	      }

	      file_out <<"// Y-rms "<<iocc<<" "<< iiter<<endl;
              for (int jk=0; jk<nstrip; jk++) {
		file_out <<ytrms[iocc][jk]<<", ";
		if ((jk+1)%8==0) file_out<<endl;
	      }

	      file_out <<"// X-str "<<iocc<<" "<< iiter<<" "<<iiterrs<<endl;
              for (int jk=0; jk<nstrip; jk++) {
		file_out <<xtoffset[iocc][jk]<<", ";
		if ((jk+1)%8==0) file_out<<endl;
	      }

	      file_out <<"// Y-str "<<iocc<<" "<< iiter<<" "<<iiterrs<<endl;

              for (int jk=0; jk<nstrip; jk++) {
		file_out <<ytoffset[iocc][jk]<<", ";
		if ((jk+1)%8==0) file_out<<endl;
	      }

	      cout <<"// X-str "<<iocc<<" "<< iiter<<" "<<iiterrs<<endl;
              for (int jk=0; jk<nstrip; jk++) {
		if (abs(xtoffset[iocc][jk])>40.0) cout <<" xtoffset "<<jk <<" "<< xtoffset[iocc][jk]<<" "<< fitmean[jk][0]<<endl;
	      }
              for (int jk=0; jk<nstrip; jk++) {
		cout <<xtoffset[iocc][jk]<<", ";
		if ((jk+1)%8==0) cout<<endl;
	      }
	      cout <<"// Y-str "<<iocc<<" "<< iiter<<" "<<iiterrs<<endl;
              for (int jk=0; jk<nstrip; jk++) {
		if (abs(ytoffset[iocc][jk])>40.0) cout <<" ytoffset "<<jk <<" "<< ytoffset[iocc][jk]<<" "<< fitmean[jk][1]<<endl;
	      }
              for (int jk=0; jk<nstrip; jk++) {
		cout <<ytoffset[iocc][jk]<<", ";
		if ((jk+1)%8==0) cout<<endl;
	      }
	    } //if (isalign >0 && iiter<nmxiter-1 && iocc<nlayer)
	  } //if (isTiming)
	} // for (int iocc = istr; iocc<iend; iocc++)
      } // for (int laye=0; laye<nlayerit; laye++)
    } // or (int ntcor=0; ntcor< (iiter<nmxiter-1) ? 2 : 1; ntcor++)

    if (isalign>0) {
      for (int ij=0; ij<nlayer; ij++) {
	shift_pos->Fill(ij, iiter, xoff[ij]);
	shift_pos->Fill(nlayer+ij, iiter, yoff[ij]);

	rms_pos->Fill(ij, iiter, xrms[ij]);
	rms_pos->Fill(nlayer+ij, iiter, yrms[ij]);

      }
    }
  } // for (int iiter=0; iiter<nmxiter; iiter++)
  if (isalign>0) {
    fileOut_align->cd();
    align->Fill();
    fileOut->cd();

    file_out <<endl;
    file_out_h<<"double xoff[nlayer] = {";
    for (int ij=0; ij<nlayer; ij++) {
      if ( ij<nlayer-1) {
	file_out_h<<xoff[ij]<<", ";
      } else {
	file_out_h<<xoff[ij]<<"};"<<endl;
      }
    }

    file_out_h<<"double yoff[nlayer] = {";
    for (int ij=0; ij<nlayer; ij++) {
      if ( ij<nlayer-1) {
	file_out_h<<yoff[ij]<<", ";
      } else {
	file_out_h<<yoff[ij]<<"};"<<endl;
      }
    }
    file_out_h<<endl;
#ifdef NEWRPCDAQ1
    file_out_h<<"double align_xstr_ydev[nlayer][nPosAlignPar] ={";
    for (int ij=0; ij<nlayer; ij++) {
      file_out_h<<"{";
      for (int ix=0; ix<nPosAlignPar; ix++) {
	if (ix<nPosAlignPar-1) {
	  file_out_h<<align_xstr_ydev[ij][ix]<<", ";
	} else {
	  file_out_h<<align_xstr_ydev[ij][ix];
	}
      }
      if (ij<nlayer-1) { file_out_h<<"},"<<endl;} else { file_out_h<<"}";}
    }
    file_out_h<<"};"<<endl;

    file_out_h<<"double align_ystr_xdev[nlayer][nPosAlignPar] ={";
    for (int ij=0; ij<nlayer; ij++) {
      file_out_h<<"{";
      for (int ix=0; ix<nPosAlignPar; ix++) {
	if (ix<nPosAlignPar-1) {
	  file_out_h<<align_ystr_xdev[ij][ix]<<", ";
	} else {
	  file_out_h<<align_ystr_xdev[ij][ix];
	}
      }
      if (ij<nlayer-1) { file_out_h<<"},"<<endl;} else { file_out_h<<"}";}
    }
    file_out_h<<"};"<<endl;

    file_out_h<<"double align_xstr_xdev[nlayer][nPosAlignPar] ={";
    for (int ij=0; ij<nlayer; ij++) {
      file_out_h<<"{";
      for (int ix=0; ix<nPosAlignPar; ix++) {
	if (ix<nPosAlignPar-1) {
	  file_out_h<<align_xstr_xdev[ij][ix]<<", ";
	} else {
	  file_out_h<<align_xstr_xdev[ij][ix];
	}
      }
      if (ij<nlayer-1) { file_out_h<<"},"<<endl;} else { file_out_h<<"}";}
    }
    file_out_h<<"};"<<endl;

    file_out_h<<"double align_ystr_ydev[nlayer][nPosAlignPar] ={";
    for (int ij=0; ij<nlayer; ij++) {
      file_out_h<<"{";
      for (int ix=0; ix<nPosAlignPar; ix++) {
	if (ix<nPosAlignPar-1) {
	  file_out_h<<align_ystr_ydev[ij][ix]<<", ";
	} else {
	   file_out_h<<align_ystr_ydev[ij][ix];
	}
      }
      if (ij<nlayer-1) { file_out_h<<"},"<<endl;} else { file_out_h<<"}";}
    }
    file_out_h<<"};"<<endl;
#endif
    file_out_h<<"double xposerrsq[nmxhits][nlayer] = {";
    for (int ix=0; ix<nmxhits; ix++) {
      file_out_h<<"{";
      for (int ij=0; ij<nlayer; ij++) {
	if ( ij<nlayer-1) {
	  file_out_h<<xposerrsq[ix][ij]<<", ";
	} else if (ix<nmxhits-1) {
	  file_out_h<<xposerrsq[ix][ij]<<"},"<<endl;
	} else {
	  file_out_h<<xposerrsq[ix][ij]<<"}};"<<endl;
	}
      }
    }

    file_out_h<<"double yposerrsq[nmxhits][nlayer] = {";
    for (int ix=0; ix<nmxhits; ix++) {
      file_out_h<<"{";
      for (int ij=0; ij<nlayer; ij++) {
	if ( ij<nlayer-1) {
	  file_out_h<<yposerrsq[ix][ij]<<", ";
	} else if (ix<nmxhits-1) {
	  file_out_h<<yposerrsq[ix][ij]<<"},"<<endl;
	} else {
	  file_out_h<<yposerrsq[ix][ij]<<"}};"<<endl;
	}
      }
    }




    if (isTiming) {
      file_out_h<<"double timeoffsetx[nlayer] = {";
      for (int ij=0; ij<nlayer; ij++) {
	if (ij<nlayer-1) {
	  file_out_h<<timeoffsetx[ij]<<", ";
	} else {
	  file_out_h<<timeoffsetx[ij]<<"};"<<endl;
	}
      }

      file_out_h<<"double timeoffsety[nlayer] = {";
      for (int ij=0; ij<nlayer; ij++) {
	if ( ij<nlayer-1) {
	  file_out_h<<timeoffsety[ij]<<", ";
	} else {
	  file_out_h<<timeoffsety[ij]<<"};"<<endl;
	}
      }
      file_out_h<<endl;

      file_out_h<<"double timeserrx2[nlayer] = {";
      for (int ij=0; ij<nlayer; ij++) {
	if ( ij<nlayer-1) {
	  file_out_h<<timeserrx2[ij]<<", ";
	} else {
	  file_out_h<<timeserrx2[ij]<<"};"<<endl;
	}
      }
      file_out_h<<endl;

      file_out_h<<"double timeserry2[nlayer] = {";
      for (int ij=0; ij<nlayer; ij++) {
	if ( ij<nlayer-1) {
	  file_out_h<<timeserry2[ij]<<", ";
	} else {
	  file_out_h<<timeserry2[ij]<<"};"<<endl;
	}
      }
      file_out_h<<endl;


      file_out_h <<"double xtoffset[nlayer][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
	file_out_h <<"// X-str "<<ij<<endl;
        for (int jk=0; jk<nstrip; jk++) {
	  if (ij<nlayer-1 || jk<nstrip-1) {
	    file_out_h <<xtoffset[ij][jk]<<", ";
	  } else {
	    file_out_h <<xtoffset[ij][jk]<<"};";
	  }
	  if ((jk+1)%8==0) file_out_h<<endl;
	}
      }
      file_out_h<<endl;

      file_out_h <<"double ytoffset[nlayer][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
	file_out_h <<"// Y-str "<<ij<<endl;
        for (int jk=0; jk<nstrip; jk++) {
	  if (ij<nlayer-1 || jk<nstrip-1) {
	    file_out_h <<ytoffset[ij][jk]<<", ";
	  } else {
	    file_out_h <<ytoffset[ij][jk]<<"};";
	  }
	  if ((jk+1)%8==0) file_out_h<<endl;
	}
      }

      file_out_h <<"double xtoffystr[nlayer][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
	file_out_h <<"// X-str "<<ij<<endl;
        for (int jk=0; jk<nstrip; jk++) {
	  if (ij<nlayer-1 || jk<nstrip-1) {
	    file_out_h <<xtoffystr[ij][jk]<<", ";
	  } else {
	    file_out_h <<xtoffystr[ij][jk]<<"};";
	  }
	  if ((jk+1)%8==0) file_out_h<<endl;
	}
      }
      file_out_h<<endl;

      file_out_h <<"double ytoffxstr[nlayer][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
	file_out_h <<"// Y-str "<<ij<<endl;
        for (int jk=0; jk<nstrip; jk++) {
	  if (ij<nlayer-1 || jk<nstrip-1) {
	    file_out_h <<ytoffxstr[ij][jk]<<", ";
	  } else {
	    file_out_h <<ytoffxstr[ij][jk]<<"};";
	  }
	  if ((jk+1)%8==0) file_out_h<<endl;
	}
      }
      file_out_h<<endl;

      file_out_h <<"double xt_slope_cor[nlayer][nstrip][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
	for (int jk=0; jk<nstrip; jk++) {
	  for (int kl=0; kl<nstrip; kl++) {
	    if (ij<nlayer-1 || jk<nstrip-1 || kl<nstrip-1) {
	      file_out_h <<xt_slope_cor[ij][jk][kl]<<", ";
	    } else {
	      file_out_h <<xt_slope_cor[ij][jk][kl]<<"};";
	    }
	  }
	  file_out_h <<"// X-str "<<ij<<" "<<jk<<endl;
	}
      }
      file_out_h<<endl;


      file_out_h <<"double yt_slope_cor[nlayer][nstrip][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
	for (int jk=0; jk<nstrip; jk++) {
	  for (int kl=0; kl<nstrip; kl++) {
	    if (ij<nlayer-1 || jk<nstrip-1 || kl<nstrip-1) {
	      file_out_h <<yt_slope_cor[ij][jk][kl]<<", ";
	    } else {
	      file_out_h <<yt_slope_cor[ij][jk][kl]<<"};";
	    }
	  }
	  file_out_h <<"// Y-str "<<ij<<" "<<jk<<endl;
	}
      }
      file_out_h<<endl;
    } // if (isTiming)
  } //  if (isalign>0)

  //  ps.Close();
  fileOut->cd();


#ifdef ISEFFICIENCY
  for (int ij=0; ij<nlayer-2; ij++) {
    for (int ix=0; ix<itset; ix++) {
      inefficiency_corx_set[ij][ix]->Divide(totalentry_set[ij][ix]);
      triggereffi_x_set[ij][ix]->Divide(totalentry_set[ij][ix]);
      triggereffi_y_set[ij][ix]->Divide(totalentry_set[ij][ix]);
    }

    for (int jk=0; jk<2*nmxiter; jk++) {
      inefficiency_corx[ij][jk]->Divide(totalentry[ij][jk]);
      inefficiency_cory[ij][jk]->Divide(totalentry[ij][jk]);

      inefficiencytrue_corx[ij][jk]->Divide(totalentry[ij][jk]);
      inefficiencytrue_cory[ij][jk]->Divide(totalentry[ij][jk]);

      inefficiency_uncx[ij][jk]->Divide(totalentry[ij][jk]);
      inefficiency_uncy[ij][jk]->Divide(totalentry[ij][jk]);

      efficiency_xpixel[ij][jk]->Divide(efficiency_xallpixel[ij][jk]);
      efficiency_ypixel[ij][jk]->Divide(efficiency_yallpixel[ij][jk]);

      triggereffi_x[ij][jk]->Divide(totalentry[ij][jk]);
      triggereffi_y[ij][jk]->Divide(totalentry[ij][jk]);

      fine_triggereffi_x[ij][jk]->Divide(fine_totalentry[ij][jk]);
      fine_triggereffi_y[ij][jk]->Divide(fine_totalentry[ij][jk]);

      for (int kl=0; kl<=totalentry[ij][jk]->GetNbinsX()+1; kl++) {
	for (int lm=0; lm <=totalentry[ij][jk]->GetNbinsY()+1; lm++) {
	  double antot = totalentry[ij][jk]->GetBinContent(kl, lm);
	  if (antot <1) { antot = 1;};
	  double xtot = triggereffi_xevt[ij][jk]->GetBinContent(kl, lm);
	  double ytot = triggereffi_yevt[ij][jk]->GetBinContent(kl, lm);
	  double effix = xtot/antot;
	  if(effix>1.) { effix=1.;}
	  if(effix<0.) { effix=0.;}
	  double erreffix = pow(effix*(1-effix)/antot, 0.5);

	  double effiy = ytot/antot;
	  if(effiy>1.) { effiy=1.;}
	  if(effiy<0.) { effiy=0.;}
	  double erreffiy = pow(effiy*(1-effiy)/antot, 0.5);

	  triggereffi_xevt[ij][jk]->SetBinContent(kl, lm, effix);
	  triggereffi_xevt[ij][jk]->SetBinError(kl, lm, erreffix);

	  triggereffi_yevt[ij][jk]->SetBinContent(kl, lm, effiy);
	  triggereffi_yevt[ij][jk]->SetBinError(kl, lm, erreffiy);
	}
      }

     //   triggereffi_xevt[ij][jk]->Divide(totalentry[ij][jk]);
     //  triggereffi_yevt[ij][jk]->Divide(totalentry[ij][jk]);

      inefficiency_xt[ij][jk]->Divide(totalentry[ij][jk]); // (total_xt[ij][jk]);
      inefficiency_yt[ij][jk]->Divide(totalentry[ij][jk]); // (total_yt[ij][jk]);

      difefficiency_uncx[ij][jk]->Add(inefficiency_uncx[ij][jk], defefficiency_uncx[ij], 1., -1.);
      difefficiency_uncy[ij][jk]->Add(inefficiency_uncy[ij][jk], defefficiency_uncy[ij], 1., -1.);

      diftriggereffi_x[ij][jk]->Add(triggereffi_x[ij][jk], deftriggereffi_x[ij], 1., -1.);
      diftriggereffi_y[ij][jk]->Add(triggereffi_y[ij][jk], deftriggereffi_y[ij], 1., -1.);

      diftriggereffi_xevt[ij][jk]->Add(triggereffi_xevt[ij][jk], deftriggereffi_x[ij], 1., -1.);
      diftriggereffi_yevt[ij][jk]->Add(triggereffi_yevt[ij][jk], deftriggereffi_y[ij], 1., -1.);

      difefficiency_xt[ij][jk]->Add(inefficiency_uncx[ij][jk], inefficiency_xt[ij][jk], 1., -1.);
      difefficiency_yt[ij][jk]->Add(inefficiency_uncy[ij][jk], inefficiency_yt[ij][jk], 1., -1.);

    }
  }
  if (isTiming) {
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<nmxiter+6; jk++) {
	xystr_xtdev[ij][jk]->Divide(nxystr_xtdev[ij][jk]);
	xystr_ytdev[ij][jk]->Divide(nxystr_ytdev[ij][jk]);
      }
    }
  }

  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_corX[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_corx[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// Xlayer="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl;
        for (int iy=0; iy<inefficiency_corx[ij][jk]->GetNbinsY(); iy++) {
	  file_out <<setprecision(4)<<inefficiency_corx[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;

    file_out<<"  double ineffi_table_corY[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_cory[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// Ylayer="<<ij<<" iter "<<jk<<" ystrip "<<ix<<endl;
        for (int iy=0; iy<inefficiency_cory[ij][jk]->GetNbinsY(); iy++) {
	  file_out <<setprecision(4)<<inefficiency_cory[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }

  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"Combined NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_uncX[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncx[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// Xlayer="<<ij<<" iter "<<jk<<"unc xstrip "<<ix<<endl;
        for (int iy=0; iy<inefficiency_uncx[ij][jk]->GetNbinsY(); iy++) {
	  file_out <<setprecision(4)<<inefficiency_uncx[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;

    file_out<<"  double ineffi_table_uncY[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncy[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// Ylayer="<<ij<<" iter "<<jk<<"unc ystrip "<<ix<<endl;
        for (int iy=0; iy<inefficiency_uncy[ij][jk]->GetNbinsY(); iy++) {
	  file_out <<setprecision(4)<<inefficiency_uncy[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }

  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"NMXITER Timing ===================="<<jk<<endl;
    file_out<<"  double ineffi_table_X_time[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_xt[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// XlayerTime ineffi="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl;
        for (int iy=0; iy<inefficiency_xt[ij][jk]->GetNbinsY(); iy++) {
	  //Remove position ininefficincy to have extra inefficiency due to timing
	  //Note that inefficiency_xt is an EFIICIENCY but inefficiency_cor/uncx is INEFFICIENCY
	  file_out <<setprecision(4)<<
	    inefficiency_xt[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	  //	    -inefficiency_corx[ij][jk]->GetBinContent(ix+1, iy+1)
	  //	    -inefficiency_uncx[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;

    file_out<<"  double ineffi_table_Y_time[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_yt[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// YlayerTime ineffi="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl;
        for (int iy=0; iy<inefficiency_yt[ij][jk]->GetNbinsY(); iy++) {
	  file_out <<setprecision(4)<<
	    inefficiency_yt[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	  //	    -inefficiency_cory[ij][jk]->GetBinContent(ix+1, iy+1)
	  //	    -inefficiency_uncy[ij][jk]->GetBinContent(ix+1, iy+1) <<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }


  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"//NMXITER Trigger========================"<<jk<<endl;
    file_out<<"  double effi_trig_X[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_x[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// Xlayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl;
        for (int iy=0; iy<triggereffi_x[ij][jk]->GetNbinsY(); iy++) {
	  file_out <<setprecision(4)<<triggereffi_x[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;

    file_out<<"  double effi_trig_Y[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_y[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// Ylayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl;
        for (int iy=0; iy<triggereffi_y[ij][jk]->GetNbinsY(); iy++) {
	  file_out <<setprecision(4)<<triggereffi_y[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }

  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"//NMXITER TriggerEVT========================"<<jk<<endl;
    file_out<<"  double effi_trig_xevt[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_xevt[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// Xlayer TrigEVT="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl;
        for (int iy=0; iy<triggereffi_xevt[ij][jk]->GetNbinsY(); iy++) {
	  file_out <<setprecision(4)<<triggereffi_xevt[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;

    file_out<<"  double effi_trig_yevt[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_yevt[ij][jk]->GetNbinsX(); ix++) {
	file_out<<"// Ylayer TrigEVT="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl;
        for (int iy=0; iy<triggereffi_yevt[ij][jk]->GetNbinsY(); iy++) {
	  file_out <<setprecision(4)<<triggereffi_yevt[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
	}
	file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }

#endif

  int nbinsx = costhe[0]->GetNbinsX();

  for (int ix=0; ix<5; ix++) {
    file_out <<"double ExpCount"<<ix<<"["<<nbinsx<<"]={";
    for (int ij=0; ij<nbinsx; ij++) {
      file_out <<setprecision(5)<<costhe[ix]->GetBinContent(ij+1)<<",";
    }
    file_out <<"};"<<endl;
  }

  int nbinsxx = costhe[5]->GetNbinsX();//10Nov
  file_out <<"double ExpCount"<<5<<"["<<nbinsxx<<"]={";
  for (int ij=0; ij<nbinsxx; ij++) {
    file_out <<setprecision(5)<<costhe[5]->GetBinContent(ij+1)<<",";
  }
  file_out <<"};"<<endl;

  if (isTiming) {
    file_out<<"delta T covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<deltatcov[ix][iy]<<", \t";
      }
      file_out<<endl;
    }
    file_out<<endl;
    file_out<<"delta T covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<deltatCount[ix][iy]<<",";
      }
      file_out<<endl;
    }
    file_out<<endl;

    file_out<<"delta Tcovariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<deltatcov[ix][iy]/TMath::Max(1,deltatCount[ix][iy])<<",";
	h_deltatcov->Fill(ix, iy, deltatcov[ix][iy]/TMath::Max(1,deltatCount[ix][iy]));
      }
      file_out<<endl;
    }
    file_out<<endl;

    file_out<<"delta T2 covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<deltatcov2[ix][iy]<<", \t";
      }
      file_out<<endl;
    }
    file_out<<endl;
    file_out<<"delta T2 covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<deltatCount2[ix][iy]<<",";
      }
      file_out<<endl;
    }
    file_out<<endl;

    file_out<<"delta T2covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<deltatcov2[ix][iy]/TMath::Max(1,deltatCount2[ix][iy])<<",";
	h_deltatcov2->Fill(ix, iy, deltatcov2[ix][iy]/TMath::Max(1,deltatCount2[ix][iy]));
      }
      file_out<<endl;
    }
    file_out<<endl;

    //Y-side----------------------------------------------
    file_out<<"delta T covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<deltatcovy[ix][iy]<<", \t";
      }
      file_out<<endl;
    }
    file_out<<endl;
    file_out<<"delta T covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<deltatCounty[ix][iy]<<",";
      }
      file_out<<endl;
    }
    file_out<<endl;

    file_out<<"delta Tcovariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<deltatcovy[ix][iy]/TMath::Max(1,deltatCounty[ix][iy])<<",";
	h_deltatcovy->Fill(ix, iy, deltatcovy[ix][iy]/TMath::Max(1,deltatCounty[ix][iy]));
      }
      file_out<<endl;
    }
    file_out<<endl;
  }

  //Covariaance of X-pos
  file_out<<"delta Xpos covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposxcov[ix][iy]<<", \t";
    }
    file_out<<endl;
  }
  file_out<<endl;

  file_out<<"delta Xpos covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposxCount[ix][iy]<<",";
    }
    file_out<<endl;
  }
  file_out<<endl;

  file_out<<"delta Xposcovariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposxcov[ix][iy]/TMath::Max(1,deltaposxCount[ix][iy])<<",";
      h_deltaposxcov->Fill(ix, iy, deltaposxcov[ix][iy]/TMath::Max(1,deltaposxCount[ix][iy]));
    }
    file_out<<endl;
  }
  file_out<<endl;

  //Covariaance of Y-pos
  file_out<<"delta Ypos covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposycov[ix][iy]<<", \t";
    }
    file_out<<endl;
  }
  file_out<<endl;
  file_out<<"delta Ypos covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposyCount[ix][iy]<<",";
    }
    file_out<<endl;
  }
  file_out<<endl;

  file_out<<"delta Yposcovariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposycov[ix][iy]/TMath::Max(1,deltaposyCount[ix][iy])<<",";
      h_deltaposycov->Fill(ix, iy, deltaposycov[ix][iy]/TMath::Max(1,deltaposyCount[ix][iy]));
    }
    file_out<<endl;
  }
  file_out<<endl;

  for (int ij=0; ij<8; ij++) {
    for (int jk=0; jk<nlayer; jk++) {
      for (int kl=0; kl<nstrip; kl++) {
	if (!posinuse[ij][jk][kl]) {file_out<<"posinuse["<<ij<<"]["<<jk<<"]["<<kl<<"] ="<<posinuse[ij][jk][kl]<<";"<<endl;}
      }
    }
  }

  for (int ij=0; ij<8; ij++) {
    for (int jk=0; jk<nlayer; jk++) {
      for (int kl=0; kl<nstrip; kl++) {
	if (!timeinuse[ij][jk][kl]) {file_out<<"timeinuse["<<ij<<"]["<<jk<<"]["<<kl<<"] ="<<timeinuse[ij][jk][kl]<<";"<<endl;}
      }
    }
  }
  if (isTiming) {
    for (int ij=0; ij<nlayer-1; ij++) {
      for (int jk=ij+1; jk<nlayer; jk++) {
	for (int kl=0; kl<=nmxiter; kl++) {
	  double ent = timex_correl[ij][jk][kl]->GetEntries();
	  if (ent>0) {
	    double mean = timex_correl[ij][jk][kl]->GetMean();
	    double rms = timex_correl[ij][jk][kl]->GetRMS();
	    time_both_cormean[kl]->Fill(ij, jk, mean);
	    time_both_corrms[kl]->Fill(ij, jk, rms);
	  }
	  ent = timey_correl[ij][jk][kl]->GetEntries();
	  if (ent>0) {
	    double mean = timey_correl[ij][jk][kl]->GetMean();
	    double rms = timey_correl[ij][jk][kl]->GetRMS();
	    time_both_cormean[kl]->Fill(jk+1, ij, mean);
	    time_both_corrms[kl]->Fill(jk+1, ij, rms);
	  }
	}
      }
    }

    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<=nmxiter; jk++) {
	double ent = timex_shift[ij][jk]->GetEntries();
	if (ent>0) {
	  double mean = timex_shift[ij][jk]->GetMean();
	  double rms = timex_shift[ij][jk]->GetRMS();
	  time_both_cormean[jk]->Fill(ij, ij, mean);
	  time_both_corrms[jk]->Fill(ij, ij, rms);
	}

	ent = timey_shift[ij][jk]->GetEntries();
	if (ent>0) {
	  double mean = timey_shift[ij][jk]->GetMean();
	  double rms = timey_shift[ij][jk]->GetRMS();
	  time_both_cormean[jk]->Fill(ij+1, ij, mean);
	  time_both_corrms[jk]->Fill(ij+1, ij, rms);
	}
      }
    }

    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<npixel+1; jk++) {
	for (int kl=0; kl<=nmxiter; kl++) {
	  double ent = timexy_correl[ij][jk][kl]->GetEntries();
	  if (ent >10) {
	    double mean = timexy_correl[ij][jk][kl]->GetMean();
	    double rms = timexy_correl[ij][jk][kl]->GetRMS();
	    timexy_cormean[jk]->Fill(ij, kl, mean);
	    timexy_corrms[jk]->Fill(ij, kl, rms);
	  }
	}
      }
    }
  } //if (isTiming)

   for(int ij=0;ij<nlayer;ij++) {
    double xx = max(1.,h_corr_tdc_ref->GetBinContent(ij+1,0));
    for(int jk=0;jk<nlayer;jk++) {
      double  yy = 100.*h_corr_tdc_ref->GetBinContent(ij+1,jk+1)/xx;
      double yer = 100.*h_corr_tdc_ref->GetBinError(ij+1,jk+1)/xx;
      h_corr_tdc_ref->SetBinContent(ij+1,jk+1,yy);
      h_corr_tdc_ref->SetBinError(ij+1,jk+1,yer);

    }
  }

  h_corr_layer_mult->Scale(1./max(1.,h_corr_layer_mult->GetBinContent(0,0)));
  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=1; jk<=nstrip; jk++) {
      double xx = max(1., h_raw_xcorstrips[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_raw_xcorstrips[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_raw_xcorstrips[ij]->GetBinError(jk, kl)/xx;
	h_raw_xcorstrips[ij]->SetBinContent(jk, kl, yy);
	h_raw_xcorstrips[ij]->SetBinError(jk, kl, yer);
      }
      xx = max(1., h_raw_ycorstrips[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_raw_ycorstrips[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_raw_ycorstrips[ij]->GetBinError(jk, kl)/xx;
	h_raw_ycorstrips[ij]->SetBinContent(jk, kl, yy);
	h_raw_ycorstrips[ij]->SetBinError(jk, kl, yer);
      }

      xx = max(1., h_raw_xstrpnhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_raw_xstrpnhits[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_raw_xstrpnhits[ij]->GetBinError(jk, kl)/xx;
	h_raw_xstrpnhits[ij]->SetBinContent(jk, kl, yy);
	h_raw_xstrpnhits[ij]->SetBinError(jk, kl, yer);
      }
      xx = max(1., h_raw_ystrpnhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_raw_ystrpnhits[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_raw_ystrpnhits[ij]->GetBinError(jk, kl)/xx;
	h_raw_ystrpnhits[ij]->SetBinContent(jk, kl, yy);
	h_raw_ystrpnhits[ij]->SetBinError(jk, kl, yer);
      }

      xx = max(1., h_raw_xystrpnhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_raw_xystrpnhits[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_raw_xystrpnhits[ij]->GetBinError(jk, kl)/xx;
	h_raw_xystrpnhits[ij]->SetBinContent(jk, kl, yy);
	h_raw_xystrpnhits[ij]->SetBinError(jk, kl, yer);
      }
      xx = max(1., h_raw_yxstrpnhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_raw_yxstrpnhits[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_raw_yxstrpnhits[ij]->GetBinError(jk, kl)/xx;
	h_raw_yxstrpnhits[ij]->SetBinContent(jk, kl, yy);
	h_raw_yxstrpnhits[ij]->SetBinError(jk, kl, yer);
      }

    }
  }


  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=1; jk<=nstrip; jk++) {
      double xx = max(1., h_xmucorstrips[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_xmucorstrips[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_xmucorstrips[ij]->GetBinError(jk, kl)/xx;
	h_xmucorstrips[ij]->SetBinContent(jk, kl, yy);
	h_xmucorstrips[ij]->SetBinError(jk, kl, yer);
      }
      xx = max(1., h_ymucorstrips[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_ymucorstrips[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_ymucorstrips[ij]->GetBinError(jk, kl)/xx;
	h_ymucorstrips[ij]->SetBinContent(jk, kl, yy);
	h_ymucorstrips[ij]->SetBinError(jk, kl, yer);
      }

       xx = max(1., h_xymucorstrips[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_xymucorstrips[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_xymucorstrips[ij]->GetBinError(jk, kl)/xx;
	h_xymucorstrips[ij]->SetBinContent(jk, kl, yy);
	h_xymucorstrips[ij]->SetBinError(jk, kl, yer);
      }
      xx = max(1., h_yxmucorstrips[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_yxmucorstrips[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_yxmucorstrips[ij]->GetBinError(jk, kl)/xx;
	h_yxmucorstrips[ij]->SetBinContent(jk, kl, yy);
	h_yxmucorstrips[ij]->SetBinError(jk, kl, yer);
      }


      xx = max(1., h_xmucornhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_xmucornhits[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_xmucornhits[ij]->GetBinError(jk, kl)/xx;
	h_xmucornhits[ij]->SetBinContent(jk, kl, yy);
	h_xmucornhits[ij]->SetBinError(jk, kl, yer);
      }
      xx = max(1., h_ymucornhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy =  h_ymucornhits[ij]->GetBinContent(jk, kl)/xx;
	double yer =  h_ymucornhits[ij]->GetBinError(jk, kl)/xx;
	 h_ymucornhits[ij]->SetBinContent(jk, kl, yy);
	 h_ymucornhits[ij]->SetBinError(jk, kl, yer);
      }

      xx = max(1., h_xymucornhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_xymucornhits[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_xymucornhits[ij]->GetBinError(jk, kl)/xx;
	h_xymucornhits[ij]->SetBinContent(jk, kl, yy);
	h_xymucornhits[ij]->SetBinError(jk, kl, yer);
      }
      xx = max(1., h_yxmucornhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
	double yy = h_yxmucornhits[ij]->GetBinContent(jk, kl)/xx;
	double yer = h_yxmucornhits[ij]->GetBinError(jk, kl)/xx;
	h_yxmucornhits[ij]->SetBinContent(jk, kl, yy);
	h_yxmucornhits[ij]->SetBinError(jk, kl, yer);
      }

    }
  }

  /*
  for(int ij=0;ij<nlayer;ij++) {
    double xx = 0.;
      for(int kl=0;kl<h_xtstrpmult_xtdchitmult[ij]->GetNbinsX();kl++) {
	double yy=0.;
	for(int lm=0;kl<h_xtstrpmult_xtdchitmult[ij]->GetNbinsY();lm++) {
	  yy += h_xtstrpmult_xtdchitmult[ij]->GetBinContent(kl+1,lm+1);
	}
	for(int lm=0;kl<h_xtstrpmult_xtdchitmult[ij]->GetNbinsY();lm++) {
	 xx =  h_xtstrpmult_xtdchitmult[ij]->GetBinContent(kl+1,lm+1)/yy;
	 h_xtstrpmult_xtdchitmult[ij]->SetBinContent(kl+1,xx);
	}
      }
     for(int kl=0;kl<h_ytstrpmult_ytdchitmult[ij]->GetNbinsX();kl++) {
	double yy=0.;
	for(int lm=0;kl<h_ytstrpmult_ytdchitmult[ij]->GetNbinsY();lm++) {
	  yy += h_ytstrpmult_ytdchitmult[ij]->GetBinContent(kl+1,lm+1);
	}
	for(int lm=0;kl<h_ytstrpmult_ytdchitmult[ij]->GetNbinsY();lm++) {
	 xx =  h_ytstrpmult_ytdchitmult[ij]->GetBinContent(kl+1,lm+1)/yy;
	 h_ytstrpmult_ytdchitmult[ij]->SetBinContent(kl+1,xx);
	}
      }

     for(int kl=0;kl<h_xtstrpmult_ytdchitmult[ij]->GetNbinsX();kl++) {
	double yy=0.;
	for(int lm=0;kl<h_xtstrpmult_ytdchitmult[ij]->GetNbinsY();lm++) {
	  yy += h_xtstrpmult_ytdchitmult[ij]->GetBinContent(kl+1,lm+1);
	}
	for(int lm=0;kl<h_xtstrpmult_ytdchitmult[ij]->GetNbinsY();lm++) {
	 xx =  h_xtstrpmult_ytdchitmult[ij]->GetBinContent(kl+1,lm+1)/yy;
	 h_xtstrpmult_ytdchitmult[ij]->SetBinContent(kl+1,xx);
	}
      }

     for(int kl=0;kl<h_ytstrpmult_xtdchitmult[ij]->GetNbinsX();kl++) {
	double yy=0.;
	for(int lm=0;kl<h_ytstrpmult_xtdchitmult[ij]->GetNbinsY();lm++) {
	  yy += h_ytstrpmult_xtdchitmult[ij]->GetBinContent(kl+1,lm+1);
	}
	for(int lm=0;kl<h_ytstrpmult_xtdchitmult[ij]->GetNbinsY();lm++) {
	 xx =  h_ytstrpmult_xtdchitmult[ij]->GetBinContent(kl+1,lm+1)/yy;
	 h_ytstrpmult_xtdchitmult[ij]->SetBinContent(kl+1,xx);
	}
      }
  }

   */
  //h_xtstrpmult_xtdchitmult[ij]





  if (!isOnlyCom) {
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=1; jk<=nstrip; jk++) {
	double xx = max(1., h_xcorstrips[ij]->GetBinContent(jk,0));
	for (int kl=1; kl<=nstrip; kl++) {
	  double yy = h_xcorstrips[ij]->GetBinContent(jk, kl)/xx;
	  double yer = h_xcorstrips[ij]->GetBinError(jk, kl)/xx;
	  h_xcorstrips[ij]->SetBinContent(jk, kl, yy);
	  h_xcorstrips[ij]->SetBinError(jk, kl, yer);
	}

	xx = max(1., h_ycorstrips[ij]->GetBinContent(jk,0));
	for (int kl=1; kl<=nstrip; kl++) {
	  double yy = h_ycorstrips[ij]->GetBinContent(jk, kl)/xx;
	  double yer = h_ycorstrips[ij]->GetBinError(jk, kl)/xx;
	  h_ycorstrips[ij]->SetBinContent(jk, kl, yy);
	  h_ycorstrips[ij]->SetBinError(jk, kl, yer);
	}

	xx = max(1., h_xycorstrips[ij]->GetBinContent(jk,0));
	for (int kl=1; kl<=nstrip; kl++) {
	  double yy = h_xycorstrips[ij]->GetBinContent(jk, kl)/xx;
	  double yer = h_xycorstrips[ij]->GetBinError(jk, kl)/xx;
	  h_xycorstrips[ij]->SetBinContent(jk, kl, yy);
	  h_xycorstrips[ij]->SetBinError(jk, kl, yer);
	}

	xx = max(1., h_yxcorstrips[ij]->GetBinContent(jk,0));
	for (int kl=1; kl<=nstrip; kl++) {
	  double yy = h_yxcorstrips[ij]->GetBinContent(jk, kl)/xx;
	  double yer = h_yxcorstrips[ij]->GetBinError(jk, kl)/xx;
	  h_yxcorstrips[ij]->SetBinContent(jk, kl, yy);
	  h_yxcorstrips[ij]->SetBinError(jk, kl, yer);
	}

	if (isTiming) {
	  xx = max(1., h_xtcorstrips[ij]->GetBinContent(jk,0));
	  for (int kl=1; kl<=nstrip; kl++) {
	    double yy = h_xtcorstrips[ij]->GetBinContent(jk, kl)/xx;
	    double yer = h_xtcorstrips[ij]->GetBinError(jk, kl)/xx;
	    h_xtcorstrips[ij]->SetBinContent(jk, kl, yy);
	    h_xtcorstrips[ij]->SetBinError(jk, kl, yer);
	  }

	  xx = max(1., h_ytcorstrips[ij]->GetBinContent(jk,0));
	  for (int kl=1; kl<=nstrip; kl++) {
	    double yy = h_ytcorstrips[ij]->GetBinContent(jk, kl)/xx;
	    double yer = h_ytcorstrips[ij]->GetBinError(jk, kl)/xx;
	    h_ytcorstrips[ij]->SetBinContent(jk, kl, yy);
	    h_ytcorstrips[ij]->SetBinError(jk, kl, yer);
	  }

	  xx = max(1., h_xytcorstrips[ij]->GetBinContent(jk,0));
	  for (int kl=1; kl<=nstrip; kl++) {
	    double yy = h_xytcorstrips[ij]->GetBinContent(jk, kl)/xx;
	    double yer = h_xytcorstrips[ij]->GetBinError(jk, kl)/xx;
	    h_xytcorstrips[ij]->SetBinContent(jk, kl, yy);
	    h_xytcorstrips[ij]->SetBinError(jk, kl, yer);
	  }

	  xx = max(1., h_yxtcorstrips[ij]->GetBinContent(jk,0));
	  for (int kl=1; kl<=nstrip; kl++) {
	    double yy = h_yxtcorstrips[ij]->GetBinContent(jk, kl)/xx;
	    double yer = h_yxtcorstrips[ij]->GetBinError(jk, kl)/xx;
	    h_yxtcorstrips[ij]->SetBinContent(jk, kl, yy);
	    h_yxtcorstrips[ij]->SetBinError(jk, kl, yer);
	  }





	} // if (isTiming)
      } //for (int jk=1; jk<=nstrip; jk++)
    } //for (int ij=0; ij<nlayer; ij++)


    if (isTiming) {
      cout<<"Problematic strips with timing, large offsets"<<endl;
      file_out<<"Problematic strips with timing, large offsets"<<endl;

      double ntentry=0;
      double nover=0;
      double nunder = 0;
      for (int ixy=0; ixy<2; ixy++) {
	int nbinxx = (ixy==0) ? time_xstrreso[0][0][0]->GetNbinsX()+1 : time_ystrreso[0][0][0]->GetNbinsX()+1;
	for (int ntt = 0; ntt<2*nmxiter; ntt++) {
	  for (int ij=0; ij<nlayer; ij++) {
	    for (int jk=0; jk<nstrip; jk++) {
	      if (ixy==0) {
		ntentry=time_xstrreso[ij][jk][ntt]->GetEntries();
		nunder =time_xstrreso[ij][jk][ntt]->GetBinContent(0);
		nover  =time_xstrreso[ij][jk][ntt]->GetBinContent(nbinxx);
	      } else {
		ntentry=time_ystrreso[ij][jk][ntt]->GetEntries();
		nunder =time_ystrreso[ij][jk][ntt]->GetBinContent(0);
		nover  =time_ystrreso[ij][jk][ntt]->GetBinContent(nbinxx);
	      }
	      if (nunder >0.25*ntentry || nover >0.25*ntentry) {
		cout<<"ixy= "<<ixy<<" il= "<< ij<<" iter= "<<ntt<<" Strip= "<<jk<<" Ents: "<<ntentry<<" "<< nunder<<" "<<nover<<" Off= "<<xtoffset[ij][jk]<<" "<<ytoffset[ij][jk]<<endl;
		file_out<<"ixy= "<<ixy<<" il= "<< ij<<" iter= "<<ntt<<" Strip= "<<jk<<" Ents: "<<ntentry<<" "<< nunder<<" "<<nover<<" Off= "<<xtoffset[ij][jk]<<" "<< ytoffset[ij][jk]<<endl;
	      }

	    } //for (int jk=1; jk<=nstrip; jk++)
	  } //for (int ij=0; ij<nlayer; ij++)
	} // (ixy)
      } //(int ntt = 0; ntt<nmxiter; ntt++)
    } // if (isTiming)

  } //   if (!isOnlyCom)

  ps.NewPage();
  gStyle->SetOptTitle(1);
  gStyle->SetOptLogy(1);
  gStyle->SetOptFit(101);
  gStyle->SetStatW(.36); //40);
  gStyle->SetStatH(.20); //30);
  gStyle->SetPadLeftMargin(0.07);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.09); //0.03);
  gStyle->SetPadRightMargin(0.01);
  TCanvas* c4c = new TCanvas("c4c", "c4c", 700, 900);
  c4c->Divide(3,4);
  gStyle->SetOptStat(1100);

  gStyle->SetStatY(0.99); gStyle->SetStatTextColor(1);

  if (nmxiter>0) {
    // if (isTiming && (!isOnlyCom)) {
//       ps.NewPage();
//       for (int ij=0; ij<nlayer; ij++) {
// 	c4c->cd(ij+1);
// 	gPad->SetLogy(0);
// 	widthx_timex[ij][0]->Draw("colz");
//       }
//       c4c->Update();

//       ps.NewPage();
//       for (int ij=0; ij<nlayer; ij++) {
// 	c4c->cd(ij+1);
// 	gPad->SetLogy(0);
// 	widthy_timey[ij][0]->Draw("colz");
//       }
//       c4c->Update();
//     }

    ps.NewPage();
    int iitr = (isOnlyCom) ? 0 : 2*nmxiter-1; //which iteration should we use

    for (int ij=0; ij<nlayer; ij++) {
      c4c->cd(ij+1);
      xlayer_exterr[ij][iitr]->GetXaxis()->SetLabelSize(0.055);
      xlayer_exterr[ij][iitr]->SetLineColor(1); xlayer_exterr[ij][iitr]->Draw();
    }
    c4c->Update();
    gStyle->SetStatY(0.79); gStyle->SetStatTextColor(2);
    for (int ij=0; ij<nlayer; ij++) {
      c4c->cd(ij+1);
      ylayer_exterr[ij][iitr]->SetLineColor(2); ylayer_exterr[ij][iitr]->Draw("sames");
    }
    c4c->Update();

    if (isTiming) {
      ps.NewPage();
      gStyle->SetStatY(0.99); gStyle->SetStatTextColor(1);
      for (int ij=0; ij<nlayer; ij++) {
	c4c->cd(ij+1);
	xtime_exterr[ij][iitr]->GetXaxis()->SetLabelSize(0.055);
	xtime_exterr[ij][iitr]->SetLineColor(1); xtime_exterr[ij][iitr]->Draw();
      }
      c4c->Update();

      gStyle->SetStatY(0.79); gStyle->SetStatTextColor(2);
      for (int ij=0; ij<nlayer; ij++) {
	c4c->cd(ij+1);
	ytime_exterr[ij][iitr]->SetLineColor(2); ytime_exterr[ij][iitr]->Draw("sames");
      }
      c4c->Update();
    }

 }




  int iitr = (isOnlyCom) ? 0 : 2*nmxiter-1; //which iteration should we use
  ps.NewPage();
  gStyle->SetOptStat(1110);
  gStyle->SetStatW(.36);
  gStyle->SetStatH(.18);
  gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
  for (int ij=0; ij<nlayer; ij++) {
    c4c->cd(ij+1);
    gPad->SetLogy(1);
    strp_xmul[ij][iitr]->GetXaxis()->SetLabelSize(0.07);
    strp_xmul[ij][iitr]->GetYaxis()->SetLabelSize(0.07);
    strp_xmul[ij][iitr]->SetLineColor(1);
    strp_xmul[ij][iitr]->ProjectionY()->Draw();
  }
  c4c->Update();
  gStyle->SetStatY(.77); gStyle->SetStatTextColor(2);
  TH1F* tmphist[nlayer]={0};
  for (int ij=0; ij<nlayer; ij++) {
    c4c->cd(ij+1);
    gPad->SetLogy(1);
    tmphist[ij] = (TH1F*)strp_ymul[ij][iitr]->ProjectionY();
    tmphist[ij]->SetLineColor(2);
    tmphist[ij]->Draw("sames");
  }

  c4c->Update();
  for (int ij=0; ij<nlayer; ij++) {
    if (tmphist[ij]) { delete tmphist[ij]; tmphist[ij]=0;}
  }

  ps.NewPage();
  gStyle->SetStatTextColor(1); gStyle->SetStatY(0.99);

  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatW(.36); //40);
  gStyle->SetStatH(.28); //30);

  TCanvas* c4a = new TCanvas("c4a", "c4a", 700, 900);
  c4a->Divide(3,4);

  //  int iitr = 2*nmxiter-1; //which iteration should we use
  TH2F* histy[nlayer]={0};
  const int nmult=6;
  TGraph* gry[nlayer][nmult]={0};
  double content[100][nmult];
  double xvl[100];
  double yvl[100];
  double xerr[100]={0.0};
  double yerr[100];

  double total[100];
  const char* namex[nmult]={"null", "ZERO", "ONE", "TWO", "THREE", "More"};
  TLegend *tleg = new TLegend(.70, .62, .90, .92, "","brNDC");
  tleg->SetFillColor(10);
  tleg->SetTextSize(0.06);
  tleg->SetBorderSize(0);
  tleg->SetTextFont(42);

  int ishft = (isOnlyCom) ? 0 : nmxiter; //Without any iteration & proper efficiency
  for (int ixx=0; ixx<nmxiter; ixx++) {
    for (int jkl=0; jkl <2; jkl++) {
      if (ixx!=0 || jkl!=0) {ps.NewPage();}
      for (int ij=0; ij<nlayer; ij++) {
	c4a->cd(ij+1);
	switch(jkl) {
	case 0 : histy[ij] = (TH2F*)strp_xmul[ij][ishft+ixx]->Clone(); break;
	case 1 : histy[ij] = (TH2F*)strp_ymul[ij][ishft+ixx]->Clone(); break;
	default : histy[ij] = (TH2F*)strp_xmul[ij][ishft+ixx]->Clone(); break;
	}
	int nbinx=2+histy[ij]->GetNbinsX();
	int nbiny=2+histy[ij]->GetNbinsY();

	for (int jk=0; jk<nbinx; jk++) {
	  total[jk]=0;
	  xvl[jk] = histy[ij]->GetXaxis()->GetBinCenter(jk);
	  for (int kl=0; kl<nbiny; kl++) {
	    content[jk][kl] = histy[ij]->GetBinContent(jk, kl);
	    total[jk] +=content[jk][kl];
	  }
	  if (total[jk]<1.0) total[jk]=1.0;
	}
	int icol=0;
	for (int kl=1; kl<nmult; kl++) {
	  for (int jk=0; jk<nbinx; jk++) {
	    yvl[jk] = content[jk][kl]/total[jk];
	    yerr[jk] =sqrt(yvl[jk]*(1-yvl[jk])/total[jk]);
	  }
	  icol++; if (icol%5==0) icol++;

	  gry[ij][kl]=new TGraphErrors(nbinx, xvl, yvl, xerr, yerr);
	  gry[ij][kl]->SetTitle(histy[ij]->GetTitle());
	  gry[ij][kl]->SetMarkerStyle(20);
	  gry[ij][kl]->SetMarkerSize(0.5);
	  gry[ij][kl]->SetMarkerColor(icol);
	  gry[ij][kl]->GetXaxis()->SetLabelSize(.064);
	  gry[ij][kl]->GetYaxis()->SetLabelSize(.054);
	  gry[ij][kl]->GetXaxis()->SetLabelOffset(.001);
	  gry[ij][kl]->GetYaxis()->SetLabelOffset(.01);
	  gry[ij][kl]->GetXaxis()->SetTitle("Position in strip");
	  gry[ij][kl]->GetXaxis()->CenterTitle();
	  gry[ij][kl]->GetXaxis()->SetTitleSize(.08);
	  gry[ij][kl]->GetXaxis()->SetTitleOffset(.8);
	  gry[ij][kl]->GetXaxis()->SetRangeUser(-0.5, 0.5);
	  gry[ij][kl]->SetMinimum(0);
	  gry[ij][kl]->SetMaximum(0.9);
	  //    gry[ij][kl]->SetNdivisions(506,"XY");
	  if (kl==1) {
	    gry[ij][kl]->Draw("AP");
	  } else {
	    gry[ij][kl]->Draw("P:same");
	  }
	  if (ixx==0 && ij==0 && jkl==0) {
	    tleg->AddEntry(gry[ij][kl],namex[kl],"p");
	  }
	  //	  if (ij==0 && jkl==1) {latex.DrawLatex(0.20, 0.98, histy[ij]->GetTitle());}
	}
	if (ij==0) tleg->Draw();
      }
      c4a->Update();
      for (int ij=0; ij<nlayer; ij++) {
	if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
	for (int jk=0; jk<6; jk++) {
	  if (gry[ij][jk]) { delete gry[ij][jk]; gry[ij][jk]=0;}
	}
      }
    } //for (int jkl=0; jkl <2; jkl++)
  } //for (int ixx=0; ixx<nmxiter; ixx++)

  //  int iitr = 2*nmxiter-1; //which iteration should we use


  ps.NewPage();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadLeftMargin(0.15);
  TH1F* histz[nlayer];
  if (nmxiter >0) {
    latex.SetTextSize(0.08);
    for (int jkl=0; jkl <4; jkl++) {
      if (isOnlyCom && jkl>=2) continue;
      for (int klm =0; klm<nmxhits; klm++) {
      //for (int klm =0; klm<nmxusedhits; klm++) {
	ps.NewPage();
	for (int ij=0; ij<nlayer; ij++) {
	  c4c->cd(ij+1);
	  switch(jkl) {
	  case 0 : histz[ij] = (TH1F*)xlayer_reso_mul[ij][2*nmxiter-2][klm]->Clone(); break;
	  case 1 : histz[ij] = (TH1F*)ylayer_reso_mul[ij][2*nmxiter-2][klm]->Clone(); break;
	  case 2 : histz[ij] = (TH1F*)xlayer_reso_mul[ij][2*nmxiter-1][klm]->Clone(); break;
	  case 3 : histz[ij] = (TH1F*)ylayer_reso_mul[ij][2*nmxiter-1][klm]->Clone(); break;
	  default : histz[ij] = (TH1F*)ylayer_reso_mul[ij][2*nmxiter-1][klm]->Clone(); break;
	  }
	  histz[ij]->GetXaxis()->SetLabelSize(.07);
	  histz[ij]->GetYaxis()->SetLabelSize(.055);
	  if (jkl<=1) {
	    histz[ij]->GetXaxis()->SetTitle("#Delta X (pitch)");
	  } else {
	    histz[ij]->GetXaxis()->SetTitle("#Delta Y (pitch)");
	  }
	  histz[ij]->GetXaxis()->CenterTitle();
	  histz[ij]->GetXaxis()->SetTitleOffset(0.8);
	  histz[ij]->GetXaxis()->SetTitleSize(0.07);

	  TFitResultPtr ptr = histz[ij]->Fit("gaus","SQ");
	  latex.DrawLatex(0.15, 0.70,Form("%g", int(1000*histz[ij]->GetMean())/1000.));
	  latex.DrawLatex(0.15, 0.56,Form("%g", int(1000*histz[ij]->GetRMS())/1000.));
	  int fitStatus = ptr;
	  if (fitStatus==0 && histz[ij]->GetEntries()>3) {
	    latex.DrawLatex(0.75, 0.70,Form("%g", int(1000*ptr->Parameter(1))/1000.));
	    latex.DrawLatex(0.75, 0.56,Form("%g", int(1000*ptr->Parameter(2))/1000.));
	  }
	}
	c4c->Update();
	for (int ij=0; ij<nlayer; ij++) {
	  if (histz[ij]) { delete histz[ij]; histz[ij]=0;}
	}
      }
    }

    if (isTiming) {
      for (int jkl=0; jkl <4; jkl++) {
	if (isOnlyCom && jkl>=2) continue;
	for (int klm =0; klm<nmxtimehit; klm++) {
	  ps.NewPage();
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    switch(jkl) {
	    case 0 : histz[ij] = (TH1F*)time_mulxreso[ij][2*nmxiter-2][klm]->Clone(); break;
	    case 1 : histz[ij] = (TH1F*)time_mulyreso[ij][2*nmxiter-2][klm]->Clone(); break;
	    case 2 : histz[ij] = (TH1F*)time_mulxreso[ij][2*nmxiter-1][klm]->Clone(); break;
	    case 3 : histz[ij] = (TH1F*)time_mulyreso[ij][2*nmxiter-1][klm]->Clone(); break;
	    default : histz[ij] = (TH1F*)time_mulxreso[ij][2*nmxiter-1][klm]->Clone(); break;
	    }
	    histz[ij]->GetXaxis()->SetLabelSize(.07);
	    histz[ij]->GetYaxis()->SetLabelSize(.055);
	    histz[ij]->GetXaxis()->SetTitle("#Delta t(ns)");
	    histz[ij]->GetXaxis()->CenterTitle();
	    histz[ij]->GetXaxis()->SetTitleOffset(0.9);
	    histz[ij]->GetXaxis()->SetTitleSize(0.06);

	    TFitResultPtr ptr = histz[ij]->Fit("gaus","SQ");
	    latex.DrawLatex(0.20, 0.76,Form("%g", int(1000*histz[ij]->GetMean())/1000.));
	    latex.DrawLatex(0.20, 0.66,Form("%g", int(1000*histz[ij]->GetRMS())/1000.));
	    int fitStatus = ptr;
	    if (fitStatus==0 && histz[ij]->GetEntries()>3) {
	      latex.DrawLatex(0.70, 0.76,Form("%g", int(1000*ptr->Parameter(1))/1000.));
	      latex.DrawLatex(0.70, 0.66,Form("%g", int(1000*ptr->Parameter(2))/1000.));
	    }
	  }
	  c4c->Update();
	  for (int ij=0; ij<nlayer; ij++) {
	    if (histz[ij]) { delete histz[ij]; histz[ij]=0;}
	  }
	}
      }
    }
    latex.SetTextSize(0.12);
  }

  ps.NewPage();
  gStyle->SetOptStat(1100);
  gStyle->SetStatH(.30);
  gStyle->SetPadLeftMargin(0.15);
  ps.NewPage();
  // gStyle->SetPadLeftMargin(0.15);
  c4c->cd(1);gPad->SetLeftMargin(0.11); h_xprob->Draw();
  c4c->cd(2); gPad->SetLeftMargin(0.11);h_reduchisqx->Draw();
  c4c->cd(3); gPad->SetLeftMargin(0.11);h_xndf->Draw();

  c4c->cd(4); gPad->SetLeftMargin(0.11);h_yprob->Draw();
  c4c->cd(5); gPad->SetLeftMargin(0.11);h_reduchisqy->Draw();
  c4c->cd(6); gPad->SetLeftMargin(0.11);h_yndf->Draw();

  if (isTiming) {
    c4c->cd(7); gPad->SetLeftMargin(0.11);h_xtprob->Draw();
    c4c->cd(8); gPad->SetLeftMargin(0.11); h_treduchisqx->Draw();
    c4c->cd(9); gPad->SetLeftMargin(0.11); h_txndf->Draw();

    c4c->cd(10); gPad->SetLeftMargin(0.11); h_ytprob->Draw();
    c4c->cd(11);gPad->SetLeftMargin(0.11); h_treduchisqy->Draw();
    c4c->cd(12); gPad->SetLeftMargin(0.11);  h_tyndf->Draw();
  }
  c4c->Update();

  ps.NewPage();

  TH1F* hist1x[2*nprob];
  TH1F* hist1y[2*nprob];
  for (int jk=0; jk<nprob; jk++) {
    hist1x[jk] = (TH1F*)h_xnprob[jk]->Clone();
    hist1y[jk] = (TH1F*)h_ynprob[jk]->Clone();

    hist1x[jk+nprob] = (TH1F*)h_xtnprob[jk]->Clone();
    hist1y[jk+nprob] = (TH1F*)h_ytnprob[jk]->Clone();
  }
  gStyle->SetStatY(0.99); gStyle->SetStatTextColor(1);
  for (int jk=0; jk<2*nprob; jk++) {
    c4c->cd(jk+1);
    gPad->SetLeftMargin(0.11);
    hist1x[jk]->SetLineColor(1); hist1x[jk]->Draw();
  }

  c4c->Update();
  gStyle->SetStatY(0.84); gStyle->SetStatTextColor(2);
  for ( int jk=0; jk<2*nprob; jk++) {
    c4c->cd(jk+1);
    // gPad->SetLeftMargin(0.11);
    hist1y[jk]->SetLineColor(2); hist1y[jk]->Draw("sames");
  }

  c4c->Update();

  for (int jk=0; jk<2*nprob; jk++) {
    if (hist1x[jk]) { delete hist1x[jk]; hist1x[jk]=0;}
    if (hist1y[jk]) { delete hist1y[jk]; hist1y[jk]=0;}
  }

  if (!isOnlyCom && plot_level>80) {
    ps.NewPage();
    gStyle->SetStatTextColor(1); gStyle->SetStatY(0.99);
    gStyle->SetOptStat(0);
    for (int iitr=0; iitr<nmxiter; iitr++) {
      for (int jkl=0; jkl <4; jkl++) {
	for (int klm =0; klm<nstr_posmx; klm++) {
	  ps.NewPage();
	  for (int ij=0; ij<nlayer; ij++) {
	    c4a->cd(ij+1);
	    switch(jkl) {
	    case 0 : histy[ij] = (TH2F*)xstr_xdev[ij][iitr][klm]->Clone(); break;
	    case 1 : histy[ij] = (TH2F*)ystr_ydev[ij][iitr][klm]->Clone(); break;
	    case 2 : histy[ij] = (TH2F*)ystr_xdev[ij][iitr][klm]->Clone(); break;
	    case 3 : histy[ij] = (TH2F*)xstr_ydev[ij][iitr][klm]->Clone(); break;
	    default : histy[ij] = (TH2F*)xstr_xdev[ij][iitr][klm]->Clone(); break;
	    }
	    histy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
	    histy[ij]->GetXaxis()->SetTitle("Strip No");
	    histy[ij]->GetXaxis()->CenterTitle();
	    histy[ij]->GetXaxis()->SetTitleOffset(0.8);
	    histy[ij]->GetXaxis()->SetTitleSize(0.07);

	    histy[ij]->GetYaxis()->SetTitle("#Delta X (Strip pitch)");
	    histy[ij]->GetYaxis()->CenterTitle();
	    histy[ij]->GetYaxis()->SetTitleOffset(0.8);
	    histy[ij]->GetYaxis()->SetTitleSize(0.07);

	    histy[ij]->Draw("colz");
	    histy[ij]->ProfileX()->Draw("sames");
	  }
	  c4a->Update();
	  for (int ij=0; ij<nlayer; ij++) {
	    if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
	  }
	}
      }
    }

    if (isTiming) {
      TProfile* str_tdev[nlayer]={0};
      for (int iitr=0; iitr<nmxiter; iitr++) {
	if (iitr>0 && iitr<nmxiter-2) continue;
	for (int jkl=0; jkl <4; jkl++) {

	  ps.NewPage();
	  for (int ij=0; ij<nlayer; ij++) {
	    c4a->cd(ij+1);
	    switch(jkl) {
	    case 0 : histy[ij] = (TH2F*)xstr_xtdev[ij][iitr]->Clone(); break;
	    case 1 : histy[ij] = (TH2F*)ystr_ytdev[ij][iitr]->Clone(); break;
	    case 2 : histy[ij] = (TH2F*)ystr_xtdev[ij][iitr]->Clone(); break;
	    case 3 : histy[ij] = (TH2F*)xstr_ytdev[ij][iitr]->Clone(); break;
	    default : histy[ij] = (TH2F*)xstr_xtdev[ij][iitr]->Clone(); break;
	    }
	    histy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
	    histy[ij]->Draw("colz");
	    sprintf(name, "xxx_%i", ij);
	    str_tdev[ij] = (TProfile*)histy[ij]->ProfileX(name);
	    str_tdev[ij]->SetMarkerSize(0.36);
	    str_tdev[ij]->SetMarkerStyle(24);
	    str_tdev[ij]->Fit("pol4","Q","same");
	  }
	  c4a->Update();
	  for (int ij=0; ij<nlayer; ij++) {
	    if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
	    if (str_tdev[ij]) { delete str_tdev[ij]; str_tdev[ij]=0;}
	  }
	}
      }
    }
  }

  if (!isOnlyCom && plot_level>70) {
    for (int lm=0; lm<nmxtimehit; lm++) { //multiplicity
      for (int jkl=0; jkl <12; jkl++) {
	if (jkl<=3 && lm >=nstr_posmx) continue;
	if (!isTiming && jkl>=4) continue;
	ps.NewPage();
	gStyle->SetOptStat(0);
	for (int ij=0; ij<nlayer; ij++) {
	  c4a->cd(ij+1);
	  switch(jkl) {

	  case 0 : histy[ij] = (TH2F*)xpos_xdev[ij][lm]->Clone(); break;
	  case 1 : histy[ij] = (TH2F*)xpos_ydev[ij][lm]->Clone(); break;
	  case 2 : histy[ij] = (TH2F*)ypos_xdev[ij][lm]->Clone(); break;
	  case 3 : histy[ij] = (TH2F*)ypos_ydev[ij][lm]->Clone(); break;

	  case 4 : histy[ij] = (TH2F*)xpos_xtdev_str[ij][lm]->Clone(); break;
	  case 5 : histy[ij] = (TH2F*)ypos_ytdev_str[ij][lm]->Clone(); break;

	  case 6 : histy[ij] = (TH2F*)xpos_xtdev[ij][lm]->Clone(); break;
	  case 7 : histy[ij] = (TH2F*)ypos_ytdev[ij][lm]->Clone(); break;

	  case 8 : histy[ij] = (TH2F*)xpos_xtdev_glb[ij][lm]->Clone(); break;
	  case 9 : histy[ij] = (TH2F*)xpos_ytdev_glb[ij][lm]->Clone(); break;
	  case 10 : histy[ij] = (TH2F*)ypos_xtdev_glb[ij][lm]->Clone(); break;
	  case 11 : histy[ij] = (TH2F*)ypos_ytdev_glb[ij][lm]->Clone(); break;


	  default : histy[ij] = (TH2F*)xpos_xdev[ij][lm]->Clone(); break;
	  }
	  if (jkl >=4 && jkl<=7 ) {
	    histy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
	  }
	  histy[ij]->GetXaxis()->SetLabelSize(.07);
	  histy[ij]->GetYaxis()->SetLabelSize(.07);

	  histy[ij]->Draw("colz");
	  histy[ij]->ProfileX()->Draw("sames");
	}
	c4a->Update();
	for (int ij=0; ij<nlayer; ij++) {
	  if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
	}
      }
    }
  }

  if (isTiming) {
    ps.NewPage();
    for (int ij=0; ij<8; ij++) {
      c4a->cd(ij+1);
      nxypos_xytdev[ij]->GetXaxis()->SetRangeUser(-1., nusedstrip);
      nxypos_xytdev[ij]->GetYaxis()->SetRangeUser(-1., nusedstrip);
      nxypos_xytdev[ij]->Draw("colz");
    }
    for (int ij=0; ij<4; ij++) {
      c4a->cd(ij+9);
      xtdev_ytdev[ij]->Draw("colz");
    }
    c4a->Update();
  }

#ifdef ISEFFICIENCY
  int iitrx = (isOnlyCom) ? 0 : 2*nmxiter-1; //which iteration should we use
  if (plot_level>100) {
    ps.NewPage();
    gStyle->SetPaintTextFormat("5.2f");
    gStyle->SetOptStat(0);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadRightMargin(0.09);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadLeftMargin(0.08);
    gStyle->SetTitleFontSize(0.04);

    TCanvas *c45=new TCanvas ("c45","Strip occupancy",500,700);
    for (int ij=0; ij<nlayer; ij++) {
      ps.NewPage();
      c45->cd();
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetLabelSize(.035);
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetLabelSize(.035);
      inefficiency_cory[ij][iitrx]->GetZaxis()->SetLabelSize(.035);
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetLabelOffset(-.001);
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetLabelOffset(.01);

      inefficiency_cory[ij][iitrx]->GetXaxis()->SetTitle("X-Strip No");
      inefficiency_cory[ij][iitrx]->GetXaxis()->CenterTitle();
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetTitleOffset(0.7);
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetTitleSize(0.04);

      inefficiency_cory[ij][iitrx]->GetYaxis()->SetTitle("Y-Strip No");
      inefficiency_cory[ij][iitrx]->GetYaxis()->CenterTitle();
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetTitleOffset(0.5);
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetTitleSize(0.04);
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_cory[ij][iitrx]->SetMarkerSize(0.62);
      inefficiency_cory[ij][iitrx]->Draw("colz");
      c45->Update();
      ps.NewPage();
      inefficiency_cory[ij][iitrx]->Draw("text25");
      c45->Update();

      ps.NewPage();
      c45->cd();
      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetLabelSize(.035);
      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetLabelSize(.035);
      inefficiency_uncx[ij][iitrx]->GetZaxis()->SetLabelSize(.035);
      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetTitle("X-Strip No");
      inefficiency_uncx[ij][iitrx]->GetXaxis()->CenterTitle();
      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetTitleOffset(0.7);
      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetTitleSize(0.04);

      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetTitle("Y-Strip No");
      inefficiency_uncx[ij][iitrx]->GetYaxis()->CenterTitle();
      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetTitleOffset(0.5);
      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetTitleSize(0.04);


      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_uncx[ij][iitrx]->SetMarkerSize(0.62);
      inefficiency_uncx[ij][iitrx]->Draw("colz");
      c45->Update();
      ps.NewPage();
      inefficiency_uncx[ij][iitrx]->Draw("text25");
      c45->Update();

      ps.NewPage();
      c45->cd();
      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetLabelSize(.035);
      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetLabelSize(.035);
      inefficiency_uncy[ij][iitrx]->GetZaxis()->SetLabelSize(.035);

      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetTitle("X-Strip No");
      inefficiency_uncy[ij][iitrx]->GetXaxis()->CenterTitle();
      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetTitleOffset(0.7);
      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetTitleSize(0.04);

      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetTitle("Y-Strip No");
      inefficiency_uncy[ij][iitrx]->GetYaxis()->CenterTitle();
      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetTitleOffset(0.5);
      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetTitleSize(0.04);

      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_uncy[ij][iitrx]->SetMarkerSize(0.62);
      inefficiency_uncy[ij][iitrx]->Draw("colz");
      c45->Update();
      ps.NewPage();
      inefficiency_uncy[ij][iitrx]->Draw("text25");
      c45->Update();
    }
    c45->Clear();
  }

  if (plot_level>90) { //GMA
    ps.NewPage();
    gStyle->SetOptStat(0);
    gStyle->SetPadTopMargin(0.13);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.07);
    gStyle->SetTitleFontSize(0.07);
    for (int ix=0; ix<itset; ix++) {
      ps.NewPage();
      for (int ij=0; ij<nlayer-2; ij++) {

	c4a->cd(ij+1);
	gPad->SetRightMargin(0.15);
	gPad->SetTopMargin(0.13);
	gPad->SetLeftMargin(0.07);

	inefficiency_corx_set[ij][ix]->GetXaxis()->SetLabelSize(.07);
	inefficiency_corx_set[ij][ix]->GetYaxis()->SetLabelSize(.07);
	inefficiency_corx_set[ij][ix]->GetZaxis()->SetLabelSize(.06);
	inefficiency_corx_set[ij][ix]->GetXaxis()->SetLabelOffset(-.001);
	inefficiency_corx_set[ij][ix]->GetYaxis()->SetLabelOffset(.01);

	inefficiency_corx_set[ij][ix]->GetXaxis()->SetTitle("X-Strip No");
	inefficiency_corx_set[ij][ix]->GetXaxis()->CenterTitle();
	inefficiency_corx_set[ij][ix]->GetXaxis()->SetTitleOffset(0.8);
	inefficiency_corx_set[ij][ix]->GetXaxis()->SetTitleSize(0.07);

	inefficiency_corx_set[ij][ix]->GetYaxis()->SetTitle("Y-Strip No");
	inefficiency_corx_set[ij][ix]->GetYaxis()->CenterTitle();
	inefficiency_corx_set[ij][ix]->GetYaxis()->SetTitleOffset(0.8);
	inefficiency_corx_set[ij][ix]->GetYaxis()->SetTitleSize(0.07);

	inefficiency_corx_set[ij][ix]->SetNdivisions(506,"XY");
	inefficiency_corx_set[ij][ix]->Draw("colz");
	cout<<"ixxxx"<<"   "<<double(get_average(inefficiency_corx_set[ij][ix]))<<endl;
	latex.DrawLatex(0.4, 0.49, Form("[%g]", int(10000*get_average(inefficiency_corx_set[ij][ix]))/100.));


      } // for (int ij=0; ij<nlayer; ij++)
      c4a->Update();
    } //for (int ix=0; ix<nset; ix++)

    for (int ix=0; ix<itset; ix++) {
      ps.NewPage();
      for (int ij=0; ij<nlayer-2; ij++) {
	c4a->cd(ij+1);
	gPad->SetRightMargin(0.15);
	gPad->SetTopMargin(0.13);
	gPad->SetLeftMargin(0.07);


	triggereffi_x_set[ij][ix]->GetXaxis()->SetLabelSize(.07);
	triggereffi_x_set[ij][ix]->GetYaxis()->SetLabelSize(.07);
	triggereffi_x_set[ij][ix]->GetZaxis()->SetLabelSize(.06);
	triggereffi_x_set[ij][ix]->GetXaxis()->SetLabelOffset(-.001);
	triggereffi_x_set[ij][ix]->GetYaxis()->SetLabelOffset(.01);

	triggereffi_x_set[ij][ix]->GetXaxis()->SetTitle("X-Strip No");
	triggereffi_x_set[ij][ix]->GetXaxis()->CenterTitle();
	triggereffi_x_set[ij][ix]->GetXaxis()->SetTitleOffset(0.8);
	triggereffi_x_set[ij][ix]->GetXaxis()->SetTitleSize(0.07);

	triggereffi_x_set[ij][ix]->GetYaxis()->SetTitle("Y-Strip No");
	triggereffi_x_set[ij][ix]->GetYaxis()->CenterTitle();
	triggereffi_x_set[ij][ix]->GetYaxis()->SetTitleOffset(0.8);
	triggereffi_x_set[ij][ix]->GetYaxis()->SetTitleSize(0.07);

	triggereffi_x_set[ij][ix]->SetNdivisions(506,"XY");
	triggereffi_x_set[ij][ix]->Draw("colz");
	latex.DrawLatex(0.4, 0.49, Form("[%g]", int(10000*get_average(triggereffi_x_set[ij][ix]))/100.));
      } // for (int ij=0; ij<nlayer; ij++)
      c4a->Update();
    } //for (int ix=0; ix<nset; ix++)

    for (int ix=0; ix<itset; ix++) {
      ps.NewPage();
      for (int ij=0; ij<nlayer-2; ij++) {
	c4a->cd(ij+1);
	gPad->SetRightMargin(0.15);
	gPad->SetTopMargin(0.13);
	gPad->SetLeftMargin(0.07);


	triggereffi_y_set[ij][ix]->GetXaxis()->SetLabelSize(.07);
	triggereffi_y_set[ij][ix]->GetYaxis()->SetLabelSize(.07);
	triggereffi_y_set[ij][ix]->GetZaxis()->SetLabelSize(.06);
	triggereffi_y_set[ij][ix]->GetXaxis()->SetLabelOffset(-.001);
	triggereffi_y_set[ij][ix]->GetYaxis()->SetLabelOffset(.01);

	triggereffi_y_set[ij][ix]->GetXaxis()->SetTitle("X-Strip No");
	triggereffi_y_set[ij][ix]->GetXaxis()->CenterTitle();
	triggereffi_y_set[ij][ix]->GetXaxis()->SetTitleOffset(0.8);
	triggereffi_y_set[ij][ix]->GetXaxis()->SetTitleSize(0.07);

	triggereffi_y_set[ij][ix]->GetYaxis()->SetTitle("Y-Strip No");
	triggereffi_y_set[ij][ix]->GetYaxis()->CenterTitle();
	triggereffi_y_set[ij][ix]->GetYaxis()->SetTitleOffset(0.8);
	triggereffi_y_set[ij][ix]->GetYaxis()->SetTitleSize(0.07);

	triggereffi_y_set[ij][ix]->SetNdivisions(506,"XY");
	triggereffi_y_set[ij][ix]->Draw("colz");
	latex.DrawLatex(0.4, 0.49, Form("[%g]", int(10000*get_average(triggereffi_x_set[ij][ix]))/100.));
      } // for (int ij=0; ij<nlayer; ij++)
      c4a->Update();
    } //for (int ix=0; ix<nset; ix++)
  } //if (plot_level>90)


  //  for (int jkl=0; jkl <26; jkl++) { //GMA
  for (int jkl=0; jkl <2; jkl++) {
    ps.NewPage();
    gStyle->SetOptStat(0);
    gStyle->SetPadTopMargin(0.13);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetTitleFontSize(0.07);
    for (int ij=0; ij<nlayer; ij++) {
      switch(jkl) {
      case 0 : histz[ij] = (TH1F*)efficiency_xpixel[ij][iitrx]->Clone(); break;
      case 1 : histz[ij] = (TH1F*)efficiency_ypixel[ij][iitrx]->Clone(); break;
      default : histz[ij] = (TH1F*)efficiency_xpixel[ij][iitrx]->Clone(); break;
      }

      c4a->cd(ij+1);
      gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.02);
      histz[ij]->GetXaxis()->SetLabelSize(.07);
      histz[ij]->GetYaxis()->SetLabelSize(.06);
      histz[ij]->GetXaxis()->SetLabelOffset(-.001);
      histz[ij]->GetYaxis()->SetLabelOffset(.01);

      histz[ij]->GetXaxis()->SetTitle("Position in strip");
      histz[ij]->GetXaxis()->CenterTitle();
      histz[ij]->GetXaxis()->SetTitleOffset(0.8);
      histz[ij]->GetXaxis()->SetTitleSize(0.07);

      histz[ij]->GetYaxis()->SetTitle("Efficiency");
      histz[ij]->GetYaxis()->CenterTitle();
      histz[ij]->GetYaxis()->SetTitleOffset(1.0);
      histz[ij]->GetYaxis()->SetTitleSize(0.07);

      histz[ij]->Draw();
    }
    c4a->Update();
    for (int ij=0; ij<nlayer; ij++) {
      if (histz[ij]) { delete histz[ij]; histz[ij]=0;}
    }
  }

  //  for (int jkl=0; jkl <26; jkl++) { //GMA
  for (int jkl=0; jkl <38; jkl++) {
    if (!isTiming && jkl>17) continue;
    if (nmxiter==0 && (jkl<=6 || jkl>=34)) continue;
    if (nmxiter<2 && jkl>=34) continue;
    if (isOnlyCom && jkl>20) continue;
    ps.NewPage();
    gStyle->SetOptStat(0);
    gStyle->SetPadTopMargin(0.13);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.07);
    gStyle->SetTitleFontSize(0.07);
    double xmx=-1000000;
    double xmn=1000000;

    for (int ij=0; ij<nlayer; ij++) {
      switch(jkl) {
      case 0 : histy[ij] = (TH2F*)totalentry[ij][iitrx]->Clone(); break;
      case 1 : histy[ij] = (TH2F*)fine_totalentry[ij][iitrx]->Clone(); break;

      case 2 : histy[ij] = (TH2F*)inefficiency_corx[ij][iitrx]->Clone(); break;
      case 3 : histy[ij] = (TH2F*)inefficiency_cory[ij][iitrx]->Clone(); break;

      case 4 : histy[ij] = (TH2F*)inefficiencytrue_corx[ij][iitrx]->Clone(); break;
      case 5 : histy[ij] = (TH2F*)inefficiencytrue_cory[ij][iitrx]->Clone(); break;

      case 6 : histy[ij] = (TH2F*)inefficiency_uncx[ij][iitrx]->Clone(); break;
      case 7 : histy[ij] = (TH2F*)inefficiency_uncy[ij][iitrx]->Clone(); break;

      case 8 : histy[ij] = (TH2F*)triggereffi_xevt[ij][iitrx]->Clone(); break;
      case 9 : histy[ij] = (TH2F*)triggereffi_yevt[ij][iitrx]->Clone(); break;

      case 10 : histy[ij] = (TH2F*)triggereffi_x[ij][iitrx]->Clone(); break;
      case 11 : histy[ij] = (TH2F*)triggereffi_y[ij][iitrx]->Clone(); break;

      case 12 : histy[ij] = (TH2F*)fine_triggereffi_x[ij][iitrx]->Clone(); break;
      case 13 : histy[ij] = (TH2F*)fine_triggereffi_y[ij][iitrx]->Clone(); break;

      case 14 : histy[ij] = (TH2F*)difefficiency_uncx[ij][iitrx]->Clone(); break;
      case 15 : histy[ij] = (TH2F*)difefficiency_uncy[ij][iitrx]->Clone(); break;

      case 16 : histy[ij] = (TH2F*)diftriggereffi_x[ij][iitrx]->Clone(); break;
      case 17 : histy[ij] = (TH2F*)diftriggereffi_y[ij][iitrx]->Clone(); break;

      case 18 : histy[ij] = (TH2F*)inefficiency_xt[ij][iitrx]->Clone(); break;
      case 19 : histy[ij] = (TH2F*)inefficiency_yt[ij][iitrx]->Clone(); break;

      case 20 : histy[ij] = (TH2F*)nxystr_xtdev[ij][nmxiter]->Clone(); break;
      case 21 : histy[ij] = (TH2F*)nxystr_ytdev[ij][nmxiter]->Clone(); break;

      case 22 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter]->Clone(); break;
      case 23 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter]->Clone(); break;

      case 24 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+1]->Clone(); break;
      case 25 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+1]->Clone(); break;

      case 26 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+2]->Clone(); break;
      case 27 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+2]->Clone(); break;

      case 28 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+3]->Clone(); break;
      case 29 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+3]->Clone(); break;
      case 30 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+4]->Clone(); break;
      case 31 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+4]->Clone(); break;

      case 32 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+5]->Clone(); break;
      case 33 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+5]->Clone(); break;

      case 34 : histy[ij] = (TH2F*)xystr_xtdev[ij][0]->Clone(); break;
      case 35 : histy[ij] = (TH2F*)xystr_ytdev[ij][0]->Clone(); break;

      case 36 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter-1]->Clone(); break;
      case 37 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter-1]->Clone(); break;

      default :histy[ij] = (TH2F*)inefficiency_uncx[ij][2*nmxiter-1]->Clone(); break;
      }
      double yy1=histy[ij]->GetMaximum();
      if (yy1>xmx) xmx=yy1;

      double yy2=histy[ij]->GetMinimum();
      if (yy2<xmn) xmn=yy2;
    }
    if (jkl>=7 && xmn<0.1) { xmn=0.1;} //Avoid no hit area

    for (int ij=0; ij<nlayer; ij++) {
      c4a->cd(ij+1);
      //      histy[ij]->SetMinimum(xmn);
      //      histy[ij]->SetMaximum(xmx);
       gPad->SetRightMargin(0.15);
       gPad->SetTopMargin(0.13);
       gPad->SetLeftMargin(0.07);
       histy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip-0.5);
      histy[ij]->GetYaxis()->SetRangeUser(-0.5, nusedstrip-0.5);

      histy[ij]->GetXaxis()->SetLabelSize(.07);
      histy[ij]->GetYaxis()->SetLabelSize(.07);
      histy[ij]->GetZaxis()->SetLabelSize(.06);
      histy[ij]->GetXaxis()->SetLabelOffset(-.001);
      histy[ij]->GetYaxis()->SetLabelOffset(.01);

      histy[ij]->GetXaxis()->SetTitle("X-Strip No");
      histy[ij]->GetXaxis()->CenterTitle();
      histy[ij]->GetXaxis()->SetTitleOffset(0.8);
      histy[ij]->GetXaxis()->SetTitleSize(0.07);

      histy[ij]->GetYaxis()->SetTitle("Y-Strip No");
      histy[ij]->GetYaxis()->CenterTitle();
      histy[ij]->GetYaxis()->SetTitleOffset(0.8);
      histy[ij]->GetYaxis()->SetTitleSize(0.07);

      //      if(jkl==47 || jkl == 48) {histy[ij]->GetZaxis()->SetRangeUser(0.,1.); }

      //      histy[ij]->GetXaxis()->SetTitle(histy[ij]->GetTitle());
      //      histy[ij]->GetXaxis()->SetTitleSize(.10);
      //      histy[ij]->GetXaxis()->SetTitleOffset(.5);
      histy[ij]->SetNdivisions(506,"XY");
      histy[ij]->Draw("colz");
      if (jkl>=2 && jkl<=13) {
	latex.DrawLatex(0.4, 0.49, Form("#scale[0.5]{%g}", int(10000*get_average1(histy[ij]))/100.));
      }
    }

    c4a->Update();
    for (int ij=0; ij<nlayer; ij++) {
      if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
    }
  }


#endif //define ISEFFICIENCY
  if (!isOnlyCom) {
    TH2F* histyy[nlayer]={0};
    for (int ixx=0; ixx<8; ixx++) {
      if (!isTiming && ixx>=4) continue;
      ps.NewPage();

      for (int ij=0; ij<nlayer; ij++) {
	c4a->cd(ij+1);
	switch(ixx) {
	case 0 : histyy[ij] = (TH2F*)h_xcorstrips[ij]->Clone(); break;
	case 1 : histyy[ij] = (TH2F*)h_ycorstrips[ij]->Clone(); break;
	case 2 : histyy[ij] = (TH2F*)h_xycorstrips[ij]->Clone(); break;
	case 3 : histyy[ij] = (TH2F*)h_yxcorstrips[ij]->Clone(); break;
	case 4 : histyy[ij] = (TH2F*)h_xtcorstrips[ij]->Clone(); break;
	case 5 : histyy[ij] = (TH2F*)h_ytcorstrips[ij]->Clone(); break;
	case 6 : histyy[ij] = (TH2F*)h_xytcorstrips[ij]->Clone(); break;
	case 7 : histyy[ij] = (TH2F*)h_yxtcorstrips[ij]->Clone(); break;
	defalt : histyy[ij] = (TH2F*)h_xcorstrips[ij]->Clone(); break;
	}
	if (histyy[ij]->GetMaximum()>histyy[ij]->GetMinimum()+0.01) {
	  histyy[ij]->SetMaximum(min(1.0, histyy[ij]->GetMaximum()));
	  histyy[ij]->SetMinimum(max(1.e-3, histyy[ij]->GetMinimum()));
	}
	histyy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
	histyy[ij]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);

	histyy[ij]->GetXaxis()->SetTitle("Extrapolated Strip");
	histyy[ij]->GetXaxis()->CenterTitle();
	histyy[ij]->GetXaxis()->SetTitleOffset(0.8);
	histyy[ij]->GetXaxis()->SetTitleSize(0.07);

	histyy[ij]->GetYaxis()->SetTitle("Observed Strip");
	histyy[ij]->GetYaxis()->CenterTitle();
	histyy[ij]->GetYaxis()->SetTitleOffset(0.7);
	histyy[ij]->GetYaxis()->SetTitleSize(0.07);

	histyy[ij]->GetZaxis()->SetLabelSize(.06);
	histyy[ij]->Draw("colz");
      }
      c4a->Update();
      for (int ij=0; ij<nlayer; ij++) {
	if (histyy[ij]) { delete histyy[ij]; histyy[ij]=0;}
      }
    }
  }

  ps.NewPage();
  gStyle->SetOptLogz(1);
  TCanvas* c11 = new TCanvas("c11","c11",700,900);
  c11->Divide(3,4);
  for(int jk=0;jk<6;jk++) {
    TH2F* histyy[nlayer]={0};
    if (jk>0) ps.NewPage();

    // if(jk<2) {
    //  gStyle->SetOptLogz(0);
    // }else {
    //gStyle->SetOptLogz(1);
     // }

    for(int ij=0;ij<nlayer;ij++) {
      c11->cd(ij+1);
      switch(jk) {
      case 0 : histyy[ij] = (TH2F*)h_raw_xcorstrips[ij]->Clone(); break;
      case 1 : histyy[ij] = (TH2F*)h_raw_ycorstrips[ij]->Clone(); break;
      case 2 : histyy[ij] = (TH2F*)h_raw_xstrpnhits[ij]->Clone(); break;
      case 3 : histyy[ij] = (TH2F*)h_raw_ystrpnhits[ij]->Clone(); break;
      case 4 : histyy[ij] = (TH2F*)h_raw_xystrpnhits[ij]->Clone(); break;
      case 5 : histyy[ij] = (TH2F*)h_raw_yxstrpnhits[ij]->Clone(); break;
      default : histyy[ij] = (TH2F*)h_raw_xcorstrips[ij]->Clone(); break;
      }

      if (histyy[ij]->GetMaximum()>histyy[ij]->GetMinimum()+0.01) {
	histyy[ij]->SetMaximum(min(1.0, histyy[ij]->GetMaximum()));
	histyy[ij]->SetMinimum(max(1.e-3, histyy[ij]->GetMinimum()));
      }
      histyy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      histyy[ij]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);

      histyy[ij]->GetXaxis()->SetTitle("Observed Strip");
      histyy[ij]->GetXaxis()->CenterTitle();
      histyy[ij]->GetXaxis()->SetTitleOffset(0.8);
      histyy[ij]->GetXaxis()->SetTitleSize(0.07);

      histyy[ij]->GetYaxis()->SetTitle("Observed Strip");
      histyy[ij]->GetYaxis()->CenterTitle();
      histyy[ij]->GetYaxis()->SetTitleOffset(0.7);
      histyy[ij]->GetYaxis()->SetTitleSize(0.07);

      histyy[ij]->GetZaxis()->SetLabelSize(.06);
      histyy[ij]->Draw("colz");
    }

    c11->Update();
    for (int ij=0; ij<nlayer; ij++) {
      if (histyy[ij]) { delete histyy[ij]; histyy[ij]=0;}
    }
  }

  for (int ij=0; ij<nMuselec; ij++) {
    for (int jk = 0; jk<nEnv; jk++) {
      ps.NewPage();
      for (int il=0; il<nlayer; il++) {
	c11->cd(il+1);
	mult2d_selecxy[il][ij][jk]->GetXaxis()->SetTitle("# of Observed X-Strip");
	mult2d_selecxy[il][ij][jk]->GetXaxis()->CenterTitle();
	mult2d_selecxy[il][ij][jk]->GetXaxis()->SetTitleOffset(0.8);
	mult2d_selecxy[il][ij][jk]->GetXaxis()->SetTitleSize(0.07);

	mult2d_selecxy[il][ij][jk]->GetYaxis()->SetTitle("# of Observed Y-Strip");
	mult2d_selecxy[il][ij][jk]->GetYaxis()->CenterTitle();
	mult2d_selecxy[il][ij][jk]->GetYaxis()->SetTitleOffset(0.7);
	mult2d_selecxy[il][ij][jk]->GetYaxis()->SetTitleSize(0.07);

	mult2d_selecxy[il][ij][jk]->GetZaxis()->SetLabelSize(.06);
	mult2d_selecxy[il][ij][jk]->Draw("colz");
      }
      c11->Update();
    }
  }

  ps.NewPage();

  for(int jk=0;jk<8;jk++) {
    TH2F* histyy[nlayer]={0};
    if (jk>0) ps.NewPage();

    // if(jk<2) {
    //  gStyle->SetOptLogz(0);
    // }else {
    //gStyle->SetOptLogz(1);
     // }

    for(int ij=0;ij<nlayer;ij++) {
      c11->cd(ij+1);
      switch(jk) {
      case 0 : histyy[ij] = (TH2F*)h_xmucorstrips[ij]->Clone(); break;
      case 1 : histyy[ij] = (TH2F*)h_ymucorstrips[ij]->Clone(); break;
      case 2 : histyy[ij] = (TH2F*)h_xymucorstrips[ij]->Clone();  ;break;
      case 3 : histyy[ij] = (TH2F*)h_yxmucorstrips[ij]->Clone();  ;break;
      case 4 : histyy[ij] = (TH2F*)h_xmucornhits[ij]->Clone(); break;
      case 5 : histyy[ij] = (TH2F*)h_ymucornhits[ij]->Clone(); break;
      case 6 : histyy[ij] = (TH2F*)h_xymucornhits[ij]->Clone(); break;
      case 7 : histyy[ij] = (TH2F*)h_yxmucornhits[ij]->Clone(); break;
      default : histyy[ij] = (TH2F*)h_xmucorstrips[ij]->Clone(); break;
      }

      if (histyy[ij]->GetMaximum()>histyy[ij]->GetMinimum()+0.01) {
	histyy[ij]->SetMaximum(min(1.0, histyy[ij]->GetMaximum()));
	histyy[ij]->SetMinimum(max(1.e-3, histyy[ij]->GetMinimum()));
      }
      histyy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      histyy[ij]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);

      histyy[ij]->GetXaxis()->SetTitle("#mu-Strip");
      histyy[ij]->GetXaxis()->CenterTitle();
      histyy[ij]->GetXaxis()->SetTitleOffset(0.8);
      histyy[ij]->GetXaxis()->SetTitleSize(0.07);

      histyy[ij]->GetYaxis()->SetTitle("Observed Strips");
      histyy[ij]->GetYaxis()->CenterTitle();
      histyy[ij]->GetYaxis()->SetTitleOffset(0.7);
      histyy[ij]->GetYaxis()->SetTitleSize(0.07);

      histyy[ij]->GetZaxis()->SetLabelSize(.06);
      histyy[ij]->Draw("colz");
    }

    c11->Update();
    for (int ij=0; ij<nlayer; ij++) {
      if (histyy[ij]) { delete histyy[ij]; histyy[ij]=0;}
    }
  }


  /*
  // if(c11) {delete c11;}
   ps.NewPage();
   for(int jk=0;jk<4;jk++) {
    TH2F* histyy[nlayer]={0};
    if (jk>0) ps.NewPage();

    if(jk<2) {
      gStyle->SetOptLogz(0);
    }else {
     gStyle->SetOptLogz(1);
    }

    for(int ij=0;ij<nlayer;ij++) {
      c11->cd(ij+1);
      switch(jk) {
      case 0 : histyy[ij] = (TH2F*)h_xtstrpmult_xtdchitmult[ij]->Clone(); break;
      case 1 : histyy[ij] = (TH2F*)h_xtstrpmult_xtdchitmult[ij]->Clone(); break;
      case 2 : histyy[ij] = (TH2F*)h_xtstrpmult_xtdchitmult[ij]->Clone(); break;
      case 3 : histyy[ij] = (TH2F*)h_xtstrpmult_xtdchitmult[ij]->Clone(); break;
      default : histyy[ij] = (TH2F*)h_xtstrpmult_xtdchitmult[ij]->Clone(); break;
      }

      if (histyy[ij]->GetMaximum()>histyy[ij]->GetMinimum()+0.01) {
	histyy[ij]->SetMaximum(min(1.0, histyy[ij]->GetMaximum()));
	histyy[ij]->SetMinimum(max(1.e-3, histyy[ij]->GetMinimum()));
      }
      histyy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      histyy[ij]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);

      histyy[ij]->GetXaxis()->SetTitle("Strip_Mult");
      histyy[ij]->GetXaxis()->CenterTitle();
      histyy[ij]->GetXaxis()->SetTitleOffset(0.8);
      histyy[ij]->GetXaxis()->SetTitleSize(0.07);

      histyy[ij]->GetYaxis()->SetTitle("TDC_Mult");
      histyy[ij]->GetYaxis()->CenterTitle();
      histyy[ij]->GetYaxis()->SetTitleOffset(0.7);
      histyy[ij]->GetYaxis()->SetTitleSize(0.07);

      histyy[ij]->GetZaxis()->SetLabelSize(.06);
      histyy[ij]->Draw("colz");
    }

    c11->Update();
    for (int ij=0; ij<nlayer; ij++) {
      if (histyy[ij]) { delete histyy[ij]; histyy[ij]=0;}
    }
  }
*/
  ps.NewPage();
  gStyle->SetOptLogz(0);
  gStyle->SetPadTopMargin(0.09); //0.01);
  gStyle->SetOptStat(0);

  TCanvas* c6 = new TCanvas("c6", "c6", 700, 900);
  c6->Divide(2,4);

  //  ps.NewPage();
  c6->cd(1); xlayer_alloccu->SetLineColor(1);    xlayer_alloccu->Draw();
  c6->cd(2); ylayer_alloccu->SetLineColor(1);    ylayer_alloccu->Draw();
  double entry = max(1.,h_xrawcorhits->GetBinContent(0,0));
  c6->cd(3);  h_xrawcorhits->Scale(1./entry);  h_xrawcorhits->Draw("colz");
  c6->cd(5);  h_yrawcorhits->Scale(1./entry);  h_yrawcorhits->Draw("colz");
  c6->cd(7);  h_xyrawcorhits->Scale(1./entry); h_xyrawcorhits->Draw("colz");

  if (isTiming) {
    c6->cd(4);  h_xtrawcorhits->Scale(1./entry);  h_xtrawcorhits->Draw("colz");
    c6->cd(6);  h_ytrawcorhits->Scale(1./entry);  h_ytrawcorhits->Draw("colz");
    c6->cd(8);  h_xytrawcorhits->Scale(1./entry); h_xytrawcorhits->Draw("colz");
  }

  c6->Update();

  ps.NewPage();
  c6->cd(1); xlayer_alloccu->SetLineColor(1);    xlayer_alloccu->Draw();
             xlayer_alloccusel->SetLineColor(2); xlayer_alloccusel->Draw("same");
  c6->cd(2); ylayer_alloccu->SetLineColor(1);    ylayer_alloccu->Draw();
             ylayer_alloccusel->SetLineColor(2); ylayer_alloccusel->Draw("same");
  c6->cd(3);  h_xcorhits->Scale(1./nTotalp);  h_xcorhits->Draw("colz");
  c6->cd(5);  h_ycorhits->Scale(1./nTotalp);  h_ycorhits->Draw("colz");
  c6->cd(7);  h_xycorhits->Scale(1./nTotalp); h_xycorhits->Draw("colz");

  if (isTiming) {
    c6->cd(4);  h_xtcorhits->Scale(1./nTotalt);  h_xtcorhits->Draw("colz");
    c6->cd(6);  h_ytcorhits->Scale(1./nTotalt);  h_ytcorhits->Draw("colz");
    c6->cd(8);  h_xytcorhits->Scale(1./nTotalt); h_xytcorhits->Draw("colz");
  }
  c6->Update();

  ps.NewPage();
  if (isalign>0) {
    c6->cd(1); shift_pos->Draw("colz");
    c6->cd(2); rms_pos->Draw("colz");
  }
  c6->cd(3); h_deltaposxcov->Draw("colz");
  c6->cd(4); h_deltaposycov->Draw("colz");
  //  if (isTiming) {
  //    c6->cd(5); time_mean_reso->Draw("colz");
  //    c6->cd(6); time_rms_reso->Draw("colz");
  //    c6->cd(7); h_deltatcov->Draw("colz");
  //    c6->cd(8); h_deltatcovy->Draw("colz");
  //    //    c6->cd(6); h_deltatcov2->Draw("colz");
  //  }
  c6->Update();


  if (isTiming) {
    ps.NewPage();
    c6->cd(1); time_mean_reso->Draw("colz");
    c6->cd(2); time_rms_reso->Draw("colz");
    c6->cd(3); time_exterr_reso->Draw("colz");
    c6->cd(4); time_corrms_reso->Draw("colz");
    c6->cd(5); h_deltatcov->Draw("colz");
    c6->cd(6); h_deltatcov2->Draw("colz");
    c6->cd(7); h_deltatcovy->Draw("colz");

    c6->Update();
  }
  TH2F* hist_x_lay = new TH2F("hist_x_lay", "hist_x_lay", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* hist_y_lay = new TH2F("hist_y_lay", "hist_y_lay", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  if (posxcoreffCount>0 && posycoreffCount>0) {
    h_posxcoreff->Scale(1./nTotallx);
    h_posycoreff->Scale(1./nTotally);
    //h_posxcoreff->
    file_out<<"Xlay_Corr_eff ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<h_posxcoreff->GetBinContent(ix+1,iy+1)<<", \t";
	// hist_x_lay->Fill(ix,iy,h_posxcoreff->GetBinContent(ix,iy)/max(1.,posxcoreffCount));
      }
      file_out<<endl;
    }
    file_out<<"}"<<endl;

    file_out<<"Ylay_Corr_eff ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
	file_out<<h_posycoreff->GetBinContent(ix+1,iy+1)<<", \t";
	// hist_y_lay->Fill(ix,iy,h_posycoreff->GetBinContent(ix,iy)/max(1.,posycoreffCount));
      }
      file_out<<endl;
    }
    file_out<<"}"<<endl;
  }

  ps.NewPage();
  gStyle->SetOptLogz(0);
  gStyle->SetPadTopMargin(0.10); //0.01);
  gStyle->SetOptStat(0);
  TCanvas* c16 = new TCanvas("c16","c16",700,900);
  c16->Divide(1,2);
  c16->cd(1);
  h_posxcoreff->Draw("colz");
  h_posxcoreff->SetTitle("Cor_lay_eff_X");
  c16->Update();
  c16->cd(2);
  h_posycoreff->Draw("colz");
  h_posycoreff->SetTitle("Cor_lay_eff_Y");
  c16->Update();

  h_posxcoreff_dif->Add(h_posxcoreff,h_posxcoreff_def,1.,-1.);
  h_posycoreff_dif->Add(h_posycoreff,h_posycoreff_def,1.,-1.);
  ps.NewPage();
  gStyle->SetOptLogz(0);
  gStyle->SetPadTopMargin(0.10); //0.01);
  gStyle->SetOptStat(0);
  TCanvas* c17 = new TCanvas("c17","c17",700,900);
  c17->Divide(1,2);
  c17->cd(1);
  h_posxcoreff_dif->Draw("colz");
  h_posxcoreff_dif->SetTitle("Cor_lay_eff_X_dif");
  c17->Update();
  c17->cd(2);
  h_posycoreff_dif->Draw("colz");
  h_posycoreff_dif->SetTitle("Cor_lay_eff_Y_dif");
  c17->Update();







  if (isTiming) {
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetPadTopMargin(0.11);
    gStyle->SetPadBottomMargin(0.07);

    if (isalign>0) {
      ps.NewPage();
      gStyle->SetPadLeftMargin(1.0);
      TCanvas* c5 = new TCanvas("c5", "c5", 700, 900);
      c5->Divide(3,4);

      for (int kl=0; kl<=nmxiter; kl++) {
	c5->cd(kl+1);
	time_both_cormean[kl]->GetXaxis()->SetLabelSize(0.08);
	time_both_cormean[kl]->GetYaxis()->SetLabelSize(0.08);
	time_both_cormean[kl]->GetZaxis()->SetLabelSize(0.065);
	time_both_cormean[kl]->SetMaximum(0.5); //min(2.,time_both_cormean[kl]->GetMaximum()));
	time_both_cormean[kl]->SetMinimum(-0.5); //max(-2.,time_both_cormean[kl]->GetMinimum()));
	time_both_cormean[kl]->GetXaxis()->SetLabelOffset(0.001);
	time_both_cormean[kl]->Draw("colz");
      }
      c5->Update();

      ps.NewPage();
      for (int kl=0; kl<=nmxiter; kl++) {
	c5->cd(kl+1);
	time_both_corrms[kl]->GetXaxis()->SetLabelSize(0.08);
	time_both_corrms[kl]->GetYaxis()->SetLabelSize(0.08);
	time_both_corrms[kl]->GetZaxis()->SetLabelSize(0.08);
	time_both_corrms[kl]->SetMaximum(2.5); //min(3.5,time_both_corrms[kl]->GetMaximum()));
	time_both_corrms[kl]->SetMinimum(0.5); //max(0.5,time_both_corrms[kl]->GetMinimum()));
	time_both_corrms[kl]->GetXaxis()->SetLabelOffset(0.001);
	time_both_corrms[kl]->Draw("colz");
      }
      c5->Update();

      ps.NewPage();
      gStyle->SetPadBottomMargin(0.07);
      TCanvas* c5b = new TCanvas("c5b", "c5b", 700, 900);
      c5b->Divide(4,5);
      for (int jk=0; jk<=npixel; jk++) { //150123
	c5b->cd(jk+1);
	timexy_cormean[jk]->GetXaxis()->SetLabelSize(0.08);
	timexy_cormean[jk]->GetYaxis()->SetLabelSize(0.08);

	timexy_cormean[jk]->GetZaxis()->SetLabelSize(0.08);
	if (timexy_cormean[jk]->GetMaximum()>timexy_cormean[jk]->GetMinimum()+0.1) {
	  timexy_cormean[jk]->SetMaximum(min(2.5,timexy_cormean[jk]->GetMaximum()));
	  timexy_cormean[jk]->SetMinimum(max(-2.5,timexy_cormean[jk]->GetMinimum()));
	}
	timexy_cormean[jk]->Draw("colz");
      }

      c5b->Update();
      ps.NewPage();
      for (int jk=0; jk<=npixel; jk++) {
	c5b->cd(jk+1);
	timexy_corrms[jk]->GetXaxis()->SetLabelSize(0.08);
	timexy_corrms[jk]->GetYaxis()->SetLabelSize(0.08);
	timexy_corrms[jk]->GetZaxis()->SetLabelSize(0.08);
	timexy_corrms[jk]->SetMaximum(2.0); //min(2.4,timexy_corrms[jk]->GetMaximum()));
	timexy_corrms[jk]->SetMinimum(0.3); //max(0.4,timexy_corrms[jk]->GetMinimum()));
	timexy_corrms[jk]->Draw("colz");
      }

      c5b->Update();
      double amx[11] ={0, 0, 0, 10.0, 0.6, 1.0, 3.0, 1.5, 6.0, 4.5, 3.5};
      double amn[11]={10000, .5, .5, -10.0, -0.6, -1.0, 0.5, -1.0, -2.0, 0.0, 0.5};
      for (int ixx=0; ixx<11; ixx++) {
	TH1F* tmphist[nmxiter];
	ps.NewPage();
	for (int jk=0; jk<nmxiter; jk++) {
	  c5->cd(jk+1);
	  switch(ixx) {
	  case 0 : tmphist[jk] = (TH1F*)time_entry[jk]->Clone(); break;
	  case 1 : tmphist[jk] = (TH1F*)time_underflow[jk]->Clone(); break;
	  case 2 : tmphist[jk] = (TH1F*)time_overflow[jk]->Clone(); break;
	  case 3 : tmphist[jk] = (TH1F*)time_offset[jk]->Clone(); break;
	  case 4 : tmphist[jk] = (TH1F*)shift_time_mnft[jk]->Clone(); break;
	  case 5 : tmphist[jk] = (TH1F*)statmean_time[jk]->Clone(); break;
	  case 6 : tmphist[jk] = (TH1F*)statrms_time[jk]->Clone(); break;
	  case 7 : tmphist[jk] = (TH1F*)statskew_time[jk]->Clone(); break;
	  case 8 : tmphist[jk] = (TH1F*)statkurt_time[jk]->Clone(); break;
	  case 9 : tmphist[jk] = (TH1F*)rms_time[jk]->Clone(); break;
	  case 10 : tmphist[jk] = (TH1F*)rms_timeused[jk]->Clone(); break;
	  }

	  if ((jk>2) && (jk>3 || tmphist[jk]->GetMaximum() >tmphist[jk]->GetMinimum()+0.1)) {
	    tmphist[jk]->SetMaximum(min(amx[ixx], tmphist[jk]->GetMaximum()));
	    tmphist[jk]->SetMinimum(max(amn[ixx], tmphist[jk]->GetMinimum()));
	  }
	  tmphist[jk]->GetXaxis()->SetLabelSize(0.08);
	  tmphist[jk]->GetYaxis()->SetLabelSize(0.08);
	  tmphist[jk]->GetXaxis()->SetLabelOffset(0.001);
	  tmphist[jk]->Draw();
	  if (ixx==9) {rms_timeused[jk]->SetLineColor(2); rms_timeused[jk]->Draw("same");}
	}
	c5->Update();
	for (int jk=0; jk<nmxiter; jk++) {
	  if (tmphist[jk]) { delete tmphist[jk]; tmphist[jk]=0;}
	}
      }

      ps.NewPage();
      gStyle->SetOptLogz(0);
      gStyle->SetTitleFontSize(0.07);
      gStyle->SetPadLeftMargin(0.09);
      gStyle->SetPadBottomMargin(0.09);
      gStyle->SetPadTopMargin(0.10); //0.03);

      TCanvas* c4b = new TCanvas("c4b", "c4b", 700, 900);
      c4b->Divide(3,4);
      for (int jkl=0; jkl<4; jkl++) {
	if (jkl>0) { ps.NewPage();}
	TH2F* tmp2hist[nlayer];
	for (int ij=0; ij<nlayer; ij++) {
	  switch(jkl) {
	  case 0 : tmp2hist[ij] = (TH2F*)correction_xtime[ij]->Clone(); break;
	  case 1 : tmp2hist[ij] = (TH2F*)fitted_rms_xtime[ij]->Clone(); break;
	  case 2 : tmp2hist[ij] = (TH2F*)correction_ytime[ij]->Clone(); break;
	  case 3 : tmp2hist[ij] = (TH2F*)fitted_rms_ytime[ij]->Clone(); break;
	  }
	  c4b->cd(ij+1);
	  tmp2hist[ij]->GetXaxis()->SetLabelSize(0.07);
	  tmp2hist[ij]->GetYaxis()->SetLabelSize(0.07);
	  tmp2hist[ij]->GetZaxis()->SetLabelSize(0.06);
	  if (jkl==0 || jkl==2) {
	    tmp2hist[ij]->SetMaximum(0.50); // 0.10); // min( 0.20, correction_xtime[ij]->GetMaximum()));
	    tmp2hist[ij]->SetMinimum(-0.50); //0.10); //max(-0.20, correction_xtime[ij]->GetMinimum()));
	  } else {
	    tmp2hist[ij]->SetMaximum(4.00); // 0.10); // min( 0.20, correction_xtime[ij]->GetMaximum()));
	    tmp2hist[ij]->SetMinimum(0.50); //0.10); //max(-0.20, correction_xtime[ij]->GetMinimum()));
	  }

	  tmp2hist[ij]->GetXaxis()->SetTitle("Strip No");
	  tmp2hist[ij]->GetXaxis()->CenterTitle();
	  tmp2hist[ij]->GetXaxis()->SetTitleOffset(0.8);
	  tmp2hist[ij]->GetXaxis()->SetTitleSize(0.08);

	  tmp2hist[ij]->GetYaxis()->SetTitle("# of iteration");
	  tmp2hist[ij]->GetYaxis()->CenterTitle();
	  tmp2hist[ij]->GetYaxis()->SetTitleOffset(0.5);
	  tmp2hist[ij]->GetYaxis()->SetTitleSize(0.08);

	  tmp2hist[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
	  tmp2hist[ij]->Draw("colz");
	}
	c4b->Update();
	for (int ij=0; ij<nlayer; ij++) {
	if (tmp2hist[ij]) { delete tmp2hist[ij]; tmp2hist[ij]=0;}
	}
      }



      ps.NewPage();
      gStyle->SetOptLogy(0);
      gStyle->SetOptStat(1110);
      gStyle->SetStatW(0.36);
      gStyle->SetStatH(0.28); //30);
      gStyle->SetPadRightMargin(0.02);
      gStyle->SetTitleFontSize(0.07);
      gStyle->SetPadTopMargin(0.10);
      gStyle->SetPadBottomMargin(0.07);

      TCanvas* c5a = new TCanvas("c5a", "c5a", 700, 900);
      c5a->Divide(3,4);

      TH1F* histx[nmxiter];
      TH1F* histy[nmxiter];

      for (int jkl=0; jkl<8; jkl++) {
	if (jkl>0) ps.NewPage();
	for (int jk=0; jk<nmxiter; jk++) {
	  switch(jkl) {
	  case 0 : histx[jk] = (TH1F*)time_offsetx[jk]->Clone();
	    histy[jk] = (TH1F*)time_offsety[jk]->Clone(); break;

	  case 1 : histx[jk] = (TH1F*)shift_time_mnftx[jk]->Clone();
	    histy[jk] = (TH1F*)shift_time_mnfty[jk]->Clone(); break;

	  case 2 : histx[jk] = (TH1F*)statmean_timex[jk]->Clone();
	    histy[jk] = (TH1F*)statmean_timey[jk]->Clone(); break;

	  case 3 : histx[jk] = (TH1F*)statrms_timex[jk]->Clone();
	    histy[jk] = (TH1F*)statrms_timey[jk]->Clone(); break;

	  case 4 : histx[jk] = (TH1F*)statskew_timex[jk]->Clone();
	    histy[jk] = (TH1F*)statskew_timey[jk]->Clone(); break;

	  case 5 : histx[jk] = (TH1F*)statkurt_timex[jk]->Clone();
	    histy[jk] = (TH1F*)statkurt_timey[jk]->Clone(); break;

	  case 6 : histx[jk] = (TH1F*)rms_timex[jk]->Clone();
	    histy[jk] = (TH1F*)rms_timey[jk]->Clone(); break;

	  case 7 : histx[jk] = (TH1F*)rms_timeusedx[jk]->Clone();
	    histy[jk] = (TH1F*)rms_timeusedy[jk]->Clone(); break;

	  default : histx[jk] = (TH1F*)statkurt_time[jk]->Clone();
	    histy[jk] = (TH1F*)statkurt_timey[jk]->Clone(); break;
	  }
	}
	gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
	for (int jk=0; jk<nmxiter; jk++) {
	  c5a->cd(jk+1);
	  histx[jk]->GetXaxis()->SetLabelSize(0.065);
	  histx[jk]->GetYaxis()->SetLabelSize(0.065);
	  histx[jk]->SetLineColor(1); histx[jk]->Draw();
	}
	c5a->Update();

	gStyle->SetStatY(0.77); gStyle->SetStatTextColor(2);
	for (int jk=0; jk<nmxiter; jk++) {
	  c5a->cd(jk+1);
	  histy[jk]->SetLineColor(2); histy[jk]->Draw("sames");
	}

	c5a->Update();

	for (int jk=0; jk<nmxiter; jk++) {
	  if (histx[jk]) { delete histx[jk]; histx[jk]=0;}
	  if (histy[jk]) { delete histy[jk]; histy[jk]=0;}
	}
      }
      ps.NewPage();






    }  //if (isalign>0)
  } // if (isTiming)

  // ps.NewPage();
  gStyle->SetOptLogy(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.24);
  gStyle->SetStatH(0.20); //30);
  gStyle->SetStatY(.99);
  gStyle->SetStatX(.99);
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetPadLeftMargin(0.07);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadBottomMargin(0.08);
  TH1F* tmphistt = new TH1F("event_time_diff_1","event_time_diff_1",event_time_diff->GetXaxis()->GetNbins(),xtbins);
  for(int ij=0;ij<event_time_diff->GetXaxis()->GetNbins();ij++) {
    double xcontwidth = event_time_diff->GetBinContent(ij+1)/(xtbins[ij+1] - xtbins[ij]);
    tmphistt->SetBinContent(ij+1,xcontwidth);
  }

  TCanvas* c5t = new TCanvas("c5t", "c5t", 700, 800);
  c5t->Divide(1,3);
  c5t->cd(1);
  //  if(event_time_diff->GetBinContent(20001)>0) { event_time_diff->SetBinContent(19990,event_time_diff->GetBinContent(20001));}
  gPad->SetLogy(0);
  gPad->SetLogx(1);
  event_time_diff->SetTitle(event_time_diff->GetTitle());
  event_time_diff->GetXaxis()->SetTitle("#DeltaT (ms)");
  event_time_diff->GetXaxis()->CenterTitle();
  event_time_diff->GetXaxis()->SetLabelSize(0.045);
  event_time_diff->GetXaxis()->SetLabelOffset(0.001);
  event_time_diff->GetYaxis()->SetLabelSize(0.045);
  event_time_diff->GetXaxis()->SetTitleSize(0.045);
  event_time_diff->GetXaxis()->SetTitleOffset(0.89);
 // event_time_diff->GetXaxis()->SetRangeUser(0.5, 199.);
 // event_time_diff->Fit("expo", "Q");
  event_time_diff->GetXaxis()->SetRangeUser(0., 200.);
  event_time_diff->Draw();

  c5t->cd(2);

  gPad->SetLogy(0);
  gPad->SetLogx(1);
  tmphistt->SetTitle(tmphistt->GetTitle());
  tmphistt->GetXaxis()->SetTitle("#DeltaT (ms)");
  tmphistt->GetXaxis()->CenterTitle();
  tmphistt->GetXaxis()->SetLabelSize(0.045);
  tmphistt->GetXaxis()->SetLabelOffset(0.001);
  tmphistt->GetYaxis()->SetLabelSize(0.045);
  tmphistt->GetXaxis()->SetTitleSize(0.045);
  tmphistt->GetXaxis()->SetTitleOffset(0.89);
  tmphistt->GetXaxis()->SetRangeUser(0., 200.);
  tmphistt->Fit("expo", "Q");
  //tmphist->GetXaxis()->SetRangeUser(0., 200.);
  tmphistt->Draw();
  c5t->Update();
   //if(tmphistt) {delete tmphistt; tmphistt =0; }
  c5t->cd(3);
  gPad->SetLogy(0);
  gPad->SetTicks(1,1);
  gPad->SetGridx(); //default is one
  gPad->SetGridy();

  h_corr_layer_mult->SetTitle(h_corr_layer_mult->GetTitle());
  h_corr_layer_mult->GetYaxis()->SetTitle("Layer (X : 12+Y)");
  h_corr_layer_mult->GetYaxis()->CenterTitle();
  h_corr_layer_mult->GetYaxis()->SetLabelSize(0.045);
  h_corr_layer_mult->GetYaxis()->SetLabelOffset(0.01);
  h_corr_layer_mult->GetYaxis()->SetTitleSize(0.045);
  h_corr_layer_mult->GetYaxis()->SetTitleOffset(0.69);

  h_corr_layer_mult->GetXaxis()->SetTitle("Layer (X : 12+ Y)");
  h_corr_layer_mult->GetXaxis()->CenterTitle();
  h_corr_layer_mult->GetXaxis()->SetLabelSize(0.045);
  h_corr_layer_mult->GetXaxis()->SetLabelOffset(0.001);
  h_corr_layer_mult->GetXaxis()->SetTitleSize(0.045);
  h_corr_layer_mult->GetXaxis()->SetTitleOffset(0.89);

  h_corr_layer_mult->GetZaxis()->SetLabelSize(0.05);
  h_corr_layer_mult->Draw("colz");
  c5t->Update();


  ps.NewPage();
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetPadLeftMargin(0.08);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.08);
  TCanvas* c5z = new TCanvas("c5z","c5z",800,900);
  c5z->Divide(2,2);
  c5z->cd(1);
  gPad->SetLogy(0);
  h_corr_tdc_ref->GetYaxis()->SetTitle("Other Layer No");
  h_corr_tdc_ref->GetYaxis()->CenterTitle();
  h_corr_tdc_ref->GetYaxis()->SetLabelSize(0.045);
  h_corr_tdc_ref->GetYaxis()->SetLabelOffset(0.01);
  h_corr_tdc_ref->GetYaxis()->SetTitleSize(0.045);
  h_corr_tdc_ref->GetYaxis()->SetTitleOffset(0.72);
  h_corr_tdc_ref->GetXaxis()->SetTitle("Ref Layer No");
  h_corr_tdc_ref->GetXaxis()->CenterTitle();
  h_corr_tdc_ref->GetXaxis()->SetLabelSize(0.045);
  h_corr_tdc_ref->GetXaxis()->SetLabelOffset(0.01);
  h_corr_tdc_ref->GetXaxis()->SetTitleSize(0.045);
  h_corr_tdc_ref->GetXaxis()->SetTitleOffset(0.8);
  h_corr_tdc_ref->GetZaxis()->SetLabelSize(0.04);
  h_corr_tdc_ref->Draw("colz");

  TH2F* h_corr_tdc_ref2 = (TH2F*)h_corr_tdc_ref->Clone();;

  for(int ij=0;ij<nlayer;ij++) {
    double xx = h_corr_tdc_ref2->GetBinContent(ij+1,0);
    for(int jk=0;jk<nlayer;jk++) {
      double  yy = h_corr_tdc_ref2->GetBinContent(ij+1,jk+1)*xx/100.;
      double yer = h_corr_tdc_ref2->GetBinError(ij+1,jk+1)*xx/100.;
      h_corr_tdc_ref2->SetBinContent(ij+1,jk+1,yy);
      h_corr_tdc_ref2->SetBinError(ij+1,jk+1,yer);

    }
  }

  c5z->cd(2);
  gPad->SetLogy(0);
  h_corr_tdc_ref2->GetYaxis()->SetTitle("Other Layer No");
  h_corr_tdc_ref2->GetYaxis()->CenterTitle();
  h_corr_tdc_ref2->GetYaxis()->SetLabelSize(0.045);
  h_corr_tdc_ref2->GetYaxis()->SetLabelOffset(0.01);
  h_corr_tdc_ref2->GetYaxis()->SetTitleSize(0.045);
  h_corr_tdc_ref2->GetYaxis()->SetTitleOffset(0.72);
  h_corr_tdc_ref2->GetXaxis()->SetTitle("Ref Layer No");
  h_corr_tdc_ref2->GetXaxis()->CenterTitle();
  h_corr_tdc_ref2->GetXaxis()->SetLabelSize(0.045);
  h_corr_tdc_ref2->GetXaxis()->SetLabelOffset(0.01);
  h_corr_tdc_ref2->GetXaxis()->SetTitleSize(0.045);
  h_corr_tdc_ref2->GetXaxis()->SetTitleOffset(0.8);
  h_corr_tdc_ref2->GetZaxis()->SetLabelSize(0.035);
  h_corr_tdc_ref2->Draw("colz");
  c5z->Update();


  c5z->cd(3);
  h_tdc_ref->GetXaxis()->SetTitle("# of Layer (TDC ref=0)");
  h_tdc_ref->GetXaxis()->CenterTitle();
  h_tdc_ref->GetXaxis()->SetLabelSize(0.045);
  h_tdc_ref->GetXaxis()->SetLabelOffset(0.01);
  h_tdc_ref->GetXaxis()->SetTitleSize(0.045);
  h_tdc_ref->GetXaxis()->SetTitleOffset(0.8);
  h_tdc_ref->GetYaxis()->SetLabelSize(0.045);
  h_tdc_ref->GetYaxis()->SetLabelOffset(0.01);
  h_tdc_ref->GetYaxis()->SetTitleSize(0.045);
  h_tdc_ref->Draw();
  c5z->Update();



  ps.Close();

  //

  fileOut->cd();
  fileOut->Write();
  fileOut_align->cd();
  fileOut_align->Write();
  file_out.close();
  file_out_h.close();
  file_outstr.close();
  fileOut->Close();
  fileOut_align->Close();
  // filecorOut->Close();

  // filecorOut->cd();
  // T3->Write();
  // filecorOut->Close();
  return 0;
  //--------------------------------------------------------------
}


/*

scp INORUN_20110514_0912.ine INORUN_20110514_0757.ine INORUN_20110514_0641.ine INORUN_20110514_0526.ine INORUN_20110514_0412.ine INORUN_20110514_0255.ine INORUN_20110514_0141.ine INORUN_20110514_0257.ine INORUN_20110514_0024.ine INORUN_20110514_0026.ine INORUN_20110513_2309.ine INORUN_20110513_2154.ine INORUN_20110513_2035.ine INORUN_20110513_2040.ine INORUN_20110513_1809.ine INORUN_20110513_1806.ine INORUN_20110513_1921.ine INORUN_20110514_0917.ine

root /data1/gobinda/ino/rpcdata/aperture/test_2479_0a_3.root
void testrandom() {

  TH1F* histx = new TH1F("histx", "1-D hist ", 100, -10., 10.);
  TH1F* histy = (TH1F*)gDirectory->Get("xlayer_reso_l0_i9");

  for (int ij=0; ij<1000000; ij++) {
    double xx = histy->GetRandom();
    histx->Fill(xx);
  }
  histx->Draw();

}
// 0c Xlayer=6 iter 19 xstrip 31
0,0,0,0,0,0,0,0,0,0,0.03571,0.08333,0.09091,0.06667,0.07407,0.04167,0.06452,0.03846,0,0,0,0,0,0.04762,0.07692,0.08333,0,0,0,0,0,0,
// od Xlayer=6 iter 19 xstrip 31
0,0,0,0,0,0,0,0,0,0,0.04167,0.09524,0.08333,0.03571,0.08333,0.04167,0.1111,0.04,0,0,0,0,0,0.05,0.08333,0.1667,0,0,0,0,0,0,


 hist_obritgap_ecalhcal_2015b_all_1.root hist_obritgap_ecalhcal_2015b_all_10.root hist_obritgap_ecalhcal_2015b_all_12.root hist_obritgap_ecalhcal_2015b_all_13.root hist_obritgap_ecalhcal_2015b_all_14.root hist_obritgap_ecalhcal_2015b_all_15.root hist_obritgap_ecalhcal_2015b_all_16.root hist_obritgap_ecalhcal_2015b_all_18.root hist_obritgap_ecalhcal_2015b_all_19.root hist_obritgap_ecalhcal_2015b_all_2.root hist_obritgap_ecalhcal_2015b_all_20.root hist_obritgap_ecalhcal_2015b_all_22.root hist_obritgap_ecalhcal_2015b_all_3.root hist_obritgap_ecalhcal_2015b_all_5.root hist_obritgap_ecalhcal_2015b_all_6.root hist_obritgap_ecalhcal_2015b_all_7.root hist_obritgap_ecalhcal_2015b_all_8.root hist_obritgap_ecalhcal_2015b_all_9.root


  double szxy=0, sz=0, sxy=0, sn=0, sz2=0;
  //  int nused = 0;
  double errsq=TimeError*TimeError;
  for(int ij = 0; ij < nhits; ij++) {
    // cout<<"ZPos = "<<trk->ClustsInTrack[ij]->GetZPos()<<", Time = "<<trk->ClustsInTrack[ij]->GetTime()<<endl;

    for (unsigned int ix=0; ix<trk->ClustsInTrack[ij]->HitsInCluster.size(); ix++) {
      cout<<" ix "<<ij<<" "<< ixe<<" "
	  <<setw(6)<<trk->ClustsInTrack[ij]->GetZPos()<<" "
	  <<setw(6)<<trk->ClustsInTrack[ij]->GetTime()<<" "
	  <<setw(6)<<trk->ClustsInTrack[ij]->HitsInCluster[ix]->GetTime()<<" "
	  <<setw(6)<<trk->ClustsInTrack[ij]->HitsInCluster[ix]->GetTrueTime()<<" "
	  <<setw(6)<<trk->ClustsInTrack[ij]->HitsInCluster[ix]->GetXTime()<<" "
	  <<setw(6)<<0.1*trk->ClustsInTrack[ij]->HitsInCluster[ix]->GetXTrueTime()<<" "
	  <<setw(6)<<trk->ClustsInTrack[ij]->HitsInCluster[ix]->GetYTime()<<" "
	  <<setw(6)<<0.1*trk->ClustsInTrack[ij]->HitsInCluster[ix]->GetYTrueTime()<<" "
	  <<setw(6)<<trk->ClustsInTrack[ij]->HitsInCluster[ix]->GetXStripNum()<<" "
	  <<setw(6)<<trk->ClustsInTrack[ij]->HitsInCluster[ix]->GetYStripNum()<<" "
	  <<endl;
    }


nxfail 0 0 96.6214 0.00295428 8 18.7671
ixx 1749 0 -998.443 1 96.326 1000
ixx 1749 1 -1001.22 1 96.326 1000
ixx 1749 2 98.8942 1 96.6214 2.27283
ixx 1749 3 97.0992 1 96.6798 0.419378
ixx 1749 4 98.5723 1 96.737 1.83537
ixx 1749 5 96.5422 1 96.7948 -0.252646
ixx 1749 6 95.0058 1 96.326 -1.32017
ixx 1749 7 97.5559 1 96.9091 0.646788
ixx 1749 8 95.3456 1 96.326 -0.980327
ixx 1749 9 94.404 1 97.0252 -2.62122
ixx 1749 10 -1004.86 1 96.326 1000
ixx 1749 11 -1008.39 1 96.326 1000

nxfail 0 0 96.896 0.00247011 11 52.1351
ixx 1747 0 97.8883 1 96.896 0.992347
ixx 1747 1 99.9424 1 96.9364 3.00604
ixx 1747 2 100.278 1 96.9757 3.30252
ixx 1747 3 98.6274 1 97.0162 1.6112
ixx 1747 4 98.4083 1 97.0558 1.3525
ixx 1747 5 97.6611 1 97.0958 0.565296
ixx 1747 6 94.6382 1 96.6489 -2.0107
ixx 1747 7 -1003.38 1 96.6489 1000
ixx 1747 8 95.8275 1 96.6489 -0.821432
ixx 1747 9 94.7365 1 97.2554 -2.51888
ixx 1747 10 94.1133 1 97.295 -3.18169
ixx 1747 11 94.3517 1 96.6489 -2.2972



c1->cd(1);
xtime_exterr_l0_i0->Draw();

c1->cd(2);
xtime_exterr_l0_i1->Draw();

c1->cd(3);
xtime_exterr_l0_i2->Draw();

c1->cd(4);
xtime_exterr_l0_i3->Draw();

c1->cd(5);
xtime_exterr_l0_i4->Draw();

c1->cd(6);
xtime_exterr_l0_i5->Draw();

c1->cd(7);
xtime_exterr_l0_i6->Draw();

c1->cd(8);
xtime_exterr_l0_i7->Draw();

c1->cd(9);
xtime_exterr_l0_i8->Draw();

c1->cd(10);
xtime_exterr_l0_i9->Draw();

c1->cd(11);
 T2->Draw("xtexter[0]");


Attaching file tst4yyy_test_pow2_5_timeres_0_5_0_ii_3.root as _file0...
Attaching file tst4yyy_test_pow2_5_timeres_0_6_0_ii_3.root as _file1...
Attaching file tst4yyy_test_pow2_5_timeres_0_7_0_ii_3.root as _file2...
Attaching file tst4yyy_test_pow2_5_timeres_0_8_0_ii_3.root as _file3...
Attaching file tst4yyy_test_pow2_5_timeres_0_9_0_ii_3.root as _file4...
Attaching file tst4yyy_test_pow2_5_timeres_1_0_1_ii_3.root as _file5...
Attaching file tst4yyy_test_pow2_5_timeres_1_0_2_ii_3.root as _file6...
Attaching file tst4yyy_test_pow2_5_timeres_1_0_3_ii_3.root as _file7...
Attaching file tst4yyy_test_pow2_5_timeres_1_0_rndm3_0_ii_3.root as _file8...
Warning in <TFile::Init>: file tst4yyy_test_pow2_5_timeres_1_0_rndm3_0_ii_3.root probably not closed, trying to recover
Info in <TFile::Recover>: tst4yyy_test_pow2_5_timeres_1_0_rndm3_0_ii_3.root, recovered key TTree:T2 at address 90291885
Warning in <TFile::Init>: successfully recovered 1 keys
Attaching file tst4yyy_test_pow2_5_timeres_1_1_0_ii_3.root as _file9...
Attaching file tst4yyy_test_pow2_5_timeres_1_2_0_ii_3.root as _file10...
Attaching file tst4yyy_test_pow2_5_timeres_1_3_0_ii_3.root as _file11...
Attaching file tst4yyy_test_pow2_5_timeres_1_4_0_ii_3.root as _file12...
Attaching file tst4yyy_test_pow2_5_timeres_true_3_ii_3.root as _file13...
Attaching file tst4yyy_test_pow2_5_timeres_true_4_ii_3.root as _file14...
Attaching file tst4yyy_test_pow2_5_timeres_true_5_ii_3.root as _file15...


P3-P8 : Missing entries in plots, but entries are in rootuple
P13-14 : Arbitray entries in raw_seloccu/raw_noiseocc plots

No entry in Delta_time_evt/event_time_diff

P51 : timex_shift_l8_l1 : Shift in timing (Not in _l2-_l5 is it due to correction term ?)
P121 : time_xreso_l3_i5 (P122 Strip 51) : IN higher iteration, not in combined one but in individual one
P141 : time_xreso_l7_i5 (P142 Strip 29)
P147 : tdc3, large offset
P148 : Strip 52, from where this lower tail ?

P699 : Correlation with xdev as a function of ystr
P703 : Correlation with ydev as a function of xstr
P707 : Last 8 strips in layer11
P710 : First layer, modulation of mean


Position
  x               x        x
-0.00099662  0.00617565  -0.00599652 -0.00407934
-0.000701789 -0.00743384 0.00857961  0.000678459
-0.000454604 -0.00110773 0.0030311   -0.000598354
0.00249214   0.00436821  -0.00673732 0.00106561
0.00138709   0.000787645 -0.00169856 -0.00149991
-0.000138803 0.00371555  -0.00171399 -0.00175985
-0.00356749  0.00166661  -0.00250379 0.000504382
-4.6832e-05  0.000837395 -0.00298925 0.000818629
0.000395917  -0.00313572 0.00476659  0.00252465
-0.00102688  -0.00383587 0.00322756  0.00183072
 0.00389311  -0.00142641 0.00234955  -0.00340903
 0.00212725  2.88e-05    -0.00021994 0.000718197

scal 8.09422e-07

 *** Break *** segmentation violation



===========================================================
There was a crash.
This is the entire stack trace of all threads:
===========================================================
#0  0x000000315bcac65e in waitpid () from /lib64/libc.so.6
#1  0x000000315bc3e609 in do_system () from /lib64/libc.so.6
#2  0x00007f3f7e6634e8 in TUnixSystem::StackTrace() () from /usr/local/lib/libCore.so
#3  0x00007f3f7e662363 in TUnixSystem::DispatchSignals(ESignals) () from /usr/local/lib/libCore.so
#4  <signal handler called>
#5  0x0000000000459478 in main () at anal_rpcdata.cc:6647
===========================================================


The lines below might hint at the cause of the crash.
If they do not help you then please submit a bug report at
http://root.cern.ch/bugs. Please post the ENTIRE stack trace
from above as an attachment in addition to anything else
that might help us fixing this issue.
===========================================================
#5  0x0000000000459478 in main () at anal_rpcdata.cc:6647
===========================================================


*/
