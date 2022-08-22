// Define CONSTANTS

#define CORRELATIONPSI4 0.8 // 0 for psi_4 uniform dist, 1 for PSI4=PSI2
#define CORRELATIONPSI3 0 // 0 for psi_3 uniform dist, 1 for PSI3=PSI2
#define NT 1 // Number of particles in the trigger distribution
#define VTWEET 0.15
//#define VTWEET 0.0001
#define VDRIET 0.08
//#define VDRIET 0.0
#define VVIERT 0.04
//#define VVIERT 0.0
#define NA 1 // Number of particles in the associate distribution
#define VTWEEA 0.2
//#define VTWEEA 0.0001
#define VDRIEA 0.1
//#define VDRIEA 0.0
#define VVIERA 0.05
//#define VVIERA 0.0
#define VTWEEV 0.09
#define VNUL 0 
#define NV 251 // Number of particles in the event-plane-determining distribution
#define VBINS 8 // Number of bins in that distribution


#include "TVectorT.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TMinuit.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TNtuple.h"
#include <fstream>


TRandom* RAND = new TRandom();
Double_t RANDVERDEEL = -TMath::Pi() / 128.;


// This TF1 distribution is the true distribution
TF1* true_dist(Double_t v2, Double_t v3, Double_t v4, Double_t psi3, Double_t psi4) {
  TF1* func = new TF1("truedist", "1 + 2*[0]*cos(2*x) + 2*[1]*cos(3*(x-[3])) + 2*[2]*cos(4*(x-[4]))", 0, 2 * TMath::Pi());
  func->SetParameter(0, v2);
  func->SetParameter(1, v3);
  func->SetParameter(2, v4);
  func->SetParameter(3, psi3);
  func->SetParameter(4, psi4);
  return func;
}


Double_t mapToTwoPi(Double_t rand) {
  while(1) {
    if(rand >= 0.0) {
      if(rand < 2*TMath::Pi())
        return rand;
      else
        rand-=2*TMath::Pi();
    }
    else
      rand+=2*TMath::Pi();
  }
}


// Get a event plane direction for psi3 or psi4 if boolean is on, currently both are from uniform distribution
Double_t randomevpl(Bool_t psi4) {
  if (psi4) {
    if(CORRELATIONPSI4==1)
      return 0;
    else if (CORRELATIONPSI4==0) {
      return RAND->Uniform(2 * TMath::Pi());
    }
    else {
      return mapToTwoPi(RAND->Gaus(0.0, (1./CORRELATIONPSI4)-1));
    }
  }
  else {
    if(CORRELATIONPSI3==1)
      return 0;
    else
      return RAND->Uniform(2 * TMath::Pi());
  }
//    RANDVERDEEL += TMath::Pi()/64.;
//    if(RANDVERDEEL > 2*TMath::Pi())
//      RANDVERDEEL = TMath::Pi()/128.;
//    return RANDVERDEEL;
}


// Returns the particles for 1 distribution (trigger, associate or v0)
TVectorD particle_dist(Int_t length, Double_t v2, Double_t v3, Double_t v4, Double_t psi3, Double_t psi4) {
  TVectorD vector = TVectorD(length);
  TF1* dist = true_dist(v2, v3, v4, psi3, psi4);
  dist->SetNpx(240);
  for(Int_t i = 0; i < length; i++) {
    vector[i] = dist->GetRandom();
  }
  delete dist;
  return vector;
}


// Puts the information of TVectorD on the position of particles in a given number of bins between 0 and 2pi.
TH1D* putinbins(TVectorD particles, Int_t nobins) {
  TH1D* hist = new TH1D("binnified", "hist", nobins, 0, 2*TMath::Pi());
  for(Int_t i = 0; i < particles.GetNoElements(); i++)
    hist->Fill(particles[i]);
  return hist;
}


// Computes the angle between phi1, phi2, expressed in interval [0,2pi)
Double_t angleBetween(Double_t phi1, Double_t phi2) {
  Double_t phi = phi1 - phi2;
  if(phi<0)
    phi += 2*TMath::Pi();
  if(phi>2*TMath::Pi())
    phi -= 2*TMath::Pi();
  return phi;
}


// Mirror [pi, 2pi) on [0, pi)
Double_t mirrorToRightHalf(Double_t phi) {
  if(phi>TMath::Pi())
    return 2*TMath::Pi() - phi;
  else
    return phi;
}


// Determine whether trigger is in, mid or out
Int_t determine_inmidout(Double_t phi) {
  phi = mirrorToRightHalf(phi);
  if((phi>=TMath::Pi()/3) && (phi<=2*TMath::Pi()/3))
    return 2;
  if((phi<=TMath::Pi()/6) || (phi>=5*TMath::Pi()/6))
    return 0;
  else
    return 1;
}


// Contains the information needed from a single event
class DiHadronEvent {
  public:
  TH1D* dihadrhist[3];
  TH2D* dihadrhistvspsi4[3];
  Int_t noTriggers[3];
  Double_t partplane;
  Double_t evplane3;
  Double_t evplane4;
  Bool_t planeSet;
  
  // constuctor
  DiHadronEvent() {
    for(Int_t i = 0; i < 3; i++) {
      dihadrhist[i] = new TH1D(Form("dihadron%d", i), Form("dihadron%d", i), 64, 0, 2*TMath::Pi());
      dihadrhist[i]->Sumw2();
      dihadrhistvspsi4[i] = new TH2D(Form("dihadron2D%d", i), Form("dihadron2D%d", i), 64, 0, 2*TMath::Pi(), 64, 0, 2*TMath::Pi());
      dihadrhistvspsi4[i]->Sumw2();
      noTriggers[i] = 0;
    }
    partplane = 0;
    evplane3 = 0;
    evplane4 = 0;
    planeSet = kFALSE;
  }

  // Computes the v_2 participant plane based on a TH1D
  void ComputeParticipantPlane(TH1D* dist, Int_t n = 2) {
    Double_t qx = 0;
    Double_t qy = 0;
    for(Int_t bin = 1; bin <= dist->GetNbinsX(); bin++) {
      qx += dist->GetBinContent(bin)*TMath::Cos(n*dist->GetBinCenter(bin));
      qy += dist->GetBinContent(bin)*TMath::Sin(n*dist->GetBinCenter(bin));
    }
    
    Double_t psi = TMath::ATan2(qy,qx)/n;
    if(n==2) {
      if(psi>=TMath::Pi()) psi-=TMath::Pi();
      if(psi<0) psi+=TMath::Pi();
      partplane = psi;
      planeSet = kTRUE;
      return;
    }
/*    if(n==3) {
      if(psi>TMath::Pi()) psi-=2*TMath::Pi()/3;
      if(psi<0) psi+=2*TMath::Pi()/3;
      partplane3 = psi;
      return;
    }
    if(n==4) {
      if(psi>TMath::Pi()) psi-=TMath::Pi()/2;
      if(psi<0) psi+=TMath::Pi()/2;
      partplane4 = psi;
      return;
    }*/
  }

  void SetEventPlanes(Double_t psi3, Double_t psi4) {
    evplane3 = psi3;
    evplane4 = psi4;
  }

  // computes the dihadron correlations based on the trigger and associate particles for in/mid/out with respect to an participant plane
  void FillCorrelation(TVectorD trig, TVectorD asso) {
    if(planeSet==kFALSE) {
      cout << "Need to compute plane first." << endl;
      return;
    }

    for(Int_t t = 0; t < trig.GetNoElements(); t++) {
      Double_t angleToPlane = angleBetween(trig[t], partplane);
      Int_t inmidout = determine_inmidout(angleToPlane);
      noTriggers[inmidout]++;
      for(Int_t a = 0; a < asso.GetNoElements(); a++) {
        dihadrhist[inmidout]->Fill(angleBetween(trig[t], asso[a]));
        dihadrhistvspsi4[inmidout]->Fill(angleBetween(trig[t], asso[a]), evplane4);
      }
    }

    return;
  }

  // copy constructor
  DiHadronEvent(DiHadronEvent* orig) {
    for(Int_t i = 0; i < 3; i++) {
      if(dihadrhist[i]) {
        delete dihadrhist[i];
        delete dihadrhistvspsi4[i];
      }
      dihadrhist[i] = new TH1D(*(orig->dihadrhist[i]));
      dihadrhistvspsi4[i] = new TH2D(*(orig->dihadrhistvspsi4[i]));
      noTriggers[i] = orig->noTriggers[i];
    }
    partplane = orig->partplane;
    evplane3 = orig->evplane3;
    evplane4 = orig->evplane4;
    planeSet = orig->planeSet;
  }

  // destructor
  ~DiHadronEvent() {
    if(dihadrhist)
      for(Int_t i = 0; i < 3; i++)
        if(dihadrhist[i]) {
          delete dihadrhist[i];
          dihadrhist[i] = 0;
        }
    if(dihadrhistvspsi4)
      for(Int_t i = 0; i < 3; i++)
        if(dihadrhistvspsi4[i]) {
          delete dihadrhistvspsi4[i];
          dihadrhistvspsi4[i] = 0;
        }
  }
};


// Simulates one event
DiHadronEvent* sim_event(Bool_t computePsi3 = kFALSE) {
  // simulate particles
  Double_t psi3 = randomevpl(kFALSE);
  Double_t psi4 = randomevpl(kTRUE);
  TVectorD leadpart = particle_dist(NT, VTWEET, VDRIET, VVIERT, psi3, psi4);
  TVectorD assopart = particle_dist(NA, VTWEEA, VDRIEA, VVIERA, psi3, psi4);
//  TVectorD evplpart = particle_dist(NV, VTWEEV, VDRIEA, VVIERA, psi3, psi4);
  TVectorD evplpart = particle_dist(NV, VTWEEV, 0, 0, psi3, psi4);

  // Put the simulated event plane in bins to simulate event plane resolution effects.
  TH1D* simulatedv0 = putinbins(evplpart, VBINS);

  DiHadronEvent* dihadronevent = new DiHadronEvent();
  dihadronevent->ComputeParticipantPlane(simulatedv0);
  dihadronevent->SetEventPlanes(psi3, psi4);
  dihadronevent->FillCorrelation(leadpart, assopart);

/*  if(computePsi3==kTRUE) {
    TH1D* leadbinned = putinbins(leadpart, 64);
    dihadronevent->ComputeParticipantPlane(leadbinned, 3);
    dihadronevent->ComputeParticipantPlane(leadbinned, 4);
  }*/

  // cleanup
  delete simulatedv0;

  return dihadronevent;
}


// A class to compute the sum of several events and compute necessary results.
class DiHadronSum {
  private:
  TH1D* dihadrondist[3];
  TH2D* dihadrondistvspsi4[3];
  TH1D* partplanedist;
  TH1D* evplanedist3;
  TH1D* evplanedist4;
  TH2D* evplanedist34;
  Double_t noTriggers[3];
  Int_t nrsum;

  public:
  // constructor
  DiHadronSum() {
    for(Int_t i = 0; i < 3; i++) {
      dihadrondist[i] = 0;
      dihadrondistvspsi4[i] = 0;
      noTriggers[i] = 0;
    }
    partplanedist = new TH1D("partplanedistribution", "Distribution of Participant Plane", 4000, 0, 2*TMath::Pi());
    evplanedist3 = new TH1D("evplanedistribution3", "Distribution of v3-Event Plane", 4000, 0, 2*TMath::Pi());
    evplanedist4 = new TH1D("evplanedistribution4", "Distribution of v4-Event Plane", 4000, 0, 2*TMath::Pi());
    evplanedist34 = new TH2D("evplanedistribution34", "Distribution of v3-v4-Event Plane", 4000, 0, 2*TMath::Pi(), 4000, 0, 2*TMath::Pi());
    nrsum = 0;
  }

  // destructor
  ~DiHadronSum() {
    for(Int_t i = 0; i < 3; i++) {
      if(dihadrondist[i]) delete dihadrondist[i];
      if(dihadrondistvspsi4[i]) delete dihadrondistvspsi4[i];
    }
    delete partplanedist;
    delete evplanedist3;
    delete evplanedist4;
    delete evplanedist34;
  }

  TH1D* partplanedistcopy()  {return new TH1D(*partplanedist); }
  TH1D* evplane3distcopy() {return new TH1D(*evplanedist3);}
  TH1D* evplane4distcopy() {return new TH1D(*evplanedist4);}
  TH2D* evplane34distcopy() {return new TH2D(*evplanedist34);}

  Double_t NbinsX() {return dihadrondist[0]->GetNbinsX();}

  // Add an extra event to the DiHadronSum
  void AddEvent(DiHadronEvent* event) {
    // make special case for first event, then set format of dihadrondist histograms as copy of other ones.
    if(nrsum == 0) {
      for(Int_t i = 0; i < 3; i++) {
        dihadrondist[i] = new TH1D(*(event->dihadrhist[i]));
        dihadrondist[i]->SetName("sumdist");
        dihadrondistvspsi4[i] = new TH2D(*(event->dihadrhistvspsi4[i]));
        dihadrondistvspsi4[i]->SetName("sumdistvspsi4");
        noTriggers[i] += event->noTriggers[i];
      }
    }
    else {
      for(Int_t i = 0; i < 3; i++) {
        dihadrondist[i]->Add(event->dihadrhist[i]);
        dihadrondistvspsi4[i]->Add(event->dihadrhistvspsi4[i]);
        noTriggers[i] += event->noTriggers[i];
      }
    }
    partplanedist->Fill(event->partplane);
    evplanedist3->Fill(event->evplane3);
    evplanedist4->Fill(event->evplane4);
    evplanedist34->Fill(event->evplane3, event->evplane4);
    nrsum++;
  }

  // Normalise and get the distribution
  TH1D* GetDist(Int_t plane) {
    TH1D* copy = new TH1D(*(dihadrondist[plane]));
    copy->Scale(1.0/(noTriggers[plane]*nrsum));
    return copy;
  }

  // Normalise and get the distribution
  TH2D* GetDistvspsi4(Int_t plane) {
    TH2D* copy = new TH2D(*(dihadrondistvspsi4[plane]));
    copy->Scale(1.0/(noTriggers[plane]*nrsum));
    return copy;
  }

  // Returns the R_n for n is even.
  Double_t ComputeR(Int_t power) {
    Double_t Rpow = 0;
    for(Int_t bin = 1; bin <= partplanedist->GetNbinsX(); bin++)
      Rpow += partplanedist->GetBinContent(bin) * TMath::Cos(power * partplanedist->GetXaxis()->GetBinCenter(bin));
    return Rpow / nrsum;
  }

  Double_t ComputeS(Int_t orderevplane, Int_t power) {
    Double_t Spow = 0;
    TH1D* dist = 0;
    if(orderevplane==3)
      dist = evplanedist3;
    else if (orderevplane==4)
      dist = evplanedist4;
    for(Int_t bin = 1; bin <= partplanedist->GetNbinsX(); bin++)
      Spow += dist->GetBinContent(bin) * TMath::Cos(power * partplanedist->GetXaxis()->GetBinCenter(bin));
    return Spow / nrsum;
  }

  Double_t ComputeC34() {
    Double_t Spow = 0;
    for(Int_t binx = 1; binx <= evplanedist34->GetNbinsX(); binx++)
      for(Int_t biny = 1; biny <= evplanedist34->GetNbinsY(); biny++)
        Spow += evplanedist34->GetBinContent(binx, biny) * TMath::Cos(3 * evplanedist34->GetXaxis()->GetBinCenter(binx) + 4 * evplanedist34->GetYaxis()->GetBinCenter(biny));
    return Spow / nrsum;
  }
};


class Fitter {
  private:
  // private constructor, because Singleton
  Fitter() {
    for(Int_t i = 0; i < 3; i++)
      dhd[i] = 0;
    r2 = 0; r4 = 0; r8s84 = 0; s4 = 0;
    binsX = 0;
    a = 0.8269933431; // = sqrt(27)/2pi
  }
  Double_t bgr; Double_t bgrerr;
  Double_t v2t; Double_t v2terr;
  Double_t v2a; Double_t v2aerr;
  Double_t v3; Double_t v3err;
  Double_t v4t; Double_t v4terr;
  Double_t v4a; Double_t v4aerr;
  Double_t finalchi2;
  Bool_t draw_with_smearing = kFALSE;

  public:
  static Fitter* instance;
  static TH1D** dhd;
  static Double_t r2;
  static Double_t r4;
  static Double_t r8s84;
  static Double_t s4;
  static Int_t binsX;
  static Double_t a;

  Double_t getv2a() {return v2a;}
  Double_t getv2aerr() {return v2aerr;}
  Double_t getv2t() {return v2t;}
  Double_t getv2terr() {return v2terr;}
  Double_t getv4a() {return v4a;}
  Double_t getv4aerr() {return v4aerr;}
  Double_t getv3() {return v3;}
  Double_t getv3err() {return v3err;}
  Double_t getchi2() {return finalchi2;}

  static Double_t flowInPlane(Int_t iPlane, Double_t v2t, Double_t v2a, Double_t v4t) {
    if(iPlane==0)
      return v2a*((v2t*(1+0.5*a*r4)+a*r2+a*v4t*r2*s4)/(1+2*a*v2t*r2+v4t*a*r4*s4));
    if(iPlane==1)
      return v2a*v2t*(1-a*r4)/(1-2*a*v4t*r4); // 0.23...=3/4pi
    if(iPlane==2)
      return v2a*((v2t*(1+0.5*a*r4)-a*r2-a*v4t*r2*s4)/(1-2*a*v2t*r2+v4t*a*r4*s4));
  }
  static Double_t v3inplane(Int_t iPlane, Double_t v3, Double_t v2t, Double_t v4t) {
    if(iPlane==0)
      return v3/(1+2*v2t*a*r2+v4t*a*r4*s4);
    if(iPlane==1)
      return v3/(1-2*a*v4t*r4);
    if(iPlane==2)
      return v3/(1-2*v2t*a*r2+v4t*a*r4*s4);
  }
  static Double_t v4inplane(Int_t iPlane, Double_t v4t, Double_t v4a, Double_t v2t) {
    if(iPlane==0)
      return v4a*(v4t*(1-0.25*a*r8s84)+0.5*a*r4*s4+v2t*a*r2*s4)/(1+2*v2t*a*r2+v4t*a*r4*s4);
    if(iPlane==1)
      return v4a*(v4t*(1+0.5*a*r8s84)-1*a*r4*s4)/(1-2*a*v4t*r4*s4);
    if(iPlane==2)
      return v4a*(v4t*(1-0.25*a*r8s84)+0.5*a*r4*s4-v2t*a*r2*s4)/(1-2*v2t*a*r2+v4t*a*r4*s4);
  }

  static Fitter* getInstance();
  static void Chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    Double_t chi2 = 0;
    Double_t delta = 0; Double_t flow;
//////    Translation to make parameters more readable
    Double_t B = par[0]; Double_t v2t = par[1]; Double_t v2a = par[2]; Double_t v3p = par[3]; Double_t v4t = par[4]; Double_t v4a = par[5];
    for(Int_t iPlane = 0; iPlane < 3; iPlane++) {
      Double_t v2 = flowInPlane(iPlane,v2t,v2a,v4t);
      Double_t v3 = v3inplane(iPlane,v3p,v2t,v4t);
      Double_t v4 = v4inplane(iPlane,v4t,v4a,v2t);
      for(Int_t pbin = 1; pbin < binsX + 1; pbin++) {
//////    Readout
        Double_t x   = dhd[iPlane]->GetXaxis()->GetBinCenter(pbin);
        Double_t val = dhd[iPlane]->GetBinContent(pbin);
        Double_t err = dhd[iPlane]->GetBinError(pbin);
//////    Central function
        Double_t flow = B*(1 +2*v2*TMath::Cos(2*x)
                             +2*v3*TMath::Cos(3*x)
                             +2*v4*TMath::Cos(4*x));
        if(err)
          delta = (flow - val)/err;
        else
          delta = 0;
        chi2 += delta*delta;
//////    Custom local parameter restrictions (position dependent)
        if(flow<0.0)
          chi2+=10000*flow*flow;
      }
    }
    if(npar&&gin&&iflag) f = chi2;
    else f = chi2;
  }

  Double_t computeChi2(Int_t plane) {
    Double_t chi2 = 0.; Double_t delta = 0.;
    Double_t p2 = flowInPlane(plane,v2t,v2a,v4t);
    Double_t p3 = v3inplane(plane,v3,v2t,v4t);
    Double_t p4 = v4inplane(plane,v4t,v4a,v2t);
    for(Int_t pbin = 1; pbin < binsX + 1; pbin++) {
      Double_t x = dhd[plane]->GetXaxis()->GetBinCenter(pbin);
      Double_t val = dhd[plane]->GetBinContent(pbin);
      Double_t err = dhd[plane]->GetBinError(pbin);
      Double_t flow = bgr*(1 +2*p2*TMath::Cos(2*x)
                             +2*p3*TMath::Cos(3*x)
                             +2*p4*TMath::Cos(4*x));
      if(err)
        delta = (flow - val)/err;
      else
        delta = 0;
      chi2 += delta*delta;
//////    Custom local parameter restrictions (position dependent)
      if(flow<0.0)
        chi2+=10000*flow*flow;
    }
    return chi2;
  }

  void SetDrawWithSmearingR() { draw_with_smearing = kTRUE; }

  void SetData(DiHadronSum* dihadsum) {
    r2 = dihadsum->ComputeR(2);
    r4 = dihadsum->ComputeR(4);
    r8s84 = dihadsum->ComputeR(8) * dihadsum->ComputeS(4,8);
    s4 = dihadsum->ComputeS(4,4);
    binsX = dihadsum->NbinsX();
    draw_with_smearing = kFALSE;
    for(Int_t i = 0; i < 3; i++)
      dhd[i] = dihadsum->GetDist(i);
  }

  // perform the simulatenous fit.
  void Fit() {
    TMinuit* myMinuit = new TMinuit(6);
    myMinuit->SetFCN(Chi2);
    myMinuit->SetPrintLevel(1);
    Double_t arglist[6];
    Int_t ierflg = 0;
    arglist[0] = 1;
    myMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    TString names[6] = {"<Background> 0", "v2t", "v2a", "v3", "v4t", "v4a"};
//    Double_t vstart[6] = {.5, .5, .5, .5, .5, .5};
    Double_t vstart[6] = {.00311, .15, .2, .008, .04, .05};
    Double_t precision = 1e-7;
    Double_t step[6] = {precision, precision, precision, precision, precision, precision};
    Double_t lolim[6] = {
    0, 0.0, 0.0, -1.0, -1.0, -1.0};
    Double_t uplim[6] = {
    0, 1.0, 1.0, 1.0,  1.0,  1.0};
    for(Int_t param = 0; param < 6; param++)
      myMinuit->mnparm(param , names[param].Data(), vstart[param], step[param], lolim[param], uplim[param],ierflg);
    if(CORRELATIONPSI4==0)
      myMinuit->FixParameter(4);
    arglist[0] = 50000;
    arglist[1] = 1.;
    Double_t edm, errdef;
    Int_t npari, nparx, istat;
    myMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    myMinuit->mnimpr();
    myMinuit->GetParameter(0,bgr,bgrerr);
    myMinuit->GetParameter(1,v2t,v2terr);
    myMinuit->GetParameter(2,v2a,v2aerr);
    myMinuit->GetParameter(3,v3,v3err);
    myMinuit->GetParameter(4,v4t,v4terr);
    myMinuit->GetParameter(5,v4a,v4aerr);
    myMinuit->mnstat(finalchi2, edm, errdef, npari, nparx, istat);
  }

  void Save(const char* mapnaam) {
    Plot(mapnaam);
    PlotPlaneSum(mapnaam);
    SavePars(mapnaam);
  }

  void Plot(const char* mapnaam) {
    TF1* funcgold = new TF1("fitdist", "[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))", 0, 2 * TMath::Pi());
    funcgold->SetLineColor(3);
    for(Int_t iPlane = 0; iPlane < 3; iPlane++) {
      funcgold->SetParameter(0,bgr);
      funcgold->SetParameter(1,flowInPlane(iPlane, VTWEET, VTWEEA, VVIERT));
      funcgold->SetParameter(2,v3inplane(iPlane, (VDRIET*VDRIEA), VTWEET, VVIERT));
      funcgold->SetParameter(3,v4inplane(iPlane, VVIERT, VVIERA, VTWEET));
      TF1* func = new TF1("fitdist", "[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))", 0, 2 * TMath::Pi());
      func->SetParameter(0,bgr);
      func->SetParameter(1,flowInPlane(iPlane, v2t, v2a, v4t));
      func->SetParameter(2,v3inplane(iPlane, v3, v2t, v4t));
      func->SetParameter(3,v4inplane(iPlane, v4t, v4a, v2t));
      TCanvas *c1 = new TCanvas("c1","c1:2-dim fit projection on Eta",800,600);
      dhd[iPlane]->Draw("E1");
      func->Draw("same");
      funcgold->DrawCopy("same");
      if(draw_with_smearing) {
        r2 = r2*1.1;
        funcgold->SetParameter(1,flowInPlane(iPlane, VTWEET, VTWEEA, VVIERT));
        funcgold->SetParameter(2,v3inplane(iPlane, (VDRIET*VDRIEA), VTWEET, VVIERT));
        funcgold->SetParameter(3,v4inplane(iPlane, VVIERT, VVIERA, VTWEET));
	funcgold->DrawCopy("same");
        r2 = r2*0.9/1.1;
        funcgold->SetParameter(1,flowInPlane(iPlane, VTWEET, VTWEEA, VVIERT));
        funcgold->SetParameter(2,v3inplane(iPlane, (VDRIET*VDRIEA), VTWEET, VVIERT));
        funcgold->SetParameter(3,v4inplane(iPlane, VVIERT, VVIERA, VTWEET));
	funcgold->DrawCopy("same");
	r2 = r2/0.9;
      }
      c1->SaveAs(Form("%s/plot%d.pdf", mapnaam, iPlane));
      delete func; delete c1;
    }
    delete funcgold;
  }

  void PlotPlaneSum(const char* mapnaam) {
    TF1* funcgold = new TF1("fitdist", "[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))", 0, 2 * TMath::Pi());
    funcgold->SetParameter(0,3*bgr);
    funcgold->SetParameter(1,VTWEET*VTWEEA);
    funcgold->SetParameter(2,VDRIET*VDRIEA);
    funcgold->SetParameter(3,VVIERT*VVIERA);
    TH1D* sum = new TH1D(*dhd[0]);
    sum->Add(dhd[1]);
    sum->Add(dhd[2]);
    TCanvas *c1 = new TCanvas("c1","c1:2-dim fit projection on Eta",800,600);
    sum->Draw("E1");
    funcgold->SetLineColor(3);
    funcgold->Draw("same");
    c1->SaveAs(Form("%s/plot_sum.pdf", mapnaam));
    delete funcgold; delete c1;
  }

  void SavePars(const char* mapnaam) {
    ofstream fout;
    fout.open(Form("%s/param.txt", mapnaam), ios::trunc);
    fout << "CORRELATIONPSI4     | " << CORRELATIONPSI4 << endl;
    fout << "CORRELATIONPSI3     | " << CORRELATIONPSI3 << endl;
    fout << "NT                  | " << NT << endl;
    fout << "VTWEET              | " << Form("%5f", VTWEET)        << " | " << Form("%5f", v2t)     << "+-" << Form("%5f", v2terr) << endl;
    fout << "VDRIET              | " << Form("%5f", VDRIET)        << endl;
    fout << "VVIERT              | " << Form("%5f", VVIERT)        << " | " << Form("%5f", v4t)     << "+-" << Form("%5f", v4terr) << endl;
    fout << "NA                  | " << NA << endl;
    fout << "VTWEEA              | " << Form("%5f", VTWEEA)        << " | " << Form("%5f", v2a)     << "+-" << Form("%5f", v2aerr) << endl;
    fout << "VDRIEA              | " << Form("%5f", VDRIEA)        << endl;
    fout << "VVIERA              | " << Form("%5f", VVIERA)        << " | " << Form("%5f", v4a)     << "+-" << Form("%5f", v4aerr) << endl;
    fout << "VDRIET*VDRIEA       | " << Form("%5f", VDRIET*VDRIEA) << " | " << Form("%5f", v3)      << "+-" << Form("%5f", v3err)  << endl;
    fout << "VVIERT*VVIERA       | " << Form("%5f", VVIERT*VVIERA) << " | " << Form("%5f", v4a*v4t) << "+-" << Form("%5f", TMath::Sqrt(v4aerr*v4aerr*v4t*v4t + v4terr*v4terr*v4a*v4a)) << endl;
    fout << "VTWEEV              | " << VTWEEV << endl;
    fout << "NV                  | " << NV << endl;
    fout << "VBINS               | " << VBINS << endl;
    fout << "BCKGR               |          | "                             << Form("%5f", bgr)     << "+-" << Form("%5f", bgrerr) << endl;
    fout << "R2                  |          | "                             << Form("%5f", r2)      << "+-" << "?" << endl;
    fout << "R4                  |          | "                             << Form("%5f", r4)      << "+-" << "?" << endl;
    fout << "R8*S84              |          | "                             << Form("%5f", r8s84)   << "+-" << "?" << endl;
    fout << "S4                  |          | "                             << Form("%5f", s4)      << "+-" << "?" << endl;
    fout << "CHI2                |          | "                             << Form("%5f", finalchi2)              << endl;
    fout.close();
  }
};
Fitter* Fitter::instance = 0;
Double_t Fitter::r2 = 0;
Double_t Fitter::r4 = 0;
Double_t Fitter::r8s84 = 0;
Double_t Fitter::s4 = 0;
Int_t Fitter::binsX = 0;
Double_t Fitter::a = 0.8269933431; // = sqrt(27)/2pi
TH1D** Fitter::dhd = (TH1D**) malloc(sizeof(TH1D*)*3);
Fitter* Fitter::getInstance() {
  if (instance == 0)
    instance = new Fitter();
  return instance;
}


// Run the simulation and fit
DiHadronSum* runPsi34(Int_t nritt = 100) {
  cout << "creating..." << endl;
  DiHadronSum* disum = new DiHadronSum();
  for(Int_t evnr = 0; evnr < nritt; evnr++) {
    DiHadronEvent* event = sim_event(kTRUE);
    disum->AddEvent(event);
    delete event;
  }
  return disum;
}


TH2D* normalizePerYbin(TH2D* hist) {
  TH2D* copy = new TH2D(*hist);
  TH1D* yproj = copy->ProjectionY();
  for(Int_t ybin = 1; ybin < copy->GetNbinsY()+1; ybin++) {
    for(Int_t xbin = 1; xbin < copy->GetNbinsY()+1; xbin++) {
      if(yproj->GetBinContent(ybin))
        copy->SetBinContent(xbin,ybin,copy->GetBinContent(xbin, ybin) / yproj->GetBinContent(ybin));
    }
  }
  return copy;
}


TH1D* diagonalProjection(TH2D* hist) {
  TH1D* dproj = hist->ProjectionX();
  for(Int_t dbin = 1; dbin < dproj->GetNbinsX()+1; dbin++) {
    dproj->SetBinContent(dbin, 0);
    for(Int_t ybin = 1; ybin < hist->GetNbinsY()+1; ybin++)
      for(Int_t xbin = 1; xbin < hist->GetNbinsY()+1; xbin++)
        if((xbin + ybin)%hist->GetNbinsY() + 1 == dbin)
          dproj->SetBinContent(dbin, dproj->GetBinContent(dbin) + hist->GetBinContent(xbin, ybin));
  }
  return dproj;
}


void run(const char* mapnaam, Int_t submap, Int_t nritt = 100, Bool_t drawWithSmearing = kFALSE) {
  cout << "creating..." << endl;
  DiHadronSum* disum = new DiHadronSum();
  gSystem->Exec(Form("mkdir %s", mapnaam));
  gSystem->Exec(Form("mkdir %s/%d", mapnaam, submap));
  for(Int_t evnr = 0; evnr < nritt; evnr++) {
    if(evnr % (nritt/1000) == 0) { cout << 1000*evnr/nritt << " promille" << endl; }
    DiHadronEvent* event = sim_event();
    disum->AddEvent(event);
    delete event;
  }
/*  cout << "outputting... vspsi" << endl;
  TFile* output = TFile::Open(Form("%s/%d/distvspsi.root",mapnaam, submap), "RECREATE");
  for(Int_t plane = 0; plane < 3; plane++) {
    output->cd();
    TH2D* plot = disum->GetDistvspsi4(plane);
    TH2D* norm = normalizePerYbin(plot);
    TH1D* diagproj = diagonalProjection(plot);
    plot->Write(Form("plane%d", plane));
    norm->Write(Form("plane%d_norm", plane));
    diagproj->Write(Form("plane%d_proj", plane));
  }
  output->Close();*/
  cout << "fitting..." << endl;
  Fitter* fitter = Fitter::getInstance();
  fitter->SetData(disum);
  delete disum;
  fitter->Fit();
  cout << "outputting..." << endl;
  if(drawWithSmearing) fitter->SetDrawWithSmearingR();
  fitter->Save(Form("%s/%d", mapnaam, submap));
}


void runErrorAnalysis(const char* mapnaam, Int_t submap, Int_t nrouter, Int_t nrinner, Bool_t setRandomToSubmap=kFALSE) {
  if(setRandomToSubmap) {
    cout << "set random seed to " << submap << endl;
    RAND->SetSeed(submap);
  }
  cout << "creating and fitting..." << endl;
  gSystem->Exec(Form("mkdir %s", mapnaam));
  gSystem->Exec(Form("mkdir %s/%d", mapnaam, submap));
//////
  TNtuple* tree = new TNtuple("uitkomsten", "bootstrap_results", "v2a:v2aerr:v2t:v2terr:v4a:v4aerr:v3:v3err:chi2:chi2in:chi2mid:chi2out");
/*  TH1D* v2a = new TH1D("v2a-dist", "distribution of fitted v2a", 100, 0.08, 0.32);
  TH1D* v2aerr = new TH1D("v2aerr-dist", "distribution of th error of the fitted v2a", 100, 0.0, 0.05);
  TH1D* v2t = new TH1D("v2t-dist", "distribution of fitted v2t", 100, 0.0, 0.37);
  TH1D* v2terr = new TH1D("v2terr-dist", "distribution of th error of the fitted v2t", 100, 0.0, 0.05);
  TH1D* v4a = new TH1D("v4a-dist", "distribution of fitted v4a", 100, 0.0, 0.20);
  TH1D* v4aerr = new TH1D("v4aerr-dist", "distribution of th error of the fitted v4a", 100, 0.0, 0.03);
  TH1D* v3 = new TH1D("v3-dist", "distribution of fitted v3a", 100, 0.0, 0.05);
  TH1D* v3err = new TH1D("v3err-dist", "distribution of th error of the fitted v3a", 100, 0.0, 0.03);
  TH1D* chi2 = new TH1D("chi2-dist", "distribution of the chi2 of the fit", 100, 0.0, 400.0);
  TH1D* chi2in = new TH1D("chi2-dist-inplane", "distribution of the chi2 in-plane of the fit", 100, 30.0, 130.0);
  TH1D* chi2mid = new TH1D("chi2-dist-midplane", "distribution of the chi2 in-plane of the fit", 100, 30.0, 130.0);
  TH1D* chi2out = new TH1D("chi2-dist-outplane", "distribution of the chi2 in-plane of the fit", 100, 30.0, 130.0);*/
  for(Int_t fitnr = 0; fitnr < nrouter; fitnr++) {
    DiHadronSum* disum = new DiHadronSum();
    for(Int_t evnr = 0; evnr < nrinner; evnr++) {
      if(evnr % (nrinner/10) == 0) { cout << 100*evnr/nrinner << " procent" << endl; }
      DiHadronEvent* event = sim_event();
      disum->AddEvent(event);
      delete event;
    }
    Fitter* fitter = Fitter::getInstance();
    fitter->SetData(disum);
    delete disum;
    fitter->Fit();
/*    v2a->Fill(fitter->getv2a());
    v2aerr->Fill(fitter->getv2aerr());
    v2t->Fill(fitter->getv2t());
    v2terr->Fill(fitter->getv2terr());
    v4a->Fill(fitter->getv4a());
    v4aerr->Fill(fitter->getv4aerr());
    v3->Fill(fitter->getv3());
    v3err->Fill(fitter->getv3err());
    chi2->Fill(fitter->getchi2());
    chi2in->Fill(fitter->computeChi2(0));
    chi2mid->Fill(fitter->computeChi2(1));
    chi2out->Fill(fitter->computeChi2(2));*/
    tree->Fill(fitter->getv2a(), fitter->getv2aerr(), fitter->getv2t(), fitter->getv2terr(), fitter->getv4a(), fitter->getv4aerr(), fitter->getv3(), fitter->getv3err(), fitter->getchi2(), fitter->computeChi2(0), fitter->computeChi2(1), fitter->computeChi2(2));
  }
  TCanvas *c1 = new TCanvas("c1","c1:v2a",800,600);
  tree->Draw("v2a");
  TH1D* v2a = (TH1D*) tree->GetHistogram()->Clone();
//  v2a->Draw();
  c1->SaveAs(Form("%s/%d/v2a.pdf", mapnaam, submap));
  delete c1;
  TCanvas *c2 = new TCanvas("c2","c2:v2aerr",800,600);
  tree->Draw("v2aerr");
//  v2aerr->Draw();
  c2->SaveAs(Form("%s/%d/v2aerr.pdf", mapnaam, submap));
  delete c2;
  TCanvas *c3 = new TCanvas("c3","c3:v2t",800,600);
  tree->Draw("v2t");
  TH1D* v2t = (TH1D*) tree->GetHistogram()->Clone();
//  v2t->Draw();
  c3->SaveAs(Form("%s/%d/v2t.pdf", mapnaam, submap));
  delete c3;
  TCanvas *c4 = new TCanvas("c4","c4:v2terr",800,600);
  tree->Draw("v2terr");
//  v2terr->Draw();
  c4->SaveAs(Form("%s/%d/v2terr.pdf", mapnaam, submap));
  delete c4;
  TCanvas *c5 = new TCanvas("c5","c5:v4a",800,600);
  tree->Draw("v4a");
  TH1D* v4a = (TH1D*) tree->GetHistogram()->Clone();
//  v4a->Draw();
  c5->SaveAs(Form("%s/%d/v4a.pdf", mapnaam, submap));
  delete c5;
  TCanvas *c6 = new TCanvas("c6","c6:v4aerr",800,600);
  tree->Draw("v4aerr");
//  v4aerr->Draw();
  c6->SaveAs(Form("%s/%d/v4aerr.pdf", mapnaam, submap));
  delete c6;
  TCanvas *c7 = new TCanvas("c7","c7:v3",800,600);
  tree->Draw("v3");
  TH1D* v3 = (TH1D*) tree->GetHistogram()->Clone();
//  v3->Draw();
  c7->SaveAs(Form("%s/%d/v3.pdf", mapnaam, submap));
  delete c7;
  TCanvas *c8 = new TCanvas("c8","c8:v3err",800,600);
  tree->Draw("v3err");
//  v3err->Draw();
  c8->SaveAs(Form("%s/%d/v3err.pdf", mapnaam, submap));
  delete c8;
  TCanvas *c9 = new TCanvas("c9","c9:chi2",800,600);
  tree->Draw("chi2");
//  chi2->Draw();
  c9->SaveAs(Form("%s/%d/chi2.pdf", mapnaam, submap));
  delete c9;
  TCanvas *ctemp = new TCanvas("ctemp","c_temporary_plots",800,600);
  tree->Draw("chi2in");
  TH1D* chi2in = (TH1D*) tree->GetHistogram()->Clone();
  tree->Draw("chi2mid");
  TH1D* chi2mid = (TH1D*) tree->GetHistogram()->Clone();
  tree->Draw("chi2out");
  TH1D* chi2out = (TH1D*) tree->GetHistogram()->Clone();
  delete ctemp;
  TFile* output = TFile::Open(Form("%s/%d/fitresults_tree_form.root", mapnaam, submap), "RECREATE");
  output->cd();
  tree->Write("fitresults");
  output->Close();
/*  TCanvas *c10 = new TCanvas("c10","c10:chi2-perplane",800,600);
  chi2in->SetLineColor(2);
  chi2in->Draw();
  chi2mid->SetLineColor(3);
  chi2mid->Draw("same");
  chi2out->SetLineColor(4);
  chi2out->Draw("same");
  TLegend* leg = new TLegend(0.6,0.7,0.98,0.9);
  leg->AddEntry(chi2in,Form("in, mean=%.4f sig=%.5f", chi2in->GetMean(), chi2in->GetStdDev()));
  leg->AddEntry(chi2mid,Form("mid, mean=%.4f sig=%.5f", chi2mid->GetMean(), chi2mid->GetStdDev()));
  leg->AddEntry(chi2out,Form("out, mean=%.4f sig=%.5f", chi2out->GetMean(), chi2out->GetStdDev()));
  leg->Draw();
  c10->SaveAs(Form("%s/%d/chi2_perplane.pdf", mapnaam, submap));
  delete c10;
  cout << "hist sig v2a: " << v2a->GetStdDev() << endl;
  cout << "hist sig v2t: " << v2t->GetStdDev() << endl;
  cout << "hist sig v4a: " << v4a->GetStdDev() << endl;
  cout << "hist sig v3: "  << v3->GetStdDev()  << endl;*/
}


// A function to create a test TF2 that should look like the delta phi - delta psi_4 dist
void drawDphiDpsi4_theory() { 
  TFile* output = TFile::Open("other_plots/2dimtheory.root", "RECREATE");
  Double_t clist[3] = {TMath::Pi()/6., TMath::Pi()/12., TMath::Pi()/6.};
  Double_t signlist[3] = {1., -1., 1.};
  for(Int_t i = 0; i < 3; i++) {
    TF2* func = new TF2("2dim_theory", "1+[4]*2*[1]*[0]*cos(4*y)+[4]*2*[2]*[0]*cos(4*(y-x))-[4]*[1]*[2]*[0]*cos(8*(y-(x/2.)))+8*[3]*cos(4*x)*[1]*[2]",0,2*TMath::Pi(),0,2*TMath::Pi());
    TH2D* hist = new TH2D(Form("2dim_theory_%d", i), Form("2D_theory_%d", i), 64, 0, 2*TMath::Pi(), 64, 0, 2*TMath::Pi());
    func->SetParameter(0,TMath::Sqrt(3)/2);
    func->SetParameter(1,VVIERT);
    func->SetParameter(2,VVIERA);
    func->SetParameter(3,clist[i]);
    func->SetParameter(4,signlist[i]);
    output->cd();
    func->SetNpy(100);
    func->SetNpx(100);
    Double_t x; Double_t y;
    for(Int_t n = 0; n < 100000000; n++) {
      func->GetRandom2(x, y);
      hist->Fill(x, y);
    }
    hist->Write(Form("f_2dim_theory_%d", i));
    TH2D* hist_norm = normalizePerYbin(hist);
    hist_norm->Write(Form("f_2dim_theory_norm_%d", i));
  }
}
