#include "TROOT.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TColor.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraph.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include <iostream>
#include <vector>


int NbinX=50;
int MinB=8;
int Ncat0=0;
int Ncat1=3;
//
//const Double_t boundary_MX_2017[5]={250., 363., 490., 584., 800};
const Double_t boundary_MVA_2017[4]={0.248, 0.45, 0.728, 1.0};
double dwid=1./NbinX;

double *bkgs = new double[NbinX];
double *sigs = new double[NbinX];

// Functions for Significance
double Func(double S, double B) {
      return S/(1.+sqrt(B));
}

TGraph *g1,*g2;
double splgraph1(double *xx, double *)
{
         return g1->Eval(xx[0]);
}
double splgraph2(double *xx, double *)
{
         return g2->Eval(xx[0]);
}

double OptFuncA (const double *xx) {
   //   const double x[5];
      double x[4];
      x[0] = xx[0];
      x[1] = xx[1];
      x[2] = xx[2];
      x[3]=0;
   //   double x = xx;
//     cout<<"dX_mx"<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
      double Svs1[3],Svs2[3],Svs3[3],Svs4[3];
      double Svb1[3],Svb2[3],Svb3[3],Svb4[3];

    TF1 *f1 = new TF1("f1",splgraph1,0.,1.);
    TF1 *f2 = new TF1("f2",splgraph2,0.,1.);
    for (int Ipr=Ncat0; Ipr<Ncat1;Ipr++){
       Svs1[Ipr]=f1->Integral(boundary_MVA_2017[0+Ipr]+x[Ipr],boundary_MVA_2017[1+Ipr]+x[Ipr+1],1.e-4)/dwid;
       Svb1[Ipr]=f2->Integral(boundary_MVA_2017[0+Ipr]+x[Ipr],boundary_MVA_2017[1+Ipr]+x[Ipr+1],1.e-4)/dwid;
    }
    for (int i=Ncat0;i<Ncat1;i++) {
//       std::cout<<i<<"=> "<<Svs1[i]<<" / " <<Svb1[i]<<std::endl;
    }
    double Sign1[3],Sign2[3],Sign3[3],Sign4[3];
    double SumOfSign=0;
    for (int i=Ncat0; i<Ncat1;i++) {
       Sign1[i]=Func(Svs1[i],Svb1[i]);
       SumOfSign+=pow(Sign1[i],2);
//       cout<<"Rs="<<Sign1[i]<<std::endl;
    }
    SumOfSign=sqrt(SumOfSign);
    double vpen=0.;
    for (int i=Ncat0; i<Ncat1;i++){
       vpen+=SumOfSign*0.1*max(0.,MinB-Svb1[i])*(MinB-Svb1[i]);
    }
//    cout<<"penalty: "<<vpen<<std::endl;
//    cout<<"Sum of Rs"<<SumOfSign<<std::endl;
//    cout<<"==================="<<endl;
    return -(SumOfSign)+vpen;
}

void AnalCatNewMVA2A() {
   gROOT->SetStyle("Plain");
   gStyle->SetTimeOffset(0);
   gStyle->SetOptStat(0);
   gStyle->SetGridStyle(2);
   gStyle->SetGridWidth(2);
   gStyle->SetGridColor(kGray);
   gStyle->SetLabelSize(0.06,"xy");

   TCanvas cc("cc","cc",1200,900);
   cc.Range(0,0,1,1);

   float mgg,mjj,MX,mbbgg,weight,iweight,mva_result;
   int catID;

   TChain* tr = new TChain("bbggSelectionTree");
   tr->Add("/afs/cern.ch/work/t/tdimova/Training/FT_leg_merge3/output_BkG.root");
   tr->SetBranchAddress("MX",&MX);
//   tr->SetBranchAddress("CMS_hgg_mass",&mgg);
//   tr->SetBranchAddress("Mjj",&mjj);
   tr->SetBranchAddress("weight",&weight);
   tr->SetBranchAddress("MVAOutputTransformed",&mva_result);
   
   float mgg1,mjj1,MX1,mbbgg1,weight1,mva_result1;
   
   TChain* tr1 = new TChain("bbggSelectionTree");
   tr1->Add("/afs/cern.ch/work/t/tdimova/Training/FT_leg_merge3/output_GluGluToHHTo2B2G_allnodes_no_unit_norm.root");
   tr1->SetBranchAddress("MX",&MX1);
//   tr1->SetBranchAddress("CMS_hgg_mass",&mgg1);
//   tr1->SetBranchAddress("Mjj",&mjj1);
   tr1->SetBranchAddress("weight",&weight1);
   tr1->SetBranchAddress("MVAOutputTransformed",&mva_result1);
   
   int Nev = tr->GetEntries();
   int Nev1 = tr1->GetEntries();
//
//Parameters
//   int NdivX=80;
//   int NdivY=4;
   int NprojY=4;
//   double lumi16_17=77.4; //35.9+41.5(fb-1)
   double lumi16_17_18=136.8; //35.9+41.5+59.4(fb-1)
   double WkoefBkg=2.9*lumi16_17_18;
   double WkoefSig=1.*lumi16_17_18;
///   
// for background
   TH1D* mva= new TH1D("mva","mva",NbinX,0.,1.);
//for signal
   TH1D* mva1= new TH1D("mva1","mva1",NbinX,0.,1.);
   
   cout<<"Number of bkg events: "<<Nev<<" ;number of signal "<<Nev1<<endl;
   float MX_;
//filling hists for background
   double Nevbkg=0,Nevsig=0;
   for(int i=0; i < Nev; i++){
	   tr->GetEntry(i);
           MX_=MX;
	   mva->Fill(mva_result,weight*WkoefBkg);
           if (MX_>250&&mva_result>0.2) Nevbkg+=weight*WkoefBkg;
   }
///
//filling hists for signal
   for(int i=0; i < Nev1; i++){
	 tr1->GetEntry(i);
         MX_=MX1;
         mva1->Fill(mva_result1,weight1*WkoefSig);
         if (MX_>250&&mva_result1>0.2) Nevsig+=weight1*WkoefSig;
   }
   cout<<"N background weighted: "<<Nevbkg<<"; N signal weighted: "<<Nevsig<<endl;
///
///
    cc.Divide(2,2);
    cc.cd(1);
    cc.cd(1)->SetLogy();
    mva->Draw();
    cc.cd(2);
    mva1->Draw();

    mva->Smooth();
    mva1->Smooth();

//
     cout<<"MVA starting positions>> ";
     for (int i=0;i<4;i++)
              cout<<boundary_MVA_2017[i]<<" ";
     cout<<endl;

    ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
    min.SetMaxFunctionCalls(1000000);
    min.SetMaxIterations(1000000);
    min.SetTolerance(0.0001);
    ROOT::Math::Functor f(&OptFuncA,3);
    double step[3] = {0.001,0.001,0.001};
    double variable[3] = {0,0,0};
    min.SetFunction(f);
// Set the free variables to be minimized!
    min.SetVariable(0,"Dmva1",variable[0], step[0]);
    min.SetVariable(1,"Dmva2",variable[1], step[1]);
    min.SetVariable(2,"Dmva3",variable[2], step[2]);
    min.SetVariableLimits(0,-0.3,0.3);
    min.SetVariableLimits(1,-0.3,0.3);
    min.SetVariableLimits(2,-0.3,0.3);
    min.SetPrintLevel(1);

   for (int i=0; i<NbinX;i++) {
        bkgs[i] = mva->GetBinContent(i+1);
        sigs[i] = mva1->GetBinContent(i+1)/Nevsig; //for Punzi
//       sigs[i] = mva1->GetBinContent(i+1);
   }
   double *mxi = new double[NbinX];
   for (int i=0; i<NbinX; i++)
       mxi[i]=0.+dwid*(i+0.5);
   g1 = new TGraph(NbinX,mxi,sigs);
   g2 = new TGraph(NbinX,mxi,bkgs);
   cc.cd(3);
   cc.cd(3)->SetLogy();
   g2->Draw();
   cc.cd(4);
   g1->Draw();
//minimizing
   Ncat0=2;
   Ncat1=3;
   min.FixVariable(0);
   min.FixVariable(1);
   min.Minimize();
   const double *xs = min.X();
   double new_b[3];
   for (int i=0;i<3;i++) {
      new_b[i]=boundary_MVA_2017[i]+xs[i];
      std::cout<<"fit res - "<<i<<" => "<<xs[i]<<" ;new B: "<<new_b[i]<<std::endl;
   }
   std::cout<<std::endl<<"==========="<<std::endl;
//
//   variable[2]=xs[2];
   Ncat0=1;
   Ncat1=2;
   min.ReleaseVariable(1);
   min.FixVariable(2);
   min.Minimize();
   xs = min.X();
   for (int i=0;i<3;i++) {
      new_b[i]=boundary_MVA_2017[i]+xs[i];
      std::cout<<"fit res - "<<i<<" => "<<xs[i]<<" ;new B: "<<new_b[i]<<std::endl;
   }
   std::cout<<std::endl<<"==========="<<std::endl;
//   
//   variable[1]=xs[1];
   Ncat0=0;
   Ncat1=1;
   min.ReleaseVariable(0);
   min.FixVariable(1);
   min.Minimize();
   xs = min.X();
   for (int i=0;i<3;i++) {
      new_b[i]=boundary_MVA_2017[i]+xs[i];
//      std::cout<<"fit res - "<<i<<" => "<<xs[i]<<" ;new B: "<<new_b[i]<<std::endl;
   }
   cout<<"==========="<<std::endl<<"Final borders for MVA:"<<std::endl;
   for (int i=0;i<3;i++) cout<<new_b[i]<<"   ";
   std::cout<<std::endl<<"==========="<<std::endl; 
   cc.Print("MY_bkg_signal_ratio_iv_mva_workA.png");
}
