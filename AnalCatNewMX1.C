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


int NbinX=40;
int NdivY=3;
double MinB=6;
// from AnalCatNewMVA2A.C
const Double_t boundary_MVA_2017[4]={0.248, 0.450, 0.728, 1.0};
/// start position for MX
const Double_t boundary_MX_2017[5]={250.,380., 510., 610., 800};

double dwid=(boundary_MX_2017[4]-boundary_MX_2017[0])/NbinX;

double *bkgs = new double[NbinX];
double *sigs = new double[NbinX];

// Punzi function
double Func(double S, double B) {
//      Punzi
      return S/(1.+sqrt(B));
}
TGraph *g1,*g2,*g3,*g4,*g5,*g6;
double splgraph1(double *xx, double *)
{
      return g1->Eval(xx[0]);
}
double splgraph2(double *xx, double *)
{
      return g2->Eval(xx[0]);
}
double splgraph3(double *xx, double *)
{
      return g3->Eval(xx[0]);
}
double splgraph4(double *xx, double *)
{
      return g4->Eval(xx[0]);
}
double splgraph5(double *xx, double *)
{
      return g5->Eval(xx[0]);
}
double splgraph6(double *xx, double *)
{
      return g6->Eval(xx[0]);
}

double OptFunc (const double *xx) {
//   const double x[5];
   double x[5];
   x[1] = xx[0];
   x[2] = xx[1];
   x[3] = xx[2];
   x[0]=x[4]=0;
//   double x = xx;
//   cout<<"dX_mx"<<x[1]<<" "<<x[2]<<" "<<x[3]<<std::endl;
   double Svs1[4],Svs2[4],Svs3[4];
   double Svb1[4],Svb2[4],Svb3[4];
   TF1 *f1 = new TF1("f1",splgraph1,boundary_MX_2017[0],boundary_MX_2017[4]);
   TF1 *f2 = new TF1("f2",splgraph2,boundary_MX_2017[0],boundary_MX_2017[4]);
   TF1 *f3 = new TF1("f3",splgraph3,boundary_MX_2017[0],boundary_MX_2017[4]);
   TF1 *f4 = new TF1("f4",splgraph4,boundary_MX_2017[0],boundary_MX_2017[4]);
   TF1 *f5 = new TF1("f5",splgraph5,boundary_MX_2017[0],boundary_MX_2017[4]);
   TF1 *f6 = new TF1("f6",splgraph6,boundary_MX_2017[0],boundary_MX_2017[4]);
   for (int Ipr=0; Ipr<4;Ipr++){
      Svs1[Ipr]=f1->Integral(boundary_MX_2017[0+Ipr]+x[Ipr],boundary_MX_2017[1+Ipr]+x[Ipr+1],1.e-4)/dwid;
      Svb1[Ipr]=f2->Integral(boundary_MX_2017[0+Ipr]+x[Ipr],boundary_MX_2017[1+Ipr]+x[Ipr+1],1.e-4)/dwid;
      Svs2[Ipr]=f3->Integral(boundary_MX_2017[0+Ipr]+x[Ipr],boundary_MX_2017[1+Ipr]+x[Ipr+1],1.e-4)/dwid;
      Svb2[Ipr]=f4->Integral(boundary_MX_2017[0+Ipr]+x[Ipr],boundary_MX_2017[1+Ipr]+x[Ipr+1],1.e-4)/dwid;
      Svs3[Ipr]=f5->Integral(boundary_MX_2017[0+Ipr]+x[Ipr],boundary_MX_2017[1+Ipr]+x[Ipr+1],1.e-4)/dwid;
      Svb3[Ipr]=f6->Integral(boundary_MX_2017[0+Ipr]+x[Ipr],boundary_MX_2017[1+Ipr]+x[Ipr+1],1.e-4)/dwid;
   }
   for (int i=0;i<4;i++) {
      std::cout<<i<<"=> "<<Svs1[i]<<" / " <<Svb1[i];
      std::cout<<" "<<Svs2[i]<<" / " <<Svb2[i];
      std::cout<<" "<<Svs3[i]<<" / " <<Svb3[i]<<std::endl;
   }
   double Sign1[4],Sign2[4],Sign3[4];
   double SumOfSign=0;
   for (int i=0; i<4;i++) {
      Sign1[i]=Func(Svs1[i],Svb1[i]);
      Sign2[i]=Func(Svs2[i],Svb2[i]);
      Sign3[i]=Func(Svs3[i],Svb3[i]);
          
     SumOfSign+=pow(Sign1[i],2)+pow(Sign2[i],2)+pow(Sign3[i],2);
//     cout<<"Rs="<<Sign1[i]<<" "<<Sign2[i]<<" "<<Sign3[i]<<std::endl;
   }
   SumOfSign=sqrt(SumOfSign);
   double vpen=0.;
   for (int i=0; i<4;i++){
      vpen+=SumOfSign*0.1*max(0.,MinB-Svb1[i])*(MinB-Svb1[i]);
      vpen+=SumOfSign*0.1*max(0.,MinB-Svb2[i])*(MinB-Svb2[i]);
      vpen+=SumOfSign*0.1*max(0.,MinB-Svb3[i])*(MinB-Svb3[i]);
   }
   cout<<"penalty: "<<vpen<<std::endl;
   cout<<"Current value: "<<SumOfSign<<std::endl;
   cout<<"==================="<<endl;
   return -(SumOfSign)+vpen;
//   return -(SumOfSign);
}


void AnalCatNewMX1() {
   gROOT->SetStyle("Plain");
   gStyle->SetTimeOffset(0);
   gStyle->SetOptStat(0);
   gStyle->SetGridStyle(2);
   gStyle->SetGridWidth(2);
   gStyle->SetGridColor(kGray);
   gStyle->SetLabelSize(0.06,"xy");
   TCanvas cc("cc","cc",1200,900);
   cc.Range(0,0,1,1);

   float mgg,mjj,MX,weight,iweight,mva_result;

   TChain* tr = new TChain("bbggSelectionTree");
//   NEW2016+2017+2018
    tr->Add("/afs/cern.ch/work/t/tdimova/Training/FT_leg_merge3/output_BkG.root");
   tr->SetBranchAddress("MX",&MX);
//   tr->SetBranchAddress("CMS_hgg_mass",&mgg);
//   tr->SetBranchAddress("Mjj",&mjj);
   tr->SetBranchAddress("weight",&weight);
   tr->SetBranchAddress("MVAOutputTransformed",&mva_result);
   
   float mgg1,mjj1,MX1,weight1,mva_result1;
   
   TChain* tr1 = new TChain("bbggSelectionTree");
//from Ivan NEW 2016+2017+2018 signal
   tr1->Add("/afs/cern.ch/work/t/tdimova/Training/FT_leg_merge3/output_GluGluToHHTo2B2G_allnodes_no_unit_norm.root");
   tr1->SetBranchAddress("MX",&MX1);
//   tr1->SetBranchAddress("CMS_hgg_mass",&mgg1);
//   tr1->SetBranchAddress("Mjj",&mjj1);
   tr1->SetBranchAddress("weight",&weight1);
   tr1->SetBranchAddress("MVAOutputTransformed",&mva_result1);
   
//   weight=1.;   
   int Nev = tr->GetEntries();
   int Nev1 = tr1->GetEntries();
//
//Parameters
   int NprojY=3;
//   double lumi16_17=77.4; //35.9+41.5(fb-1)
   double lumi16_17_18=136.8; //35.9+41.5+59.4(fb-1)
   double WkoefBkg=2.9*lumi16_17_18;
   double WkoefSig=1.*lumi16_17_18;
///   
// for background
   TH1D* mva= new TH1D("mva","mva",NdivY,boundary_MVA_2017);
   TH1D* mx = new TH1D("mx","mx",NbinX,250.,800.);
   TH2D* hax = new TH2D("hax","bkg",NbinX,250.,800.,3,boundary_MVA_2017);

//for signal
   TH1D* mva1= new TH1D("mva1","mva1",NdivY,boundary_MVA_2017);
   TH1D* mx1 = new TH1D("mx1","mx1",NbinX,250.,800.);
   TH2D* hax1 = new TH2D("hax1"," ",NbinX,250.,800.,3,boundary_MVA_2017);


   cout<<"Number of bkg events: "<<Nev<<" ;number of signal "<<Nev1<<endl;
   float MX_;
//filling hists for background
   double Nevbkg=0,Nevsig=0;
   for(int i=0; i < Nev; i++){
	tr->GetEntry(i);
        MX_=MX;
//	mva->Fill(mva_result,weight);
//	mx->Fill(MX_,weight);
        hax->Fill(MX_,mva_result,weight*WkoefBkg);
           if (MX_>250 && mva_result>0.25) Nevbkg+=weight*WkoefBkg;
   }
///
//filling hists for signal
   for(int i=0; i < Nev1; i++){
	tr1->GetEntry(i);
          MX_=MX1;
//	   mva1->Fill(mva_result1,weight1);
//	   mx1->Fill(MX_,weight1);
           hax1->Fill(MX_,mva_result1,weight1*WkoefSig);
         if (MX_>250&& mva_result1>0.25) Nevsig+=weight1*WkoefSig;
   }
///
///
   cout<<"Number of weighted bkg events: "<<Nevbkg<<" ;number of weighted signal "<<Nevsig<<endl;

    cc.Divide(2,5);

    TH1D * hax_X[NdivY];
    TH1D * hax1_X[NdivY];
    for (int ip=0; ip<NdivY; ip++) {   
       hax_X[ip] = hax->ProjectionX(Form("pxB%d",ip+1),ip+1,ip+1);
       hax1_X[ip] = hax1->ProjectionX(Form("pxS%d",ip+1),ip+1,ip+1);
    }
   
     cout<<"MVA values>> ";
     for (int i=0;i<4;i++)
              cout<<boundary_MVA_2017[i]<<" ";
     cout<<endl;
     cout<<"MX starting positions>> ";
     for (int i=0;i<5;i++)
              cout<<boundary_MX_2017[i]<<" ";
     cout<<endl;
	
      ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
      min.SetMaxFunctionCalls(1000000);
      min.SetMaxIterations(1000000);
      min.SetTolerance(0.001);
      ROOT::Math::Functor f(&OptFunc,3);
      double step[3] = {0.1,0.1,0.1};
      double variable[3] = {0,0,0};
      min.SetFunction(f);
      // Set the free variables to be minimized!
      min.SetVariable(0,"Dmx1",variable[0], step[0]);
      min.SetVariable(1,"Dmx2",variable[1], step[1]);
      min.SetVariable(2,"Dmx3",variable[2], step[2]);
      min.SetVariableLimits(0,-40.,40.);
      min.SetVariableLimits(1,-40.,40.);
      min.SetVariableLimits(2,-40.,40.);
      min.SetPrintLevel(1);

     for (int ipr=0;ipr<NprojY;ipr++) {
         hax_X[ipr]->Smooth();
         hax1_X[ipr]->Smooth();
         for (int i=0; i<NbinX;i++) {
           bkgs[i] = hax_X[ipr]->GetBinContent(i+1);
           sigs[i] = hax1_X[ipr]->GetBinContent(i+1)/Nevsig; //For Punzi formulae eff-cy
         }
         double *mxi = new double[NbinX];
         for (int i=0; i<NbinX; i++)
               mxi[i]=250.+dwid*(i+0.5);
         if (ipr==0) {
             g1 = new TGraph(NbinX,mxi,sigs);
             g2 = new TGraph(NbinX,mxi,bkgs);
         } else if (ipr==1) {
             g3 = new TGraph(NbinX,mxi,sigs);
             g4 = new TGraph(NbinX,mxi,bkgs);
         } else {
             g5 = new TGraph(NbinX,mxi,sigs);
             g6 = new TGraph(NbinX,mxi,bkgs);
         }
      }
      cc.cd(4);
      g1->Draw();
      cc.cd(5);
      g2->Draw();
      cc.cd(6);
      g3->Draw();
      cc.cd(7);
      g4->Draw();
      cc.cd(8);
      g5->Draw();
      cc.cd(9);
      g6->Draw();
      min.Minimize();

      const double *xs = min.X();
      double new_b[3];
      for (int i=0;i<3;i++) {
            new_b[i]=boundary_MX_2017[1+i]+xs[i];
	   std::cout<<"fit res - "<<i<<" => "<<xs[i]<<" ;new B: "<<new_b[i]<<std::endl;
      }
      std::cout<<std::endl<<"==========="<<std::endl;


 cc.Print("MY_bkg_signal_ratio_iv_mx_work.png");

}
