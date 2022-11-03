//////////////////////// EBL and TAU fits //////////////////////////////////////
//
// Authors:  Antonio Stamerra                       ////////////////////
//
//
//  Fits on EBL and Tau curves
//
// Provide input file of Tau (matrix of Z(first line) and E(first coloumn))
////////////////////////////////////////////////////////////////////////////

//############################
//  CHECK here for python module for EBL:
//   https://github.com/me-manu/ebltable
//###############################

//----- COMPILATION WITH following commands:
//
// gcc -O -Wall -fPIC -pthread -I$ROOTSYS/include -c ./SEDestimation.cc
// g++ -O SEDestimation.o -L$ROOTSYS/lib -lCore -lCint -lHist -lGraf -lGPad -lRint -lPostscript -lMatrix -lPhysics -pthread -lm -ldl -rdynamc -o SEDestimator.exe
//----------------------
//
#include <TROOT.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TPad.h>
#include <TNamed.h>
#include <TMultiGraph.h>
#include <TLatex.h>
#include <TClass.h>
#include <TGClient.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>
#include <math.h>
#include <TMath.h>
#include <stdio.h>
#include <stdlib.h>
#include <TGaxis.h>

using namespace std;

TString version="Version April-2012";

Bool_t DEBUG = kTRUE;
//Bool_t DEBUG = kFALSE;

//user function

Double_t fun1(Double_t *x,Double_t *par) {
    Double_t arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    Double_t fitval = par[0] + par[1]*x[0]+par[2]*pow(x[0],2)+par[3]*pow(x[0],3);
    return fitval;
}

Double_t fun2(Double_t *x,Double_t *par) {
    Double_t arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    Double_t fitval = par[0] + par[1]*x[0]+par[2]*pow(x[0],2)+par[3]*pow(x[0],3);
    return fitval;
}

Double_t fun3(Double_t *x,Double_t *par) {
    Double_t arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    Double_t fitval = par[0] + par[1]*x[0]+par[2]*pow(x[0],2)+par[3]*pow(x[0],3)+par[4]*pow(x[0],4)+par[5]*pow(x[0],5)+par[6]*pow(x[0],6)+par[7]*pow(x[0],7)+par[8]*pow(x[0],8);
    return fitval;
}

//TString filename="";


///---- function to compute EBL absorption and interpolate in energy and redshift
Double_t ebl(Double_t z, Double_t e, Int_t EBLm)
{
    Double_t highene=0, lowene=0;
    Double_t ebl11=0,ebl12=0,ebl21=0,ebl22=0;
    TString filename="";
    switch (EBLm) {
        case 1:
            filename += "./tau_franceschini08.out";
            break;
        case 2:
            filename += "./tau_gilmore09.out";
            break;
        case 3:
            filename += "./tau_kneiske04bf.out";
            break;
        case 4:
            filename += "./tau_kneiske04lowSFR.out";
            break;
        case 5:
            filename += "./tau_franceschini_Highz.out";
            break;
        case 6:
            filename += "./tau_dominguez10.out";
            break;
    }
    TString row;
    Double_t highz=0, lowz=0;
    Int_t col=0,len;
    e/=1000; //convertion to TeV; table energy in TeV but energy input in GeV
    
    //-----------------------
    // opening the file with absorptions values
    // File format:
    // first 2 lines are skipped (comments)
    //  zvalues ->  columns
    //  energy values: rows
    //
    ifstream fil (filename, ifstream::in);
    //TFile *fil = new TFile(filename);
    if(!fil) return 0.;
    
    // skip first 2 lines
    row.ReadLine(fil);
    row.ReadLine(fil);
    
    // read first line with readshift values
    row.ReadLine(fil);
    while (highz < z)
    {
        lowz=highz;
        sscanf(row.Data(),"%lf %n",&highz, &len);
        row.Remove(0,len);
        col++;
    }
    
    while (1)
    {
        // Read columns with EBL values and store the values corresponding
        // to the read z  (variables ebl11, ebl12, ebl21, ebl22)
        //       |  z1      z2         z1 < z < z2
        //      ---------------      ene1 < e < ene2
        //  ene1 | ebl11  ebl12
        //  ene2 | ebl21  ebl22
        //-----------------------//
        
        // read second line and read energy value
        // read next row, up to energy value exceeding ene2 (highene)
        row.ReadLine(fil);
        sscanf(row.Data(),"%lf %n",&highene, &len);
        row.Remove(0,len);
        if (highene < e){
            lowene=highene;
            // read colums up to column corresponding to z2 (highz)
            for (Int_t i=1;i<=col;i++)
            {
                ebl11=ebl12;
                sscanf(row.Data(),"%lf %n",&ebl12, &len);
                row.Remove(0,len);
            }
        }
        else{
            // read colums up to column corresponding to z2 (highz)
            for (Int_t i=1;i<=col;i++)
            {
                ebl21=ebl22;
                sscanf(row.Data(),"%lf %n",&ebl22, &len);
                row.Remove(0,len);
            }
            break; // table_energy > input energy
        }
    } // while(1) loop
    
    
    Double_t zfrac=(z-lowz)/(highz-lowz);
    Double_t efrac=(e-lowene)/(highene-lowene);
    Double_t zinterp1=ebl11;
    Double_t zinterp2=ebl22;
    zinterp1 = ebl11+zfrac*(ebl12-ebl11);
    zinterp2 = ebl21+zfrac*(ebl22-ebl21);
    
    Double_t ebl_interp = zinterp1+efrac*(zinterp2-zinterp1);
    
    /*
     cout << "\t  " <<  lowz << "\t   " << highz << endl;
     cout << lowene << " " << ebl11 << "\t "<< ebl12 << endl;
     cout << highene << " " << ebl21 << "\t "<< ebl22 << endl;
     cout << "             \t       z:" << z<<",E(GeV):"<<e*1000 <<" -> tau=" << ebl_interp <<endl;
     */
    
    return ebl_interp;
}


//int main(int argc, char **argv)
// void SEDestimation_lin(Double_t flux=4.7e-8, Double_t slope=2.2, Double_t enebreak=200, Double_t slope2=3.5, Double_t z=0.1, Int_t EBLmodel=1, TString source="...", Float_t timeh=50, Int_t ismidzd =0, Int_t funlaw=0)
void SEDestimation_lin(Double_t z=0.1, Int_t EBLmodel=1)
{
    //Double_t flux=4.7e-8.,Double_t slope=2.2, Double_t enebreak=200, Double_t slope2=3.5, Double_t z=0.1, Int_t EBLmodel=1, TString source="...
    
    int ikl=0;
    
    gStyle->SetFillColor(0);
    gStyle->SetOptStat(0);
    
    //------ VARIABLES To BE ADJUSTED -------////
    
    // FLUX IN UNITS PH/cm2/s
    // energy in GeV
    
    // energy limits for computation and plot
    //    const Double_t startEne=0.1; // GeV
    const Double_t startEne=0.1; // GeV
    const Double_t startEneLog=TMath::Log10(startEne);
    const Double_t stopEne=20000; // GeV
    const Double_t stopEneLog=TMath::Log10(stopEne);
    const Int_t npoints=60; // number of points for sampling
    const Double_t LogStep=TMath::Log10(stopEne/startEne)/npoints; // Log(E[GeV])
    
    // define limits for histogram
    const Double_t GevToErg = 1.6022e-3;
    Double_t LogFluxMinGeV = -10.5;
    Double_t LogFluxMaxGeV = -6.0;
    Double_t LogFluxMinErg = LogFluxMinGeV + log10(GevToErg);
    Double_t LogFluxMaxErg = LogFluxMaxGeV + log10(GevToErg);
    
    
    // variables
    Double_t tau;  // opacity
    Double_t a;  // normalizazion for Fermi flux (ph/cm2/s)
    Double_t absorption[npoints];
    Double_t SpectrumEbl[npoints], logSpect[npoints];
    Double_t SpectrumNoEbl[npoints], logSpectNoEbl[npoints];
    Double_t enebin, logE[npoints], EneLin[npoints], FluxEBL[npoints];
    
   /*
    // string variables for output and plot legend
    TString FermiFlux = Form("Fermi flux E>100MeV: %.2e ph/cm^2/s",flux);
    //cout << FermiFlux << endl;
    TString phidx = Form("HE photon index: %.2f",slope);
    TString eb = "  E_{break}:";
    eb += enebreak;
    eb+= " GeV";
    TString ph2idx = Form("ph.idx2: %.2f",slope2);
    //cout << phidx << eb <<  ph2idx << endl;
    TString redsh = Form("z = %.3f",z);
    //cout << redsh << endl;
    
    TString outfile="./SEDplots/SED_";
    outfile+=source;
    TString outpng=outfile+".png";
    TString outgif=outfile+".gif";
    TString outpdf=outfile+".pdf";
    TString outroot=outfile+".root";
    TString outeps=outfile+".eps";
    */
    
    
    //------------- Plot tau and make a fit ------------//
    
    
    TGraph * grTau = new TGraph(npoints);
    TGraph * grTran = new TGraph(npoints);
    
    cout << "Redshift: " << z << "  Model: ("<< EBLmodel << ")"<< endl;
    
    for(Int_t i=0;i<npoints;i++)
    {
        logE[i] = log10(startEne) + i*LogStep;
        EneLin[i] = pow(10,logE[i]);
        tau=0;
        enebin = pow(10,logE[i]);
        
        absorption[i] = 1;
        if (z!=0) {
            tau=ebl(z,enebin, EBLmodel); // call function to read values from file
            absorption[i]=exp(-1*tau);
        }
        grTau->SetPoint(i,enebin,tau);
        grTran->SetPoint(i,logE[i],absorption[i]);
        
        //output on screen
        cout << "E(GeV): " << EneLin[i] << " tau:"  << tau << "  Absorption:" <<  absorption[i] << endl;

    }
    
    
    //---------- fill histograms
    TCanvas *cTau = new TCanvas("cTau","Tau and Abs");
    cTau->SetFillColor(0);
    cTau->SetBorderMode(0);
    cTau->SetBorderSize(2);
    cTau->SetFrameBorderMode(1);
    
    
    grTran->Draw("LA");

    TF1 *f1 = new TF1("f1",fun1,1.8,3.,5);
    
    TF1 *f2 = new TF1("f2",fun2,2,3.5,5);

    TF1 *f3 = new TF1("f3",fun3,1.8,3.5,8);
    
    f3->SetParameters(1,-1,1,1,1,1,1,1,1,1);
//    f->Draw("same");
//    grTran->Fit(f1,"R");
  //  grTran->Fit(f2,"R");
    grTran->Fit(f3,"R");

    f1->Draw("same");

    
    return;
   
}



