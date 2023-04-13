#include <TCanvas.h>
#include <TMath.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
using namespace std;
Double_t pmt_response(Double_t *x, Double_t *par);
Double_t pmt_response_S0real(Double_t *x, Double_t *par);
Double_t pmt_response_Sn1real(Double_t *x, Double_t *par);
Double_t pmt_response_Snreal(Double_t *x, Double_t *par);

void PMT_analysis(){
    //Load data 
    TH1F* h1 = new TH1F("hist", "PhSpec; Channel; Events", 1000,0, 1000);
    ifstream inp; 
    // Fill histogram
    double x;
    inp.open("PMT_data/data/ch00/run_18390_ch00_low.dat");
    for (int i =1; i <= 4097; i++){
        inp.ignore(10, '\n');
        inp >> x;
        h1->SetBinContent(i,x);
    }

    //Fit data
    //First gaussian to pedestal first
    int pedpos = 270;
    TF1* f2 = new TF1("Pedestal_fit", "gaus", pedpos-10, pedpos+3);
    f2->SetParNames("Amplitude", "Mean", "Sigma");
    f2->SetParameter(0, h1->GetBinContent(h1->FindBin(pedpos)));
    f2->SetParameter(1, pedpos);
    f2->SetParameter(2,1);
    TFitResultPtr r1 = h1->Fit("Pedestal_fit", "WORS");
    double pedsigma = f2->GetParameter(2);
    double pedsigma_error = f2->GetParError(2);
    double pedpeak = f2->GetParameter(1);
    double pedpeak_error = f2->GetParError(1);


    int pmin = pedpos-10;  //pedpos-10
    int pmax = 1000;
    int npar = 8;

    TF1 *f1 = new TF1("fitfunc", pmt_response, pmin, pmax, npar);
    f1->SetParNames("N", "w", "#alpha", "#mu", "Q_{0}", "#sigma_{0}", "Q_{1}", "#sigma_{1}");
    f1->SetParameter(0, 500*h1->GetEntries()); //Nev
    //Set parameters limits
    f1->SetParLimits(1, 0, 1); //w
    f1->SetParLimits(2, 0, 1); //alpha
    f1->SetParLimits(3, 0, 0.7); //mu
    f1->SetParLimits(4, 0, 1000); //Q0
    f1->SetParLimits(5, 0, 100); //sigma0
    f1->SetParLimits(6, 0, 300); //Q1
    f1->SetParLimits(7, 0, 200); //sigma1
    //Set Parameter
    f1->SetParameter(1, 0.0); //w
    f1->SetParameter(2, 0.0); //alpha
    f1->SetParameter(3, 0.0); //mu
    f1->SetParameter(4, pedpeak); // Q0 (pedestal peak)
    f1->SetParameter(5, pedsigma); // sigma0 (pedestal sigma)
    f1->SetParameter(6, 30); // Q1
    f1->SetParameter(7, 12); //sigma1

    //Canvas
    TCanvas* c = new TCanvas("c0", "Npe fit", 800,600);
    c->SetLogy();
    h1->Draw("hist");
    //h1->GetXaxis()->SetRangeUser(200,550);
    TFitResultPtr r0 = h1->Fit("fitfunc", "LEORS");
    //f2->Draw("LFsame");
    //f1->Draw("LFsame");
  //Save fit info into a txt file
   fstream fs;
        fs.open("fit_info_ET795_1300V_200ns_18450.txt", ios::out);
        string lin;
// Make a backup of stream buffer
        streambuf* sb_cout = cout.rdbuf();
        streambuf* sb_cin = cin.rdbuf();
// Get the file stream buffer
        streambuf* sb_file = fs.rdbuf();
// Now cout will point to file
        cout.rdbuf(sb_file);
// everything printed after this line will go to the fs
        r0->Print("V");
        cout << " Q_0 " << r0->Parameter(4) << " error " << r0->Error(4) << endl;
        cout << " sigma_0 " << r0->Parameter(5) << " error " << r0->Error(5) << endl;
        cout << " w " << r0->Parameter(1) << " error " << r0->Error(1) << endl;
        cout << " alpha " << r0->Parameter(2) << " error " << r0->Error(2) << endl;
        cout << " mu " << r0->Parameter(3) << " error " << r0->Error(3) << endl;
        cout << " Q_1 " << r0->Parameter(6) << " error " << r0->Error(6) << endl;
        cout << " sigma_1 " << r0->Parameter(7) << " error " << r0->Error(7) << endl;
        cout << " Chisquare test " << h1->Chisquare(f1, "R") << endl;

// get the previous buffer from backup
        cout.rdbuf(sb_cout);
// everything printed after this line will go to whereever it went before
        fs.close();



}
 Double_t pmt_response(Double_t *x, Double_t *par){
        double par2[8];
        for(int i=1; i<8;i++) par2[i-1] = par[i];

        double val = 0;
        for(int n=0; n<10; n++){
            par2[7] = n;
            val += pmt_response_Snreal(x, par2);
        }
        return par[0]*val;
    }

    Double_t pmt_response_S0real(Double_t *xx, Double_t *par){
        double x =xx[0];

        double w = par[0];
        double alpha = par[1];
        double mu = par[2];
        double Q0 = par[3];
        double sigma0 = par[4];

        double pn = TMath::Power(mu, 0)* TMath::Exp(-mu)/TMath::Factorial(0);
        double t1 = -TMath::Exp(-TMath::Power(Q0 - x, 2)/(2. * sigma0*sigma0)); //Define sigma0
        double t2 = (w - 1.)/(TMath::Sqrt(2.*TMath::Pi())*sigma0);
        double t3 = x - Q0 >=0 ? alpha * TMath::Exp((Q0 - x) * alpha) : 0;

        return pn * (t1*t2 + t3);
    }
    Double_t pmt_response_Sn1real(Double_t *xx, Double_t *par){
        double x = xx[0];

        double w = par[0];
        double alpha = par[1];
        double mu = par[2];
        double Q0 = par[3];
        double sigma0 = par[4];
        double Q1 = par[5];
        double sigma1 = par[6];
        int n = (int)(par[7]+0.1);

        double pn = TMath::Power(mu, n)*TMath::Exp(-mu)/TMath::Factorial(n);
        double t1 = -TMath::Exp(-TMath::Power(Q0 + n* Q1 - x,2)/(2.  * (sigma0*sigma0 + n * sigma1*sigma1)));
        double t2 = (w - 1.)/(TMath::Sqrt(2.*TMath::Pi())*TMath::Sqrt(sigma0*sigma0 + n * sigma1*sigma1));
        double t3 = 0.5 * w * alpha * TMath::Exp(0.5 * alpha * (2. * Q0 + 2. * n * Q1 - 2. * x + n * alpha * sigma1*sigma1));
        double t4 = TMath::Erfc((Q0 - x + n * (Q1 + alpha * sigma1*sigma1))/(TMath::Sqrt(2. * n) * sigma1) );
        return pn * (t1*t2+t3*t4);
    }

    Double_t pmt_response_Snreal(Double_t *x, Double_t *par){
        int n = (int)(par[7]+0.1);

        double val = 0;
        if(n==0) val = pmt_response_S0real(x, par);
        else val = pmt_response_Sn1real(x, par);

        return val;
    }