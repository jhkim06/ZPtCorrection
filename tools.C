class Chebyshev {
public:
   Chebyshev(int n, double xmin, double xmax) :
      fA(xmin), fB(xmax),
      fT(std::vector<double>(n) )  {}

   double operator() (const double * xx, const double *p) {
      double x = (2.0 * xx[0] - fA -fB)/(fB-fA);
      int order = fT.size();
      if (order == 1) return p[0];
      if (order == 2) return p[0] + x*p[1];
      // build the polynomials
      fT[0] = 1;
      fT[1] = x;
      for (int i = 1; i< order; ++i) {
         fT[i+1] =  2 *x * fT[i] - fT[i-1];
      }
      double sum = p[0]*fT[0];
      for (int i = 1; i<= order; ++i) {
         sum += p[i] * fT[i];
      }
      return sum;
   }

private:
   double fA;
   double fB;
   std::vector<double> fT; // polynomial
   std::vector<double> fC; // coefficients
};

TF1* fit_cheb(int n, double min, double max, TH1* h)
{
    Chebyshev * cheb = new Chebyshev(n, min, max);
    TF1 * f1 = new TF1("f1",cheb,min,max,n+1,"Chebyshev");
    for (int i = 0; i <=n; ++i) f1->SetParameter(i,1);
    h->Fit(f1);

    return f1;
}

TH1* fit_ZptWeightHist(TH1* hZptWeight, bool isDetUnfold=true)
{

    TFile outfile("ZptWeight_final.root","UPDATE"); 
    // extract histogram in each mass bin
    TFile fhist("hists/unfold_input.root"); // to get binning definition
    TUnfoldBinning* pt_binning_Rec = NULL; 
    TString PtRec_binName = "Detector/Pt_FineCoarse/Rec_Pt";
    if(!isDetUnfold)
        PtRec_binName = "Detector/Pt_FineCoarse/Gen_Pt";
    pt_binning_Rec = (TUnfoldBinning*)fhist.Get(PtRec_binName);

    double mass[6] = {51., 65., 82., 102., 201., 321.};
    TString iter_;
    TH1* hout = pt_binning_Rec->CreateHistogram("fit_to_bin"); 
    hout->Sumw2();
    for(int ibin = 1; ibin <=hout->GetXaxis()->GetNbins(); ibin++) 
    {
        hout->SetBinContent(ibin, 1.);
    }

    for(int i = 0; i < 6; i++)
    {
        iter_.Form("%d",i);

        TH1* hist_to_fit = pt_binning_Rec->ExtractHistogram("histo_fit_test", hZptWeight, 0, true, "pt[];mass[UOC" + iter_ + "]");
        int fit_degree = 10;
        TF1* fit_func = fit_cheb(fit_degree, 0., 100., hist_to_fit);

        for(int ibin = 1; ibin <=hist_to_fit->GetXaxis()->GetNbins(); ibin++)    
        {
            double pt = hist_to_fit->GetBinCenter(ibin);
            int bin = pt_binning_Rec->GetGlobalBinNumber(pt, mass[i]);
            hout->SetBinContent(bin, fit_func->Eval(pt));
            cout << "pt: " << pt << " bin: " << bin << " eval: " << fit_func->Eval(pt) << endl;
        }
        delete hist_to_fit;
        delete fit_func;
    }
 
    outfile.cd();
    hout->SetName("ZptWeight_for_next_iter0");
    hout->Write();

    //hZptWeight->SetName("hZptWeight_iter1");
    hZptWeight->Write();
    outfile.Write();

    return hout;
}

void draw_plots()
{

    //TString fin_string = "./ZptWeight_electron_current_result.root";
    TString fin_string = "/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/output/2016/electron_detector_dressedDRp1_extended/DetUNFOLD_electron_2016.root";
    TFile fin(fin_string);
    TString h1_name = "unfolded/Mass/histo_Data_UnfoldModel";
    TString h2_name = "unfolded/Mass/histo_Data";
    TString h3_name = "unfolded/Mass/histo_DY";

    TH1 * h1 = (TH1*) fin.Get(h1_name);
    TH1 * h2 = (TH1*) fin.Get(h2_name);
    TH1 * h3 = (TH1*) fin.Get(h3_name);

    TFile fhist("hists/unfold_input.root");
    TUnfoldBinning* pt_binning_Rec = NULL; 
    TString PtRec_binName = "Detector/Mass_FineCoarse/Gen_Mass";
    pt_binning_Rec = (TUnfoldBinning*)fhist.Get(PtRec_binName);

    TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,800,480);
    //h1 = (TH1*) pt_binning_Rec->ExtractHistogram("histo1", h1, 0, false, "pt[O];mass[UO]");
    //h2 = (TH1*) pt_binning_Rec->ExtractHistogram("histo2", h2, 0, false, "pt[O];mass[UO]");
    //h3 = (TH1*) pt_binning_Rec->ExtractHistogram("histo3", h3, 0, false, "pt[O];mass[UO]");

    h1 = (TH1*) pt_binning_Rec->ExtractHistogram("histo1", h1, 0, false, "mass[UO];pt[O]");
    h2 = (TH1*) pt_binning_Rec->ExtractHistogram("histo2", h2, 0, false, "mass[UO];pt[O]");
    h3 = (TH1*) pt_binning_Rec->ExtractHistogram("histo3", h3, 0, false, "mass[UO];pt[O]");

    h1->SetMinimum(0.9); 
    h1->SetMaximum(1.1); 
    h1->SetMarkerStyle(20);
    h1->SetLineColor(kBlack);
    h1->SetStats(false);
    h1->SetTitle("");
    //h1->GetXaxis()->SetTitle("p_{T} bin index");
    h1->GetXaxis()->SetTitle("Mass bin index");
    //h1->GetYaxis()->SetTitle("Z p_{T} weights");
    h1->GetYaxis()->SetTitle("ratio");
    h1->Divide(h2);
    h1->Draw();

    h3->Divide(h2);

    //h2->SetLineColor(kRed);
    //h2->Draw("hist same");
    //
    h3->SetLineColor(kRed);
    //h3->Draw("hist same");

    double npt_bin = 9;
    TLine *l1 = new TLine(npt_bin+0.5,  0.9, npt_bin+0.5,  1.1);
    TLine *l2 = new TLine(npt_bin*2+0.5,0.9,npt_bin*2+0.5,1.1);
    TLine *l3 = new TLine(npt_bin*3+0.5,0.9,npt_bin*3+0.5,1.1);
    TLine *l4 = new TLine(npt_bin*4+0.5,0.9,npt_bin*4+0.5,1.1);
    TLine *l5 = new TLine(npt_bin*5+0.5,0.9,npt_bin*5+0.5,1.1);
    //l1->Draw();
    //l2->Draw();
    //l3->Draw();
    //l4->Draw();
    //l5->Draw();
    l1->SetLineColor(kBlue);
    l2->SetLineColor(kBlue);
    l3->SetLineColor(kBlue);
    l4->SetLineColor(kBlue);
    l5->SetLineColor(kBlue);

    c1->SaveAs("test.pdf");
}

void draw_plots_mass()
{

    TString fin_string = "./ZptWeight_electron_current_result.root";
    TFile fin(fin_string);
    TH1 * h1 = (TH1*) fin.Get("hDYMassReweighted_iter1");
    TH1 * h2 = (TH1*) fin.Get("histo_DoubleEG_Mass");
    TH1 * h3 = (TH1*) fin.Get("histo_DYJets_Mass");

    TFile fhist("hists/unfold_input.root");
    TUnfoldBinning* pt_binning_Rec = NULL; 
    TString MassRec_binName = "Detector/Mass_FineCoarse/Rec_Mass";
    pt_binning_Rec = (TUnfoldBinning*)fhist.Get(MassRec_binName);

    TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,800,480);
    h1 = (TH1*) pt_binning_Rec->ExtractHistogram("histo1", h1, 0, false, "mass[UO];pt[O]");
    h2 = (TH1*) pt_binning_Rec->ExtractHistogram("histo2", h2, 0, false, "mass[UO];pt[O]");
    h3 = (TH1*) pt_binning_Rec->ExtractHistogram("histo3", h3, 0, false, "mass[UO];pt[O]");

    h1->SetMinimum(0.); 
    h1->SetMaximum(2.); 
    h1->SetMarkerStyle(20);
    h1->SetLineColor(kBlack);
    h1->SetStats(false);
    h1->SetTitle("");
    h1->GetXaxis()->SetTitle("Mass bin index");
    //h1->GetYaxis()->SetTitle("Z p_{T} weights");
    h1->GetYaxis()->SetTitle("ratio");
    h1->Divide(h2);
    h1->Draw();

    h3->Divide(h2);

    //h2->SetLineColor(kRed);
    //h2->Draw("hist same");
    //
    h3->SetLineColor(kRed);
    h3->Draw("hist same");

    double npt_bin = 18;
    TLine *l1 = new TLine(18+0.5,0., 18+0.5,2.);
    TLine *l2 = new TLine(18*2+0.5,0,18*2+0.5,2.);
    TLine *l3 = new TLine(18*3+0.5,0,18*3+0.5,2.);
    TLine *l4 = new TLine(18*4+0.5,0,18*4+0.5,2.);
    TLine *l5 = new TLine(18*5+0.5,0,18*5+0.5,2.);
    //l1->Draw();
    //l2->Draw();
    //l3->Draw();
    //l4->Draw();
    //l5->Draw();
    l1->SetLineColor(kBlue);
    l2->SetLineColor(kBlue);
    l3->SetLineColor(kBlue);
    l4->SetLineColor(kBlue);
    l5->SetLineColor(kBlue);

    c1->SaveAs("test.pdf");
}


double get_piecewise_parameters(double knot, int par_index, double prev_p0, double prev_p1, double prev_p2, double prev_p3, double current_p3)
{
    if(par_index == 0)
    {
        return prev_p0 + pow(knot, 3) * (prev_p3 - current_p3);
    }

    if(par_index == 1)
    {
        return prev_p1 - 3 * pow(knot, 2) * (prev_p3 - current_p3);
    }

    if(par_index == 2)
    {
        return prev_p2 + 3 * knot * (prev_p3 - current_p3);
    }

    return 0.;
}


double piecewise_cubic(double *x, double *par)
{

    // bin edges 0., 2., 4., 6., 8., 10., 12., 14., 18., 23, 28., 34., 40., 47.5, 55., 65., 75., 87.5, 100.
    // bin edges 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 16., 18., 20.5, 23, 25.5, 28., 31., 34., 37., 40., 43.75, 47.5, 51.25, 55., 60., 65., 70., 75., 81.25, 87.5, 93.75, 100.
    // 6., 20., 40., 70
    double knot1 = 6.; // 8, 35, 65, 100
    double knot2 = 12.;
    double knot3 = 23.;
    double knot4 = 40.;
    double knot5 = 65.;

    double second_cubic_par0 = get_piecewise_parameters(knot1, 0, par[0], par[1], par[2], par[3], par[4]);
    double second_cubic_par1 = get_piecewise_parameters(knot1, 1, par[0], par[1], par[2], par[3], par[4]);
    double second_cubic_par2 = get_piecewise_parameters(knot1, 2, par[0], par[1], par[2], par[3], par[4]);

    double third_cubic_par0 = get_piecewise_parameters(knot2, 0, second_cubic_par0, second_cubic_par1, second_cubic_par2, par[4], par[5]);
    double third_cubic_par1 = get_piecewise_parameters(knot2, 1, second_cubic_par0, second_cubic_par1, second_cubic_par2, par[4], par[5]);
    double third_cubic_par2 = get_piecewise_parameters(knot2, 2, second_cubic_par0, second_cubic_par1, second_cubic_par2, par[4], par[5]);

    double fourth_cubic_par0 = get_piecewise_parameters(knot3, 0, third_cubic_par0, third_cubic_par1, third_cubic_par2, par[5], par[6]);
    double fourth_cubic_par1 = get_piecewise_parameters(knot3, 1, third_cubic_par0, third_cubic_par1, third_cubic_par2, par[5], par[6]);
    double fourth_cubic_par2 = get_piecewise_parameters(knot3, 2, third_cubic_par0, third_cubic_par1, third_cubic_par2, par[5], par[6]);

    double fifth_cubic_par0 = get_piecewise_parameters(knot4, 0, fourth_cubic_par0, fourth_cubic_par1, fourth_cubic_par2, par[6], par[7]);
    double fifth_cubic_par1 = get_piecewise_parameters(knot4, 1, fourth_cubic_par0, fourth_cubic_par1, fourth_cubic_par2, par[6], par[7]);
    double fifth_cubic_par2 = get_piecewise_parameters(knot4, 2, fourth_cubic_par0, fourth_cubic_par1, fourth_cubic_par2, par[6], par[7]);

    double sixth_cubic_par0 = get_piecewise_parameters(knot5, 0, fifth_cubic_par0, fifth_cubic_par1, fifth_cubic_par2, par[7], par[8]);
    double sixth_cubic_par1 = get_piecewise_parameters(knot5, 1, fifth_cubic_par0, fifth_cubic_par1, fifth_cubic_par2, par[7], par[8]);
    double sixth_cubic_par2 = get_piecewise_parameters(knot5, 2, fifth_cubic_par0, fifth_cubic_par1, fifth_cubic_par2, par[7], par[8]);

    return (x[0] < knot1) *                                (par[0]+x[0]*par[1]+x[0]*x[0]*par[2]+x[0]*x[0]*x[0]*par[3])
        +  (x[0] >= knot1 && x[0] < knot2)   *             (second_cubic_par0 + x[0] * second_cubic_par1 + x[0] * x[0] * second_cubic_par2  + x[0]*x[0]*x[0]*par[4])
        +  (x[0] >= knot2 && x[0] < knot3)   *             (third_cubic_par0 +  x[0] * third_cubic_par1 +  x[0] * x[0] * third_cubic_par2  +  x[0]*x[0]*x[0]*par[5])
        +  (x[0] >= knot3 && x[0] < knot4)   *             (fourth_cubic_par0 + x[0] * fourth_cubic_par1 + x[0] * x[0] * fourth_cubic_par2  + x[0]*x[0]*x[0]*par[6])
        +  (x[0] >= knot4 && x[0] < knot5)   *             (fifth_cubic_par0 +  x[0] * fifth_cubic_par1 +  x[0] * x[0] * fifth_cubic_par2  +  x[0]*x[0]*x[0]*par[7])
        +  (x[0] >= knot5)   *                             (sixth_cubic_par0 +  x[0] * sixth_cubic_par1 +  x[0] * x[0] * sixth_cubic_par2  +  x[0]*x[0]*x[0]*par[8]);
        //+  (x[0] >= knot4)   *                           (fourth_cubic_par0 - fourth_cubic_par2 * pow(knot4, 2) - 2 * par[6] * pow(knot4, 3) + (fourth_cubic_par1 + 2 * fourth_cubic_par2 * knot4 + 3 * par[6] * pow(knot4, 2))*x[0] ); // Linear function
        //+  (x[0] >= knot4)   *                             (fourth_cubic_par0 + knot4 * fourth_cubic_par1 + pow(knot4,2) * fourth_cubic_par2  + pow(knot4,3)*par[6]); // constant
}

void draw(TH1* nominator, TH1* denominator, int color, TString option)
{
    nominator->Divide(denominator);
    nominator->SetMinimum(0.7);
    nominator->SetMaximum(1.3);
    nominator->SetMarkerStyle(20);
    nominator->SetMarkerColor(color);
    nominator->SetLineColor(color);
    nominator->Draw(option);
}

