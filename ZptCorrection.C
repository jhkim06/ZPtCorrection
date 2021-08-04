#include "TH1.h"
#include "TFile.h"
#include "TUnfoldBinning.h"
#include "TTree.h"
#include "TChain.h"
#include "TVectorD.h" 
#include <fstream> 

#include <iostream>
#include <cstdio>
#include <string>
#include <regex>

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
    double knot1 = 6.;
    double knot2 = 20.;
    double knot3 = 40.;
    double knot4 = 70.;

    double second_cubic_par0 = get_piecewise_parameters(knot1, 0, par[0], par[1], par[2], par[3], par[4]);
    double second_cubic_par1 = get_piecewise_parameters(knot1, 1, par[0], par[1], par[2], par[3], par[4]);
    double second_cubic_par2 = get_piecewise_parameters(knot1, 2, par[0], par[1], par[2], par[3], par[4]);

    double third_cubic_par0 = get_piecewise_parameters(knot2, 0, second_cubic_par0, second_cubic_par1, second_cubic_par2, par[4], par[5]);
    double third_cubic_par1 = get_piecewise_parameters(knot2, 1, second_cubic_par0, second_cubic_par1, second_cubic_par2, par[4], par[5]);
    double third_cubic_par2 = get_piecewise_parameters(knot2, 2, second_cubic_par0, second_cubic_par1, second_cubic_par2, par[4], par[5]);

    double fourth_cubic_par0 = get_piecewise_parameters(knot3, 0, third_cubic_par0, third_cubic_par1, third_cubic_par2, par[5], par[6]);
    double fourth_cubic_par1 = get_piecewise_parameters(knot3, 1, third_cubic_par0, third_cubic_par1, third_cubic_par2, par[5], par[6]);
    double fourth_cubic_par2 = get_piecewise_parameters(knot3, 2, third_cubic_par0, third_cubic_par1, third_cubic_par2, par[5], par[6]);

    double fifth_cubic_par0 = fourth_cubic_par0 + pow(knot4,3) * (par[6]-par[7]);
    double fifth_cubic_par1 = fourth_cubic_par1 - 3 * pow(knot4,2) * (par[6]-par[7]);
    double fifth_cubic_par2 = fourth_cubic_par2 + 3 * knot4 * (par[6]-par[7]);

    return (x[0] < knot1) *                                (par[0]+x[0]*par[1]+x[0]*x[0]*par[2]+x[0]*x[0]*x[0]*par[3]) 
        +  (x[0] >= knot1 && x[0] < knot2)   *             (second_cubic_par0 + x[0] * second_cubic_par1 + x[0] * x[0] * second_cubic_par2  + x[0]*x[0]*x[0]*par[4])
        +  (x[0] >= knot2 && x[0] < knot3)   *             (third_cubic_par0 +  x[0] * third_cubic_par1 +  x[0] * x[0] * third_cubic_par2  +  x[0]*x[0]*x[0]*par[5])
        +  (x[0] >= knot3 && x[0] < knot4)   *             (fourth_cubic_par0 + x[0] * fourth_cubic_par1 + x[0] * x[0] * fourth_cubic_par2  + x[0]*x[0]*x[0]*par[6])
        //+  (x[0] >= knot4)   *                           (fifth_cubic_par0 +  x[0] * fifth_cubic_par1 +  x[0] * x[0] * fifth_cubic_par2  +  x[0]*x[0]*x[0]*par[7]);
        //+  (x[0] >= knot4)   *                           (fourth_cubic_par0 - fourth_cubic_par2 * pow(knot4, 2) - 2 * par[6] * pow(knot4, 3) + (fourth_cubic_par1 + 2 * fourth_cubic_par2 * knot4 + 3 * par[6] * pow(knot4, 2))*x[0] );
        +  (x[0] >= knot4)   *                             (fourth_cubic_par0 + knot4 * fourth_cubic_par1 + pow(knot4,2) * fourth_cubic_par2  + pow(knot4,3)*par[6]); // constant
}


void GetZptReweight(bool isElectron = true, int massBin = 2, TString histFile = "hists/ISR_detector_plots_electron_new.root")
{
    gROOT->SetBatch();
    TH1::AddDirectory(kFALSE);

    TF1 *fitFcn = NULL;

    // Input histograms for Data and Backgrounds at reconstruction level
    if(!isElectron)
    {
        histFile = "hists/ISR_detector_plots_muon_new.root";
    }

    TFile fhist(histFile);
    // Output file
    TString output_postfix = "electron";
    TString massBinStr = "";
    TUnfoldBinning* pt_binning_Rec = NULL;
    TUnfoldBinning* pt_binning_RecAnalysis = NULL;

    if(massBin ==0)
    {
        massBinStr = "_m50to64";
        if(!isElectron)
            massBinStr = "_m40to64";
    }
    else if(massBin ==1)
    {
        massBinStr = "_m64to81";
        if(!isElectron)
            massBinStr = "_m60to81";
    }
    else if(massBin ==2)
    {
        massBinStr = "_m81to101";
    }
    // m100to320
    else if(massBin ==3)
    {
        massBinStr = "_m101to200";
    }
    else if(massBin ==4)
    {
        massBinStr = "_m200to320";
    }
    else if(massBin ==5)
    {
        massBinStr = "_m101to320";
    }
    else
    {
        // Get
        TString Rec_binName = "Detector/Pt_FineCoarse/Rec_Pt";
        TString RecAnalysis_binName = "Detector/Pt_FineCoarse/Rec_Pt";
        pt_binning_Rec = (TUnfoldBinning*)fhist.Get(Rec_binName);
        pt_binning_RecAnalysis = (TUnfoldBinning*)fhist.Get(RecAnalysis_binName);
    }

    if(!isElectron)
    {
        output_postfix = "muon";
    }
    TFile outfile("ZptWeight_"+output_postfix+massBinStr+".root","UPDATE");

    TH1* hDYPtReweighted = NULL;

    int stop = 0;
    int iter = 0;

    TH1* hDataPt = NULL;
    TH1* hBkgPt = NULL;
    TH1* hZptWeight = NULL;
    TH1* hZptWeight_temp = NULL;
    TH1* hDYPtBfReweighted = NULL;
    TH1* hDYPtReweighted_previous = NULL;

    TH1* hZptWeight_temp_m81to101 = NULL;

    int previous_iter = 0;
    TString previous_iter_;
    do
    {
        previous_iter_.Form("%d", previous_iter + 1);
        hZptWeight_temp = (TH1*)outfile.Get("ZptWeight_iter"+previous_iter_);
        if(hZptWeight_temp != NULL) cout << previous_iter_ << " th iteration found!" << endl;
        previous_iter++;
    }
    while(hZptWeight_temp != NULL);

    iter = previous_iter - 1;
    previous_iter_.Form("%d", iter);

    if(iter != 0)
    {
        TString current_iter_;
        current_iter_.Form("%d", iter + 1);
        hZptWeight_temp = (TH1*)outfile.Get("ZptWeightInputFor_iter"+current_iter_);
    }

    TString dir = "Detector"+massBinStr;
    TString histDir = "/Pt_FineCoarse/";
    //TString dir = "Detector_m81to101";
    TString dataName = "DoubleEG";
    TString channelName = "EE";
    if(isElectron == false)
    {
        dataName = "DoubleMuon";
        channelName = "MuMu";
    }

    // Get Data and Background histograms
    hDataPt=(TH1*)fhist.Get(dir+histDir+"histo_"+dataName);
    hBkgPt =(TH1*)fhist.Get(dir+histDir+"histo_ZZ_pythia");
    hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_WZ_pythia"));
    hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_WW_pythia"));
    hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_TTLL_powheg"));
    hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_DYJets10to50ToTauTau"));
    hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_DYJetsToTauTau"));

    TH1* hDataPt_m81to101 = pt_binning_Rec->ExtractHistogram("histo_"+dataName + "_m81to101", hDataPt, 0, true, "pt[];mass[UOC2]"); 
    TH1* hBkgPt_m81to101 = pt_binning_Rec->ExtractHistogram("histo_Bkg_m81to101", hBkgPt, 0, true, "pt[];mass[UOC2]");

    hDYPtBfReweighted=(TH1*)fhist.Get(dir+histDir+"histo_DYJetsTo"+channelName);
    hDYPtBfReweighted->Add((TH1*)fhist.Get(dir+histDir+"histo_DYJets10to50To"+channelName));
    hDYPtReweighted=(TH1*)fhist.Get(dir+histDir+"histo_DYJetsTo"+channelName);
    hDYPtReweighted->Add((TH1*)fhist.Get(dir+histDir+"histo_DYJets10to50To"+channelName));

    TH1* hDYPtBfReweighted_m81to101 = pt_binning_Rec->ExtractHistogram("histo_DYPtBfReweighted_m81to101", hDYPtBfReweighted, 0, true, "pt[];mass[UOC2]");
    TH1* hDYPtReweighted_m81to101 = pt_binning_Rec->ExtractHistogram("histo_DYPtReweighted_m81to101", hDYPtReweighted, 0, true, "pt[];mass[UOC2]");

    if(hZptWeight_temp != NULL)
    {
        cout << "Get " + previous_iter_ + " th reweighted DY histogram..." << endl;
        hDYPtReweighted_previous=(TH1*)outfile.Get("hDYReweighted_iter" + previous_iter_);
    }
    else
    {
        hDYPtReweighted_previous=(TH1*)fhist.Get(dir+histDir+"histo_DYJetsTo"+channelName);
        hDYPtReweighted_previous->Add((TH1*)fhist.Get(dir+histDir+"histo_DYJets10to50To"+channelName));
    }

    if(hZptWeight_temp == NULL)
    {
        cout << "No previous ZptWeight exits..." << endl;
        hZptWeight_temp=(TH1*)fhist.Get(dir+histDir+"/histo_"+dataName);
        hZptWeight_temp_m81to101 = pt_binning_Rec->ExtractHistogram("histo_"+dataName + "_temp_m81to101", hZptWeight_temp, 0, true, "pt[];mass[UOC2]");
        hZptWeight_temp_m81to101->Add(hBkgPt_m81to101, -1);                   // Subtract bkg from data

        hZptWeight_temp_m81to101->Scale(hDYPtBfReweighted_m81to101->Integral()/hZptWeight_temp_m81to101->Integral());
        hZptWeight_temp_m81to101->Divide(hDYPtBfReweighted_m81to101);
    
        // Piece-wise fit
        fitFcn = new TF1("fitFcn", piecewise_cubic, 0, 100, 7);
        hZptWeight_temp_m81to101->Fit(fitFcn);
        //cout<<"Eval fn(99.) = " << fitFcn->Eval(99.0) << endl;

        outfile.cd();
        hZptWeight_temp_m81to101->Write();
    }


    // DY MC tree
    // Reconstruction level variables
    bool evt_tag_analysisevnt_sel_rec;
    bool evt_tag_dielectron_rec;
    bool evt_tag_dimuon_rec;

    Double_t evt_weight_total_rec;
    Double_t evt_weight_recoSF_rec;
    Double_t evt_weight_idSF_rec;
    Double_t evt_weight_isoSF_rec;
    Double_t evt_weight_trigSF_rec;
    Double_t evt_weight_trigSFDZ_rec;

    Double_t dilep_pt_rec;
    Double_t dilep_mass_rec;

    // Generator level variables
    bool evt_tag_ditau_hardprocess;
    bool evt_tag_dielectron_hardprocess;
    bool evt_tag_dimuon_hardprocess;
    bool evt_tag_dimuon_promptfinal;
    bool evt_tag_dielectron_promptfinal;

    Double_t dilep_pt_FSRgamma_gen_ispromptfinal;
    Double_t dilep_mass_FSRgamma_gen_ispromptfinal;
    Double_t evt_weight_total_gen;

    Long64_t nentries;

    TTree* tsignal;

    TChain * chain = new TChain("recoTree/SKFlat","ZptCorrection");
    string line;
    ifstream myfile("/home/jhkim/ISR_Run2/ZptCorrection/2016ISR/ntuples/list.txt");
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            if( line.length() != 0 )
            {
                cout << line << '\n';
                char file_path[300];
                strcpy(file_path, line.c_str());
                chain->Add(file_path);
            }
        }
        myfile.close();
    }
    tsignal = chain;

    tsignal->SetBranchAddress("evt_tag_analysisevnt_sel_rec_Nominal",&evt_tag_analysisevnt_sel_rec);
    tsignal->SetBranchAddress("evt_tag_dielectron_rec_Nominal",&evt_tag_dielectron_rec);
    tsignal->SetBranchAddress("evt_tag_dimuon_rec_Nominal",&evt_tag_dimuon_rec);
    tsignal->SetBranchAddress("evt_weight_total_rec",&evt_weight_total_rec);
    tsignal->SetBranchAddress("evt_weight_recoSF_rec_Nominal",&evt_weight_recoSF_rec);
    tsignal->SetBranchAddress("evt_weight_idSF_rec_Nominal",&evt_weight_idSF_rec);
    tsignal->SetBranchAddress("evt_weight_isoSF_rec_Nominal",&evt_weight_isoSF_rec);
    tsignal->SetBranchAddress("evt_weight_trigSF_rec_Nominal",&evt_weight_trigSF_rec);
    tsignal->SetBranchAddress("evt_weight_trigSFDZ_rec_Nominal",&evt_weight_trigSFDZ_rec);
    tsignal->SetBranchAddress("dilep_pt_rec_Nominal",&dilep_pt_rec);
    tsignal->SetBranchAddress("dilep_mass_rec_Nominal",&dilep_mass_rec);

    tsignal->SetBranchAddress("evt_tag_ditau_hardprocess",&evt_tag_ditau_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dielectron_hardprocess",&evt_tag_dielectron_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dimuon_hardprocess",&evt_tag_dimuon_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dimuon_promptfinal",&evt_tag_dimuon_promptfinal);
    tsignal->SetBranchAddress("evt_tag_dielectron_promptfinal",&evt_tag_dielectron_promptfinal);
    tsignal->SetBranchAddress("dilep_pt_FSRgamma_gen_ispromptfinal",&dilep_pt_FSRgamma_gen_ispromptfinal);
    tsignal->SetBranchAddress("dilep_mass_FSRgamma_gen_ispromptfinal",&dilep_mass_FSRgamma_gen_ispromptfinal);
    tsignal->SetBranchAddress("evt_weight_total_gen",&evt_weight_total_gen);

    nentries=tsignal->GetEntries();

    do
    {
        cout << "Start " << iter+1 << " th iteration..." << endl;

        // Get Z pT weight for each mass bin seperately
        TH1* hDYPtReweighted_temp;
        if(massBin < 6)
        {
            hDYPtReweighted_temp = (TH1*)hDYPtReweighted->Clone("hDYReweighted_temp");
            hDYPtReweighted_temp->Reset(); // Reset before apply the current Z pt reweight in this iteration
        }
        else
        {
            cout << "Histogram created using TUnfold package!" << endl;
            hDYPtReweighted_temp = pt_binning_Rec->CreateHistogram("hDYReweighted_temp");
        }

        cout<<"Start signal loop"<<endl;
        for(Long64_t i=0;i<nentries;i++)
        //for(int i=0;i<100000;i++)
        {
            tsignal->GetEntry(i);
            if(i%10000000==0) cout<<i<< " /" << nentries << " (" << (double)(100.*i/nentries) << "%)" << endl;

            if(evt_tag_ditau_hardprocess == false)
            {
                Double_t GenZpt = 0., GenZmass = 0.;          
                if(evt_tag_ditau_hardprocess) continue; 
                if(isElectron)
                    if(evt_tag_dimuon_hardprocess) continue;
                if(!isElectron)
                    if(evt_tag_dielectron_hardprocess) continue;
                // No tautau event saved in the input tree

                GenZpt   = dilep_pt_FSRgamma_gen_ispromptfinal;
                GenZmass = dilep_mass_FSRgamma_gen_ispromptfinal;

                if(GenZpt < 0 || GenZmass < 0 ) continue;

                // TODO Get mass also!

                bool recChannel = evt_tag_dielectron_rec;
                Double_t totWeight = evt_weight_total_gen*evt_weight_total_rec*evt_weight_recoSF_rec*evt_weight_idSF_rec*evt_weight_trigSF_rec*evt_weight_trigSFDZ_rec;
                if(!isElectron)
                {
                    recChannel = evt_tag_dimuon_rec;
                    totWeight = evt_weight_total_gen*evt_weight_total_rec*evt_weight_isoSF_rec*evt_weight_idSF_rec*evt_weight_trigSF_rec;
                }

                if(evt_tag_analysisevnt_sel_rec && recChannel)      // TODO require generator cuts 
                {
                    double low_mass = 15.;
                    double high_mass = 3000.;


                    if(dilep_mass_rec > low_mass && dilep_mass_rec < high_mass && dilep_pt_rec < 3000.)
                    {
                        //     
                        if( GenZpt > 100.) GenZpt = 99.9;
            
                        // Use Z pt weights from the 81 < mass < 101 for the first iteration  
                        //if( GenZmass > high_mass) GenZmass = 319.5;
                        //if( GenZmass < low_mass) GenZmass = low_mass + 0.5;

                        int bin;
                        Double_t zpt_weight = 1.;
                        if(massBin < 6)
                        {
                            bin = hZptWeight_temp->FindBin(GenZpt);
                            zpt_weight = hZptWeight_temp->GetBinContent(bin);
                            hDYPtReweighted_temp->Fill(dilep_pt_rec, totWeight*zpt_weight);
                        }
                        else
                        {
                            double zpt_weight = fitFcn->Eval(GenZpt); // get Z pt weight from the previous fit
                            //bin = hZptWeight_temp_m81to101->FindBin(GenZpt);
                            //zpt_weight = hZptWeight_temp_m81to101->GetBinContent(bin);
                            if(zpt_weight == 0 ) zpt_weight = 1.;

                            //if(fabs(zpt_weight_-zpt_weight)/ zpt_weight > 0.2)
                            //{
                            //    cout << "fit : " << zpt_weight_ << " bin: " << zpt_weight << endl;
                            //    cout << "mass: " << dilep_mass_rec << " pt: " << dilep_pt_rec << "gen pt: " << GenZpt << endl;
                            //    // fit : 94895.3 bin: 1
                            //    // mass: 69.0003 pt: 35.6966
                            //}
                            //cout << "gen pt, mass: " << GenZpt << " " << GenZmass << " bin: " << bin << " weight " << zpt_weight << endl;
                            //cout << "rec pt, mass: " << dilep_pt_rec << " " << dilep_mass_rec << endl;

                            int recBin = pt_binning_Rec->GetGlobalBinNumber(dilep_pt_rec, dilep_mass_rec);
                            hDYPtReweighted_temp->Fill(recBin, totWeight*zpt_weight);
                            //hDYPtReweighted_temp->Fill(recBin, totWeight);
                        }

                    }// Z peak
                }
            }
        }// DY loop
        iter++;
        TString iter_; iter_.Form("%d",iter);


        const TVectorD* tvecd = pt_binning_Rec->GetDistributionBinning(1);
        int nMassBins = tvecd->GetNrows() - 1; // [50., 64., ] 
        if(pt_binning_Rec->HasUnderflow(1)) nMassBins++;
        if(pt_binning_Rec->HasOverflow(1)) nMassBins++;

        int startBin = 1;
        int endBin = 1;
        double nominator_sum = 0.;
        double denominator_sum = 0.;
        // Normalize for each mass region seperately
        // pt_binning_Rec->ExtractHistogram("histo_"+dataName + "_m81to101", hDataPt, 0, true, "pt[];mass[UOC2]");
        // print out normalisation factor
        double norm_;
        norm_ = hDYPtBfReweighted->Integral()/ hDYPtReweighted_temp->Integral();
        cout << "Integral bf reweight: " << hDYPtBfReweighted->Integral() << endl;
        cout << "Integral af reweight: " << hDYPtReweighted_temp->Integral() << endl;
        cout << "norm: " << norm_ << endl;
        hDYPtReweighted_temp->SetName("hDYReweighted_iter" + iter_);

        TH1* hNormalisation_temp = pt_binning_Rec->CreateHistogram("hNormalisation_temp");
        hNormalisation_temp->Sumw2();
        //hZptWeight_temp->Scale(norm_);  // Normalisae the fit function not the weight histogram
        //cout << "bin n: " << hDYPtReweighted_temp->GetXaxis()->GetNbins() << endl;
        for(int ibin = 1; ibin <= hDYPtReweighted_temp->GetXaxis()->GetNbins(); ibin++)
        {
            int bin_number;
            // #[^\s]+ \(Rec_Pt:mass\[[^\s]+]:pt\[[^\s]+]\) 
            string temp_cuurent_bin_name(pt_binning_Rec->GetBinName(ibin));
            string current_mass_bin_name = temp_cuurent_bin_name.substr(temp_cuurent_bin_name.find("mass[")+5, temp_cuurent_bin_name.find("]:pt[")-temp_cuurent_bin_name.find("mass[")-5);
            string temp_next_bin_name(pt_binning_Rec->GetBinName(ibin+1));
            string next_mass_bin_name = temp_next_bin_name.substr(temp_next_bin_name.find("mass[")+5, temp_next_bin_name.find("]:pt[")-temp_next_bin_name.find("mass[")-5);

            //cout << pt_binning_Rec->GetBinName(ibin) << " mass bin name: " << current_mass_bin_name << endl;
            //cout << mass_bin << " " << pt_bin << endl;
            if(next_mass_bin_name.compare(current_mass_bin_name) != 0)
            {

                denominator_sum += hDYPtReweighted_temp->GetBinContent(ibin);
                nominator_sum += hDYPtBfReweighted->GetBinContent(ibin);
                cout << "current mass bin name: " << current_mass_bin_name << " normalization: " << nominator_sum/ denominator_sum << endl; 
                
                endBin = ibin;
                for(int jbin = startBin; jbin <= endBin; jbin++)
                {
                    hNormalisation_temp->SetBinContent(jbin, nominator_sum/ denominator_sum);                   
                }
                startBin=ibin+1;

                denominator_sum = 0;
                nominator_sum = 0; 

            }
            else
            {
                nominator_sum += hDYPtBfReweighted->GetBinContent(ibin);
                denominator_sum += hDYPtReweighted_temp->GetBinContent(ibin);
            }
        }

        hDYPtReweighted_temp->Multiply(hNormalisation_temp);
/*
        // Update reweight histogram
        TH1* hDataPt_temp;
        TH1* hBkgPt_temp;

        hDataPt_temp=(TH1*)hDataPt->Clone("hDataPt_temp");
        hBkgPt_temp=(TH1*)hBkgPt->Clone("hBkgPt_temp");

        hDataPt_temp->Add(hBkgPt_temp, -1);
        hDataPt_temp->Divide(hDYPtReweighted_temp); // Next weight before NORMALISATION applied

        // Save weight and normalisation information, this can be used for reweighting
        hZptWeight = (TH1*)hZptWeight_temp->Clone("ZptWeight_iter" + iter_);

        for(int bin = 1; bin < hZptWeight_temp->GetNbinsX()+1; bin++)
        {
            float current_weight = hZptWeight_temp->GetBinContent(bin);
            hZptWeight_temp->SetBinContent(bin, current_weight*(hDataPt_temp->GetBinContent(bin))); // for next iteration
            cout << bin << " th bin.." << " next weight: " << hZptWeight_temp->GetBinContent(bin) << " current weight: " << current_weight << endl;
        }
        TString next_iter_; next_iter_.Form("%d",iter + 1);
        TH1* hZptWeightInputForNextIter = (TH1*)hZptWeight_temp->Clone("ZptWeightInputFor_iter" + next_iter_);

        // TODO save this comparison
        TH1* hDYPt_comparison = (TH1*)hDYPtReweighted_previous->Clone("hDYcomparison_iter"+ iter_);
        hDYPt_comparison->Divide(hDYPtReweighted_temp);
        hDYPtReweighted_previous->Reset(); // Reset for next iteration
        hDYPtReweighted_previous->Add(hDYPtReweighted_temp); // Save current ZPt reweighted histogram

        cout << iter << " th iteration finished " << endl;
*/
        outfile.cd();
        //hZptWeight->Write();
        hDYPtReweighted_temp->Write();
        hNormalisation_temp->Write();
/*
        hZptWeightInputForNextIter->Write();
        hDYPt_comparison->Write();

        std::cout << "initial N: " << hDYPtBfReweighted->Integral() << " final N: " << hDYPtReweighted_temp->Integral() << std::endl;
        delete hDataPt_temp;
        delete hBkgPt_temp;
        delete hDYPt_comparison;
        delete hZptWeightInputForNextIter;
*/
        delete hZptWeight;
        delete hDYPtReweighted_temp;

        //delete tsignal;
    }
    while(iter != 1); // TODO add condition to exit
    delete chain;

    outfile.cd();
    outfile.Close();
}

