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


void GetZptReweight(bool isElectron = true, bool isDetUnfold = true, TString histFile = "hists/ISR_detector_plots_electron_new.root", bool onlyNormalisation = false)
{
    gROOT->SetBatch();
    TH1::AddDirectory(kFALSE);

    // Vector of fit functions for each iteration step
    vector<TF1*> fitFcn;

    if(!isElectron)
    {
        histFile = "hists/ISR_detector_plots_muon_new.root";
    }

    TFile fhist(histFile);
    TString output_postfix = "electron";
    TUnfoldBinning* pt_binning = NULL;
    TUnfoldBinning* mass_binning = NULL;

    TFile* fhist2 = NULL; 
    if(isDetUnfold) 
    {
        fhist2 = new TFile("hists/DY.root", "Read");
    }
    else
    {
        fhist2 = new TFile("/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/inFiles/2016/electron_detector_dressedDRp1_extended/efficiency/acceptance_dRp1.root", "Read");
    }

    // Get
    TString Pt_binName;
    TString Mass_binName;
    if(isDetUnfold)
    {
        Pt_binName = "Detector/Pt_FineCoarse/Rec_Pt";
        Mass_binName = "Detector/Mass_FineCoarse/Rec_Mass";
    }
    else
    {
        Pt_binName = "acceptance/Pt/Gen_Pt";
        Mass_binName = "acceptance/Mass/Gen_Mass";
    }
    pt_binning = (TUnfoldBinning*)fhist.Get(Pt_binName);
    mass_binning = (TUnfoldBinning*)fhist.Get(Mass_binName);

    // Output file
    if(!isElectron)
    {
        output_postfix = "muon";
    }
    TFile outfile("ZptWeight_"+output_postfix+".root","UPDATE");

    int stop = 0;
    int iter = 0;

    TH1* hDataPt = NULL;
    TH1* hBkgPt = NULL;
    TH1* hZptWeight = NULL;
    TH1* hZptWeight_for_next_iteration = NULL;
    TH1* hDYPtBfReweighted = NULL;
    TH1* hDYPtBfReweighted_from_tree = NULL;

    TH1* hDataMass = NULL;
    TH1* hBkgMass = NULL;
    TH1* hDYMassBfReweighted = NULL;
    TH1* hDYMassBfReweighted_from_tree = NULL;

    int previous_iter = 0;
    TString previous_iter_;
    do
    {
        previous_iter_.Form("%d", previous_iter);
        TH1* hZptWeight_temp = (TH1*)outfile.Get("hZptWeight_iter"+previous_iter_);
        if(hZptWeight_temp != NULL)
        {
            cout << previous_iter_ << " th iteration found!" << endl;
            delete hZptWeight_temp;
            previous_iter++;
        }
        else
        {
            if(previous_iter==0) break;
            else
            {
                previous_iter--;
                previous_iter_.Form("%d", previous_iter);
                break;
            }
        }
    }
    while(true);

    iter = previous_iter;
    TString next_iter_;
    next_iter_.Form("%d", iter+1);

    if(iter != 0)
    {
        hZptWeight_for_next_iteration = (TH1*)outfile.Get("ZptWeight_for_next_iter" + next_iter_);
    }

    if(onlyNormalisation)
    {
        // Z pt weight input after fit result
        TFile zpt_external("ZptWeight_final.root","READ");
        hZptWeight_for_next_iteration = (TH1*) zpt_external.Get("ZptWeight_for_next_iter0");
        zpt_external.Close();
    }

    TString dir = "Detector"; // Dressed_DRp1
    TString histDir = "/Pt_FineCoarse/";
    TString dataName = "DoubleEG";
    TString channelName = "EE";
    if(isElectron == false)
    {
        dataName = "DoubleMuon";
        channelName = "MuMu";
    }
    if(!isDetUnfold)
    {
        dir = "acceptance";
        histDir = "Pt";
    }

    // Get Data and Background histograms

    if(isDetUnfold)
    {
        hDataPt=(TH1*)fhist.Get(dir+histDir+"histo_"+dataName);

        hBkgPt =(TH1*)fhist.Get(dir+histDir+"histo_ZZ_pythia");
        hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_WZ_pythia"));
        hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_WW_pythia"));
        hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_TTLL_powheg"));
        hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_DYJets10to50ToTauTau"));
        hBkgPt->Add((TH1*)fhist.Get(dir+histDir+"histo_DYJetsToTauTau"));
        hBkgPt->Add((TH1*)fhist2->Get("Detector_DY_Fake/Pt_FineCoarse/histo_DYJets"));
        hBkgPt->Add((TH1*)fhist2->Get("Detector_DY_Fake/Pt_FineCoarse/histo_DYJets10to50"));

        hDYPtBfReweighted=(TH1*)fhist2->Get("Detector_Dressed_DRp1_Fiducial/Pt_FineCoarse/histo_DYJets");
        hDYPtBfReweighted->Add((TH1*)fhist2->Get("Detector_Dressed_DRp1_Fiducial/Pt_FineCoarse/histo_DYJets10to50"));
    }
    else
    {
        hDataPt=(TH1*)fhist.Get("acceptance/Pt/histo_Data");

        hBkgPt =(TH1*)fhist2->Get("Dressed_DRp1_DY_Fake/Pt_FineCoarse/histo_DYJets");
        hBkgPt->Add((TH1*)fhist2->Get("Dressed_DRp1_DY_Fake/Pt_FineCoarse/histo_DYJets10to50"));

        hDYPtBfReweighted=(TH1*)fhist2->Get("Acceptance/Pt_FineCoarse/histo_DYJets");
        hDYPtBfReweighted->Add((TH1*)fhist2->Get("Acceptance/Pt_FineCoarse/histo_DYJets10to50"));
    }

    hDataPt->Add(hBkgPt, -1);

    // Mass
    histDir = "/Mass_FineCoarse/";

    if(isDetUnfold)
    {
        hDataMass=(TH1*)fhist.Get(dir+histDir+"histo_"+dataName);

        hBkgMass =(TH1*)fhist.Get(dir+histDir+"histo_ZZ_pythia");
        hBkgMass->Add((TH1*)fhist.Get(dir+histDir+"histo_WZ_pythia"));
        hBkgMass->Add((TH1*)fhist.Get(dir+histDir+"histo_WW_pythia"));
        hBkgMass->Add((TH1*)fhist.Get(dir+histDir+"histo_TTLL_powheg"));
        hBkgMass->Add((TH1*)fhist.Get(dir+histDir+"histo_DYJets10to50ToTauTau"));
        hBkgMass->Add((TH1*)fhist.Get(dir+histDir+"histo_DYJetsToTauTau"));

        hBkgMass->Add((TH1*)fhist2->Get("Detector_DY_Fake/Mass_FineCoarse/histo_DYJets"));
        hBkgMass->Add((TH1*)fhist2->Get("Detector_DY_Fake/Mass_FineCoarse/histo_DYJets10to50"));

        hDYMassBfReweighted=(TH1*)fhist2->Get("Detector_Dressed_DRp1_Fiducial/Mass_FineCoarse/histo_DYJets");
        hDYMassBfReweighted->Add((TH1*)fhist2->Get("Detector_Dressed_DRp1_Fiducial/Mass_FineCoarse/histo_DYJets10to50"));
    }
    else
    {
        hDataMass=(TH1*)fhist.Get("acceptance/Mass/histo_Data");

        hBkgMass =(TH1*)fhist2->Get("Dressed_DRp1_DY_Fake/Mass_FineCoarse/histo_DYJets");
        hBkgMass->Add((TH1*)fhist2->Get("Dressed_DRp1_DY_Fake/Mass_FineCoarse/histo_DYJets10to50"));

        hDYMassBfReweighted=(TH1*)fhist2->Get("Acceptance/Mass_FineCoarse/histo_DYJets");
        hDYMassBfReweighted->Add((TH1*)fhist2->Get("Acceptance/Mass_FineCoarse/histo_DYJets10to50"));
    }
    hDataMass->Add(hBkgMass, -1);


    //
    if(hZptWeight_for_next_iteration == NULL)
    {
        cout << "No previous ZptWeight exits..." << endl;
        hZptWeight_for_next_iteration = (TH1*) hDataPt->Clone("hZptWeight_for_next_iteration");

        // TODO make a function for the following normalisation procedure
        int startBin = 1;  
        int endBin = 1;  
        double nominator_sum = 0.;
        double denominator_sum = 0.; 
        TH1* hNormalisation_temp = pt_binning->CreateHistogram("hNormalisation_temp");
        hNormalisation_temp->Sumw2();
        for(int ibin = 1; ibin <= hZptWeight_for_next_iteration->GetXaxis()->GetNbins(); ibin++)
        {
            // Regular expression, #[^\s]+ \(Rec_Pt:mass\[[^\s]+]:pt\[[^\s]+]\) 
            string temp_cuurent_bin_name(pt_binning->GetBinName(ibin));
            string current_mass_bin_name = temp_cuurent_bin_name.substr(temp_cuurent_bin_name.find("mass[")+5, temp_cuurent_bin_name.find("]:pt[")-temp_cuurent_bin_name.find("mass[")-5);
       
            string temp_next_bin_name;
            string next_mass_bin_name; 
            if(ibin < hZptWeight_for_next_iteration->GetXaxis()->GetNbins())
            {
                temp_next_bin_name = (pt_binning->GetBinName(ibin+1));
                next_mass_bin_name = temp_next_bin_name.substr(temp_next_bin_name.find("mass[")+5, temp_next_bin_name.find("]:pt[")-temp_next_bin_name.find("mass[")-5);
            }

            //cout << pt_binning->GetBinName(ibin) << " mass bin name: " << current_mass_bin_name << endl;
            //cout << mass_bin << " " << pt_bin << endl;
            if(ibin == hZptWeight_for_next_iteration->GetXaxis()->GetNbins() || next_mass_bin_name.compare(current_mass_bin_name) != 0)
            {

                denominator_sum += hZptWeight_for_next_iteration->GetBinContent(ibin);
                nominator_sum += hDYPtBfReweighted->GetBinContent(ibin);
                cout << "current mass bin name: " << current_mass_bin_name << " norminator: " << nominator_sum << " denominator: " <<  denominator_sum << endl; 
                
                endBin = ibin;
                for(int jbin = startBin; jbin <= endBin; jbin++)
                {
                    hNormalisation_temp->SetBinContent(jbin, nominator_sum/ denominator_sum);  // 
                }
                startBin=ibin+1;

                denominator_sum = 0;
                nominator_sum = 0; 

            }
            else
            {
                nominator_sum += hDYPtBfReweighted->GetBinContent(ibin);
                denominator_sum += hZptWeight_for_next_iteration->GetBinContent(ibin);
            }
        }

        hZptWeight_for_next_iteration->Multiply(hNormalisation_temp);
        delete hNormalisation_temp;

        hZptWeight_for_next_iteration->Divide(hDYPtBfReweighted);

        for(int ibin = 1; ibin <= hZptWeight_for_next_iteration->GetXaxis()->GetNbins(); ibin++)  
        {
            string temp_cuurent_bin_name(pt_binning->GetBinName(ibin));
            string current_mass_bin_name = temp_cuurent_bin_name.substr(temp_cuurent_bin_name.find("mass[")+5, temp_cuurent_bin_name.find("]:pt[")-temp_cuurent_bin_name.find("mass[")-5);
            string current_pt_bin_name = temp_cuurent_bin_name.substr(temp_cuurent_bin_name.find(":pt[")+4, temp_cuurent_bin_name.find("])")-temp_cuurent_bin_name.find(":pt[")-4);
            cout << "current_mass_bin_name: " << current_mass_bin_name << " current_pt_bin_name: " << current_pt_bin_name << endl;

            if(current_mass_bin_name=="ufl" or current_mass_bin_name=="ofl" or current_pt_bin_name =="ofl")
            {
                hZptWeight_for_next_iteration->SetBinContent(ibin, 1.); // Will not apply Z pt correction for the underflow/overflow bins 
            }
        }
    
        outfile.cd();
        hZptWeight_for_next_iteration->SetName("hZptWeight_iter0");
        hZptWeight_for_next_iteration->SetTitle("Data(Bkg subtracted)/DY MC iteration 0");
        hZptWeight_for_next_iteration->Write();
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

    bool evt_tag_pass_kinematic_cut_el_FSRgammaDRp1_gen;
    bool evt_tag_pass_kinematic_cut_mu_FSRgammaDRp1_gen;
    bool evt_tag_pass_kinematic_cut_el_FSRgamma_gen;
    bool evt_tag_pass_kinematic_cut_mu_FSRgamma_gen;

    Double_t dilep_pt_FSRgammaDRp1_gen_ispromptfinal;
    Double_t dilep_mass_FSRgammaDRp1_gen_ispromptfinal;
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

    tsignal->SetBranchAddress("pass_kinematic_cut_el_FSRgammaDRp1_gen",&evt_tag_pass_kinematic_cut_el_FSRgammaDRp1_gen);
    tsignal->SetBranchAddress("pass_kinematic_cut_mu_FSRgammaDRp1_gen",&evt_tag_pass_kinematic_cut_mu_FSRgammaDRp1_gen);
    tsignal->SetBranchAddress("pass_kinematic_cut_el_FSRgamma_gen",&evt_tag_pass_kinematic_cut_el_FSRgamma_gen);
    tsignal->SetBranchAddress("pass_kinematic_cut_mu_FSRgamma_gen",&evt_tag_pass_kinematic_cut_mu_FSRgamma_gen);
    tsignal->SetBranchAddress("evt_tag_ditau_hardprocess",&evt_tag_ditau_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dielectron_hardprocess",&evt_tag_dielectron_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dimuon_hardprocess",&evt_tag_dimuon_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dimuon_promptfinal",&evt_tag_dimuon_promptfinal);
    tsignal->SetBranchAddress("evt_tag_dielectron_promptfinal",&evt_tag_dielectron_promptfinal);
    tsignal->SetBranchAddress("dilep_pt_FSRgammaDRp1_gen_ispromptfinal",&dilep_pt_FSRgammaDRp1_gen_ispromptfinal);
    tsignal->SetBranchAddress("dilep_mass_FSRgammaDRp1_gen_ispromptfinal",&dilep_mass_FSRgammaDRp1_gen_ispromptfinal);
    tsignal->SetBranchAddress("dilep_pt_FSRgamma_gen_ispromptfinal",&dilep_pt_FSRgamma_gen_ispromptfinal);
    tsignal->SetBranchAddress("dilep_mass_FSRgamma_gen_ispromptfinal",&dilep_mass_FSRgamma_gen_ispromptfinal);
    tsignal->SetBranchAddress("evt_weight_total_gen",&evt_weight_total_gen);

    nentries=tsignal->GetEntries();

    TString iter_;
    do
    {
        cout << "Start " << iter+1 << " th iteration..." << endl;

        TH1* hDYPtReweighted_temp;
        TH1* hDYMassReweighted_temp;
        cout << "Histogram created using TUnfold package!" << endl;
        hDYPtReweighted_temp = pt_binning->CreateHistogram("hDYReweighted_temp");
        hDYMassReweighted_temp = mass_binning->CreateHistogram("hDYMassReweighted_temp");
        TH1* hZptWeight_temp = NULL;
        if(iter > 0)
        {
            hZptWeight_temp = (TH1*)outfile.Get("hZptWeight_iter"+previous_iter_);
        }
        if(iter==0)
        {
            hDYPtBfReweighted_from_tree = pt_binning->CreateHistogram("hDYPtBfReweighted_from_tree");
            hDYMassBfReweighted_from_tree = mass_binning->CreateHistogram("hDYMassBfReweighted_from_tree");
        }

        cout<<"Start signal loop"<<endl;
        for(Long64_t i=0;i<nentries;i++)
        //for(int i=0;i<100000;i++)
        {
            tsignal->GetEntry(i);
            if(i%10000000==0) cout<<i<< " /" << nentries << " (" << (double)(100.*i/nentries) << "%)" << endl;

            Double_t truthZpt = 0., truthZmass = 0.;          
            Double_t measuredZpt = 0., measuredZmass = 0.;
            if(isElectron && !evt_tag_dielectron_promptfinal) continue;
            if(!isElectron && !evt_tag_dimuon_promptfinal) continue;
            // No tautau event saved in the input tree
            //
            measuredZmass = dilep_mass_rec;
            measuredZpt = dilep_pt_rec;

            truthZmass = dilep_mass_FSRgammaDRp1_gen_ispromptfinal;
            truthZpt   = dilep_pt_FSRgammaDRp1_gen_ispromptfinal;

            if(!isDetUnfold)  
            {
                measuredZmass = dilep_mass_FSRgammaDRp1_gen_ispromptfinal;
                measuredZpt = dilep_pt_FSRgammaDRp1_gen_ispromptfinal;

                truthZmass = dilep_mass_FSRgamma_gen_ispromptfinal;
                truthZpt   = dilep_pt_FSRgamma_gen_ispromptfinal;
            }

            if(truthZpt < 0 || truthZmass < 0 ) continue;

            bool event_selection = evt_tag_analysisevnt_sel_rec && evt_tag_dielectron_rec && evt_tag_pass_kinematic_cut_el_FSRgammaDRp1_gen;
            Double_t totWeight = evt_weight_total_gen*evt_weight_total_rec*evt_weight_recoSF_rec*evt_weight_idSF_rec*evt_weight_trigSF_rec*evt_weight_trigSFDZ_rec;

            if(!isDetUnfold)
            {
                event_selection = evt_tag_pass_kinematic_cut_el_FSRgammaDRp1_gen && evt_tag_pass_kinematic_cut_el_FSRgamma_gen;    
                totWeight = evt_weight_total_gen; 
            } 

            if(!isElectron)
            {
                event_selection = evt_tag_analysisevnt_sel_rec && evt_tag_dimuon_rec && evt_tag_pass_kinematic_cut_mu_FSRgammaDRp1_gen;
                totWeight = evt_weight_total_gen*evt_weight_total_rec*evt_weight_isoSF_rec*evt_weight_idSF_rec*evt_weight_trigSF_rec;
            }

            if(event_selection)       
            {
                double low_mass = 15.;
                double high_mass = 3000.;

                if(measuredZmass > low_mass && measuredZmass < high_mass && measuredZpt < 3000. && 
                    truthZmass > low_mass && truthZmass < high_mass && truthZpt < 3000.) // TODO add generator requirement 
                {
                    //     
                    //if( truthZpt > 100.) truthZpt = 99.9;
            
                    int bin;
                    //double zpt_weight = fitFcn->Eval(truthZpt); // get Z pt weight from the previous fit
                    double zpt_weight = 1.;
                    bin = pt_binning->GetGlobalBinNumber(truthZpt, truthZmass);
                    zpt_weight = hZptWeight_for_next_iteration->GetBinContent(bin);
                    if(zpt_weight == 0 ) zpt_weight = 1.;

                    int recBin = pt_binning->GetGlobalBinNumber(measuredZpt, measuredZmass);
                    int recBin_mass = mass_binning->GetGlobalBinNumber(measuredZmass, measuredZpt);
                    hDYPtReweighted_temp->Fill(recBin, totWeight*zpt_weight);

                    if(iter > 0)
                    {
                        double zpt_weight_ = hZptWeight_temp->GetBinContent(bin);
                        if(zpt_weight_ == 0 ) zpt_weight_ = 1.;
                        hDYMassReweighted_temp->Fill(recBin_mass, totWeight*zpt_weight_);
                    }
                    if(iter==0)
                    {
                        hDYPtBfReweighted_from_tree->Fill(recBin, totWeight);
                        hDYMassBfReweighted_from_tree->Fill(recBin_mass, totWeight);
                    }
                    //hDYPtReweighted_temp->Fill(recBin, totWeight);
                }// Z peak
            }
        }// DY loop
        if(iter > 0 )
            hDYMassReweighted_temp->SetName("hDYMassReweighted_iter" + previous_iter_); // Reweighted DY MC at the current iteration

        iter++;
        iter_.Form("%d",iter);
        previous_iter_.Form("%d",iter);

        const TVectorD* tvecd = pt_binning->GetDistributionBinning(1);
        int nMassBins = tvecd->GetNrows() - 1; // [50., 64., ] 
        if(pt_binning->HasUnderflow(1)) nMassBins++;
        if(pt_binning->HasOverflow(1)) nMassBins++;

        int startBin = 1;
        int endBin = 1;
        double nominator_sum = 0.;
        double denominator_sum = 0.;

        // Normalize for each mass region seperately
        double norm_;
        norm_ = hDYPtBfReweighted->Integral()/ hDYPtReweighted_temp->Integral();
        cout << "Integral bf reweight: " << hDYPtBfReweighted->Integral() << endl;
        cout << "Integral af reweight: " << hDYPtReweighted_temp->Integral() << endl;
        cout << "norm: " << norm_ << endl;
        hDYPtReweighted_temp->SetName("hDYReweighted_iter" + iter_); // Reweighted DY MC at the current iteration

        TH1* hNormalisation_temp = pt_binning->CreateHistogram("hNormalisation_temp");
        hNormalisation_temp->Sumw2();
        hNormalisation_temp->SetName("hNormalization_iter" + iter_);
        //cout << "bin n: " << hDYPtReweighted_temp->GetXaxis()->GetNbins() << endl;
        for(int ibin = 1; ibin <= hDYPtReweighted_temp->GetXaxis()->GetNbins(); ibin++)
        {
            int bin_number;
            // Regular expression, #[^\s]+ \(Rec_Pt:mass\[[^\s]+]:pt\[[^\s]+]\) 
            string temp_cuurent_bin_name(pt_binning->GetBinName(ibin));
            string current_mass_bin_name = temp_cuurent_bin_name.substr(temp_cuurent_bin_name.find("mass[")+5, temp_cuurent_bin_name.find("]:pt[")-temp_cuurent_bin_name.find("mass[")-5);

            string temp_next_bin_name;
            string next_mass_bin_name; 
            if(ibin < hZptWeight_for_next_iteration->GetXaxis()->GetNbins())
            {
                temp_next_bin_name = (pt_binning->GetBinName(ibin+1));
                next_mass_bin_name = temp_next_bin_name.substr(temp_next_bin_name.find("mass[")+5, temp_next_bin_name.find("]:pt[")-temp_next_bin_name.find("mass[")-5);
            }

            //cout << pt_binning->GetBinName(ibin) << " mass bin name: " << current_mass_bin_name << endl;
            //cout << mass_bin << " " << pt_bin << endl;
            if(ibin == hDYPtReweighted_temp->GetXaxis()->GetNbins() || next_mass_bin_name.compare(current_mass_bin_name) != 0)
            {

                denominator_sum += hDYPtReweighted_temp->GetBinContent(ibin);
                nominator_sum += hDYPtBfReweighted->GetBinContent(ibin);
                
                cout << "current mass bin name: " << current_mass_bin_name << " norminator: " << nominator_sum << " denominator: " <<  denominator_sum << endl; 
                
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

        hZptWeight_for_next_iteration->Multiply(hNormalisation_temp);   
        hDYPtReweighted_temp->Multiply(hNormalisation_temp);
        // Actuall end of iteration process

        outfile.cd();
        hZptWeight_for_next_iteration->SetName("hZptWeight_iter" + iter_);
        hZptWeight_for_next_iteration->SetTitle("Data(Bkg subtracted)/DY MC iteration " + iter_);
        hZptWeight_for_next_iteration->Write();

        // Save normalized fit function
        // 1. Use Z pt weight from 81 < mass < 101 region only

        // Update reweight histogram
        TH1* hZptWeight_update_temp=(TH1*) hDataPt->Clone("hZptWeight_update_temp_iter" + iter_);  ;

        TH1* hNormalisation_temp_ = pt_binning->CreateHistogram("hNormalisation_temp_");
        hNormalisation_temp_->Sumw2();
        startBin = 1; 
        endBin = 1;
        nominator_sum = 0.;
        denominator_sum = 0.;
        for(int ibin = 1; ibin <= hZptWeight_update_temp->GetXaxis()->GetNbins(); ibin++)
        {
            int bin_number;
            // Regular expression, #[^\s]+ \(Rec_Pt:mass\[[^\s]+]:pt\[[^\s]+]\) 
            string temp_cuurent_bin_name(pt_binning->GetBinName(ibin));
            string current_mass_bin_name = temp_cuurent_bin_name.substr(temp_cuurent_bin_name.find("mass[")+5, temp_cuurent_bin_name.find("]:pt[")-temp_cuurent_bin_name.find("mass[")-5);

            string temp_next_bin_name;
            string next_mass_bin_name; 
            if(ibin < hZptWeight_update_temp->GetXaxis()->GetNbins())
            {
                temp_next_bin_name = (pt_binning->GetBinName(ibin+1));
                next_mass_bin_name = temp_next_bin_name.substr(temp_next_bin_name.find("mass[")+5, temp_next_bin_name.find("]:pt[")-temp_next_bin_name.find("mass[")-5);
            }
            //cout << pt_binning->GetBinName(ibin) << " mass bin name: " << current_mass_bin_name << endl;
            //cout << mass_bin << " " << pt_bin << endl;
            if(ibin == hZptWeight_update_temp->GetXaxis()->GetNbins() || next_mass_bin_name.compare(current_mass_bin_name) != 0)
            {

                denominator_sum += hZptWeight_update_temp->GetBinContent(ibin);
                nominator_sum += hDYPtBfReweighted->GetBinContent(ibin);

                cout << "current mass bin name: " << current_mass_bin_name << " norminator: " << nominator_sum << " denominator: " <<  denominator_sum << endl; 
                
                endBin = ibin;
                for(int jbin = startBin; jbin <= endBin; jbin++)
                {
                    if(current_mass_bin_name=="ufl" or current_mass_bin_name=="ofl")
                        hNormalisation_temp_->SetBinContent(jbin, 1.);                   
                    else
                        hNormalisation_temp_->SetBinContent(jbin, nominator_sum/ denominator_sum);                   
                }
                startBin=ibin+1;

                denominator_sum = 0;
                nominator_sum = 0; 

            }
            else
            {
                nominator_sum += hDYPtBfReweighted->GetBinContent(ibin);
                denominator_sum += hZptWeight_update_temp->GetBinContent(ibin);
            }
        }

        hZptWeight_update_temp->Multiply(hNormalisation_temp_);
        delete hNormalisation_temp_;
        hZptWeight_update_temp->Divide(hDYPtReweighted_temp); 

        hZptWeight_for_next_iteration->Multiply(hZptWeight_update_temp);
        delete hZptWeight_update_temp;

        for(int ibin = 1; ibin <= hZptWeight_for_next_iteration->GetXaxis()->GetNbins(); ibin++)  
        {
            string temp_cuurent_bin_name(pt_binning->GetBinName(ibin));
            string current_mass_bin_name = temp_cuurent_bin_name.substr(temp_cuurent_bin_name.find("mass[")+5, temp_cuurent_bin_name.find("]:pt[")-temp_cuurent_bin_name.find("mass[")-5);
            string current_pt_bin_name = temp_cuurent_bin_name.substr(temp_cuurent_bin_name.find(":pt[")+4, temp_cuurent_bin_name.find("])")-temp_cuurent_bin_name.find(":pt[")-4);
            cout << "current_mass_bin_name: " << current_mass_bin_name << " current_pt_bin_name: " << current_pt_bin_name << endl;

            if(current_mass_bin_name=="ufl" or current_mass_bin_name=="ofl" or current_pt_bin_name =="ofl")
            {
                hZptWeight_for_next_iteration->SetBinContent(ibin, 1.); // Will not apply Z pt correction for the underflow/overflow bins 
            }
        }

        cout << iter << " th iteration finished " << endl;

        outfile.cd();
        // TODO Save only for the first iteration
        if(iter==1)
        {
            hDataPt->Write(); 
            hDYPtBfReweighted->Write();
            hDYPtBfReweighted_from_tree->Write();

            hDataMass->SetName("histo_DoubleEG_Mass");
            hDYMassBfReweighted->SetName("histo_DYJets_Mass");
            hDataMass->Write(); 
            hDYMassBfReweighted->Write();
            hDYMassBfReweighted_from_tree->Write();
        }
        hDYPtReweighted_temp->Write();
        if(iter>1) hDYMassReweighted_temp->Write();
        hNormalisation_temp->Write();
        delete hZptWeight;
        delete hDYPtReweighted_temp;

        //delete tsignal;
    }
    while(iter != 1); // TODO add condition to exit

    outfile.cd();
    iter_.Form("%d",iter+1);
    hZptWeight_for_next_iteration->SetName("ZptWeight_for_next_iter" + iter_);    
    hZptWeight_for_next_iteration->Write();
    delete chain;

    outfile.cd();
    outfile.Close();

}
