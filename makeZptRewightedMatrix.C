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

void makeZptRewightedMatrix(TString inFile_Zpt, TString inFile_binDef, TString outFile="ZptWeight_matrix.root") 
{
    gROOT->SetBatch();
    TH1::AddDirectory(kFALSE);

    bool isElectron = true;

    TFile fhist_Zpt(inFile_Zpt);
    TH1* hZptWeight = (TH1*)fhist_Zpt.Get("hZptWeight_iter1");

    TFile fhist_binDef(inFile_binDef);
    TString ptRec_binName = "Detector/Pt_FineCoarse/Rec_Pt";
    TString ptGen_binName = "Detector/Pt_FineCoarse/Gen_Pt";
    TString massRec_binName = "Detector/Mass_FineCoarse/Rec_Mass";
    TString massGen_binName = "Detector/Mass_FineCoarse/Gen_Mass";

    TUnfoldBinning* pt_binning_Rec = (TUnfoldBinning*)fhist_binDef.Get(ptRec_binName);
    TUnfoldBinning* pt_binning_Gen = (TUnfoldBinning*)fhist_binDef.Get(ptGen_binName);
    TUnfoldBinning* mass_binning_Rec = (TUnfoldBinning*)fhist_binDef.Get(massRec_binName);
    TUnfoldBinning* mass_binning_Gen = (TUnfoldBinning*)fhist_binDef.Get(massGen_binName);

    // Create matrix
    TH2* hmatrix_pt = TUnfoldBinning::CreateHistogramOfMigrations(pt_binning_Gen, pt_binning_Rec, "reweighted_pt_matrix");
    TH2* hmatrix_mass = TUnfoldBinning::CreateHistogramOfMigrations(mass_binning_Gen, mass_binning_Rec, "reweighted_mass_matrix");
    hmatrix_pt->Sumw2();
    hmatrix_mass->Sumw2();

    TFile outfile(outFile,"UPDATE"); 

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

    Double_t dilep_pt_FSRgammaDRp1_gen_ispromptfinal;
    Double_t dilep_mass_FSRgammaDRp1_gen_ispromptfinal;
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
    tsignal->SetBranchAddress("evt_tag_ditau_hardprocess",&evt_tag_ditau_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dielectron_hardprocess",&evt_tag_dielectron_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dimuon_hardprocess",&evt_tag_dimuon_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dimuon_promptfinal",&evt_tag_dimuon_promptfinal);
    tsignal->SetBranchAddress("evt_tag_dielectron_promptfinal",&evt_tag_dielectron_promptfinal);
    tsignal->SetBranchAddress("dilep_pt_FSRgammaDRp1_gen_ispromptfinal",&dilep_pt_FSRgammaDRp1_gen_ispromptfinal);
    tsignal->SetBranchAddress("dilep_mass_FSRgammaDRp1_gen_ispromptfinal",&dilep_mass_FSRgammaDRp1_gen_ispromptfinal);
    tsignal->SetBranchAddress("evt_weight_total_gen",&evt_weight_total_gen);

    nentries=tsignal->GetEntries();

    cout<<"Start signal loop"<<endl;
    for(Long64_t i=0;i<nentries;i++)
    //for(int i=0;i<100000;i++)
    {
        tsignal->GetEntry(i);
        if(i%10000000==0) cout<<i<< " /" << nentries << " (" << (double)(100.*i/nentries) << "%)" << endl;

        Double_t genZpt = 0., genZmass = 0.;
        if(isElectron && !evt_tag_dielectron_promptfinal) continue;
        if(!isElectron && !evt_tag_dimuon_promptfinal) continue;
        // No tautau event saved in the input tree

        genZpt   = dilep_pt_FSRgammaDRp1_gen_ispromptfinal;
        genZmass = dilep_mass_FSRgammaDRp1_gen_ispromptfinal;

        if(genZpt < 0 || genZmass < 0 ) continue;

        bool event_selection = evt_tag_analysisevnt_sel_rec && evt_tag_dielectron_rec && evt_tag_pass_kinematic_cut_el_FSRgammaDRp1_gen;
        Double_t genWeight = evt_weight_total_gen;
        Double_t recWeight = evt_weight_total_rec*evt_weight_recoSF_rec*evt_weight_idSF_rec*evt_weight_trigSF_rec*evt_weight_trigSFDZ_rec;
        if(!isElectron)
        {
            event_selection = evt_tag_analysisevnt_sel_rec && evt_tag_dimuon_rec && evt_tag_pass_kinematic_cut_mu_FSRgammaDRp1_gen;
            recWeight = evt_weight_total_rec*evt_weight_isoSF_rec*evt_weight_idSF_rec*evt_weight_trigSF_rec;
        }
        Double_t totWeight = genWeight * recWeight; 

        if(event_selection)
        {
            double low_mass = 15.;
            double high_mass = 3000.;
            if(dilep_mass_rec > low_mass && dilep_mass_rec < high_mass && dilep_pt_rec < 3000. &&
                genZmass > low_mass && genZmass < high_mass && genZpt < 3000.) // TODO add generator requirement
            {
                //
                int recBin_pt_index = pt_binning_Rec->GetGlobalBinNumber(dilep_pt_rec, dilep_mass_rec);
                int genBin_pt_index = pt_binning_Gen->GetGlobalBinNumber(genZpt, genZmass);
                int recBin_mass_index = mass_binning_Rec->GetGlobalBinNumber(dilep_mass_rec, dilep_pt_rec);
                int genBin_mass_index = mass_binning_Gen->GetGlobalBinNumber(genZmass, genZpt);

                int recBin_index_for_zpt_weight = pt_binning_Rec->GetGlobalBinNumber(genZpt, genZmass);
                double zpt_weight = hZptWeight->GetBinContent(recBin_index_for_zpt_weight);

                hmatrix_pt->Fill(genBin_pt_index, recBin_pt_index, totWeight * zpt_weight);
                hmatrix_mass->Fill(genBin_mass_index, recBin_mass_index, totWeight * zpt_weight);

                // bin zero
                hmatrix_pt->Fill(genBin_pt_index, 0., zpt_weight * evt_weight_total_gen * (1.- recWeight)); 
                hmatrix_mass->Fill(genBin_mass_index, 0., zpt_weight * evt_weight_total_gen * (1.- recWeight)); 
            }
        }
    }

    outfile.cd();
    hmatrix_pt->Write();
    hmatrix_mass->Write();
    outfile.Write();
}
