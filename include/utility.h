#ifndef UTILITY_H
#define UTILITY_H

#include "TString.h"
#include <string>
#include <vector>


//Author Ezekiel Wertz
//A location to hold useful functions that do a task

namespace utility{

  //hard-coded output directory
  string out_dir_temp = "/volatile/halla/sbs/ktevans/KateJackSBSAnalysis";
  TString output_directory = TString(out_dir_temp);

  //function to determine if what is found was a number or a string
  bool check_number(const char *myChar);

  //convert from degrees to radians
  double DegToRad(double myNum);

  //convert from radians to degrees
  double RadToDeg(double datRad);

  //convert an int to a TString. Just a helper function
  TString intToTString(int datInt);

  //convert a double to a TString. Just a helper function
  TString doubToTString(double datDoub, int prec);

  //This is a helper function that will tell you the default directory for all analysis output
  TString getOutputDir();

  //Helper function to make output file name for parsed data root file
  TString makeOutputFileNameParse(TString exp, TString pass, TString Kin, int SBS_field,TString target);

  //Helper function to make output file name for parse mc root file
  TString makeOutputFileName_MCParse(TString exp, TString Kin, int SBS_field,TString target);

  //Helper function to make output file name for data and mc comparison 
  TString makeOutputFileName_DataMCComp(TString exp, TString pass, TString Kin, int SBS_field,TString target);  

  //Helper function to make report file name for data and mc comparison
  TString makeReportFileName_DataMCComp(TString exp, TString pass, TString Kin, int SBS_field,TString target);

  //Helper function to make output file name for stability studies
  TString makeOutputFileName_Stability(TString exp, TString pass, TString Kin, int SBS_field,TString target, int slice_mode, TString left_right);

  //Helper function to make report file name for stability studies
  TString makeReportFileName_Stability(TString exp, TString pass, TString Kin, int SBS_field,TString target);
  
  //Helper function to make output file name for yield or ratio information
  TString makeYieldReportFileName(TString exp, TString pass, TString Kin, int SBS_field,TString target);

  //Helper function to make output file name for HCal Efficiency from data
  TString makeOutputFileNameHCalEffParse(TString exp, TString pass, TString Kin, int SBS_field,TString target);

  //Helper function to make output file for HCal Efficiency Bootstrap
  TString makeOutputFileName_HCalBS(TString exp, TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event);

  //Helper function for making file name for HCal Efficiency report 
  TString makeReportFileName_HCalEff(TString exp, TString pass, TString Kin, int SBS_field,TString target);
 
  //Helper function for making file name for HCal Efficiency report with boot strap
  TString makeReportFileName_HCalEffBS(TString exp, TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event);

  //Helper function for output file name for scale field study for MC
  TString makeOutputFileName_MCsf(TString exp, TString Kin, int SBS_field,double sf,TString target);
 
  //Helper function to make output file name for HCal Efficiency from MC
  TString makeOutputFileName_HCalEffMC(TString exp, TString Kin);

  //Helper function to make output file name for HCal Efficiency for uniformity
  TString makeOutputFileNameHCalEffUniformity(TString exp,TString pass, TString Kin, int SBS_field,TString target);

  TString makeOutputFileNameHCalEffMapComp(TString exp,TString pass1, TString pass2, TString Kin1, TString Kin2, int SBS_field1, int SBS_field2,TString target1,TString target2);

  //Helper function to make output file name for HCal Efficiency map analysis
  TString makeOutputFileName_HCalEffMap(TString exp,TString pass, TString Kin,TString target);

  //Helper function to find MC histogram files. Right now this only supports John Boyd simulation files. But in future might support others.
  vector<string> findHistFiles(TString replay_type,TString histDirectory,TString partialName);

  //Helper function to match all the MC root files with the already found MC hist files. This function is not ideal, it modifies both of the vectors by the time the function is done. Though it does not have a return type
  void matchMCFiles(TString replay_type,vector<string>& histFiles,vector<string>& rootFiles, TString rootDirectory);

  //Function that searches through two vectors of strings and removes entries without a matching pair
  //since we are sending references to vectors it should modifiy the original reference
  void syncJobNumbers(vector<string>& proton_vec,vector<string>& neutron_vec);

  //Used with regular sim files form JLab-HPC, not from J. Boyd.
  //Function that searches through two vectors of pair of string and float vector and removes entries without a matching pair
  //since we are sending references to vectors it should modifiy the original reference
  //overloaded
  void syncJobNumbers(vector<pair<string,vector<float>>>& vec1,vector<pair<string,vector<float>>>& vec2);

  //function to parse simc hist files replayed by jboyd. Need to verify that this is still valid.
  double searchSimcHistFile(const TString &pattern, const TString &filename);

  //Used with regular sim files, not from J. Boyd.
  //dir1 is path to .csv , dir2 is path to .root , partialName is the search word, vec1 stores the root file absolute paths, csvData is a vector to store the CSV info and the root file path
  void SyncFilesCSV(TString dir1, TString dir2,TString partialName,vector<string>& vec1,vector<pair<string,vector<float>>>& csvData);

  //Function for stability analysis initially. Calculate the mean of the vector of doubles
  double calculateMean(vector<double> myVec);

  //Function for stability analysis initially. Calculate the Std Dev of the vector of doubles
  double calculateStDev(vector<double> myVec);

  //Function for stability analysis initially. Calculate the weighted mean of the vector of doubles
  double calculateWeightMean(vector<double> myVec, vector<double> myVec_uncert);

  //Function for stability analysis initially. Calculate the weighted Std Dev of the vector of doubles
  double calculateWeightStDev(vector<double> myVec, vector<double> myVec_uncert);

  //Function for stability analysis initially. Calculate the pull of the distribution of the vector of doubles
  double calculatePull(vector<double> myVec, vector<double> myVec_uncert);

  //Function to customize TGraphErrors
  void customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor, double markerSize,string graphTitle, string xAxisLabel, string yAxisLabel, double TitleOffsetX, double TitleOffsetY, double LabelOffsetX, double LabelOffsetY);

  //Function that puts TH2D minimally on Canvas
  TCanvas* printTH2D(TH2D* myHisto, TString title);

  //Function to set a time from a string
  TDatime SetTime(string time_str);
}
#endif
