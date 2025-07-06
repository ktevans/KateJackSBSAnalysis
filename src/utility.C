#include "../include/utility.h"
#include <iostream>
#include <filesystem>
#include <utility>
#include "TMath.h"
#include <regex>


//Author Ezekiel Wertz
//Implementation for useful task functions

namespace utility{
  //function to determine if what is found was a number or a string
  bool check_number(const char *myChar){
   bool mybool=false;
  	for(int i=0; i< strlen(myChar); i++){
		//So we ever find not a digit. Immediately return mybool=false. So we have a string
		//cout << mybool << endl;
		if(isdigit(myChar[i])==false){
		return mybool;
		}
  	}
  //If we get through the for loop we never found not a digit. Make it true and return
  mybool=true;
  return mybool;
  }//end check_number

  //convert from degrees to radians
  double DegToRad(double myNum){
  return (myNum * (TMath::Pi()/180.0));
  }
  //convert from radians to degrees
  double RadToDeg(double datRad){
  return ((datRad * 180.0)/(TMath::Pi()));
  }

  //convert an int to a TString. Just a helper function
  TString intToTString(int datInt){
  stringstream ss;
  ss << datInt;
  std::string datInt_temp;
  ss >> datInt_temp;
  TString datInt_string(datInt_temp.c_str());
  //cout << datInt_string << endl;
  return datInt_string;
  }

  //convert a double to a TString. Just a helper function
  TString doubToTString(double datDoub, int prec){
  stringstream ss;
  ss << setprecision(prec) << datDoub;
  std:string datDoub_temp;
  ss >> datDoub_temp;
  TString datDoub_string(datDoub_temp.c_str());
  return datDoub_string;
  }

  //This is a helper function that will tell you the default directory for all analysis output
  TString getOutputDir(){
  return output_directory;
  }

  //Helper function to make output file name for parsed data root file
  TString makeOutputFileNameParse(TString exp,TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/KJ_parsed_%s_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function to make output file name for parse mc root file
  TString makeOutputFileName_MCParse(TString exp, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/KJ_MC_pn_parsed_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),Kin.Data(),target.Data(),SBS_field);
  return outfile;
  }

  //Helper function to make output file name for data and mc comparison
  TString makeOutputFileName_DataMCComp(TString exp, TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/KJ_DataMC_Compare_%s_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  return outfile;
  }

  //Helper function to make report file name for data and mc comparison
  TString makeReportFileName_DataMCComp(TString exp, TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/KJ_DataMC_Compare_%s_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  return outfile;
  }

  //Helper function to make output file name for for stability studies
  TString makeOutputFileName_Stability(TString exp, TString pass, TString Kin, int SBS_field,TString target, int slice_mode, TString left_right){

  //We only want left_right for mode 1

  TString outfile;
  if(slice_mode == 1){
  outfile = Form("%s/KJ_Stability_%s_%s_%s_%s_%i_mode_%i_%s.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field,slice_mode,left_right.Data());
  }else{
  outfile = Form("%s/KJ_Stability_%s_%s_%s_%s_%i_mode_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field,slice_mode);
  }
  return outfile;
  }

  //Helper function to make report file name for for stability studies
  TString makeReportFileName_Stability(TString exp, TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/KJ_Stability_%s_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  return outfile;
  }

  //Helper function to make output file name for yield or ratio information
  TString makeYieldReportFileName(TString exp,TString pass,  TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/yieldRep_KJ_%s_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(), pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function to make output file name for HCal Efficiency from data
  TString makeOutputFileNameHCalEffParse(TString exp,TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/KJ_HCalEff_parsed_%s_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function to make output file name for HCal Efficiency from MC
  TString makeOutputFileName_HCalEffMC(TString exp, TString Kin){
  TString outfile = Form("%s/KJ_HCalEff_MC_%s_%s",(getOutputDir()).Data(),exp.Data(),Kin.Data());

  return outfile;
  }

  //Helper function to make output file name for HCal Efficiency for uniformity
  TString makeOutputFileNameHCalEffUniformity(TString exp,TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/KJ_HCalEff_uniformity_%s_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  TString makeOutputFileNameHCalEffMapComp(TString exp,TString pass1, TString pass2, TString Kin1, TString Kin2, int SBS_field1, int SBS_field2,TString target1,TString target2){
  TString outfile;

  if(SBS_field1 < 0 || SBS_field2 < 0){
  outfile = Form("%s/KJ_HCalEff_map_Compare_%s_%s_%s_%s_%s_%s_%s.root",(getOutputDir()).Data(),exp.Data(),pass1.Data(),pass2.Data(),Kin1.Data(),Kin2.Data(),target1.Data(),target2.Data());
  }else{
  outfile = Form("%s/KJ_HCalEff_map_Compare_%s_%s_%s_%s_%s_%s_%s_%i_%i.root",(getOutputDir()).Data(),exp.Data(),pass1.Data(),pass2.Data(),Kin1.Data(),Kin2.Data(),target1.Data(),target2.Data(),SBS_field1,SBS_field2);
  }
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function to make output file name for HCal Efficiency map analysis
  TString makeOutputFileName_HCalEffMap(TString exp,TString pass, TString Kin,TString target){
  TString outfile = Form("%s/KJ_HCalEff_map_%s_%s_%s_%s.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data());
  //cout << outfile << endl;
  return outfile;

  }

  //Helper function to make output file for HCal Efficiency Bootstrap
  TString makeOutputFileName_HCalBS(TString exp,TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event){
  TString outfile = Form("%s/KJ_HCalEff_bootstrap_%s_%s_%s_%s_%i_event_%i_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field,first_event,last_event);
  //cout << outfile << endl;
  return outfile;
  }
  //Helper function for making file name for HCal Efficiency report
  TString makeReportFileName_HCalEff(TString exp,TString pass,  TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/efficiencyRep_KJ_%s_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(), pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function for making file name for HCal Efficiency report with boot strap
  TString makeReportFileName_HCalEffBS(TString exp,TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event){
  TString outfile = Form("%s/efficiencyRep_KJ_%s_%s_%s_%s_%i_event_%i_%i.txt",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field,first_event,last_event);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function for output file name for scale field study for MC
  TString makeOutputFileName_MCsf(TString exp, TString Kin, int SBS_field,double sf,TString target){
  TString outfile = Form("%s/MC_KJ_%s_%s_%s_%i_sf%.1f.root",(getOutputDir()).Data(),exp.Data(),Kin.Data(),target.Data(),SBS_field,sf);
  //cout << outfile << endl;
   return outfile;
  }

  //Used for J. Boyd Simulation files
  //Helper function to find MC histogram files. Right now this only supports John Boyd simulation files. But in future might support others.
  vector<string> findHistFiles(TString replay_type,TString histDirectory,TString partialName){
  vector<string> myFiles;
  const string& histDir_str = histDirectory.Data();
  const string& partialName_str = partialName.Data();
	if(replay_type == "jboyd"){
		//populate the vector with all matching hist files
		for(const auto& entry : filesystem::directory_iterator(histDir_str)){
			//check that the file is a regular file and that the partial name is found somewhere in the filename
			if(entry.is_regular_file() && (entry.path().filename().string().find(partialName_str) != string::npos) ){
			//if we find one that should be a match put it on the vector
			myFiles.push_back(entry.path().string());
			}
		}
	}else{
	cout << "Error: MC file type not supported. Replay Type: " << replay_type << " This needs fixing!" << endl;
	}
  return myFiles;
  }

  //Used for J. Boyd Sim files
  //Helper function to match all the MC root files with the already found MC hist files. This function is not ideal, it modifies both of the vectors by the time the function is done. Though it does not have a return type
  void matchMCFiles(TString replay_type,vector<string>& histFiles,vector<string>& rootFiles, TString rootDirectory){
	//vector to track file indices that don't match
	vector<int> unmatchedIndices;
	int matches = 0;
	const string& rootDir_str = rootDirectory.Data();
	if(replay_type == "jboyd"){
		//loop over the values in the hist file vector
		for(size_t i = 0; i < histFiles.size(); ++i){
		cout << "Looking for hist to root dir  matches, " << i << "/" << histFiles.size() << " total matches " << matches << "\r";
       		cout.flush();
		filesystem::path histPath = histFiles[i];
		string hist_filename = histPath.filename().string();

		//modfiy the filename for the expect root file name
		string expect_root_filename = "replayed_digitized_" + hist_filename;
		expect_root_filename = expect_root_filename.substr(0, expect_root_filename.find_last_of('.')) + ".root";

		//boolean to track if we find match
		bool matchFound = false;

			//iterate over the root file directory and see if there is a match
			for(const auto& entry : filesystem::directory_iterator(rootDir_str)){
				//check if there is a match
				if(entry.is_regular_file() && (entry.path().filename() == expect_root_filename)){
				//if we did put the file name on the root file vector and update other parameters
				rootFiles.push_back(entry.path().string());
				matchFound = true;
				matches++;
				break; //Since we found the match we can stop searching for this iteration
				}//conditional
			}//inner for loop

			if(!matchFound){
			//if match not found indicate that index needs removal
			unmatchedIndices.push_back(i);
			}
		}// end outer for loop

		//Remove any unmatched indices
		//Need to iterate in reverse because of implementation of vectors and to get correct indices
		for(auto j = unmatchedIndices.rbegin(); j != unmatchedIndices.rend() ; ++j){

		histFiles.erase(histFiles.begin() + *j);
		}

	}else{
        cout << "Error: MC file type not supported during matching. Replay Type: " << replay_type << " This needs fixing!" << endl;
        }
  //There is no return here so this is a little confusing but when the function terminates the root file vector should be populate and the non matching hist vector elements should be removed. So this function modifies both vectors to have matching information, since we are sending references to the original vector.
  }//end matching function

  //Used for J. Boyd Sim files
  //Function that searches through two vectors of strings and removes entries without a matching pair
  //since we are sending references to vectors it should modifiy the original reference
  void syncJobNumbers(vector<string>& proton_vec,vector<string>& neutron_vec){
  regex jobNumberRegex("job([0-9]+)");
  std::smatch match;
  std::unordered_set<string> jobNumbers_pro, jobNumbers_neu;

  	//Get job numbers for proton files
  	for(const auto& str: proton_vec){
  		//search for files with the job numbers
		if(std::regex_search(str,match,jobNumberRegex)){
		//store references if we find one
		jobNumbers_pro.insert(match[0]);
		}
  	}
  	//Get job numbers for nuetron files
  	for(const auto& str: neutron_vec){
        	//search for files with the job numbers
        	if(std::regex_search(str,match,jobNumberRegex)){
        	//store references if we find one
        	jobNumbers_neu.insert(match[0]);
        	}
  	}
	//remove unpaired entries from proton vec
	//this codes a little odd, unfamiliar library.
	//I suspect it is doing something like searching from the begin to end of the vector and removing the string. But the string is made on the fly and is some how checking if one string found in the first vector is also in the second vector
	proton_vec.erase(std::remove_if(proton_vec.begin(),proton_vec.end(),[&](const string& str){
			if(std::regex_search(str,match,jobNumberRegex) && jobNumbers_neu.find(match[0]) == jobNumbers_neu.end() ){
			return true;
			}
		return false;
		}),proton_vec.end());
	//remove unpaired entries from neutron vec
        neutron_vec.erase(std::remove_if(neutron_vec.begin(),neutron_vec.end(),[&](const string& str){
                        if(std::regex_search(str,match,jobNumberRegex) && jobNumbers_pro.find(match[0]) == jobNumbers_pro.end() ){
                        return true;
                        }
                return false;
                }),neutron_vec.end());

  }// end sync job function

  //Used with regular sim files form JLab-HPC, not from J. Boyd.
  //Function that searches through two vectors of pair of string and float vector and removes entries without a matching pair
  //since we are sending references to vectors it should modifiy the original reference
  //overloaded
  void syncJobNumbers(vector<pair<string,vector<float>>>& vec1,vector<pair<string,vector<float>>>& vec2){
  //Extract Job Ids from vec1
  std::unordered_set<float> jobIds1;
	for(const auto& item: vec1){
		if(!item.second.empty()){
		jobIds1.insert(item.second[0]);
		}
	}
  //Extract JobIds from vec2
  std::unordered_set<float> jobIds2;
  	for(const auto& item: vec2){
                if(!item.second.empty()){
                jobIds2.insert(item.second[0]);
                }
        }
  //Remove unpaired entries from vec1
  vec1.erase(std::remove_if(vec1.begin(),vec1.end(), [&](const pair<string,vector<float>>& item){return item.second.empty() || jobIds2.find(item.second[0]) == jobIds2.end();}), vec1.end());

  //Remove unpaired entries from vec2
  vec2.erase(std::remove_if(vec2.begin(),vec2.end(), [&](const pair<string,vector<float>>& item){return item.second.empty() || jobIds1.find(item.second[0]) == jobIds1.end();}), vec2.end());

  }//end sync Jobs function

  //Used for J.Boyd Sim files
  //function to parse simc hist files replayed by jboyd. Need to verify that this is still valid.
  double searchSimcHistFile(const TString &pattern, const TString &filename){
	ifstream inputFile(filename.Data());//Try to open the file with the corresponding name of interest

  	//make sure the file exists. If it does not return
	if(inputFile.fail()){
	cerr << "Error: There was a problem with the input hist file!" << endl;
	cerr << "Filename: " << filename << endl;
	return -99; //This value does not matter as long as it is nonsensical for the other parameters.
	}

  TString myLine;
  double foundValue = -1; //Default to a value that would not make sense physically for any parameter
  TString previousLine; //Keep track of the previous line for exceptions
  	//search through each line of the file for the designated pattern
	while(myLine.ReadLine(inputFile)){
		//If we found the pattern in the line
		if(myLine.Contains(pattern)){
		//handle the case of luminosity or other patterns with scientific notation
			if(myLine.EndsWith("ub^-1") || myLine.EndsWith("geV^2")){
			//Extract the numerical value
			TString delim("=");
			TObjArray* tokens = myLine.Tokenize(delim);
				//check that we found enough elements
				if(tokens->GetEntries() >= 2){
				TObjString* value = dynamic_cast<TObjString*>(tokens->At(1));
					if(value){
					TString valueString = value->String();
					//convert to double
					foundValue = valueString.Atof();
					delete tokens; //cleanup
					return foundValue;
					}
				}
			delete tokens; //cleanup if for some reason we don't find a value this way
			}else{
			//For other nonspecific patterns
			TString valueString;
			TString delim("=");
			TObjArray* tokens = myLine.Tokenize(delim);
				if(tokens->GetEntries() >= 2){
				TObjString* value = dynamic_cast<TObjString*>(tokens->At(1));
					if(value){
					valueString = value->String();
					}
				}
				if(valueString.IsFloat()){
				//convert to double
				foundValue = valueString.Atof();
                                delete tokens; //cleanup
                                return foundValue;
				}
			delete tokens; //cleanup if for some reason we don't find a value this way
			}
		}//end conditional on pattern line check
	previousLine = myLine;
	}//end of while loop
  //Handle luminosity and other execptions when an appropriate line is not found
  	if((pattern == "luminosity" && previousLine.EndsWith("ub^-1")) || (pattern == "genvol")){
	TString delim("=");
	TObjArray* tokens = myLine.Tokenize(delim);
		if(tokens->GetEntries() >= 2){
                TObjString* value = dynamic_cast<TObjString*>(tokens->At(1));
			if(value){
                        TString valueString = value->String();
                       	foundValue = valueString.Atof();
			}

		}
	delete tokens; //cleanup
	}
  inputFile.close();

  return foundValue;
  }//end of hist search function

  //Used with regular sim files form JLab-HPC, not from J. Boyd.
  //dir1 is path to .csv , dir2 is path to .root , partialName is the search word, vec1 stores the root file absolute paths, csvData is a vector to store the CSV info and the root file path
  void SyncFilesCSV(TString dir1, TString dir2,TString partialName,vector<string>& vec1,vector<pair<string,vector<float>>>& csvData){
  const string& dir1_str = dir1.Data();
  const string& dir2_str = dir2.Data();
  const string& partialName_str = partialName.Data();

  //cout << dir1_str << " " << dir2_str << " " << partialName_str << endl;

  //populate vec1 with digitized and replayed root files
  string replay_partialName = "replayed_" + partialName_str;
  //cout << replay_partialName << endl;
  	//loop over every entry in dir2
  	for(const auto& entry : filesystem::directory_iterator(dir2_str)){
		//check that the entry is a regular file, that the file contains the partialName info, and is a root file
		if(entry.is_regular_file() && entry.path().filename().string().find(replay_partialName) != std::string::npos && entry.path().extension() == ".root"){
		//cout << entry.path().string() << endl;
		//if it meets the conditional keep it
		vec1.push_back(entry.path().string());
		}
	}//end loop
  //Read CSV file
  filesystem::path csvPath;
  bool csvFound = false;

  	//loop over every entry in dir1 to try and find the csv file
  	for(const auto& entry : filesystem::directory_iterator(dir1_str)){
		//check that the entry is a regular file, that the file contains the partialName info, and is a csv file
		if(entry.is_regular_file() && entry.path().filename().string().find(partialName_str) != std::string::npos && entry.path().filename().string().find("_summary.csv") != std::string::npos){
		//if we find it save the path and break the loop
		//cout << entry.path().string() << endl;
		csvPath = entry.path();
		csvFound = true;
		break;
		}

	}//end loop
  	//Send an error if we did not find any CSV file
  	if(!csvFound){
	cerr << "CSV file matching " << partialName_str << " not found in " << dir1_str << endl;
	}
  //open an ifstream so we can read info from the CSV file
  ifstream csvFile(csvPath);
  string line;
  	//loop over the lines in the file to find the information we care about
  	while(getline(csvFile,line)){
		//skip the header or any empty lines
		if(line.empty() || line[0] == 'j'){
		continue;
		}
	std::istringstream ss(line);
	float jobid;
	//First column is jobid, stored as float
	ss >> jobid;
	char comma;
	//Initialize a vector with jobid as the first element
	vector<float> rowData = {jobid};
	float data;
		//loop over the line and get all the information stored between the commas
		while(ss >> comma >> data){
		rowData.push_back(data);
		}
	//store data in csvData vector. Path is a placeholder for now
	csvData.push_back({"",rowData});
	}//end outer loop

	//Attempt to match .root files in vec1 with CSV data based on job ID
	for(auto& dataPair : csvData){
	regex jobIdPattern("job_([0-9]+)\\.root$");
		for(const auto& path : vec1 ){
		std::smatch matches;
			if(std::regex_search(path,matches,jobIdPattern) && matches.size()>1){
			int jobId = stoi(matches[1].str());
				if(jobId == static_cast<int>(dataPair.second[0])){
				//match found update path in csvData
				dataPair.first = path;
				//cout << path << endl;
				//stop searching after a match
				break;
				}
			}
		}
	}
  //conditionally remove entries from csvData if no match is found
  csvData.erase(std::remove_if(csvData.begin(),csvData.end(),[](const auto& pair){return pair.first.empty(); }),csvData.end());

  //Everything is stored in references to the orignal structures so we do not have to return anything.
  }//End CSV function

  //Function for stability analysis initially. Calculate the mean of the vector of doubles
  double calculateMean(vector<double> myVec){
  double sum = 0;

  //loop over the entries of the vector
  int vec_length = myVec.size();
  for(int i=0; i<vec_length; i++){
	sum += myVec[i];
	}
  double mean = sum/vec_length;
  return mean;
  }

  //Function for stability analysis intitially. Calculate the Std Dev of the vector of doubles
  double calculateStDev(vector<double> myVec){
  double mean = calculateMean(myVec);
  double numerator = 0;
  int vec_length = myVec.size();

  	//loop over the entries of the vector
  	for(int i=0; i<vec_length; i++){
        numerator += pow((myVec[i] - mean),2);
        }

  double std_dev = sqrt(numerator/(vec_length-1));
  return std_dev;
  }

  //Function for stability analysis initially. Calculate the weight mean of the vector of doubles
  double calculateWeightMean(vector<double> myVec, vector<double> myVec_uncert){
  double numerator = 0;
  double sum_weights = 0;
  int vec_length = myVec.size();
  int vec_uncert_length = myVec_uncert.size();

  	if(vec_length != vec_uncert_length){
	cout << "The lengths of the vectors while calculating weighted mean are not the same. Figure it out!" << endl;
	}

	//loop over the entries of the vector
	for(int i=0; i< vec_length; i++){
	double value = myVec[i];
	double weight = 1/pow(myVec_uncert[i],2);

	numerator += weight*value;
	sum_weights += weight;
	}
  double weight_mean = numerator/sum_weights;
  return weight_mean;
  }


  //Function for stability analysis initially. Calculate the weighted Std Dev of the vector of doubles
  double calculateWeightStDev(vector<double> myVec, vector<double> myVec_uncert){
  double numerator = 0;
  double sum_weights = 0;
  int vec_length = myVec.size();
  int vec_uncert_length = myVec_uncert.size();
  double weight_mean = calculateWeightMean(myVec,myVec_uncert);

        if(vec_length != vec_uncert_length){
        cout << "The lengths of the vectors while calculating weighted mean are not the same. Figure it out!" << endl;
        }

  	//loop over the entries of the vector
  	for(int i=0; i<vec_length; i++){
	double value = myVec[i];
        double weight = 1/pow(myVec_uncert[i],2);
	numerator += weight*pow((value-weight_mean),2);
	sum_weights += weight;
	}

  double std_dev = sqrt(numerator/sum_weights);
  return std_dev;
  }

  //Function for stability analysis initially. Calculate the pull of the distribution of the vector of doubles
  double calculatePull(vector<double> myVec, vector<double> myVec_uncert){
  int vec_length = myVec.size();
  int vec_uncert_length = myVec_uncert.size();
  double mean = calculateMean(myVec);
  double pull_sum = 0;

  	//loop over the entries
	for(int i=0; i<vec_length;i++){
	double pull = (mean - myVec[i])/myVec_uncert[i];
	pull_sum += pull;
	}
  double pull_avg = pull_sum/vec_length;
  return pull_avg;
  }

  void customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor, double markerSize,
                                 string graphTitle ="", string xAxisLabel="", string yAxisLabel="",
                                 double TitleOffsetX = 1.4, double TitleOffsetY = 2,
                                 double LabelOffsetX = 0.01, double LabelOffsetY = 0.01){
  // Set marker style, color, and size
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(markerColor);
  graph->SetMarkerSize(markerSize);

  // Set graph title and axis labels
  graph->SetTitle(graphTitle.c_str());
  graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
  graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

  // Adjust axis title offsets to provide more space
  graph->GetXaxis()->SetTitleOffset(TitleOffsetX); // Adjust as needed
  graph->GetYaxis()->SetTitleOffset(TitleOffsetY); // Adjust as needed


  // Adjust axis label offsets to provide more space
  graph->GetXaxis()->SetLabelOffset(LabelOffsetX); // Adjust as needed
  graph->GetYaxis()->SetLabelOffset(LabelOffsetY); // Adjust as needed

  }

  TCanvas* printTH2D(TH2D* myHisto, TString title){

  TCanvas* daCanvas = new TCanvas(myHisto->GetName(),myHisto->GetName(),1600,1200);

  TH2D* myHisto_Clone = (TH2D*) myHisto->Clone();

  TString Xaxis_name = (myHisto_Clone->GetXaxis())->GetTitle();
  TString Yaxis_name = (myHisto_Clone->GetYaxis())->GetTitle();
  TString new_Title = title + "_" + Yaxis_name + "__" + Xaxis_name;
  myHisto_Clone->SetName(new_Title.Data());
  myHisto_Clone->SetTitle(new_Title.Data());
  myHisto_Clone->SetStats(0);
  myHisto_Clone->Draw("COLZ");

  daCanvas->Update();

  return daCanvas;
  }

  TDatime SetTime(string time_str){

    TDatime time(time_str.c_str());
    time.Set(time.GetYear(),time.GetMonth(),time.GetDay(),time.GetHour(),time.GetMinute(),0);
      
    return time;
  }

} //end namespace
