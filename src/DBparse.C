#include "../include/DBparse.h"

namespace DBparse {

  void DB_load(DBInfo &request){

    for(int ivar = 0; ivar < request.var_req.size(); ivar++){
      DBrequest var_info = request.var_req[ivar];
      bool found_var = false;

      TString DIR = DB_dir;

      auto it = DBFileMap.find(var_info.var_names);
      TString file = it->second;
      fstream file_csv; file_csv.open(DIR + file);
      
      // Here we check that the DB files are there
      if (it == DBFileMap.end()) { // This is not in the list of DB files
        if(!var_info.mandatory){
	  cout<<"DB variable NOT found but NOT mandatory: "<<var_info.var_names<<endl;
	  continue;
	}
	if(var_info.mandatory){
	  cout<<"!!!!WARNING, ERROR, DBparse::DB_load: DB variable NOT found and IS mandatory: "<<var_info.var_names<<endl;
	  exit(0);
	}
      }
      else if (!file_csv.is_open()) {// This is in the list but the file is not there
	if(!var_info.mandatory){
	  cout<<"DB file NOT found but NOT mandatory: "<<file<<endl;
	  continue;
	}
	if(var_info.mandatory){
	  cout<<"!!!!WARNING, ERROR, DBparse::DB_load: File "<<file<<" is not found"<<endl;
	  exit(0);
	}
      }
      else {// The file is found and everything is fine
	cout<<"DB file found: "<<file<<endl;
      }
      
      string line;
      int iline = 0;
      
      while (getline(file_csv, line)) {
	// Create a stringstream to parse the line
	stringstream ss(line);
	string cell;
	iline++;
	
	vector<string> val;

	// Split the line into cells using a comma as a delimiter
	while (getline(ss, cell, ',')) {
	  val.push_back(cell);  // Put one line into vectros
	}
	if(it->first == "He3 Polarization"){// He3 polarization DB
	  if(iline == 1) continue;
	  request.He3Pol.insert(std::make_pair(utility::SetTime(val[0]),stod(val[1])));
	}// end He3 Pol
      }      
    }
  }

}
