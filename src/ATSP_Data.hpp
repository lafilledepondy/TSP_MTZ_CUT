#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

class ATSPDataC{
public:
  int size;
  vector<vector<int>> distances;

  ATSPDataC(string filename){
    ifstream inf(filename);
    if (!inf){
      size = -1;
      return;
    }

    string temp;
    do{
      inf >> temp >> ws;
      if (temp == "DIMENSION:")
        inf >> size >> ws;
    } while (temp != "EDGE_WEIGHT_SECTION");
    distances.resize(size, vector<int>(size));
    for (size_t i = 0; i < size; ++i){
      for (size_t j = 0; j < size; ++j){
        inf >> distances[i][j] >> ws;
      }
    }
    inf.close();
  }

  void printData(){
    cout << "==================================" << endl;
    cout << "Size : " << size << endl;
    cout << "Distances :" << endl;
    for (size_t i = 0; i < size; ++i){
      for (size_t j = 0; j < size; ++j){
        cout << distances[i][j] << " ";
      }
      cout << endl;
    }
    cout << "==================================" << endl;
  }
};
