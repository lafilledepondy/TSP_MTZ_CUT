#include <string>
#include <sstream>
#include <memory>
#include "gurobi_c++.h"

using namespace std;

#include "ATSP_Data.hpp"

#pragma once

class ATSP_MTZ {

private:
    ATSPDataC data;
    std::unique_ptr<GRBEnv> env;
    std::unique_ptr<GRBModel> model;
    int status;

    vector<vector<GRBVar>> x; // x[i][j] == var decision arc i->j

public:
    // Constructeur
    ATSP_MTZ(ATSPDataC data);

    // Setters & Getters
    void setterX(vector<vector<GRBVar>> &x) { this->x = x; }
    void setterStatus(int status) { this->status = status; }
    GRBModel *getterModel() { return this->model.get(); }
    vector<vector<GRBVar>> &getterX() { return this->x; }
    int getterStatus() { return this->status; }

    void solve(); // build + solve model
    void printSolution();  // affiche sol
};
