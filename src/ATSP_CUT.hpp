#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include "gurobi_c++.h"

using namespace std;

#include "ATSP_Data.hpp"

#pragma once

class ATSP_CUT
{
private:
    ATSPDataC data;

public:
    ATSP_CUT(ATSPDataC data);

    void solve();
};