#pragma once
#include <vector>
#include "gurobi_c++.h"
#include "ATSP_Data.hpp"

class ATSP_CUT : public GRBCallback
{
private:
    ATSPDataC data;

    void findSubtour(const std::vector<std::vector<int>> &sol, std::vector<int> &S);

public:
    ATSP_CUT(ATSPDataC data);
    int solve();
    void printSolution(int status);

protected:
    void callback() override;
};
