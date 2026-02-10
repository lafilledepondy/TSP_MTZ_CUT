#pragma once
#include <memory>
#include <algorithm>
#include "gurobi_c++.h"
#include "ATSP_Data.hpp"

class ATSP_CUT : public GRBCallback
{
private:
    ATSPDataC data;
    std::unique_ptr<GRBEnv> env;
    std::unique_ptr<GRBModel> model;
    vector<vector<GRBVar>> x;
    int status;

    void findSubtour_S(const int sol, std::vector<int> &S);

public:
    void setterX(vector<vector<GRBVar>> &x) { this->x = x; }
    void setterStatus(int status) { this->status = status; }
    GRBModel *getterModel() { return this->model.get(); }
    vector<vector<GRBVar>> &getterX() { return this->x; }
    int getterStatus() { return this->status; }

    ATSP_CUT(ATSPDataC data);
    void solve();
    void printSolution();
};

class ATSP_CUT_Callback : public GRBCallback
{
private:
    int n;
    vector<vector<GRBVar>> &x;

public:
    ATSP_CUT_Callback(int n, vector<vector<GRBVar>> &x)
        : n(n), x(x) {}

protected:
    void callback()
    {
        try
        {
            // Séparation uniquement sur solutions entières
            if (where == GRB_CB_MIPSOL)
            {

                vector<int> succ(n, -1);

                // Reconstruction du successeur de chaque sommet
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        if (i != j && getSolution(x[i][j]) > 0.5)
                            succ[i] = j;

                vector<bool> visited(n, false);

                // Recherche de sous-tours
                for (int start = 0; start < n; ++start)
                {
                    if (visited[start])
                        continue;

                    vector<int> S;
                    int cur = start;

                    while (!visited[cur])
                    {
                        visited[cur] = true;
                        S.push_back(cur);
                        cur = succ[cur];
                    }

                    // Sous-tour strict -> contrainte violée
                    if ((int)S.size() < n)
                    {
                        GRBLinExpr cut = 0;
                        for (int i = 0; i < n; ++i)
                            for (int j : S)
                                if (find(S.begin(), S.end(), i) == S.end())
                                    cut += x[i][j];

                        addLazy(cut >= 1);
                        return; // une coupe suffit
                    }
                }
            }
        }
        catch (GRBException e)
        {
            cout << "Erreur callback : " << e.getMessage() << endl;
        }
    }
};
