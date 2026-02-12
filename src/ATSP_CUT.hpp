#pragma once
#include <memory>
#include <algorithm>
#include "gurobi_c++.h"
#include "ATSP_Data.hpp"
#include "hi_pr.hpp"

bool findSubtour_S(const std::vector<std::vector<double>> &sol, std::vector<int> &S);

class ATSP_CUT : public GRBCallback
{
private:
    ATSPDataC data;
    std::unique_ptr<GRBEnv> env;
    std::unique_ptr<GRBModel> model;
    vector<vector<GRBVar>> x;
    int status;
    int lazyCuts;
    int userCuts;

public:
    enum class SolveMode
    {
        IntegerMIP,
        FractionalLP
    };

private:
    SolveMode mode;

public:
    void setterX(vector<vector<GRBVar>> &x) { this->x = x; }
    void setterStatus(int status) { this->status = status; }
    GRBModel *getterModel() { return this->model.get(); }
    vector<vector<GRBVar>> &getterX() { return this->x; }
    int getterStatus() { return this->status; }
    int getLazyCuts() const { return lazyCuts; }
    int getUserCuts() const { return userCuts; }
    int getTotalCuts() const { return lazyCuts + userCuts; }
    SolveMode getMode() const { return mode; }

    ATSP_CUT(ATSPDataC data, SolveMode mode = SolveMode::IntegerMIP);
    void solve();
    void printSolution();
};

class ATSP_CUT_Callback : public GRBCallback
{
private:
    int n;
    vector<vector<GRBVar>> &x;
    int *lazyCuts;
    int *userCuts;

public:
    ATSP_CUT_Callback(int n, vector<vector<GRBVar>> &x, int *lazyCuts, int *userCuts)
        : n(n), x(x), lazyCuts(lazyCuts), userCuts(userCuts) {}

protected:
    void callback()
    {
        try
        {
            // Séparation uniquement sur solutions entières
            if (where == GRB_CB_MIPSOL)
            {
                vector<vector<double>> sol(n, vector<double>(n, 0.0));
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        if (i != j)
                            sol[i][j] = getSolution(x[i][j]);

                vector<int> S;
                if (findSubtour_S(sol, S))
                {
                    vector<bool> inS(n, false);
                    for (int v : S)
                        inS[v] = true;

                    GRBLinExpr cut = 0;
                    for (int i = 0; i < n; ++i)
                        if (!inS[i])
                            for (int j : S)
                                cut += x[i][j];

                    addLazy(cut >= 1);
                    if (lazyCuts)
                        (*lazyCuts)++;
                    return; // une coupe suffit
                }
            }

            if (where == GRB_CB_MIPNODE)
            {
                if (getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL)
                    return;

                vector<vector<double>> sol(n, vector<double>(n, 0.0));
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        if (i != j)
                            sol[i][j] = getNodeRel(x[i][j]);

                double **cap = new double *[n];
                for (int i = 0; i < n; ++i)
                {
                    cap[i] = new double[n];
                    for (int j = 0; j < n; ++j)
                        cap[i][j] = (i == j) ? 0.0 : sol[i][j];
                }

                for (int sink = 1; sink < n; ++sink)
                {
                    long *dist = new long[n];
                    double val = 0.0;
                    directed_min_cut(cap, n, 0, sink, val, dist);

                    if (val < 1.0 - 1e-6)
                    {
                        vector<int> S;
                        S.reserve(n);
                        for (int v = 0; v < n; ++v)
                            if (dist[v] <= n - 1)
                                S.push_back(v);

                        if (S.empty() || static_cast<int>(S.size()) == n)
                        {
                            delete[] dist;
                            continue;
                        }

                        vector<bool> inS(n, false);
                        for (int v : S)
                            inS[v] = true;

                        GRBLinExpr cut = 0;
                        for (int i = 0; i < n; ++i)
                            if (!inS[i])
                                for (int j : S)
                                    cut += x[i][j];

                        addCut(cut >= 1);
                        if (userCuts)
                            (*userCuts)++;
                        delete[] dist;
                        break;
                    }

                    delete[] dist;
                }

                for (int i = 0; i < n; ++i)
                    delete[] cap[i];
                delete[] cap;
            }
        }
        catch (GRBException e)
        {
            cout << "Erreur callback : " << e.getMessage() << endl;
        }
    }
};
