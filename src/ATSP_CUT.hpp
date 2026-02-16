#pragma once
#include <memory>
#include <algorithm>
#include "gurobi_c++.h"
#include "ATSP_Data.hpp"
#include "hi_pr.hpp"

// detecte 1 sous tour dans sol => remplit S si trouve
bool findSubtour_S(const std::vector<std::vector<double>> &sol, std::vector<int> &S);

// ======================================================================
// ============== CLASS ATSP_CUT :: GRBCALLBACK =========================
// ======================================================================
class ATSP_CUT : public GRBCallback{
private:
    ATSPDataC data;
    std::unique_ptr<GRBEnv> env;
    std::unique_ptr<GRBModel> model;
    int status;

    vector<vector<GRBVar>> x; // x[i][j] == var decision arc i->j

    int lazyCuts; // nb lazy cuts ajoutees
    int userCuts; // "  user   "     "

public:
    enum class SolveMode{
        IntegerMIP, // solve entier 
        FractionalLP // solve frac (relax LP)
    };

private:
    SolveMode mode; // mode courant (int or double)

public:
    // Setters & Getters
    void setterX(vector<vector<GRBVar>> &x) { this->x = x; }
    void setterStatus(int status) { this->status = status; }
    GRBModel *getterModel() { return this->model.get(); }
    vector<vector<GRBVar>> &getterX() { return this->x; }
    int getterStatus() { return this->status; }
    int getLazyCuts() const { return lazyCuts; }
    int getUserCuts() const { return userCuts; }
    int getTotalCuts() const { return lazyCuts + userCuts; } // total cuts (user + lazy)
    SolveMode getMode() const { return mode; }

    // Constructeur
    ATSP_CUT(ATSPDataC data, SolveMode mode = SolveMode::IntegerMIP);

    void solve(); // build + solve model
    void printSolution();  // affiche sol
};

// ======================================================================
// ============== CLASS ATSP_CUT_CALLBACK :: GRBCALLBACK ================
// ======================================================================
class ATSP_CUT_Callback : public GRBCallback{
private:
    int n; // taille instance
    vector<vector<GRBVar>> &x; // ref vars x
    
    int *lazyCuts; // ptr comptaur lazy
    int *userCuts; //  "     "     user

public:
    ATSP_CUT_Callback(int n, vector<vector<GRBVar>> &x, int *lazyCuts, int *userCuts)
        : n(n), x(x), lazyCuts(lazyCuts), userCuts(userCuts) {}

protected:
    void callback(){
        try{
        // ================= QUESTION 3 =================
        // sep contraintes (11) sol int            
            // si sol entiere trouvee
            if (where == GRB_CB_MIPSOL){
                
                // reconstruit matrice sol[i][j]
                vector<vector<double>> sol(n, vector<double>(n, 0.0));
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        if (i != j)
                            sol[i][j] = getSolution(x[i][j]);

                vector<int> S;

                // cherche sous tour
                if (findSubtour_S(sol, S))
                {
                    // construit indicateur inS
                    vector<bool> inS(n, false);
                    for (int v : S)
                        inS[v] = true;

                    GRBLinExpr cut = 0;

                    // cut == sum i notin S j in S x[i][j]
                    for (int i = 0; i < n; ++i)
                        if (!inS[i])
                            for (int j : S)
                                cut += x[i][j];

                    addLazy(cut >= 1); // ajoute lazy cut (contrainte (11))
                    if (lazyCuts)
                        (*lazyCuts)++; // +1 compteur
                    return;  // 1 coupe suffit
                }
            }

            if (where == GRB_CB_MIPSOL){
                cout << "SKIPPED CALLBACK ON MIPSOL" << endl;
            }

        // ================= QUESTION 5 =================
        // sep contraintes (11) sol frac via min cut
            if (where == GRB_CB_MIPNODE){
                // seulement si relax optimale
                if (getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL)
                    {return;} // sinon stop
            
                // reconstruit sol frac x[i][j]
                vector<vector<double>> sol(n, vector<double>(n, 0.0));
                for (int i = 0; i < n; ++i)
                    {for (int j = 0; j < n; ++j)
                        {if (i != j)
                            {sol[i][j] = getNodeRel(x[i][j]);}}}

                // cap[i][j] == sol[i][j]                            
                double **cap = new double *[n];
                for (int i = 0; i < n; ++i){
                    cap[i] = new double[n];
                    {for (int j = 0; j < n; ++j)
                        {cap[i][j] = (i == j) ? 0.0 : sol[i][j];}}
                }

                // test min cut 0 -> sink
                for (int sink = 1; sink < n; ++sink){
                    long *dist = new long[n]; // labels coupe
                    double val = 0.0; // valeur min cut

                    directed_min_cut(cap, n, 0, sink, val, dist); // calcule min cut

                    // si val < 1 => contrainte 11 violee
                    if (val < 1.0 - 1e-6){
                        vector<int> S; // ensemble cote source
                        S.reserve(n);

                        // construit S depuis dist
                        for (int v = 0; v < n; ++v)
                            {if (dist[v] <= n - 1)
                                {S.push_back(v);}}
                        
                        // si trivial => ignore
                        if (S.empty() || static_cast<int>(S.size()) == n){
                            delete[] dist; // sanitize
                            continue;
                        }

                        vector<bool> inS(n, false); // indicateur S
                        for (int v : S)
                            {inS[v] = true;}

                        // cut == sum i notin S j in S x[i][j]
                        GRBLinExpr cut = 0;
                        for (int i = 0; i < n; ++i)
                            {if (!inS[i])
                                {for (int j : S)
                                    {cut += x[i][j];}}}

                        addCut(cut >= 1); // ajoute user cut
                        if (userCuts)
                            {(*userCuts)++;} // add +1 to user cut compteur
                        delete[] dist; // sanitize
                        break; // 1 coupe suffit
                    }

                    delete[] dist; // sanitize
                }

                for (int i = 0; i < n; ++i)
                    {delete[] cap[i];} // sanitize
 
                delete[] cap; // sanitize
            }
        }
        catch (GRBException e)
        {
            cout << "Erreur callback : " << e.getMessage() << endl;
        }
    }
};
