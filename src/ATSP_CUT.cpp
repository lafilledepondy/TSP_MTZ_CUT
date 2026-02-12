#include "ATSP_CUT.hpp"

using namespace std;

ATSP_CUT::ATSP_CUT(ATSPDataC data) : data(data) {}

void ATSP_CUT::solve()
{
    try
    {
        env = std::make_unique<GRBEnv>(true);
        env->set("LogFile", "atsp_cut.log");
        env->start();

        model = std::make_unique<GRBModel>(*env);
        GRBModel &modelRef = *model;

        // Variables
        x.assign(data.size, vector<GRBVar>(data.size));
        vector<vector<GRBVar>> &xRef = x;
        vector<GRBVar> u(data.size);

        for (int i = 0; i < data.size; ++i)
        {
            // 1 <= u_i <= n-1 for i in N\{0}
            if (i != 0)
            {
                u[i] = modelRef.addVar(1.0, data.size - 1, 0.0, GRB_INTEGER, "u(" + to_string(i) + ")");
            }
            else
            {
                u[i] = modelRef.addVar(0.0, 0.0, 0.0, GRB_INTEGER, "u(" + to_string(i) + ")");
            }
            for (int j = 0; j < data.size; ++j)
            {
                if (i != j)
                {
                    x[i][j] = modelRef.addVar(0.0, 1.0, data.distances[i][j], GRB_BINARY, "x(" + to_string(i) + "," + to_string(j) + ")");
                }
            }
        }

        // Fonction objective
        GRBLinExpr obj = 0;
        for (int i = 0; i < data.size; ++i)
        {
            for (int j = 0; j < data.size; ++j)
            {
                if (i != j)
                {
                    obj += data.distances[i][j] * x[i][j];
                }
            }
        }
        modelRef.setObjective(obj, GRB_MINIMIZE);

        // Contraintes
        // sum_j j!=i x[i][j] == 1 for all i in N
        // sum_j j!=i x[j][i] == 1 for all i in N
        for (int i = 0; i < data.size; ++i)
        {
            GRBLinExpr out = 0, in = 0;
            for (int j = 0; j < data.size; ++j)
                if (i != j)
                {
                    out += x[i][j];
                    in += x[j][i];
                }
            modelRef.addConstr(out == 1);
            modelRef.addConstr(in == 1);
        }

        modelRef.set(GRB_IntParam_LazyConstraints, 1);

        std::unique_ptr<ATSP_CUT_Callback> cb;
        cb = std::unique_ptr<ATSP_CUT_Callback>(new ATSP_CUT_Callback(data.size, x));
        modelRef.setCallback(cb.get());

        modelRef.set(GRB_DoubleParam_TimeLimit, 180.0);
        modelRef.set(GRB_IntParam_Threads, 1);

        modelRef.write("model.lp");
        modelRef.optimize();
        setterStatus(modelRef.get(GRB_IntAttr_Status));
    }
    catch (GRBException e)
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...)
    {
        cout << "Exception during optimization" << endl;
    }
}

bool findSubtour_S(const vector<vector<double>> &sol, vector<int> &S)
{
    int n = static_cast<int>(sol.size());
    vector<bool> visited(n, false);

    for (int start = 0; start < n; ++start)
    {
        if (visited[start])
            continue;

        vector<int> cycle;
        int current = start;

        while (!visited[current])
        {
            visited[current] = true;
            cycle.push_back(current);

            // find successor of current
            bool foundNext = false;
            for (int j = 0; j < n; ++j)
            {
                if (current != j && sol[current][j] > 0.5)
                {
                    current = j;
                    foundNext = true;
                    break;
                }
            }

            if (!foundNext)
                break;
        }

        // If we closed a cycle and it’s not the full tour → subtour
        if (current == start && cycle.size() < n)
        {
            S = cycle;
            return true;
        }
    }

    S.clear();
    return false;
}

void ATSP_CUT::printSolution()
{
    if (!model)
    {
        cerr << "Fail! (Model not available; run solve() first)" << endl;
        return;
    }

    GRBModel &modelRef = *model;
    vector<vector<GRBVar>> &xRef = x;
    int status = getterStatus();

    if (status == GRB_OPTIMAL || (status == GRB_TIME_LIMIT && modelRef.get(GRB_IntAttr_SolCount) > 0))
    {
        // le solveur a calculé la solution optimale ou une solution réalisable
        //  (le temps limite a été atteint avant de pouvoir prouver l'optimalité)
        cout << "Succes! (Status: " << status << ")" << endl; //< (cf. documentation)
        // Affiche le temps de résolution
        cout << "Runtime : " << modelRef.get(GRB_DoubleAttr_Runtime) << " seconds" << endl;
        // Ecris la solution dans le fichier donné en paramètre (optionnel)
        modelRef.write("solution.sol");
        // Affiche la valeur de la solution
        cout << "Objective value = " << modelRef.get(GRB_DoubleAttr_ObjVal) << endl;
        for (size_t i = 0; i < data.size; ++i)
        {
            for (size_t j = 0; j < data.size; ++j)
            {
                if (i != j && xRef[i][j].get(GRB_DoubleAttr_X) > 0.5)
                {
                    cout << "x(" << i << ", " << j << ") = " << xRef[i][j].get(GRB_DoubleAttr_X) << endl;
                }
            }
        }
    }
    else
    {
        // le modèle est irréalisable (ou faux)
        // ou bien aucune solution n'a pu être calculé durant le temps limite imparti
        cerr << "Fail! (Status: " << status << ")" << endl; //< (cf. documentation)
    }
}
