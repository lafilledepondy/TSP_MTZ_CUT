#include "ATSP_CUT.hpp"

using namespace std;

ATSP_CUT::ATSP_CUT(ATSPDataC data) : data(data) {}

void ATSP_CUT::solve()
{
    SolveOptions options;
    options.timeLimitSec = 600.0;
    options.threads = 1;
    options.logFile = "atsp_cut.log";
    solveWithOptions(options);
}

SolveMetrics ATSP_CUT::solveWithOptions(const SolveOptions &options)
{
    SolveMetrics metrics;
    try
    {
        env = std::make_unique<GRBEnv>(true);
        if (!options.logFile.empty())
        {
            env->set("LogFile", options.logFile);
        }
        env->start();

        model = std::make_unique<GRBModel>(*env);
        GRBModel &modelRef = *model;

        char xType = options.relax ? GRB_CONTINUOUS : GRB_BINARY;
        char uType = options.relax ? GRB_CONTINUOUS : GRB_INTEGER;

        // Variables
        x.assign(data.size, vector<GRBVar>(data.size));
        vector<vector<GRBVar>> &xRef = x;
        vector<GRBVar> u(data.size);

        for (int i = 0; i < data.size; ++i)
        {
            // 1 <= u_i <= n-1 for i in N\{0}
            if (i != 0)
            {
                u[i] = modelRef.addVar(1.0, data.size - 1, 0.0, uType, "u(" + to_string(i) + ")");
            }
            else
            {
                u[i] = modelRef.addVar(0.0, 0.0, 0.0, uType, "u(" + to_string(i) + ")");
            }
            for (int j = 0; j < data.size; ++j)
            {
                if (i != j)
                {
                    x[i][j] = modelRef.addVar(0.0, 1.0, data.distances[i][j], xType, "x(" + to_string(i) + "," + to_string(j) + ")");
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

        if (!options.relax)
        {
            modelRef.set(GRB_IntParam_LazyConstraints, 1);
        }

        std::unique_ptr<ATSP_CUT_Callback> cb;
        if (!options.relax)
        {
            cb = std::unique_ptr<ATSP_CUT_Callback>(new ATSP_CUT_Callback(data.size, x));
            modelRef.setCallback(cb.get());
        }

        modelRef.set(GRB_DoubleParam_TimeLimit, options.timeLimitSec);
        modelRef.set(GRB_IntParam_Threads, options.threads);
        if (options.writeModel)
        {
            modelRef.write("model.lp");
        }
        modelRef.optimize();
        setterStatus(modelRef.get(GRB_IntAttr_Status));

        metrics.status = modelRef.get(GRB_IntAttr_Status);
        metrics.solCount = modelRef.get(GRB_IntAttr_SolCount);
        metrics.runtime = modelRef.get(GRB_DoubleAttr_Runtime);
        metrics.nodeCount = modelRef.get(GRB_DoubleAttr_NodeCount);
        metrics.hasSolution = (metrics.solCount > 0);
        if (metrics.hasSolution)
        {
            metrics.objVal = modelRef.get(GRB_DoubleAttr_ObjVal);
        }
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

    return metrics;
}

void ATSP_CUT::findSubtour_S(const int sol, vector<int> &S)
{
    int n = data.size;
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
                if (current != j && x[current][j].get(GRB_DoubleAttr_X) > 0.5)
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
            return;
        }
    }

    S.clear();
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
