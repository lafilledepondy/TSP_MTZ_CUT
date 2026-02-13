#include "ATSP_CUT.hpp"
#include <chrono>

using namespace std;

ATSP_CUT::ATSP_CUT(ATSPDataC data, SolveMode mode)
    : data(data), status(0), lazyCuts(0), userCuts(0), mode(mode) {}

static bool findFractionalCutS(const vector<vector<double>> &sol, vector<int> &S)
{
    int n = static_cast<int>(sol.size());

    // Construction du graphe capacite a partir de la solution fractionnaire x.
    // La capacite de l'arc (i,j) vaut x_ij, et 0 sur les boucles.
    double **cap = new double *[n];
    for (int i = 0; i < n; ++i)
    {
        cap[i] = new double[n];
        for (int j = 0; j < n; ++j)
            cap[i][j] = (i == j) ? 0.0 : sol[i][j];
    }

    // Separation de coupes de type "cut" sur solutions fractionnaires :
    // on cherche un min-cut dirige separant la source 0 d'un puits s.
    // Si la valeur du min-cut est < 1, la contrainte sum_{i notin S, j in S} x_ij >= 1 est violee.
    for (int sink = 1; sink < n; ++sink)
    {
        long *dist = new long[n];
        double val = 0.0;
        directed_min_cut(cap, n, 0, sink, val, dist);

        if (val < 1.0 - 1e-6)
        {
            // Reconstruction de l'ensemble S du cote source du min-cut.
            S.clear();
            S.reserve(n);
            for (int v = 0; v < n; ++v)
                if (dist[v] <= n - 1)
                    S.push_back(v);

            delete[] dist;
            for (int i = 0; i < n; ++i)
                delete[] cap[i];
            delete[] cap;

            // Si S est un sous-ensemble propre, on a trouve une coupe violee.
            if (!S.empty() && static_cast<int>(S.size()) < n)
                return true;

            continue;
        }

        delete[] dist;
    }

    for (int i = 0; i < n; ++i)
        delete[] cap[i];
    delete[] cap;

    S.clear();
    return false;
}

void ATSP_CUT::solve()
{
    try
    {
        // Reinitialise les compteurs de coupes.
        lazyCuts = 0;
        userCuts = 0;

        env = std::make_unique<GRBEnv>(true);
        env->set("LogFile", "atsp_cut.log");
        env->start();

        model = std::make_unique<GRBModel>(*env);
        GRBModel &modelRef = *model;

        // Variables : x_ij (arcs) et u_i (ordre MTZ).
        x.assign(data.size, vector<GRBVar>(data.size));
        vector<vector<GRBVar>> &xRef = x;
        vector<GRBVar> u(data.size);
        const char xType = (mode == SolveMode::FractionalLP) ? GRB_CONTINUOUS : GRB_BINARY;
        const char uType = (mode == SolveMode::FractionalLP) ? GRB_CONTINUOUS : GRB_INTEGER;

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

        // Fonction objective : minimiser la somme des couts des arcs choisis.
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

        // Contraintes de degre :
        // - un arc sortant par sommet
        // - un arc entrant par sommet
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

        modelRef.set(GRB_DoubleParam_TimeLimit, 180.0);
        modelRef.set(GRB_IntParam_Threads, 1);

        if (mode == SolveMode::IntegerMIP)
        {
            // Separation de sous-tours en mode MIP entier via coupes "lazy".
            modelRef.set(GRB_IntParam_LazyConstraints, 1);
            std::unique_ptr<ATSP_CUT_Callback> cb;
            cb = std::unique_ptr<ATSP_CUT_Callback>(new ATSP_CUT_Callback(data.size, x, &lazyCuts, &userCuts));
            modelRef.setCallback(cb.get());

            modelRef.write("model.lp");
            modelRef.optimize();
            setterStatus(modelRef.get(GRB_IntAttr_Status));
        }
        else
        {
            // Mode LP fractionnaire : boucle "solve -> separation -> ajout de coupe".
            auto start = std::chrono::steady_clock::now();
            double timeLimit = 180.0;

            while (true)
            {
                auto now = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now - start).count();
                double remaining = timeLimit - elapsed;
                if (remaining <= 0.0)
                    break;

                // Resolution du LP courant dans le temps restant.
                modelRef.set(GRB_DoubleParam_TimeLimit, remaining);
                modelRef.optimize();
                setterStatus(modelRef.get(GRB_IntAttr_Status));

                int status = getterStatus();
                if (status != GRB_OPTIMAL && status != GRB_TIME_LIMIT)
                    break;

                if (modelRef.get(GRB_IntAttr_SolCount) == 0)
                    break;

                // Extraction de la solution fractionnaire x pour la separation.
                vector<vector<double>> sol(data.size, vector<double>(data.size, 0.0));
                for (int i = 0; i < data.size; ++i)
                    for (int j = 0; j < data.size; ++j)
                        if (i != j)
                            sol[i][j] = x[i][j].get(GRB_DoubleAttr_X);

                vector<int> S;
                // Si aucune coupe violee n'est trouvee, on arrete la separation.
                if (!findFractionalCutS(sol, S))
                    break;

                vector<bool> inS(data.size, false);
                for (int v : S)
                    inS[v] = true;

                // Contrainte de coupe : au moins un arc entre V\S et S.
                GRBLinExpr cut = 0;
                for (int i = 0; i < data.size; ++i)
                    if (!inS[i])
                        for (int j : S)
                            cut += x[i][j];

                modelRef.addConstr(cut >= 1);
                userCuts++;

                if (getterStatus() == GRB_TIME_LIMIT)
                    break;
            }
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
}

bool findSubtour_S(const vector<vector<double>> &sol, vector<int> &S)
{
    int n = static_cast<int>(sol.size());
    vector<bool> visited(n, false);

    // Separation sur solution entiere : on cherche un cycle strict (sous-tour).
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

            // Trouver le successeur (arc sortant choisi depuis current).
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

        // Si on ferme un cycle strict (taille < n), on a un sous-tour.
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
