#include "ATSP_MTZ.hpp"

ATSP_MTZ::ATSP_MTZ(ATSPDataC data) : data(data) {}

void ATSP_MTZ::solve()
{
    try
    {
        env = std::make_unique<GRBEnv>(true);
        env->set("LogFile", "atsp_mtz.log");
        env->start();

        model = std::make_unique<GRBModel>(*env);
        GRBModel &modelRef = *model;

        // Variables
        x = vector<vector<GRBVar>>(data.size, vector<GRBVar>(data.size));
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
                    xRef[i][j] = modelRef.addVar(0.0, 1.0, data.distances[i][j], GRB_BINARY, "x(" + to_string(i) + "," + to_string(j) + ")");
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
                    obj += data.distances[i][j] * xRef[i][j];
                }
            }
        }
        modelRef.setObjective(obj, GRB_MINIMIZE);

        // Contraintes
        // sum_j j!=i x[i][j] == 1 for all i in N
        for (int i = 0; i < data.size; ++i)
        {
            GRBLinExpr gaucheMembre = 0; //
            for (int j = 0; j < data.size; ++j)
            {
                if (i != j)
                {
                    gaucheMembre += xRef[i][j];
                }
            }
            modelRef.addConstr(gaucheMembre == 1, "out(" + to_string(i) + ")");
        }

        // sum_j j!=i x[j][i] == 1 for all i in N
        for (int i = 0; i < data.size; ++i)
        {
            GRBLinExpr gaucheMembre = 0; //
            for (int j = 0; j < data.size; ++j)
            {
                if (i != j)
                {
                    gaucheMembre += xRef[j][i];
                }
            }
            modelRef.addConstr(gaucheMembre == 1, "in(" + to_string(i) + ")");
        }

        // u_0 = 0
        modelRef.addConstr(u[0] == 0, "u(0)");

        // u_j >= u_i + 1 - (n-1)(1-x[i][j]) for all i,j in N, i!=j, j!=0
        for (int i = 1; i < data.size; ++i)
        {
            for (int j = 1; j < data.size; ++j)
            {
                if (i != j && j != 0)
                {
                    modelRef.addConstr(u[j] >= u[i] + 1 - (data.size - 1) * (1 - xRef[i][j]), "subtour(" + to_string(i) + "," + to_string(j) + ")");
                }
            }
        }

        modelRef.set(GRB_DoubleParam_TimeLimit, 600.0); //< définition du temps limite (en secondes)
        modelRef.set(GRB_IntParam_Threads, 1);          //< définition du nombre de threads pouvant être utilisé
        modelRef.write("model.lp");                     //< écriture du modèle PLNE dans le fichier donné en paramètre (optionnel)
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

void ATSP_MTZ::printSolution()
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
                    cout << "x(" << i << "," << j << ") = " << xRef[i][j].get(GRB_DoubleAttr_X) << endl;
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