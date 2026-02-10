#include "ATSP_CUT.hpp"

using namespace std;

ATSP_CUT::ATSP_CUT(ATSPDataC data) : data(data) {}

int ATSP_CUT::solve()
{
    try
    {
        GRBEnv env;
        GRBModel model(env);

        // Variables
        vector<vector<GRBVar>> x(data.size, vector<GRBVar>(data.size));
        vector<GRBVar> u(data.size);

        for (int i = 0; i < data.size; ++i)
        {
            // 1 <= u_i <= n-1 for i in N\{0}
            if (i != 0)
            {
                u[i] = model.addVar(1.0, data.size - 1, 0.0, GRB_INTEGER, "u(" + to_string(i) + ")");
            }
            else
            {
                u[i] = model.addVar(0.0, 0.0, 0.0, GRB_INTEGER, "u(" + to_string(i) + ")");
            }
            for (int j = 0; j < data.size; ++j)
            {
                if (i != j)
                {
                    x[i][j] = model.addVar(0.0, 1.0, data.distances[i][j], GRB_BINARY, "x(" + to_string(i) + "," + to_string(j) + ")");
                }
            }
        }

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j)
                    x[i][j] = model.addVar(0.0, 1.0, data.distances[i][j], GRB_BINARY);

        // Fonction objective
        GRBLinExpr obj = 0;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    obj += data.distances[i][j] * x[i][j];
                }
            }
        }
        model.setObjective(obj, GRB_MINIMIZE);

        // Contraintes
        // sum_j j!=i x[i][j] == 1 for all i in N
        // sum_j j!=i x[j][i] == 1 for all i in N
        for (int i = 0; i < n; ++i)
        {
            GRBLinExpr out = 0, in = 0;
            for (int j = 0; j < n; ++j)
                if (i != j)
                {
                    out += x[i][j];
                    in += x[j][i];
                }
            model.addConstr(out == 1);
            model.addConstr(in == 1);
        }

        model.set(GRB_IntParam_LazyConstraints, 1);
        model.set(GRB_DoubleParam_TimeLimit, 180.0);

        model.setCallback(this);
        model.optimize();

        return model.get(GRB_IntAttr_Status);
    }
    catch (GRBException e)
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return -1;
    }
    catch (...)
    {
        cout << "Exception during optimization" << endl;
        return -1;
    }
}

void ATSP_CUT::findSubtour(const vector<vector<int>> &sol, vector<int> &S)
{
    int n = data.size;
    vector<int> visited(n, 0);
    int bestSize = n + 1;

    for (int start = 0; start < n; ++start)
    {
        if (visited[start])
            continue;

        vector<int> tour;
        int i = start;
        while (!visited[i])
        {
            visited[i] = 1;
            tour.push_back(i);
            for (int j = 0; j < n; ++j)
                if (sol[i][j])
                {
                    i = j;
                    break;
                }
        }

        if (tour.size() < bestSize)
        {
            bestSize = tour.size();
            S = tour;
        }
    }
}

void ATSP_CUT::callback()
{
    int n = data.size;

    if (where == GRB_CB_MIPSOL)
    {
        vector<vector<int>> sol(n, vector<int>(n, 0));

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j && getSolution(x[i][j]) > 0.5)
                    sol[i][j] = 1;

        vector<int> S;
        findSubtour(sol, S);

        if ((int)S.size() == n)
            return;

        vector<int> inS(n, 0);
        for (int i : S)
            inS[i] = 1;

        GRBLinExpr cut = 0;
        for (int i = 0; i < n; ++i)
            if (!inS[i])
                for (int j : S)
                    cut += x[i][j];

        addLazy(cut >= 1);
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

        vector<int> S;
        vector<int> visited(n, 0);
        vector<int> stack = {0};
        visited[0] = 1;

        while (!stack.empty())
        {
            int i = stack.back();
            stack.pop_back();
            for (int j = 0; j < n; ++j)
            {
                if (!visited[j] && sol[i][j] > 1e-6)
                {
                    visited[j] = 1;
                    stack.push_back(j);
                }
            }
        }

        for (int i = 0; i < n; ++i)
            if (visited[i])
                S.push_back(i);

        if ((int)S.size() == n)
            return;

        vector<int> inS(n, 0);
        for (int i : S)
            inS[i] = 1;

        GRBLinExpr cut = 0;
        double capacity = 0.0;

        for (int i = 0; i < n; ++i)
            if (!inS[i])
                for (int j : S)
                {
                    cut += x[i][j];
                    capacity += sol[i][j];
                }

        if (capacity < 1.0 - 1e-6)
            addCut(cut >= 1);
    }
}

void ATSP_CUT::printSolution(int status)
{
    if (status == GRB_OPTIMAL || (status == GRB_TIME_LIMIT && model.get(GRB_IntAttr_SolCount) > 0))
    {
        cout << "Succes! (Status: " << status << ")" << endl;
        cout << "Runtime : " << model.get(GRB_DoubleAttr_Runtime) << " seconds" << endl;
        model.write("solution.sol");
        cout << "Objective value = " << model.get(GRB_DoubleAttr_ObjVal) << endl;
        for (int i = 0; i < data.size; ++i)
        {
            for (int j = 0; j < data.size; ++j)
            {
                if (i != j && x[i][j].get(GRB_DoubleAttr_X) > 0.5)
                {
                    cout << "x(" << i << "," << j << ") = " << x[i][j].get(GRB_DoubleAttr_X) << endl;
                }
            }
        }
    }
    else
    {
        cout << "No solution found. (Status: " << status << ")" << endl;
    }
}
