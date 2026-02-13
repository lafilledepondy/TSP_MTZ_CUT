#include <iostream>
#include <iomanip>
#include "ATSP_MTZ.hpp"
#include "ATSP_CUT.hpp"

using namespace std;


// --------------- PRINT pour le tableau de synthèse --------------- 
// Extrait le nom de fichier depuis un chemin complet
static string baseName(const string &path) {
  size_t pos = path.find_last_of("/\\");
  
  if (pos == string::npos){
    return path;
  }

  return path.substr(pos + 1);
}

// Convertit le statut Gurobi en libelle court pour parametre dans le terminal
static string statusToString(int status) {
  switch (status) {
  case GRB_OPTIMAL:
    return "OPT";
  case GRB_TIME_LIMIT:
    return "TL";
  case GRB_INFEASIBLE:
    return "INF";
  case GRB_UNBOUNDED:
    return "UNB";
  case GRB_INF_OR_UNBD:
    return "INF_OR_UNB";
  default:
    return to_string(status);
  }
}

static bool tryGetDoubleAttr(GRBModel &model, GRB_DoubleAttr attr, double &out) {
  try {
    out = model.get(attr);
    return true;
  }
  catch (GRBException &) {
    return false;
  }
}

// Affiche main
static void printSummary(const string &instance, const string &mode, GRBModel &model, int status, int cuts) {
  int solCount = 0;

  try {
    solCount = model.get(GRB_IntAttr_SolCount);
  }
  catch (GRBException &) {
  }

  double obj = 0.0;
  bool hasObj = (solCount > 0) && tryGetDoubleAttr(model, GRB_DoubleAttr_ObjVal, obj);

  double bound = 0.0;
  bool hasBound = tryGetDoubleAttr(model, GRB_DoubleAttr_ObjBound, bound);

  double runtime = 0.0;
  bool hasTime = tryGetDoubleAttr(model, GRB_DoubleAttr_Runtime, runtime);

  double nodes = 0.0;
  bool hasNodes = tryGetDoubleAttr(model, GRB_DoubleAttr_NodeCount, nodes);

  cout << "RESULT instance=" << instance
       << " mode=" << mode
       << " obj=" << (hasObj ? to_string(obj) : string("NA"))
       << " bound=" << (hasBound ? to_string(bound) : string("NA"))
       << " nodes=" << (hasNodes ? to_string(static_cast<long long>(nodes)) : string("NA"))
       << " cuts=" << cuts
       << " status=" << statusToString(status)
       << " time=" << (hasTime ? to_string(runtime) : string("NA"))
       << endl;
}
// --------------- END--------------- 

int main(int argc, char **argv) {
  if (argc < 2) {
    cout << "usage : " << argv[0] << " ATSPFilename [MTZ|CUT|CUT_LP] [--summary]" << endl;
    return 0;
  }

  string mode = "MTZ";
  bool summary = false;
  for (int i = 2; i < argc; ++i) {  string arg = argv[i];
    if (arg == "--summary"){
      summary = true;
    }
    else {
      mode = arg;
    }
  }

  // Mode CUT sol entier
  if (mode == "CUT" || mode == "CUT_INT") {
    ATSPDataC data(argv[1]);
    ATSP_CUT solver(data, ATSP_CUT::SolveMode::IntegerMIP);
    solver.solve();

    if (summary && solver.getterModel())
      printSummary(baseName(argv[1]), "CUT", *solver.getterModel(), solver.getterStatus(), solver.getTotalCuts());

    return 0;
  }

  // Mode CUT sol frac
  if (mode == "CUT_LP" || mode == "CUT_Q") {
    ATSPDataC data(argv[1]);
    ATSP_CUT solver(data, ATSP_CUT::SolveMode::FractionalLP);
    solver.solve();

    if (summary && solver.getterModel()) {
      printSummary(baseName(argv[1]), "CUT_LP", *solver.getterModel(), solver.getterStatus(), solver.getTotalCuts());
} 

    return 0;
  }

  // MTZ 
  ATSPDataC data(argv[1]);
  ATSP_MTZ solver(data);
  solver.solve();

  if (summary && solver.getterModel()) {
    printSummary(baseName(argv[1]), "MTZ", *solver.getterModel(), solver.getterStatus(), 0);
  }
  else {
    solver.printSolution();
  }

  return 0;
}
