#include <iostream>
#include "ATSP_MTZ.hpp"
#include "ATSP_CUT.hpp"

using namespace std;

int main(int argc, char **argv)
{
  if (argc < 2)
  {
    cout << "usage : " << argv[0] << " ATSPFilename" << endl;
    return 0;
  }

  if (argc < 3)
  {
    cout << "No solver specified, using MTZ by default." << endl;
  }

  if (argc >= 3 && string(argv[2]) == "CUT")
  {
    ATSPDataC data(argv[1]);
    // data.printData();

    ATSP_CUT solver(data);
    solver.solve();

    return 0;
  }

  ATSPDataC data(argv[1]);
  // data.printData();

  ATSP_MTZ solver(data);
  solver.solve();

  return 0;
}
