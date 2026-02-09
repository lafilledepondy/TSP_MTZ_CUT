#include <iostream>
#include "ATSP_MTZ.cpp"

using namespace std;

int main(int argc, char **argv)
{
  if (argc < 2)
  {
    cout << "usage : " << argv[0] << " ATSPFilename" << endl;
    return 0;
  }

  ATSPDataC data(argv[1]);
  data.printData();

  ATSP_MTZ solver(data);
  solver.solve();

  // ===> Test code for reading data already implemented in ATSPDataC constructor and printData() method.
  // cout << "Size : " << data.size << endl;
  // cout << "Distances :" << endl;
  // for (size_t i = 0; i < data.size; ++i){
  //   for (size_t j = 0; j < data.size; ++j){
  //     cout << data.distances[i][j] << " ";
  //   }
  //   cout << endl;
  // }

  return 0;
}
