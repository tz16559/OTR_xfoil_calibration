#include <fstream>
#include <iostream>
#include <vector>
typedef unsigned char BYTE;

std::vector<BYTE> readFile(const char* filename)
{
    // open the file:
    std::streampos fileSize;
    std::ifstream file(filename, std::ios::binary);

    // get its size:
    file.seekg(0, std::ios::end);
    fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    // read the data:
    std::vector<BYTE> fileData(fileSize);
    file.read((char*) &fileData[0], fileSize);
    return fileData;
}


int main(){

    std::vector<BYTE> fileData = readFile("./test_pgm/filament1_calib/V_00000000.pgm");

    double *sum = new double[fileData.size()];

    for(int i=0; i<fileData.size(); i++){
        sum[i] = (double)fileData[i];
    }
    std::vector<BYTE> fileData_1 = readFile("./test_pgm/filament1_calib/V_00000001.pgm");
     for(int i=0; i<fileData_1.size(); i++){
        sum[i] += (double)fileData[i];
    }

}
