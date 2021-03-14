#include <vector> 
#include <cmath>

class DiffPlot{
    
public:
    std::vector<std::vector<double>> holesX;
    std::vector<std::vector<double>> holesY;
    std::vector<std::vector<int>> order;
    std::vector<double> meanHoleX;
    std::vector<double> meanHoleY;

    int nImages, nHoles;
    
public:
    
    DiffPlot(int noI, int noH);
    virtual ~DiffPlot();
    void IDHoles();
    void ID(int picID);
    void pltDevFromAvg();
    void calcAvgPos();
    void ingest(std::vector<std::vector<double>> &X, std::vector<std::vector<double>> &Y);
    void pltDevRelative();
    void pltAbsDev();
    void pltAbsPos();
    void pltCentres();

};
