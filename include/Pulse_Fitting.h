#ifndef PULSE_FITTING_H
#define PULSE_FITTING_H

#include <tuple>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include "File_Loader.h" // For EventList

struct PDFParams {
    double ratio1, ratio2, ratio3;
    double scale1, scale2, scale3;
    double loc;
};

extern PDFParams pdfParams_;
extern std::vector<double> log_fact_table;
extern std::vector<double> log_lambda_table;
std::vector<double> makeLogFactorialTable(int max_k);
std::vector<double> makeLogLambdaTable();
double getLogLambda(double lam);

class Pulse_Fitting {
    public:
        Pulse_Fitting(const EventList& events, double binWidth = 1.0, double minGap = 10.0);

        void setWindow(double start_us, double stop_us);
        void setBackgroundWindow(double start_us);
        void analyze();

        const std::vector<std::tuple<double, double, int, double, bool>>& getSignalPulses() const { return signalPulses_; }
        const std::vector<std::tuple<double, double, int, double, bool>>& getBackgroundPulses() const { return backgroundPulses_; }

    private:
        double binWidth_;
        double fineBinWidth_ = 0.25;
        double minGap_;
        double startAfterUs_;
        double stopAfterUs_;
        double backgroundAfterUs_;
        std::vector<double> peTimes_;

        std::map<std::pair<int, double>, std::vector<std::vector<double>>> pdfCache_;

        std::vector<std::tuple<double, double, int, double, bool>> signalPulses_;
        std::vector<std::tuple<double, double, int, double, bool>> backgroundPulses_;
        double peBackgroundRate_;
        double eventBackgroundRate_;

        // === HELPER METHODS === //

        void extractTimes(const EventList& events);
        
        std::vector<double> applyTimeWindow(const std::vector<double>& times, double start, double end);
        
        std::tuple<double, int, double, double> movingWindow(const std::vector<double>& times, int startIdx);
        
        bool makeHistogram(const std::vector<double>& times, int i, double binWidth,
                        double& windowWidth, int& j, double& startTime, double& endTime,
                        std::vector<int>& hist, std::vector<double>& xCenters);
        
        void fitRegion(const std::vector<double>& data_us,
                    std::vector<std::tuple<double, double, int, double, bool>>& output);

        std::vector<double> analyticPDF(const std::vector<double>& x, int shift = 0);
        
        std::vector<std::vector<double>> generatePDFLookup(const std::vector<double>& xCenters);

        double poissonLogLikelihood(const std::vector<int>& observed,
                                    const std::vector<double>& expected);

        double negLogLikelihood(const std::vector<double>& params,
                                const std::vector<int>& observed,
                                const std::vector<std::vector<double>>& pdfLookup,
                                int nPulses);

        std::vector<int> findGradientPeaks(const std::vector<int>& hist, double threshold, int ignoreIdx);
                                
        bool fitPulses(const std::vector<int>& hist, const std::vector<double>& xCenters,
                    const std::vector<std::vector<double>>& pdfLookup,
                    std::vector<double>& fittedPEs, std::vector<double>& fittedDTs);
};

#endif // PULSE_FITTING_H