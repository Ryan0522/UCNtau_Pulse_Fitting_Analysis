#ifndef PULSE_FITTING_H
#define PULSE_FITTING_H

#include <tuple>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include "File_Loader.h" // For EventList

struct PDFParams {
    // parameters for the PDF model of PE response from the PMTs
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
        // events: raw PE hits (list of 'event'); binWidth: coarse hist bin (us); minGap: break windows (us)
        Pulse_Fitting(const EventList& events, double binWidth = 1.0, double minGap = 10.0);

        void setWindow(double start_us, double stop_us); // signal window [start, stop) in us
        void setBackgroundWindow(double start_us); // background window [start, start+60s)
        void analyze(); // build windows, fit pulses, fill outputs

        const std::vector<std::tuple<double, double, int, double, bool>>& getSignalPulses() const { return signalPulses_; }
        const std::vector<std::tuple<double, double, int, double, bool>>& getBackgroundPulses() const { return backgroundPulses_; }

    private:
        double binWidth_; // primary histogram bin (us)
        double fineBinWidth_ = 0.25; // fallback giner bining (us)
        double minGap_; // max inter-hit gap before closing a window (us)
        double startAfterUs_; // signal start time (us)
        double stopAfterUs_; // signal stop time (us)
        double backgroundAfterUs_; // bg start (us), bg end = start + 60s
        std::vector<double> peTimes_; // all PE times (us)

        std::map<std::pair<int, double>, std::vector<std::vector<double>>> pdfCache_; // keyed by (nbins, binWidth)

        // each tuple: (pulse_time_us, PE, window_index, window_width_us, is_pileup)
        std::vector<std::tuple<double, double, int, double, bool>> signalPulses_;
        std::vector<std::tuple<double, double, int, double, bool>> backgroundPulses_;
        double peBackgroundRate_;
        double eventBackgroundRate_;

        // === HELPER METHODS === //

        void extractTimes(const EventList& events); // copy realtime to peTimes_
        
        std::vector<double> applyTimeWindow(const std::vector<double>& times, double start, double end); // [start, end)
        
        std::tuple<double, int, double, double> movingWindow(const std::vector<double>& times, int startIdx); // grow window by minGap_
        
        bool makeHistogram(const std::vector<double>& times, int i, double binWidth,
                        double& windowWidth, int& j, double& startTime, double& endTime,
                        std::vector<int>& hist, std::vector<double>& xCenters); // build window hist from times[i...j)
        
        void fitRegion(const std::vector<double>& data_us,
                    std::vector<std::tuple<double, double, int, double, bool>>& output);

        std::vector<double> analyticPDF(const std::vector<double>& x, int shift = 0); // tri-exp mixture over bins (normalized)
        
        std::vector<std::vector<double>> generatePDFLookup(const std::vector<double>& xCenters); // cached shifted PDFs

        double poissonLogLikelihood(const std::vector<int>& observed,
                                    const std::vector<double>& expected);

        double negLogLikelihood(const std::vector<double>& params,
                                const std::vector<int>& observed,
                                const std::vector<std::vector<double>>& pdfLookup,
                                int nPulses); // seed candidates

        std::vector<int> findGradientPeaks(const std::vector<int>& hist, double threshold, int ignoreIdx); // NLOpt fit over PE, DT per pulse
                                
        bool fitPulses(const std::vector<int>& hist, const std::vector<double>& xCenters,
                    const std::vector<std::vector<double>>& pdfLookup,
                    std::vector<double>& fittedPEs, std::vector<double>& fittedDTs);
};

#endif // PULSE_FITTING_H