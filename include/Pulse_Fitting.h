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
        Pulse_Fitting(const EventList& events);
        void setBinWidths(double coarse_us, double fine_us) {
            binWidth_ = coarse_us > 0 ? coarse_us : binWidth_; 
            fineBinWidth_ = fine_us > 0 ? fine_us : fineBinWidth_; 
        }
        void setMinGapUs(double mg_us) { if (mg_us >= 0) minGap_ = mg_us; }
        void setSegmentId(int seg) { segmentId_ = seg; }
        void setConfigKnobs(double shift_us, double seed_pe_default, double gradient_threshold, int guard_bin) {
            shiftUs_ = shift_us; seedPE_ = seed_pe_default; gradThr_ = gradient_threshold; guardBin_ = std::max(0, guard_bin);
        }

        template <class Cfg>
        void initFromConfig(const Cfg& c) {
            setBinWidths(c.bin_width_us, c.fine_bin_width_us);
            setMinGapUs(c.min_gap_us);
            setConfigKnobs(c.shift_us, c.seed_pe_default, c.gradient_threshold, c.guard_bin);
            pileupMinPulses_ = c.pileup_min_pulses;
            capturePlots_ = c.plot_fits;
            preBins_ = (int)std::llround(shiftUs_ / std::max(1e-12, binWidth_));
            finePreBins_ = (int)std::llround(shiftUs_ / std::max(1e-12, fineBinWidth_));
        }

        double shift_us() const { return shiftUs_; }

        // NEW (2025/8/15): plotting capture controls
        void enablePlotCapture(bool enable, int pileupMinPulses) {
            capturePlots_ = enable;
            pileupMinPulses_ = pileupMinPulses;
        }
        const std::vector<int>& getSingleHist() const { return single_hist_; }
        const std::vector<double>& getSingleX() const { return single_x_; }
        const std::vector<double>& getSingleTotal() const { return single_total_; }
        const std::vector<std::vector<double>>& getSingleComps() const { return single_components_; }

        const std::vector<int>& getPileHist() const { return pile_hist_; }
        const std::vector<double>& getPileX() const { return pile_x_; }
        const std::vector<double>& getPileTotal() const { return pile_total_; }
        const std::vector<std::vector<double>>& getPileComps() const { return pile_components_; } 

        void setWindow(double start_us, double stop_us); // signal window [start, stop) in us
        void setBackgroundWindow(double start_us); // background window [start, start+60s)
        void analyze(); // build windows, fit pulses, fill outputs

        const std::vector<std::tuple<double, double, int, double, bool, bool>>& getSignalPulses() const { return signalPulses_; }
        const std::vector<std::tuple<double, double, int, double, bool, bool>>& getBackgroundPulses() const { return backgroundPulses_; }

        const double getBackgroundRate() const { return peBackgroundRate_; }

    private:
        double binWidth_ = 1.0; // primary histogram bin (us)
        double fineBinWidth_ = 0.25; // fallback giner bining (us)

        double minGap_ = 0.0; // max inter-hit gap before closing a window (us)
        double startAfterUs_ = 0; // signal start time (us)
        double stopAfterUs_ = 1e12; // signal stop time (us)
        double backgroundAfterUs_ = 0; // bg start (us), bg end = start + 60s
        int segmentId_ = 12; // default
        std::vector<double> peTimes_; // all PE times (us)

        // --- NEW configurable knobs (defaults if config omits them)
        double shiftUs_ = 5.0;
        double seedPE_ = 20.0;
        double gradThr_ = 2.0;
        int guardBin_ = 1;
        int preBins_;
        int finePreBins_;

        std::map<std::pair<int, double>, std::vector<std::vector<double>>> pdfCache_; // keyed by (nbins, binWidth)

        // each tuple: (pulse_time_us, PE, window_index, window_width_us, is_pileup, is_fineGrid)
        std::vector<std::tuple<double, double, int, double, bool, bool>> signalPulses_;
        std::vector<std::tuple<double, double, int, double, bool, bool>> backgroundPulses_;
        double peBackgroundRate_;
        double eventBackgroundRate_;

        // NEW (August 26, 2025): for tail subtraction
        std::vector<std::pair<double, double>> carryPulses_; // stores (abs_time_us, PE) for pulses whose tails may leak into future windows
        std::vector<double> buildCarryExpected(double winStartUs, double binWidth, const std::vector<double>& xCenters) const;
        // --- end NEW ---

        // NEW (2025/8/15): capture buffers
        bool capturePlots_ = false;
        int pileupMinPulses_ = 2;
        bool haveFirstSingle_ = false;
        bool haveFirstPileup_ = false;
        
        std::vector<int> single_hist_;
        std::vector<double> single_x_, single_total_;
        std::vector<std::vector<double>> single_components_;

        std::vector<int> pile_hist_;
        std::vector<double> pile_x_, pile_total_;
        std::vector<std::vector<double>> pile_components_;

        // === HELPER METHODS === //

        void extractTimes(const EventList& events); // copy realtime to peTimes_
        
        std::vector<double> applyTimeWindow(const std::vector<double>& times, double start, double end); // [start, end)
        
        std::tuple<double, int, double, double> movingWindow(const std::vector<double>& times, int startIdx); // grow window by minGap_
        
        bool makeHistogram(const std::vector<double>& times, int i, double binWidth,
                        double& windowWidth, int& j, double& startTime, double& endTime,
                        std::vector<int>& hist, std::vector<double>& xCenters); // build window hist from times[i...j)
        
        void fitRegion(const std::vector<double>& data_us,
                    std::vector<std::tuple<double, double, int, double, bool, bool>>& output);
        
        std::vector<std::vector<double>> generatePDFLookup(const std::vector<double>& xCenters); // cached shifted PDFs

        double poissonLogLikelihood(const std::vector<int>& observed,
                                    const std::vector<double>& expected);

        double negLogLikelihood(const std::vector<double>& params,
                                const std::vector<int>& observed,
                                const std::vector<std::vector<double>>& pdfLookup,
                                int nPulses,
                                const std::vector<double>* fixedExpected); // seed candidates

        std::vector<int> findGradientPeaks(const std::vector<int>& hist, double threshold);
                                
        bool fitPulses(const std::vector<int>& hist, const std::vector<double>& xCenters,
                    const std::vector<std::vector<double>>& pdfLookup,
                    std::vector<double>& fittedPEs, std::vector<double>& fittedDTs, 
                    const double binWidth, const std::vector<double>* fixedExpected = nullptr); // NLOpt fit over PE, DT per pulse
};

#endif // PULSE_FITTING_H