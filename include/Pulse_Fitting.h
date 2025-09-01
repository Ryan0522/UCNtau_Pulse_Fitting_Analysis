#ifndef PULSE_FITTING_H
#define PULSE_FITTING_H

#include <tuple>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include "File_Loader.h" // For EventList

struct WindowStat {
    int windowIndex;
    double startTimeUs;
    double binWidthUs;
    int nPulsesChosen;
    double logL;
    int nObserved; // PE his sum count
    double nExpected; // PDF fit sum count
};

struct OutlierRecord {
    int windowIndex;
    double startTimeUs;
    double binWidthUs;
    int nObserved;
    double nExpected;
    double ratioExpOverObs;
    std::vector<int> hist;
    std::vector<double> totalExpected;
};

class Pulse_Fitting {
    public:
        Pulse_Fitting(const EventList& events);
        void setBinWidths(double coarse_us, double fine_us) {
            binWidth_us_ = coarse_us > 0 ? coarse_us : binWidth_us_; 
            fineBinWidth_us_ = fine_us > 0 ? fine_us : fineBinWidth_us_; 
        }
        void setMinGapUs(double mg_us) { if (mg_us >= 0) minGap_us_ = mg_us; }
        void setSegmentId(int seg) { segmentId_ = seg; }
        void setConfigKnobs(double shift_us, double seed_pe_default, double gradient_threshold, int guard_bin) {
            shiftUs_ = shift_us; seedPE_ = seed_pe_default; gradThr_ = gradient_threshold; guardBin_ = std::max(0, guard_bin);
        }

        template <class Cfg>
        void initFromConfig(const Cfg& c) {
            setBinWidths(c.bin_width_us, c.fine_bin_width_us);
            setMinGapUs(c.min_gap_us);
            setConfigKnobs(c.shift_us, c.seed_pe_default, c.gradient_threshold, c.guard_bin);
            enablePlotCapture(c.plot_fits, c.pileup_min_pulses);
            enableOutlierCapture(c.plot_outliers, c.outlier_min_obs, c.outlier_ratio_low);
            
            preBins_ = (int)std::llround(shiftUs_ / std::max(1e-12, binWidth_us_));
            finePreBins_ = (int)std::llround(shiftUs_ / std::max(1e-12, fineBinWidth_us_));
            use_coinc_ = c.use_coinc;
            coinc_win_us_ = c.coinc_win_us;

            seeding_window_ = c.seeding_window;
            pe_min_thresh_ = c.pe_min_thresh;
        }

        double shift_us() const { return shiftUs_; }

        // NEW (2025/8/15): plotting capture controls
        void enablePlotCapture(bool enable, int pileupMinPulses) {
            capturePlots_ = enable;
            pileupMinPulses_ = pileupMinPulses;
        }

        // NEW (Sept 1, 2025): plotting outlier controls
        void enableOutlierCapture(bool enable, int min_obs = 60, double ratio_low = 0.80) {
            captureOutliers_ = enable;
            outlierMinObs_ = std::max(0, min_obs);
            outlierRatioLow_ = std::max(0.0, std::min(1.0, ratio_low));
        }
        const std::vector<OutlierRecord>& getOutliers() const { return outliers_; }

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
        const std::vector<WindowStat>& getWindowStats() const { return windowStats_; }

    private:
        const EventList events_;

        double binWidth_us_ = 1.0; // primary histogram bin (us)
        double fineBinWidth_us_ = 0.25; // fallback giner bining (us)

        double minGap_us_ = 0.0; // max inter-hit gap before closing a window (us)
        double startAfterUs_ = 0; // signal start time (us)
        double stopAfterUs_ = 1e12; // signal stop time (us)
        double backgroundAfterUs_ = 0; // bg start (us), bg end = start + 60s
        int segmentId_ = 12; // default
        std::vector<double> peTimes_us_; // all PE times (us)

        // --- NEW configurable knobs (defaults if config omits them)
        double shiftUs_ = 5.0;
        double seedPE_ = 20.0;
        double gradThr_ = 2.0;
        int guardBin_ = 1;
        int preBins_;
        int finePreBins_;

        int seeding_window_ = 5;
        double pe_min_thresh_ = 5.0;

        bool captureOutliers_ = false;
        int outlierMinObs_ = 60;
        double outlierRatioLow_ = 0.8;
        std::vector<OutlierRecord> outliers_;
        // == end NEW ==

        std::map<std::pair<int, double>, std::vector<std::vector<double>>> pdfCache_; // keyed by (nbins, binWidth)

        // each tuple: (pulse_time_us, PE, window_index, window_width_us, is_pileup, is_fineGrid)
        std::vector<std::tuple<double, double, int, double, bool, bool>> signalPulses_;
        std::vector<std::tuple<double, double, int, double, bool, bool>> backgroundPulses_;
        double peBackgroundRate_;
        double eventBackgroundRate_;
        bool fitting_bg_;

        std::vector<std::pair<double, double>> carryPulses_; // stores (abs_time_us, PE) for pulses whose tails may leak into future windows
        std::vector<double> buildCarryExpected(double winStartUs, double binWidth, const std::vector<double>& xCenters) const;

        //NEW (August 29, 2025): window stat
        std::vector<WindowStat> windowStats_;
        int curWindowIndex_ = -1;
        double curWindowStartUs_ = 0.0;
        double curWindowBinWidth_us_ = 0.0;
        // -- end NEW ---

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

        bool use_coinc_;
        double coinc_win_us_;

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

        std::vector<int> findGradientPeaks(const std::vector<int>& hist, double threshold);
                                
        bool fitPulses(const std::vector<int>& hist, const std::vector<double>& xCenters,
                    const std::vector<std::vector<double>>& pdfLookup,
                    std::vector<double>& fittedPEs, std::vector<double>& fittedDTs, 
                    const double binWidth, const std::vector<double>* fixedExpected = nullptr); // NLOpt fit over PE, DT per pulse
};

#endif // PULSE_FITTING_H