#ifndef PDF_LOOKUP_H
#define PDF_LOOKUP_H

#include <string>
#include <vector>
#include <unordered_map>

class PDF_Lookup {
    public:
        struct BasePDF {
            double binWidth = 1; // in microseconds (coarse)
            double fineBinWidth = 0.25; // in microseconds (fine)
            std::vector<double> p; // pdf for binWidth (0 us offset)
            std::vector<double> p_f; // pdf for fine bin (0 us offset)
        };

        PDF_Lookup(double coarse_bw, double fine_bw, const std::string& csv_path);

        bool load_csv(const std::string& path);
        bool has_segment(int seg_id) const { return seg_.count(seg_id) > 0; }
        const BasePDF* get(int seg_id) const {
            auto it = seg_.find(seg_id);
            return (it == seg_.end()) ? nullptr : &it->second;
        }

    private:
        double coarseBinWidth_;
        double fineBinWidth_;
        std::string csv_path_;
        std::unordered_map<int, BasePDF> seg_; // seg_id -> bases
        
        static int parse_segment_name(const std::string& header); // "Segment_12" -> 12
        static void normalize_series(std::vector<double>& v, double dt);
        static void downsample(const std::vector<double>& src, double dt_src,
                                    std::vector<double>& dst, double dt_dst);
        static bool approx_int_multiple(double value, double step, double rel_eps = 1e-9);
};

#endif // PDF_LOOKUP_H