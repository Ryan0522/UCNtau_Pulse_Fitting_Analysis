#include "File_Loader.h"
#include "PDF_Lookup.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <iostream>
#include <vector>

static TH1D* make_hist(const std::vector<double>& y, double dx, const char* name, const char* title) {
    int n = (int)y.size();
    TH1D* h = new TH1D(name, title, n, 0.0, n*dx);
    for (int i=0;i<n;++i) h->SetBinContent(i+1, y[i]);
    h->GetXaxis()->SetTitle("time (#mu s)");
    h->GetYaxis()->SetTitle("PDF");
    return h;
}

int main(int argc, char** argv) {
    Config cfg;
    try {
        cfg = load_config(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Config error: " << e.what() << "\n";
        return 1;
    }

    PDF_Lookup lookup(cfg.bin_width_us, cfg.fine_bin_width_us, cfg.pdf_csv_path);

    std::vector<int> segments = {12, 34, 56, 78};

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c","Segment PDFs",1200,800);
    c->Divide(2,2);

    for (size_t i=0;i<segments.size() && i<4;i++) {
        std::cout << "Segment " << segments[i] << "; Plotting..." << std::endl;
        int seg = segments[i];
        const auto* base = lookup.get(seg);
        c->cd((int)i+1);
        if (!base) {
            TLatex tl; tl.SetNDC(); tl.DrawLatex(0.2,0.5,Form("Segment %d not found", seg));
            continue;
        }
        TH1D* hf = make_hist(base->p_f, base->fineBinWidth, Form("hf%d",seg), Form("Segment %d (fine)",seg));
        TH1D* hc = make_hist(base->p,   base->binWidth,     Form("hc%d",seg), Form("Segment %d (coarse)",seg));
        hf->SetLineColor(kBlue);  hf->SetLineWidth(2);
        hc->SetLineColor(kRed+1); hc->SetLineWidth(2);

        double ymax = std::max(hf->GetMaximum(), hc->GetMaximum()) * 1.2;
        hf->SetMaximum(ymax); hf->SetTitle(Form("Segment %d PDF (fine vs coarse)", seg));
        hf->Draw("HIST");
        hc->Draw("HIST SAME");

        auto leg = new TLegend(0.60,0.70,0.88,0.88);
        leg->AddEntry(hf, "fine", "l");
        leg->AddEntry(hc, "coarse", "l");
        leg->Draw();
        gPad->SetGrid();
    }

    c->SaveAs((cfg.output_folder + "/graphs/segment_pdfs.png").c_str());
    std::cout << "Wrote " << cfg.output_folder << "/graphs/segment_pdfs.png\n";
    return 0;
}
