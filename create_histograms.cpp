
#include <cassert>
#include <iostream>
#include <numeric>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"

void create_histograms() {

  std::vector<std::unique_ptr<TF1> > functions;

  functions.push_back(std::make_unique<TF1>("function_0", "gaus(0) + pol1(3)", 0, 100));
  functions.back()->SetParameter(0, 40);
  functions.back()->SetParameter(1, 10);
  functions.back()->SetParameter(2, 2);
  functions.back()->SetParameter(3, 20);
  functions.back()->SetParameter(4, -20. / 100.);

  functions.push_back(std::make_unique<TF1>("function_1", "gaus(0) + gaus(3) + pol1(6)", 0, 100));
  functions.back()->SetParameter(0, 30);
  functions.back()->SetParameter(1, 15);
  functions.back()->SetParameter(2, 4);
  functions.back()->SetParameter(3, 40);
  functions.back()->SetParameter(4, 80);
  functions.back()->SetParameter(5, 2);
  functions.back()->SetParameter(6, 10);
  functions.back()->SetParameter(7, -10. / 100.);

  functions.push_back(std::make_unique<TF1>("function_2", "gaus(0) + expo(3) + pol0(5)", 0, 100));
  functions.back()->SetParameter(0, 40);
  functions.back()->SetParameter(1, 85);
  functions.back()->SetParameter(2, 2);
  functions.back()->SetParameter(3, 0.1);
  functions.back()->SetParameter(4, -0.1);
  functions.back()->SetParameter(5, 1);

  functions.push_back(std::make_unique<TF1>("function_3", "gaus(0) + gaus(3) + gaus(6)", 0, 100));
  functions.back()->SetParameter(0, 10);
  functions.back()->SetParameter(1, 12.5);
  functions.back()->SetParameter(2, 3);
  functions.back()->SetParameter(3, 10);
  functions.back()->SetParameter(4, 50);
  functions.back()->SetParameter(5, 4);
  functions.back()->SetParameter(6, 20);
  functions.back()->SetParameter(7, 90);
  functions.back()->SetParameter(8, 2);

  functions.push_back(std::make_unique<TF1>("function_4", "gaus(0) + gaus(3) + pol1(6)", 0, 100));
  functions.back()->SetParameter(0, 10);
  functions.back()->SetParameter(1, 50);
  functions.back()->SetParameter(2, 4);
  functions.back()->SetParameter(3, 30);
  functions.back()->SetParameter(4, 82.5);
  functions.back()->SetParameter(5, 3);
  functions.back()->SetParameter(6, 0);
  functions.back()->SetParameter(7, 10. / 100.);

  const std::vector<double> scaling_factors = {0.15, 0.25, 0.10, 0.15, 0.35};
  assert(std::accumulate(scaling_factors.cbegin(),
			 scaling_factors.cend(),
			 decltype(scaling_factors)::value_type(0)) == 1.0);
  assert(functions.size() == scaling_factors.size());

  auto output = std::make_unique<TFile>("histograms.root", "RECREATE");

  std::vector<std::unique_ptr<TH1D> > response_functions;
  int number_of_fills = 1e3;
  int scaling_factor_index = 0;
  auto data = std::make_unique<TH1D>("data", "", 100, 0, 100);
  for (const auto& function : functions) {
    double scaling_factor = scaling_factors[scaling_factor_index++];
    std::string name = function->GetName();
    name = "response_function" + name.substr(name.find('_'));
    response_functions.push_back(std::make_unique<TH1D>(name.c_str(), std::to_string(scaling_factor).c_str(),
							data->GetXaxis()->GetNbins(),
							data->GetXaxis()->GetXmin(),
							data->GetXaxis()->GetXmax()));
    response_functions.back()->FillRandom(function->GetName(), number_of_fills);
    response_functions.back()->Write();
    response_functions.back()->Scale(scaling_factor);
    data->Add(response_functions.back().get());
  }
  data->Write();

  auto canvas = std::make_unique<TCanvas>("canvas", "canvas");
  data->SetLineColor(1);
  data->Draw("hist same");

  int line_color = 2;
  for (const auto& response_function : response_functions) {
    response_function->SetLineColor(line_color++);
    response_function->Draw("hist same");
  }

  canvas->WaitPrimitive();

}
