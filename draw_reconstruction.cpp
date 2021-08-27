
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"

std::vector<double> get_parameters(const int& deconvolution_number,
				   const std::string& parameter_type) {

  std::vector<double> output;

  const std::string file_name = "output_" + std::to_string(deconvolution_number) + ".txt";
  std::ifstream file(file_name.c_str());
  if (!file) {
    throw std::runtime_error("Could not open " + file_name);
  }

  std::string line;
  std::string identifier;
  std::string response_function;
  double parameter = 0;
  double error = 0;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    if (iss >> identifier >> response_function >> parameter >> error) {
      if (identifier == parameter_type) {
	output.push_back(parameter);
      }
    }
  }

  if (parameter_type == "normalized") {
    double sum = std::accumulate(output.cbegin(),
				 output.cend(),
				 decltype(output)::value_type(0));
    if (std::abs(1 - sum) > 1e-10) {
      throw std::runtime_error("Normalized parameters sum to " + std::to_string(sum) + " instead of 1.0");
    }
  }

  return output;

}

void draw_reconstruction(int deconvolution_number,
			 std::string parameter_type = "normalized") {

  auto input = std::make_unique<TFile>("histograms.root");
  if (!input) {
    throw std::runtime_error("Could not open ROOT file");
  }

  TH1D* data = static_cast<TH1D*>(input->Get("data"));
  if (!data) {
    throw std::runtime_error("Could not find data histogram");
  }

  std::vector<TH1D*> response_matrix;
  for (const auto& object : *input->GetListOfKeys()) {
    auto key = static_cast<TKey*>(object);
    if (std::string(key->GetClassName()) == "TH1D") {
      std::string name = object->GetName();
      if (name.find("response_function_") != std::string::npos) {
	TH1D* response_function = static_cast<TH1D*>(input->Get(name.c_str()));
	response_matrix.push_back(response_function);
      }
    }
  }
  assert(response_matrix.size() > 0);

  std::vector<double> parameters = get_parameters(deconvolution_number, parameter_type);
  assert(response_matrix.size() == parameters.size());

  auto reconstruction = std::make_unique<TH1D>("reconstruction", "",
					       data->GetXaxis()->GetNbins(),
					       data->GetXaxis()->GetXmin(),
					       data->GetXaxis()->GetXmax());
  for (std::size_t i = 0; i < response_matrix.size(); ++i) {
    reconstruction->Add(response_matrix[i], parameters[i]);
  }

  if (parameter_type == "normalized") {
    reconstruction->Scale(data->Integral() / reconstruction->Integral());
  }

  auto canvas = std::make_unique<TCanvas>("canvas", "canvas");
  data->SetLineColor(1);
  data->SetLineStyle(1);
  data->SetLineWidth(1);
  data->Draw("hist same");
  reconstruction->SetLineColor(2);
  reconstruction->SetLineStyle(2);
  reconstruction->SetLineWidth(1);
  reconstruction->Draw("hist same");
  canvas->WaitPrimitive();

}
