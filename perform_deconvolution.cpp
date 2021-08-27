
#include <cassert>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

TH1D* data;
std::vector<TH1D*> response_matrix;

double global_chi_squared;
double global_num_bins;

void normalize(const double* params,
	       const double* errors,
	       std::vector<double>& normalized_params,
	       std::vector<double>& normalized_errors) {

  double params_sum = 0;
  double params_sum_error = 0;
  for (std::size_t i = 0; i < response_matrix.size(); ++i) {
    params_sum += params[i];
    params_sum_error += std::pow(errors[i], 2);
  }
  params_sum_error = std::sqrt(params_sum_error);

  for (std::size_t i = 0; i < response_matrix.size(); ++i) {
    normalized_params.push_back(params[i] / params_sum);
    normalized_errors.push_back(normalized_params[i] *
				std::sqrt(std::pow(errors[i] / params[i], 2) +
					  std::pow(params_sum_error / params_sum, 2)));
  }
}

std::string minimizer_status(const int& minimizer_status_value) {

  // Minimizer status
  // https://root.cern.ch/doc/master/classROOT_1_1Minuit2_1_1Minuit2Minimizer.html

  if (minimizer_status_value == 0) {
    return "Minimization converged";
  }
  else if (minimizer_status_value == 1) {
    return "Covariance was made positive definite";
  }
  else if (minimizer_status_value == 2) {
    return "Hesse is invalid";
  }
  else if (minimizer_status_value == 3) {
    return "Estimated distance to minimum is above max";
  }
  else if (minimizer_status_value == 4) {
    return "Reached call limit";
  }
  else if (minimizer_status_value == 5) {
    return "Any other failure";
  }
  throw std::runtime_error("Unknown minimizer status value: "
			   + std::to_string(minimizer_status_value));

}

double calculate_chi_squared(TH1D* data,
			     std::vector<TH1D*> response_matrix,
			     const double* parameters,
			     int& bin_count) {

  double chi_squared = 0;
  bin_count = 0;

  // [1, NbinsX] inclusive to avoid underflow/overflow bins
  int number_of_bins = data->GetXaxis()->GetNbins();
  for (int bin = 1; bin <= number_of_bins; ++bin) {
    double observed = data->GetBinContent(bin);
    double expected = 0;
    if (observed > 0) {
      for (std::size_t index = 0; index < response_matrix.size(); ++index) {
	expected += parameters[index] * response_matrix[index]->GetBinContent(bin);
      }
      chi_squared += std::pow(observed - expected, 2) / observed;
      //chi_squared += 2 * (expected - observed + observed * std::log(observed / expected));
      ++bin_count;
    }
  }

  return chi_squared;

}

double minimization_function(const double* parameters) {

  int bins = 0;
  double chi_squared = calculate_chi_squared(data,
					     response_matrix,
					     parameters,
					     bins);

  global_num_bins = bins;
  global_chi_squared  = chi_squared;

  return global_chi_squared;

}

std::string time_stamp(std::chrono::time_point<std::chrono::system_clock> time_point) {

  std::time_t time = std::chrono::system_clock::to_time_t(time_point);
  std::ostringstream buffer;
  buffer << std::put_time(localtime(&time), "%Y-%m-%d %H:%M:%S");
  return buffer.str();

}

void perform_deconvolution(int deconvolution_number,
			   std::string minimizer_name = "Minuit",
			   std::string algorithm_name = "Migrad",
			   int max_iterations = 1e5,
			   int max_function_calls = 1e6,
			   double tolerance = 0.01,
			   double error_definition = 1.0,
			   int print_level = 0) {

  auto start = std::chrono::system_clock::now();
  const std::string start_time_stamp = time_stamp(start);
  std::cout << start_time_stamp << std::endl;

  auto input = std::make_unique<TFile>("histograms.root");
  if (!input) {
    throw std::runtime_error("Could not open ROOT file");
  }

  data = static_cast<TH1D*>(input->Get("data"));
  if (!data) {
    throw std::runtime_error("Could not find data histogram");
  }

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

  for (const auto& response_function : response_matrix) {
    assert(data->GetXaxis()->GetNbins() == response_function->GetXaxis()->GetNbins());
    assert(data->GetXaxis()->GetXmin() == response_function->GetXaxis()->GetXmin());
    assert(data->GetXaxis()->GetXmax() == response_function->GetXaxis()->GetXmax());
  }

  ROOT::Math::Minimizer* minimizer =
    ROOT::Math::Factory::CreateMinimizer(minimizer_name, algorithm_name);
  ROOT::Math::Functor minimizer_function(&minimization_function, response_matrix.size());

  minimizer->SetPrintLevel(print_level);
  minimizer->SetMaxIterations(max_iterations);
  minimizer->SetMaxFunctionCalls(max_function_calls);
  minimizer->SetTolerance(tolerance);
  minimizer->SetErrorDef(error_definition);
  minimizer->SetFunction(minimizer_function);

  std::vector<std::string> parameter_names;
  std::vector<double> start_values;
  std::vector<double> step_sizes;
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;

  for (std::size_t i = 0; i < response_matrix.size(); ++i) {
    parameter_names.push_back(response_matrix[i]->GetName());
    start_values.push_back(1.0 / response_matrix.size());
    step_sizes.push_back(0.001);
    lower_bounds.push_back(0.0);
    upper_bounds.push_back(10.0);
  }
  /*
  for (std::size_t i = 0; i < response_matrix.size(); ++i) {
    std::cout << parameter_names[i] << "\t"
	      << start_values[i] << "\t"
	      << step_sizes[i] << "\t"
	      << lower_bounds[i] << "\t"
	      << upper_bounds[i] << std::endl;
  }
  */

  for (std::size_t i = 0; i < response_matrix.size(); ++i) {
    minimizer->SetLimitedVariable(i, parameter_names[i].c_str(),
				  start_values[i], step_sizes[i],
				  lower_bounds[i], upper_bounds[i]);
  }

  minimizer->Minimize();

  auto end = std::chrono::system_clock::now();
  const std::string end_time_stamp = time_stamp(end);
  std::cout << end_time_stamp << std::endl;

  const std::string output_file_name = "output_" + std::to_string(deconvolution_number) + ".txt";
  std::ofstream output_file(output_file_name.c_str());

  output_file << "Start Timestamp:                 " << start_time_stamp << std::endl;
  output_file << "End Timestamp:                   " << end_time_stamp << std::endl;
  output_file << "Elapsed Time:                    " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
  output_file << "Number of response functions:    " << response_matrix.size() << std::endl;
  output_file << "Free and constrained variables:  " << minimizer->NDim() << std::endl;
  output_file << "Free variables:                  " << minimizer->NFree() << std::endl;

  output_file << "Minimizer:                       " << minimizer_name << std::endl;
  output_file << "Algorithm:                       " << algorithm_name << std::endl;
  output_file << "Max Iterations:                  " << minimizer->MaxIterations() << std::endl;
  output_file << "Max Function Calls:              " << minimizer->MaxFunctionCalls() << std::endl;
  output_file << "Actual Number of Function Calls: " << minimizer->NCalls() << std::endl;
  output_file << "Tolerance:                       " << minimizer->Tolerance() << std::endl;
  output_file << "Precision:                       " << minimizer->Precision() << std::endl;
  output_file << "Error Definition:                " << minimizer->ErrorDef() << std::endl;
  output_file << "Estimated Distance To Minimum:   " << minimizer->Edm() << std::endl;
  output_file << "Print Level:                     " << minimizer->PrintLevel() << std::endl;

  output_file << "Status:                          " << minimizer->Status() << " (" << minimizer_status(minimizer->Status()) << ")" <<std::endl;
  output_file << "Minimum Value:                   " << minimizer->MinValue() << std::endl;

  output_file << "global_chi_squared:              " << global_chi_squared << std::endl;
  output_file << "global_num_bins:                 " << global_num_bins << std::endl;

  output_file << "Reduced global chi squared:      " << global_chi_squared / (global_num_bins - minimizer->NFree()) << std::endl;

  const double* minimum_parameters = minimizer->X();
  const double* minimum_errors = minimizer->Errors();

  std::vector<double> normalized_minimum_parameters;
  std::vector<double> normalized_minimum_errors;
  normalize(minimum_parameters,
	    minimum_errors,
	    normalized_minimum_parameters,
	    normalized_minimum_errors);

  for (std::size_t i = 0; i < response_matrix.size(); ++i) {
    output_file << std::left
		<< "unnormalized "
		<< std::setw(parameter_names[i].size() + 5) << parameter_names[i]
		<< std::setw(20) << minimum_parameters[i]
		<< std::setw(20) << minimum_errors[i] << std::endl;
  }

  for (std::size_t i = 0; i < response_matrix.size(); ++i) {
    output_file << std::left
		<< "normalized "
		<< std::setw(parameter_names[i].size() + 5) << parameter_names[i]
		<< std::setw(35) << std::setprecision(25) << normalized_minimum_parameters[i]
		<< std::setw(35) << std::setprecision(25) << normalized_minimum_errors[i] << std::endl;
  }

  minimizer->Clear();
  output_file.close();

}
