// #pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "powell_nr3.hpp"

struct leastSquaredDistance_QGauss {
  double* xvalue;
  double* yvalue;
  size_t nbValues;
  double minValue = -100.0;

  double qgauss(double x, double q, double beta, double normalizer) {
    double term = 1.0 - (1.0 - q) * beta * x * x;
    if (term < 0.0) {
      return 0.0;
    }
    return normalizer * std::pow(term, 1.0 / (1.0 - q));
  }

  double qgauss_integrate(double a, double b, double q, double beta, double normalizer, int n = 1000) {
    double h = (b - a) / n;
    double x = a;
    double sum = 0.0;

    // Using the Trapezoidal rule
    for (int i = 0; i < n; i++) {
      sum += ((qgauss(x, q, beta, normalizer) + qgauss(x + h, q, beta, normalizer)) * 0.5) * h;
      x += h;
    }

    return sum;
  }

  double operator()(std::vector<double>& X) {
    double q = X[0];
    double beta = X[1];
    double normalizer = X[2];

    double diff2_sum = 0.0;
    for (size_t i = 0; i < nbValues; i++) {
      double diff = qgauss_integrate(minValue, xvalue[i], q, beta, normalizer, 1000) - yvalue[i];
      diff2_sum += diff * diff;
    }
    return diff2_sum;
  }
};

struct CDF {
  std::vector<double> value;  // variable value
  std::vector<double> cdens;  // cumulated density value

  double q;
  double beta;
  double normalizer;

  CDF(std::vector<double>& raw, int nb, bool sym = false) {
    std::vector<double> raw_copy(raw.begin(), raw.end());
    if (sym == true){
      size_t nn = raw_copy.size();
      for (size_t i = 0 ; i < nn ; i++) {
        raw_copy.push_back(-raw_copy[i]);
      }  
    }
    std::sort(raw_copy.begin(), raw_copy.end());

    // Downsample the sorted copy
    if (raw_copy.empty() || nb <= 0 || (size_t)nb > raw_copy.size()) {
        std::cerr << "Invalid input for downsampling" << std::endl;
        return;
    }

    std::cout << "initial number of values = " << raw_copy.size() << std::endl;
    size_t nbValues = raw_copy.size();
    size_t step = (nbValues) / nb; // This ensures that the last element is included in the result

    try {
        for (size_t i = 0; i < nbValues; i++) {
            if (i % step == 0) {  // Keep the element if the index is a multiple of the step
                value.push_back(raw_copy[i]);
                cdens.push_back(static_cast<double>(i) / nbValues);
            }
        }
    } catch (const std::bad_alloc& e) {
        std::cerr << "Allocation error: " << e.what() << std::endl;
        return;
    }

    std::cout << "final number of values = " << value.size() << std::endl;
  }

  double findQGaussParameters(double q_ini = 1.2, double beta_ini = 1.0, double normalizer_ini = 1.0) {
    std::vector<double> X{q_ini, beta_ini, normalizer_ini};
    std::vector<double> dX{0.01, 0.01, 0.01};
    leastSquaredDistance_QGauss func;
    func.xvalue = &value[0];
    func.yvalue = &cdens[0];
    func.nbValues = value.size();
    func.minValue = std::min(value[0], -value.back());
    Powell<leastSquaredDistance_QGauss> powell(func);
    X = powell.minimize(X, dX); 

    q = X[0];
    beta = X[1];
    normalizer = X[2];
    return powell.fret;
  }
};

#if 0

#include <fstream>
#include <iostream>
#include <random>

std::vector<double> generate_gaussian_values(int num_values, double mean, double stddev) {
  // Create a random number engine
  std::random_device rd;
  std::mt19937 gen(rd());

  // Create a normal distribution with the specified mean and standard deviation
  std::normal_distribution<> distr(mean, stddev);

  // Generate the random values and store them in a vector
  std::vector<double> values;
  for (int i = 0; i < num_values; i++) {
    values.push_back(distr(gen));
  }

  return values;
}

std::vector<double> generate_weibull_values(int num_values, double k, double lambda) {
  // Create a random number engine
  std::random_device rd;
  std::mt19937 gen(rd());

  // Create a uniform distribution between 0 and 1
  std::uniform_real_distribution<> distr(0.0, 1.0);

  // Generate the random values and store them in a vector
  std::vector<double> values;
  for (int i = 0; i < num_values; i++) {
    double u = distr(gen);
    double x = lambda * std::pow(-std::log(1.0 - u), 1.0/k);
    values.push_back(x);
  }

  return values;
}




int main(int argc, char const* argv[]) {

  //std::vector<double> raw_values = generate_gaussian_values(100000, 0.0, 1.0);
  std::vector<double> raw_values =generate_weibull_values(100000, 1.5, 1.0);
 

  CDF cdf(raw_values, 1000);
  std::ofstream file("cdf2.txt");
  for (size_t i = 0; i < cdf.value.size(); i++) {
    file << cdf.value[i] << ' ' << cdf.cdens[i] << '\n';
  }

  return 0;
}

#endif