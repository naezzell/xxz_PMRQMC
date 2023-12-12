//
// This program implements Permutation Matrix Representation Quantum Monte Carlo for arbitrary spin-1/2 Hamiltonians.
//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//
// ExExFloat datatype and calculation of divided differences are described in the paper:
// L. Gupta, L. Barash, I. Hen, Calculating the divided differences of the exponential function by addition and removal of inputs, Computer Physics Communications 254, 107385 (2020)
//

#include"mainqmc.hpp"
#include <fstream>
#include <string>
#include <sstream>

double get_cpu_time(){ return (double)clock() / CLOCKS_PER_SEC;}

int main(int argc, char* argv[]){
	if(steps < Nbins*stepsPerMeasurement){
		std::cout<<"Error: steps cannot be smaller than Nbins*stepsPerMeasurement."<<std::endl;
		exit(1);
	}
	double start_time = get_cpu_time();
	double Rsum[N_all_observables] = {0}; double sgn_sum = 0;
	double over_bins_sum[N_all_observables] = {0}; double over_bins_sum_sgn = 0;
	double over_bins_sum_cov[N_all_observables] = {0}; double mean[N_all_observables]; double stdev[N_all_observables];
	int i,k,o=0; divdiff_init(); divdiff dd(q+4,500); d=&dd; init();
  // performs QMC simulations while saving bitstrings every measurement step
  std::ostringstream oss_bs;
  oss_bs << "data/XXZ_" << "n_" << argv[1] << "_Jp_" << argv[2] << "_delta_" << argv[3] << "_bias_" << argv[4] << "_seed_" << rng_seed << "_bitstrings.txt";
  std::string sum_bs_fname = oss_bs.str();
  std::ofstream mybsfile (sum_bs_fname);
  if (mybsfile.is_open()) {
	for(step=0;step<Tsteps;step++) update();
	for(measurement_step=0;measurement_step<measurements;measurement_step++){
		for(step=0;step<stepsPerMeasurement;step++) update(); measure();
      mybsfile << lattice << std::endl;
	}
  }
  else std::cout << "Unable to open file";
  // bins sign as default measurement
	for(i=0;i<Nbins;i++) sgn_sum += bin_mean_sgn[i]; sgn_sum /= Nbins;
	for(i=0;i<Nbins;i++) over_bins_sum_sgn += (bin_mean_sgn[i] - sgn_sum)*(bin_mean_sgn[i] - sgn_sum); over_bins_sum_sgn /= (Nbins*(Nbins-1));
  // computing average/ std of all observables and printing to screen
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
      //std::cout << "Observable #" << ++o << ": "<< name_of_observable(k) << std::endl;
		for(i=0;i<Nbins;i++) Rsum[k] += bin_mean[k][i]; Rsum[k] /= Nbins;
		for(i=0;i<Nbins;i++) over_bins_sum[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean[k][i] - Rsum[k]); over_bins_sum[k] /= (Nbins*(Nbins-1));
		for(i=0;i<Nbins;i++) over_bins_sum_cov[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean_sgn[i] - sgn_sum); over_bins_sum_cov[k] /= (Nbins*(Nbins-1));
		mean[k] = Rsum[k]/sgn_sum*(1 + over_bins_sum_sgn/sgn_sum/sgn_sum) - over_bins_sum_cov[k]/sgn_sum/sgn_sum;
		stdev[k] = fabs(Rsum[k]/sgn_sum)*sqrt(over_bins_sum[k]/Rsum[k]/Rsum[k] + over_bins_sum_sgn/sgn_sum/sgn_sum - 2*over_bins_sum_cov[k]/Rsum[k]/sgn_sum);
		//std::cout << "mean(O) = " << mean[k] << std::endl;
		//std::cout << "std.dev.(O) = " << stdev[k] << std::endl;
	}
	divdiff_clear_up();

  // save things to file
  std::ostringstream oss;
  oss << "data/XXZ_" << "n_" << argv[1] << "_Jp_" << argv[2] << "_delta_" << argv[3] << "_bias_" << argv[4] << "_seed_" << rng_seed << "_summary.txt";
  std::string sum_fname = oss.str();
  std::ofstream myfile (sum_fname);
  if (myfile.is_open())
    {
      myfile << std::setprecision(9);
      // input data for simulations
      myfile << "RNG seed = " << rng_seed << std::endl;
      myfile << "Parameters: beta = " << beta << ", Tsteps = " << Tsteps << ", steps = " << steps << std::endl;
      myfile << "XXZ params: n = " << argv[1] << ", Jp = " << argv[2] << ", delta = " << argv[3] << ", bias = " << argv[4] <<  std::endl;
      myfile << "mean(sgn(W)) = " << sgn_sum << std::endl;
      myfile << "std.dev.(sgn(W)) = " << sqrt(over_bins_sum_sgn) << std::endl;
      // warning that qmax was needed
      if(qmax_achieved) std::cout << std::endl << "Warning: qmax = " << qmax << " was achieved. The results may be incorrect. The qmax parameter should be increased." << std::endl;
      // warning that a fundamental cycle was not used
      for(i=0;i<Ncycles;i++) if(!cycles_used[i]) myfile << "Warning: cycle No. " << i << " of length " << cycle_len[i] << " was not used" << std::endl;
      // information about average and max q used
      myfile << "mean(q) = " << meanq / measurements << std::endl;
      myfile << "max(q) = "<< maxq << std::endl;
      // add expectation values of valid observables measured
      for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
          myfile << "Observable #" << ++o << ": "<< name_of_observable(k) << std::endl;
          myfile << "mean(O) = " << mean[k] << std::endl;
          myfile << "std.dev.(O) = " << stdev[k] << std::endl;
        }
      myfile << "wall-clock cpu time = " << get_cpu_time()-start_time << " seconds" << std::endl;
      myfile.close();
    }
  else std::cout << "Unable to open file";

	return 0;
}
