//
// Copyright (c) 2014-23, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/SteadyStateSolver.h>
#include <Kripke.h>
#include <Kripke/Core/Comm.h>
#include <Kripke/Kernel.h>
#include <Kripke/ParallelComm.h>
#include <Kripke/Core/PartitionSpace.h>
#include <Kripke/Timing.h>
#include <Kripke/SweepSolver.h>
#include <vector>
#include <stdio.h>

#ifdef KRIPKE_USE_CALIPER
#include <caliper/cali.h>
#endif

using namespace Kripke::Core;

/**
  Run solver iterations.
*/
int Kripke::SteadyStateSolver (Kripke::Core::DataStore &data_store, size_t max_iter, bool block_jacobi)
{
  KRIPKE_TIMER(data_store, Solve);

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");

  Kripke::Core::Comm const &comm = data_store.getVariable<Kripke::Core::Comm>("comm");
  if(comm.rank() == 0){
    printf("\n");
    printf("Steady State Solve\n");
    printf("==================\n\n");
  }

  CALI_MARK_COMM_REGION_BEGIN("sweep_comm");

  // Intialize unknowns
  Kripke::Kernel::kConst(data_store.getVariable<Kripke::Field_Flux>("psi"), 0.0);

  // Loop over iterations
  double part_last = 0.0;
#ifdef KRIPKE_USE_CALIPER

  CALI_CXX_MARK_LOOP_BEGIN(mainloop_annotation, "solve");

#endif
  for(size_t iter = 0;iter < max_iter;++ iter){
#ifdef KRIPKE_USE_CALIPER
    CALI_CXX_MARK_LOOP_ITERATION(mainloop_annotation, static_cast<int>(iter));
#endif

    /*
     * Compute the RHS:  rhs = LPlus*S*L*psi + Q
     */


    // Discrete to Moments transformation (phi = L*psi)
    Kripke::Kernel::kConst(data_store.getVariable<Field_Moments>("phi"), 0.0);
    Kripke::Kernel::LTimes(data_store);

    // Compute Scattering Source Term (psi_out = S*phi)
    Kripke::Kernel::kConst(data_store.getVariable<Kripke::Field_Moments>("phi_out"), 0.0);
    Kripke::Kernel::scattering(data_store);


    // Compute External Source Term (psi_out = psi_out + Q)
    Kripke::Kernel::source(data_store);


    // Moments to Discrete transformation (rhs = LPlus*psi_out)
    Kripke::Kernel::kConst(data_store.getVariable<Kripke::Field_Flux>("rhs"), 0.0);
    Kripke::Kernel::LPlusTimes(data_store);


    /*
     * Sweep (psi = Hinv*rhs)
     */
     {
      // Create a list of all groups
      int num_subdomains = pspace.getNumSubdomains(SPACE_PQR);
      std::vector<SdomId> sdom_list(num_subdomains);
      for(SdomId i{0};i < num_subdomains;++ i){
        sdom_list[*i] = i;
      }

      // Sweep everything
      Kripke::SweepSolver(data_store, sdom_list, block_jacobi);
    }



    /*
     * Population edit and convergence test
     */
    double part = Kripke::Kernel::population(data_store);
    if(comm.rank() == 0){
      printf("  iter %d: particle count=%e, change=%e\n", (int)iter, part, (part-part_last)/part);
      fflush(stdout);
    }
    part_last = part;
  }
#ifdef KRIPKE_USE_CALIPER
  CALI_CXX_MARK_LOOP_END(mainloop_annotation);
#endif

  if(comm.rank() == 0){
    printf("  Solver terminated\n");
  }
  CALI_MARK_COMM_REGION_END("sweep_comm");
  return(0);
}




