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
  printf("I entered this loop beginning");
  KRIPKE_TIMER(data_store, Solve);

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");
    printf("after variable declaration");

  Kripke::Core::Comm const &comm = data_store.getVariable<Kripke::Core::Comm>("comm");
  if(comm.rank() == 0){
	  printf("Inside the loop to print comm rank");
    printf("\n");
    printf("Steady State Solve\n");
    printf("==================\n\n");
  }

   printf("I entered this loop mid loop");
  // Intialize unknowns
  Kripke::Kernel::kConst(data_store.getVariable<Kripke::Field_Flux>("psi"), 0.0);

  // Loop over iterations
  double part_last = 0.0;
#ifdef KRIPKE_USE_CALIPER
    printf("before cali begin");

  CALI_CXX_MARK_LOOP_BEGIN(mainloop_annotation, "solve");
     printf("after cali begin");

#endif
  for(size_t iter = 0;iter < max_iter;++ iter){
#ifdef KRIPKE_USE_CALIPER
	    printf("before cali loop iteration");
    CALI_CXX_MARK_LOOP_ITERATION(mainloop_annotation, static_cast<int>(iter));
      printf("after cali loop iteration");

#endif

    /*
     * Compute the RHS:  rhs = LPlus*S*L*psi + Q
     */


    // Discrete to Moments transformation (phi = L*psi)
    Kripke::Kernel::kConst(data_store.getVariable<Field_Moments>("phi"), 0.0);
    Kripke::Kernel::LTimes(data_store);

  printf("getting discrete transformation");


    // Compute Scattering Source Term (psi_out = S*phi)
    Kripke::Kernel::kConst(data_store.getVariable<Kripke::Field_Moments>("phi_out"), 0.0);
    Kripke::Kernel::scattering(data_store);
 printf("computing source term");


    // Compute External Source Term (psi_out = psi_out + Q)
    Kripke::Kernel::source(data_store);

printf("compute external source");

    // Moments to Discrete transformation (rhs = LPlus*psi_out)
    Kripke::Kernel::kConst(data_store.getVariable<Kripke::Field_Flux>("rhs"), 0.0);
    Kripke::Kernel::LPlusTimes(data_store);

printf("discrete transformation");





    /*
     * Sweep (psi = Hinv*rhs)
     */
    {printf("top of sweepsolver");
      // Create a list of all groups
      int num_subdomains = pspace.getNumSubdomains(SPACE_PQR);
      std::vector<SdomId> sdom_list(num_subdomains);
      for(SdomId i{0};i < num_subdomains;++ i){
        sdom_list[*i] = i;
      }

      // Sweep everything
      Kripke::SweepSolver(data_store, sdom_list, block_jacobi);
     printf("after sweep solver");
    }



printf("do we hang after sweep solver");
    /*
     * Population edit and convergence test
     */
    double part = Kripke::Kernel::population(data_store);
    if(comm.rank() == 0){
	    printf("convergent test");
      printf("  iter %d: particle count=%e, change=%e\n", (int)iter, part, (part-part_last)/part);
      fflush(stdout);
      printf("convergent test passed");
    }
    part_last = part;

    printf("did we hang in the test loop?");

  }
#ifdef KRIPKE_USE_CALIPER
  printf("before cali end");
  CALI_CXX_MARK_LOOP_END(mainloop_annotation);
  printf("after cali end");
#endif

  if(comm.rank() == 0){
    printf("  Solver terminated\n");
  }
printf("at the end of everything");
  return(0);
}




