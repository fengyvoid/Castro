#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const int coord_type = dgeom.Coord();
    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    if (problem::num_vortices > problem::max_num_vortices) {
        amrex::Error("num_vortices too large, please increase max_num_vortices and the size of xloc_vortices");
    }

    //set problem center
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        if (coord_type == 1 || coord_type == 2) {
            problem::center[n] = 0.0_rt;
        }
        else {
            problem::center[n] = 0.5_rt * (problo[n] + probhi[n]); //why
        }
    }

    // Read initial model
    read_model_file(problem::model_name);

    if (ParallelDescriptor::IOProcessor()) {
        for (int i = 0; i < model::npts; i++) {
            std::cout << i << " " << model::profile(0).r(i) << " " << model::profile(0).state(i,model::idens) << std::endl;
        }
    }

    // velocity perturbation stuff

    Real offset = (probhi[0] - problo[0]) / problem::num_vortices;

    for (int i = 0; i < problem::num_vortices; i++) {
        problem::xloc_vortices[i] = (static_cast<Real>(i) + 0.5_rt) * offset + problo[0];
    }

    // Set up Castro data logs for this problem

    if (castro::sum_interval > 0 && amrex::ParallelDescriptor::IOProcessor()) {

        Castro::problem_data_logs.resize(1);

        Castro::problem_data_logs[0].reset(new std::fstream);
        Castro::problem_data_logs[0]->open("problem_diag.out", std::ios::out | std::ios::app);
        if (!Castro::problem_data_logs[0]->good()) {
            amrex::FileOpenFailed("problem_diag.out");
        }

    }

}

#endif
