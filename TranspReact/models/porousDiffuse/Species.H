#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#define SOLID_ID 0
#define TEMP_ID 1
#define NUM_SPECIES 2

namespace tr_species
{
    extern amrex::Vector<std::string> specnames;
    void init();
    void close();
    int find_id(std::string specname);
}

#endif
