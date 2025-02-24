#pragma once

#include "fvCFD.H"
#include <utility>

namespace Foam
{

namespace pyroModel
{

class irreversibleArrheniusReaction;

struct pyroModel
{

    const fvMesh& m_mesh;
    IOdictionary m_dict;
    volScalarField m_porosity;
    label m_nSubTimeSteps;
    List<word> m_speciesName;
    PtrList<volScalarField> m_species;
    PtrList<volScalarField> m_reactionRates;
    scalarField m_rho;
    scalarField m_cp;
    scalarField m_molWeight;
    scalarField m_kappa;
    List<std::pair<bool,label>> m_is_gas; 
    PtrList<irreversibleArrheniusReaction>  m_reactions;

    pyroModel(const fvMesh& mesh);
    ~pyroModel();

    void evolve();
};

}
}
