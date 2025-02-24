#include "pyroModel.h"
#include "irreversibleArrheniusReaction.h"

using namespace Foam;
using namespace pyroModel;

pyroModel::pyroModel(const fvMesh& mesh)
:
    m_mesh(mesh),
    m_dict (
        IOobject
        (
            "pyrolisisDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    m_porosity
    (
        IOobject
        (
            "porosity",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        m_mesh
    )
    m_nSubTimeSteps(m_dict.lookupOrDefault<label>("nSubTimeSteps",1))
{
    Info << "Reading pyrolisis model" << endl;

    List<dictionary> speciesList(m_dict.lookup("species"));

    Info << "Found species: " << speciesList << endl;

    m_speciesName.resize(speciesList.size());
    m_species.resize(speciesList.size());
    m_reactionRates.resize(speciesList.size());

    m_rho.resize(speciesList.size());
    m_cp.resize(speciesList.size());
    m_kappa.resize(speciesList.size());
    m_molarWeight.resize(speciesList.size());
    m_is_gas.resize(speciesList.size());

    forAll(speciesList, specieI)
    {
        m_speciesName[specieI] = speciesList[specieI].dictName();
        m_species.set (
            specieI,
            new volScalarField
            (
                IOobject
                (
                    speciesList[specieI].dictName() + ".solid",
                    m_mesh.time().timeName(),
                    m_mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                m_mesh,
                dimensionedScalar("Y",dimless,0.)
            )
        );

        m_rho[specieI] = scalar(readScalar(specieslList[specieI].lookup("rho")));
        m_cp[specieI] = scalar(readScalar(specieslList[specieI].lookup("cp")));
        m_kappa[specieI] = scalar(readScalar(specieslList[specieI].lookup("kappa")));
        m_molarWeight[specieI] = scalar(readScalar(specieslList[specieI].lookup("molarWeight")));
        m_is_gas[specieI].first = speciesList[specieI].lookupOrDefault<bool>("gas",false);

        if (m_is_gas[specieI])
        {
            m_reactionRates.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        speciesList[specieI].dictName() + ".RR",
                        m_mesh.time().timeName(),
                        m_mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    m_mesh,
                    dimensionedScalar("Rdot",dimDensity/dimTime,0.)
                )            
            );

            m_is_gas[specieI].second = m_reactionRates.size() - 1;
        }
    }

    List<dictionary> reactionList(m_dict.lookup("reactions"));
    
    m_reactions.resize(reactionList.size());

    forAll(reactionList, reactI)
    {
        m_reactions.set (
            reactI,
            new irreversibleArrheniusReaction
            (
                reactionList[reactI],
                speciesList
            )
        );        
    }
}

pyroModel::~pyroModel()
{}

void evolve()
{
    Info << "Updating solid composition" << endl;

    const volScalarField& T = mesh.lookupObject<volScalarField>("T");

    scalarField Y(m_species.size(),0.);

    /* Reset reaction rate */
    forAll(m_species, specieI)
    {
        m_reactionRates[specieI] *= 0.;
    }

    /* Go cell-by-cell */
    forAll(m_mesh.C(), cellI)
    {
        // Skip if outside the solid
        if ( m_porosity[cellI] > 0.999 )
        {
            continue;
        }

        scalar dt = m_mesh.time().deltaT().value();

        scalar sub_dt = dt / m_nSubTimeSteps;

        forAll(m_species, specieI)
        {
            //- Convert from [vol_specie/m3] to [mol_specie/m3] by multiplying by density and dividing
            //  by molar weight
            Y[specieI] =  m_species[specieI][cellI] * ( m_rho[specieI] / m_molarWeight[specieI] );
        }

        /* Use Euler marching */
        for (label ts = 0; ts < m_nSubTimeSteps; ts++)
        {
            scalarField ndot(m_species.size(),0.);

            forAll(m_reactions, reactI)
            {
                ndot += m_reactions[reactI].computeMolarSources(T[cellI], Y);
            }

            Y += sub_dt * ndot;

            Y = max(Y,0.);
        }

        /* Update reaction rate and concentration*/
        forAll(m_species, specieI)
        {
            scalar Y0 = m_species[specieI][cellI] * m_rho[specieI];
            m_species[specieI][cellI] = Y[specieI] * m_molarWeight[specieI] / m_rho[specieI];
            
            if( m_is_gas[specieI].first )
            {
                m_reactionRates[m_is_gas[specieI].second][cellI] = ( (Y *  m_molarWeight[specieI] ) - Y_0 )  / dt;
                m_species[specieI][cellI] = 0.;             
            }
        }

        /* Update local porosity */
        scalar vol_species(0.);
        forAll(m_species, specieI)
        {
            vol_species += m_species[specieI][cellI];
        }
        
        m_porosity[cellI] = vol_species;
    }
}