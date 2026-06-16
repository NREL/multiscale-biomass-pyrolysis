/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "HybridErgunBiomass.H"
#include "Ergun.H"
#include "Ganser.H"
#include "HaiderLevenspiel.H"
#include "HolzerSommerfeld.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(HybridErgunBiomass, 0);
    addToRunTimeSelectionTable(dragModel, HybridErgunBiomass, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::HybridErgunBiomass::HybridErgunBiomass
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject),
    Ergun_(dict, interface, false)
{
    word bioDrag(dict.lookup("biomassDrag"));
    if (bioDrag == "Ganser")
    {
        Biomass_.set(new Ganser(dict, interface, false));
    }
    else if (bioDrag == "HaiderLevenspiel")
    {
        Biomass_.set(new HaiderLevenspiel(dict, interface, false));
    }
    else if (bioDrag == "HolzerSommerfeld")
    {
        Biomass_.set(new HolzerSommerfeld(dict, interface, false));
    }
    else
    {
        FatalErrorInFunction<< "Available biomassDrag laws are:\n"
                            << "    Ganser\n"
                            << "    HaiderLevenspiel\n"
                            << "    HolzerSommerfeld\n"
                            << abort(FatalError);
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::HybridErgunBiomass::~HybridErgunBiomass()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::HybridErgunBiomass::CdRe() const
{
    return
        pos0(interface_.continuous() - 0.8)*Biomass_->CdRe()
      + neg(interface_.continuous() - 0.8)*Ergun_.CdRe();
}


// ************************************************************************* //
