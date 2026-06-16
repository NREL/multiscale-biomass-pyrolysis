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

#include "Ganser.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Ganser, 0);
    addToRunTimeSelectionTable(dragModel, Ganser, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Ganser::Ganser
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject),
    sphericity_(readScalar(dict.lookup("sphericity")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Ganser::~Ganser()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::Ganser::CdRe() const
{
    const volScalarField alpha2
    (
        max(interface_.continuous(), interface_.continuous().residualAlpha())
    );

    const volScalarField Res(alpha2*interface_.Re());

    const scalar K1 = 3.0/(1.0 + 2.0/sqrt(sphericity_));

    const scalar K2 =
        pow
        (
            10.0,
            1.8148*pow(-log10(sphericity_), 0.5743)
        );

    const volScalarField CdRe
    (
        (24.0/K1)
       *(1.0 + 0.1118*pow(Res*K1*K2, 0.6567))
      + (0.4305*sqr(Res*K2)*K1)/(Res*K1*K2 + 3305.0)
    );

    return CdRe*pow(alpha2, -2.65);

}


// ************************************************************************* //
