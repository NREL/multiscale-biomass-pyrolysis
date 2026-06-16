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

#include "HaiderLevenspiel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(HaiderLevenspiel, 0);
    addToRunTimeSelectionTable(dragModel, HaiderLevenspiel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::HaiderLevenspiel::HaiderLevenspiel
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

Foam::dragModels::HaiderLevenspiel::~HaiderLevenspiel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::HaiderLevenspiel::CdRe() const
{
    const volScalarField alpha2
    (
        max(interface_.continuous(), interface_.continuous().residualAlpha())
    );

    const volScalarField Res(alpha2*interface_.Re());

    const scalar A =
    exp(2.3288 - 6.4581*sphericity_ + 2.4486*sqr(sphericity_));

    const scalar B =
        0.0964 + 0.5565*sphericity_;

    const scalar C =
        exp
        (
            4.905
          - 13.8944*sphericity_
          + 18.4222*pow(sphericity_,2)
          - 10.2599*pow(sphericity_,3)
        );

    const scalar D =
        exp
        (
            1.4681
          + 12.2584*sphericity_
          - 20.7322*pow(sphericity_,2)
          + 15.8855*pow(sphericity_,3)
        );

    const volScalarField CdRe
    (
        24.0*(1.0 + A*pow(Res, B))
      + (C*Res)/(1.0 + D/Res)
    );

    return CdRe*pow(alpha2, -2.65);

}


// ************************************************************************* //
