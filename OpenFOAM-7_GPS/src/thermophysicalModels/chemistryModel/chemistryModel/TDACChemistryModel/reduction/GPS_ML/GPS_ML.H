/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

Class
    Foam::chemistryReductionMethods::GPS_ML

Description
    Machine Learning-Aided Global Pathway Selection

SourceFiles
    GPS_ML.C

\*---------------------------------------------------------------------------*/

#ifndef GPS_ML_H
#define GPS_ML_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
//#include "CloudSubModelBase.H"
#include "vector.H"
#include "TimeFunction1.H"
//#include "c_api_experimental.h"
#include "tf_utils.H" // Loads in tensorflow functionality

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryReductionMethods
{

/*---------------------------------------------------------------------------*\
                           Class ode Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class GPS_ML
:
    public chemistryReductionMethod<CompType, ThermoType>
{
    // Private data

        //- List of label for the search initiating set
        labelList searchInitSet_;
protected:

	TF_Graph* graph;
        TF_Output input_ph;
        TF_Output output;  


public:

    //- Runtime type information
    TypeName("GPS_ML");


    // Constructors

        //- Construct from components
        GPS_ML
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    //- Destructor
    virtual ~GPS_ML();


    // Member Functions

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p
        );
        
        virtual void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p,
            const scalarField &data,
            const labelField &always_inds,
            const labelField &never_inds,
            const labelField &vari_inds
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryReductionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GPS_ML.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
