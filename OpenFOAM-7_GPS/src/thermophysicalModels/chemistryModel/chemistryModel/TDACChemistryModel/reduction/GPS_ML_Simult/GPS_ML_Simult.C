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

\*---------------------------------------------------------------------------*/
#include "scope_guard.hpp"

#include "GPS_ML_Simult.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryReductionMethods::GPS_ML_Simult<CompType, ThermoType>::GPS_ML_Simult
(
    const IOdictionary& dict,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    chemistryReductionMethod<CompType, ThermoType>(dict, chemistry),
    searchInitSet_(this->coeffsDict_.subDict("initialSet").size())
{
    label j=0;

    dictionary initSet = this->coeffsDict_.subDict("initialSet");

    for (label i=0; i<chemistry.nSpecie(); i++)
    {
        if (initSet.found(chemistry.Y()[i].member()))
        {
            searchInitSet_[j++] = i;
        }
    }

    if (j<searchInitSet_.size())
    {
        FatalErrorInFunction
            << searchInitSet_.size()-j
            << " species in the initial set is not in the mechanism "
            << initSet
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryReductionMethods::GPS_ML_Simult<CompType, ThermoType>::~GPS_ML_Simult()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::chemistryReductionMethods::GPS_ML_Simult<CompType, ThermoType>::reduceMechanism
(
    const scalarField &c,
    const scalar T,
    const scalar p
)
{
}

template<class CompType, class ThermoType>
void Foam::chemistryReductionMethods::GPS_ML_Simult<CompType, ThermoType>::reduceMechanism
(
    const scalarField &c,
    const scalar T,
    const scalar p,
    const scalarField &data,
    const labelField &always_inds,
    const labelField &never_inds,
    const labelField &vari_inds
)
{   
    const int n_always_inds = always_inds.size();
    const int n_vari_inds = vari_inds.size();
//   Info<<"n_vari_inds: "<<n_vari_inds<<endl;
    
    label speciesNumber = 0;

/*
    Info<<"OUTPUT: ";
    for (label i=0; i<n_vari_inds; i++)
    {
        Info<<data[i]<<", ";
    }
    Info<<endl;
*/
    for (label i=0; i<this->nSpecie_; i++)
    {
        this->activeSpecies_[i] = false;
    }
    
    for (label i=0; i<n_always_inds; i++)
    {
        this->activeSpecies_[always_inds[i]] = true;
        speciesNumber++;
    }
    
    for (label i=0; i<n_vari_inds; i++) {
        if (data[i] > 0.5 && !this->activeSpecies_[vari_inds[i]])
        {
            this->activeSpecies_[vari_inds[i]] = true;
            speciesNumber++;
        }
    }
/*    
    Info<<"Species kept:, ";
    for (label i=0; i<53; i++) {
        if (this->activeSpecies_[i]) {
            Info<<1<<", ";
        }
        
        else {
            Info<<0<<", ";
        }
    }
    Info<<endl;
Info<<"End of ML"<<endl;
*/
//---------------------------Code update ends here ---------------------------    
//




    // Put a flag on the reactions containing at least one removed species
    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];
        this->chemistry_.reactionsDisabled()[i] = false;

        forAll(R.lhs(), s)
        {
            label ss = R.lhs()[s].index;
            if (!this->activeSpecies_[ss])
            {
                this->chemistry_.reactionsDisabled()[i] = true;
                break;
            }
        }
        if (!this->chemistry_.reactionsDisabled()[i])
        {
            forAll(R.rhs(), s)
            {
                label ss = R.rhs()[s].index;
                if (!this->activeSpecies_[ss])
                {
                    this->chemistry_.reactionsDisabled()[i] = true;
                    break;
                }
            }
        }
    }

    this->NsSimp_ = speciesNumber;
    scalarField& simplifiedC(this->chemistry_.simplifiedC());
    simplifiedC.setSize(this->NsSimp_+2);
    DynamicList<label>& s2c(this->chemistry_.simplifiedToCompleteIndex());
    s2c.setSize(this->NsSimp_);
    Field<label>& c2s(this->chemistry_.completeToSimplifiedIndex());

    label j = 0;
    for (label i=0; i<this->nSpecie_; i++)
    {
        if (this->activeSpecies_[i])
        {
            s2c[j] = i;
            simplifiedC[j] = c[i];
            c2s[i] = j++;
            if (!this->chemistry_.active(i))
            {
                this->chemistry_.setActive(i);
            }
        }
        else
        {
            c2s[i] = -1;
        }
    }
    simplifiedC[this->NsSimp_] = T;
    simplifiedC[this->NsSimp_+1] = p;
    this->chemistry_.setNsDAC(this->NsSimp_);
    // change temporary Ns in chemistryModel
    // to make the function nEqns working
    this->chemistry_.setNSpecie(this->NsSimp_);

//Info<<"End of GPS"<<endl;
}


// ************************************************************************* //
