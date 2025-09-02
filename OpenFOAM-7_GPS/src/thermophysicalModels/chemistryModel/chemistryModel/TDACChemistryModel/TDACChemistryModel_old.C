/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "TDACChemistryModel.H"
#include "UniformField.H"
#include "localEulerDdtScheme.H"
#include "clockTime.H"
#include "labelField.H"
#include "scope_guard.hpp"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::TDACChemistryModel<ReactionThermo, ThermoType>::TDACChemistryModel
(
    ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
    graph
    (
        tf_utils::LoadGraph("./SandiaD.pb")
    ),
    input_ph
    (
        {TF_GraphOperationByName(graph, "input_1"), 0}
    ),
    output
    (
        {TF_GraphOperationByName(graph, "dense_1/Sigmoid"), 0}
    ),
    variableTimeStep_
    (
        this->mesh().time().controlDict().lookupOrDefault
        (
            "adjustTimeStep",
            false
        )
     || fv::localEulerDdt::enabled(this->mesh())
    ),
    timeSteps_(0),
    NsDAC_(this->nSpecie_),
    completeC_(this->nSpecie_, 0),
    reactionsDisabled_(this->reactions_.size(), false),
    specieComp_(this->nSpecie_),
    completeToSimplifiedIndex_(this->nSpecie_, -1),
    simplifiedToCompleteIndex_(this->nSpecie_),
    tabulationResults_
    (
        IOobject
        (
            thermo.phasePropertyName("TabulationResults"),
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        scalar(0)
    )
{
    basicSpecieMixture& composition = this->thermo().composition();

    speciesTable speciesTab = composition.species();

    const HashTable<List<specieElement>>& specComp =
        dynamicCast<const reactingMixture<ThermoType>&>(this->thermo())
       .specieComposition();

    forAll(specieComp_, i)
    {
        specieComp_[i] = specComp[this->Y()[i].member()];
    }

    mechRed_ = chemistryReductionMethod<ReactionThermo, ThermoType>::New
    (
        *this,
        *this
    );

    if (mechRed_->active())
    {
        forAll(this->Y(), i)
        {
            IOobject header
            (
                this->Y()[i].name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ
            );

            if (!header.typeHeaderOk<volScalarField>(true))
            {
                composition.setInactive(i);
                this->Y()[i].writeOpt() = IOobject::NO_WRITE;
            }
        }
    }

    tabulation_ = chemistryTabulationMethod<ReactionThermo, ThermoType>::New
    (
        *this,
        *this
    );

    if (mechRed_->log())
    {
        cpuReduceFile_ = logFile("cpu_reduce.out");
        nActiveSpeciesFile_ = logFile("nActiveSpecies.out");
    }

    if (tabulation_->log())
    {
        cpuAddFile_ = logFile("cpu_add.out");
        cpuGrowFile_ = logFile("cpu_grow.out");
        cpuRetrieveFile_ = logFile("cpu_retrieve.out");
    }

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::TDACChemistryModel<ReactionThermo, ThermoType>::~TDACChemistryModel()
{
    tf_utils::DeleteGraph(graph);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void Foam::TDACChemistryModel<ReactionThermo, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalarField& dcdt
) const
{
    const bool reduced = mechRed_->active();

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    dcdt = Zero;

    forAll(this->reactions_, i)
    {
        if (!reactionsDisabled_[i])
        {
            const Reaction<ThermoType>& R = this->reactions_[i];

            scalar omegai = R.omega
            (
                p, T, c, pf, cf, lRef, pr, cr, rRef
            );

            forAll(R.lhs(), s)
            {
                label si = R.lhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sl = R.lhs()[s].stoichCoeff;
                dcdt[si] -= sl*omegai;
            }
            forAll(R.rhs(), s)
            {
                label si = R.rhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sr = R.rhs()[s].stoichCoeff;
                dcdt[si] += sr*omegai;
            }
        }
    }
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::TDACChemistryModel<ReactionThermo, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const scalar kf = R.kf(p, T, c);
    const scalar kr = R.kr(kf, p, T, c);

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s=1; s<Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(c[lRef], 0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(c[si], 0), exp);
        }
    }
    cf = max(c[lRef], 0);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1)
        {
            if (cf > small)
            {
                pf *= pow(cf, exp - 1);
            }
            else
            {
                pf = 0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    pr = kr;
    for (label s=1; s<Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(c[rRef], 0), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(c[si], 0), exp);
        }
    }
    cr = max(c[rRef], 0);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1)
        {
            if (cr > small)
            {
                pr *= pow(cr, exp - 1);
            }
            else
            {
                pr = 0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1);
        }
    }

    return pf*cf - pr*cr;
}


template<class ReactionThermo, class ThermoType>
void Foam::TDACChemistryModel<ReactionThermo, ThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    scalarField& dcdt
) const
{
    const bool reduced = mechRed_->active();

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    if (reduced)
    {
        this->c_ = completeC_;

        for (label i=0; i<NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(c[i], 0);
        }
    }
    else
    {
        for (label i=0; i<this->nSpecie(); i++)
        {
            this->c_[i] = max(c[i], 0);
        }
    }

    omega(this->c_, T, p, dcdt);

    scalar rho = 0;
    for (label i=0; i<this->c_.size(); i++)
    {
        const scalar W = this->specieThermo_[i].W();
        rho += W*this->c_[i];
    }

    scalar cp = 0;
    for (label i=0; i<this->c_.size(); i++)
    {
        cp += this->c_[i]*this->specieThermo_[i].cp(p, T);
    }
    cp /= rho;

    scalar dT = 0;
    for (label i=0; i<this->nSpecie_; i++)
    {
        label si;
        if (reduced)
        {
            si = simplifiedToCompleteIndex_[i];
        }
        else
        {
            si = i;
        }

        const scalar hi = this->specieThermo_[si].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    dcdt[this->nSpecie_] = -dT;

    dcdt[this->nSpecie_ + 1] = 0;
}


template<class ReactionThermo, class ThermoType>
void Foam::TDACChemistryModel<ReactionThermo, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    const bool reduced = mechRed_->active();

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    if (reduced)
    {
        this->c_ = completeC_;
        for (label i=0; i<NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(c[i], 0);
        }
    }
    else
    {
        forAll(this->c_, i)
        {
            this->c_[i] = max(c[i], 0);
        }
    }

    J = Zero;
    dcdt = Zero;
    scalarField hi(this->c_.size());
    scalarField cpi(this->c_.size());
    forAll(hi, i)
    {
        hi[i] = this->specieThermo_[i].ha(p, T);
        cpi[i] = this->specieThermo_[i].cp(p, T);
    }

    scalar omegaI = 0;

    forAll(this->reactions_, ri)
    {
        if (!reactionsDisabled_[ri])
        {
            const Reaction<ThermoType>& R = this->reactions_[ri];
            scalar kfwd, kbwd;
            R.dwdc
            (
                p,
                T,
                this->c_,
                J,
                dcdt,
                omegaI,
                kfwd,
                kbwd,
                reduced,
                completeToSimplifiedIndex_
            );
            R.dwdT
            (
                p,
                T,
                this->c_,
                omegaI,
                kfwd,
                kbwd,
                J,
                reduced,
                completeToSimplifiedIndex_,
                this->nSpecie_
            );
        }
    }

    scalar cpMean = 0;
    scalar dcpdTMean = 0;
    forAll(this->c_, i)
    {
        cpMean += this->c_[i]*cpi[i];
        dcpdTMean += this->c_[i]*this->specieThermo_[i].dcpdT(p, T);
    }

    scalar dTdt = 0;
    forAll(hi, i)
    {
        if (reduced)
        {
            const label si = completeToSimplifiedIndex_[i];
            if (si != -1)
            {
                dTdt += hi[i]*dcdt[si];
            }
        }
        else
        {
            dTdt += hi[i]*dcdt[i];
        }
    }
    dTdt /= -cpMean;
    dcdt[this->nSpecie_] = dTdt;

    for (label i = 0; i < this->nSpecie_; i++)
    {
        J(this->nSpecie_, i) = 0;
        for (label j = 0; j < this->nSpecie_; j++)
        {
            const label sj = reduced ? simplifiedToCompleteIndex_[j] : j;
            J(this->nSpecie_, i) += hi[sj]*J(j, i);
        }
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        J(this->nSpecie_, i) += cpi[si]*dTdt;
        J(this->nSpecie_, i) /= -cpMean;
    }

    J(this->nSpecie_, this->nSpecie_) = 0;
    for (label i = 0; i < this->nSpecie_; i++)
    {
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        J(this->nSpecie_, this->nSpecie_) +=
            cpi[si]*dcdt[i]
          + hi[si]*J(i, this->nSpecie_);
    }
    J(this->nSpecie_, this->nSpecie_) += dTdt*dcpdTMean;
    J(this->nSpecie_, this->nSpecie_) /= -cpMean;
    J(this->nSpecie_, this->nSpecie_) += dTdt/T;
}


template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::TDACChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    timeSteps_++;

    const bool reduced = mechRed_->active();

    label nAdditionalEqn = (tabulation_->variableTimeStep() ? 1 : 0);

    basicSpecieMixture& composition = this->thermo().composition();

    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();
    scalar reduceMechCpuTime_ = 0;
    scalar addNewLeafCpuTime_ = 0;
    scalar growCpuTime_ = 0;
    scalar solveChemistryCpuTime_ = 0;
    scalar searchISATCpuTime_ = 0;

    this->resetTabulationResults();

    scalar nActiveSpecies = 0;
    scalar nAvg = 0;

    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);

    scalarField phiq(this->nEqns() + nAdditionalEqn);
    scalarField Rphiq(this->nEqns() + nAdditionalEqn);
    
    label num_chem_cells = 0;
    forAll(rho, celli)
    {
        num_chem_cells++;
    }
    
    const int n_always_inds = 8;
    int always_inds_arr[n_always_inds] = {1, 3, 4, 5, 12, 13, 15, 47};
    labelField always_inds(n_always_inds);
    for (label i=0; i<n_always_inds; i++)
    {  
        always_inds[i] = always_inds_arr[i];
    }
    
    const int n_never_inds = 6;
    int never_inds_arr[n_never_inds] = {37, 39, 41, 42, 43,48};
    labelField never_inds(n_never_inds);
    for (label i=0; i<n_never_inds; i++)
    {  
        never_inds[i] = never_inds_arr[i];
    }
   
    const int n_vari_inds = 39;
    int vari_inds_arr[n_vari_inds] = {0, 2, 6, 7, 8, 9, 10, 11, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 40, 44, 45, 46, 49, 50, 51, 52};
    labelField vari_inds(n_vari_inds);
    for (label i=0; i<n_vari_inds; i++)
    {  
        vari_inds[i] = vari_inds_arr[i];
    }
    
    const int n_input_inds = 22;
    int input_inds[n_input_inds-2] = {13,	5,	4,	1,	14,	3,	15,	2,	12,	9,	0,	47,	6,	7,	10,	11,	16,	17,	19,	20};
   
    float input_maxes[n_input_inds] = {3.45E+03,	1.00E+00,	1.00E+00,	1.00E+00,	1.00E+00,	1.00E+00,	1.00E+00,	1.00E+00,	1.00E+00,	1.00E+00,	5.98E-04,	1.00E+00,	1.00E+00,	1.00E+00,	1.00E+00,	1.00E+00,	2.59E-01,	1.00E+00,	1.00E+00,	1.00E+00,	2.83E-04,	1.00E+00};
    float input_mins[n_input_inds] = {801.282148162797,	1.00E+00,	1.98E-27,	7.28E-27,	7.55E-27,	1.86E-20,	8.71E-21,	1.34E-06,	-1.68E-16,	-2.43E-15,	-7.40E-16,	-8.43E-19,	-1.26E-16,	-2.65E-17,	-3.30E-14,	-4.01E-16,	6.18E-33,	1.10E-07,	2.66E-39,	-9.15E-15,	-2.56E-17,	-1.42E-15};

    float input_vals[num_chem_cells][n_input_inds];

    label k = 0;
    scalar mol_frac = 0.0;

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        scalar pi = p[celli];
        scalar Ti = T[celli];

        for (label i=0; i<this->nSpecie_; i++)
        {
            c[i] = rhoi*this->Y_[i][celli]/this->specieThermo_[i].W();
            c0[i] = c[i];
            phiq[i] = this->Y()[i][celli];
        }
        
        input_vals[celli][0] = (Ti-input_mins[0])/(input_maxes[0]-input_mins[0]);
        input_vals[celli][1] = 0.0;
            
        scalar csum = 0.0;
        for (label i=0; i<this->nSpecie_; i++) 
        {
            csum += c[i];
        }

        for (label i=2; i < n_input_inds; i++) 
        {
            mol_frac = c[input_inds[i-2]]/csum;
            input_vals[celli][i] = (mol_frac - input_mins[i])/(input_maxes[i] - input_mins[i]);
        }
            
        k++;
    }

    auto status = TF_NewStatus();
    SCOPE_EXIT{ TF_DeleteStatus(status); };
    auto options = TF_NewSessionOptions();
    SCOPE_EXIT{ TF_DeleteSessionOptions(options); };
    auto sess = TF_NewSession(graph, options, status);
    SCOPE_EXIT{ tf_utils::DeleteSession(sess); };

    int num_cells = num_chem_cells;
    const std::vector<std::int64_t> input_dims = {num_cells, n_input_inds};
    
    auto input_tensor = tf_utils::CreateTensor(TF_FLOAT,
                                           input_dims.data(), input_dims.size(),
                                          &input_vals, num_cells*n_input_inds*sizeof(float));
 
    // Initialize output_tensor to nullptr. DO NOT allocate it.
    // TF_SessionRun will allocate it.
    TF_Tensor* output_tensor = nullptr;
    
    // Schedule the deletion of the tensors.
    SCOPE_EXIT{ tf_utils::DeleteTensor(input_tensor); };
    SCOPE_EXIT{ if (output_tensor) tf_utils::DeleteTensor(output_tensor); };

    TF_SessionRun(sess,
                nullptr, // Run options.
                &input_ph, &input_tensor, 1, // Input tensor ops, input tensor values, number of inputs.
                &output, &output_tensor, 1, // Output tensor ops, output tensor values, number of outputs.
                nullptr, 0, // Target operations, number of targets.
                nullptr, // Run metadata.
                status // Output status.
                );

    const auto data = static_cast<float*>(TF_TensorData(output_tensor));

    k = 0;
    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        scalar pi = p[celli];
        scalar Ti = T[celli];

        for (label i=0; i<this->nSpecie_; i++)
        {
            c[i] = rhoi*this->Y_[i][celli]/this->specieThermo_[i].W();
            c0[i] = c[i];
            phiq[i] = this->Y()[i][celli];
        }
        phiq[this->nSpecie()]=Ti;
        phiq[this->nSpecie() + 1]=pi;
        if (tabulation_->variableTimeStep())
        {
            phiq[this->nSpecie() + 2] = deltaT[celli];
        }

        scalar timeLeft = deltaT[celli];
        Rphiq = Zero;
        clockTime_.timeIncrement();

        if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
        {
            for (label i=0; i<this->nSpecie(); i++)
            {
                c[i] = rhoi*Rphiq[i]/this->specieThermo_[i].W();
            }

            searchISATCpuTime_ += clockTime_.timeIncrement();
        }
        else
        {
            scalar timeTmp = clockTime_.timeIncrement();

            if (reduced)
            {
                Info<<"celli: "<<celli<<endl;
                Info<<"INPUT: ";
                for (label i=0; i < n_input_inds; i++)
                {
                    Info<<input_vals[celli][i]<<", ";
                }
                Info<<endl;
                scalarField data_cell(n_vari_inds);
                for (label i=0; i<n_vari_inds; i++)
                {
                    data_cell[i] = data[n_vari_inds*celli+i];
                }
                mechRed_->reduceMechanism(c, Ti, pi, data_cell, always_inds, never_inds, vari_inds);
                nActiveSpecies += mechRed_->NsSimp();
                nAvg++;
                scalar timeIncr = clockTime_.timeIncrement();
                reduceMechCpuTime_ += timeIncr;
                timeTmp += timeIncr;
                k++;
 		        data_cell.clear();
            }

            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                if (reduced)
                {
                    completeC_ = c;
                    this->solve
                    (
                        simplifiedC_, Ti, pi, dt, this->deltaTChem_[celli]
                    );

                    for (label i=0; i<NsDAC_; i++)
                    {
                        c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                    }
                }
                else
                {
                    this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
                }
                timeLeft -= dt;
            }

            {
                scalar timeIncr = clockTime_.timeIncrement();
                solveChemistryCpuTime_ += timeIncr;
                timeTmp += timeIncr;
            }

            if (tabulation_->active())
            {
                forAll(c, i)
                {
                    Rphiq[i] = c[i]/rhoi*this->specieThermo_[i].W();
                }
                if (tabulation_->variableTimeStep())
                {
                    Rphiq[Rphiq.size()-3] = Ti;
                    Rphiq[Rphiq.size()-2] = pi;
                    Rphiq[Rphiq.size()-1] = deltaT[celli];
                }
                else
                {
                    Rphiq[Rphiq.size()-2] = Ti;
                    Rphiq[Rphiq.size()-1] = pi;
                }
                label growOrAdd =
                    tabulation_->add(phiq, Rphiq, rhoi, deltaT[celli]);
                if (growOrAdd)
                {
                    this->setTabulationResultsAdd(celli);
                    addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
                }
                else
                {
                    this->setTabulationResultsGrow(celli);
                    growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
                }
            }

            if (reduced)
            {
                this->nSpecie_ = mechRed_->nSpecie();
            }
            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);
        }

        for (label i=0; i<this->nSpecie_; i++)
        {
            this->RR_[i][celli] =
                (c[i] - c0[i])*this->specieThermo_[i].W()/deltaT[celli];
        }
    }

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_()
            << this->time().timeOutputValue()
            << "    " << solveChemistryCpuTime_ << endl;
    }

    if (mechRed_->log())
    {
        cpuReduceFile_()
            << this->time().timeOutputValue()
            << "    " << reduceMechCpuTime_ << endl;
    }

    if (tabulation_->active())
    {
        tabulation_->update();
        tabulation_->writePerformance();

        if (tabulation_->log())
        {
            cpuRetrieveFile_()
                << this->time().timeOutputValue()
                << "    " << searchISATCpuTime_ << endl;

            cpuGrowFile_()
                << this->time().timeOutputValue()
                << "    " << growCpuTime_ << endl;

            cpuAddFile_()
                << this->time().timeOutputValue()
                << "    " << addNewLeafCpuTime_ << endl;
        }
    }

    if (reduced && nAvg && mechRed_->log())
    {
        nActiveSpeciesFile_()
            << this->time().timeOutputValue()
            << "    " << nActiveSpecies/nAvg << endl;
    }

    if (Pstream::parRun())
    {
        List<bool> active(composition.active());
        Pstream::listCombineGather(active, orEqOp<bool>());
        Pstream::listCombineScatter(active);

        forAll(active, i)
        {
            if (active[i])
            {
                composition.setActive(i);
            }
        }
    }

    forAll(this->Y(), i)
    {
        if (composition.active(i))
        {
            this->Y()[i].writeOpt() = IOobject::AUTO_WRITE;
        }
    }

    return deltaTMin;
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::TDACChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::TDACChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ReactionThermo, class ThermoType>
void Foam::TDACChemistryModel<ReactionThermo, ThermoType>::
setTabulationResultsAdd
(
    const label celli
)
{
    tabulationResults_[celli] = 0;
}


template<class ReactionThermo, class ThermoType>
void Foam::TDACChemistryModel<ReactionThermo, ThermoType>::
setTabulationResultsGrow(const label celli)
{
    tabulationResults_[celli] = 1;
}


template<class ReactionThermo, class ThermoType>
void Foam::TDACChemistryModel<ReactionThermo, ThermoType>::
setTabulationResultsRetrieve
(
    const label celli
)
{
    tabulationResults_[celli] = 2;
}


// ************************************************************************* //
