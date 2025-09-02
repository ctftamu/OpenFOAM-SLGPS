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
        //tf_utils::LoadGraph("./sup_model.pb")
        //tf_utils::LoadGraph("./model_c_500_a_0.001_l_16_dt_200.pb")
        //tf_utils::LoadGraph("./model_TDAC_a_0.001_l_16.pb")
        //tf_utils::LoadGraph("./model_LOW_TEMP_a_0.001_l_16.pb")
        //tf_utils::LoadGraph("./model_LOW_TEMP_a_0.005_l_16.pb")
        //tf_utils::LoadGraph("./model_LOW_TEMP_a_0.000000000001_l_16.pb")
        //tf_utils::LoadGraph("./model_EXTRA_PROD_a_0.000000000001_l_16.pb")
        //tf_utils::LoadGraph("./model_EXTRA_PROD_a_0.0000000000001_l_16.pb")
        //tf_utils::LoadGraph("./model_EXTRA_PROD_a_0.0000000000005_l_16.pb")
        //tf_utils::LoadGraph("./model_EXTRA_PROD_a_0.000000000001_v2_l_16.pb")
        tf_utils::LoadGraph("./model_EXTRA_PROD_a_0.00000000000075_l_16.pb")
    ),

    // Defines an input operation - this must correspond to the name specified in python API
         input_ph
             (
                     //{TF_GraphOperationByName(graph, "dense_5_input"), 0} // the operation would look like "input_placeholder:0" in model.inputs[0] (keras)
                     //{TF_GraphOperationByName(graph, "dense_9_input"), 0} // the operation would look like "input_placeholder:0" in model.inputs[0] (keras)
                     //{TF_GraphOperationByName(graph, "dense_60_input"), 0} // the operation would look like "input_placeholder:0" in model.inputs[0] (keras)
                     //{TF_GraphOperationByName(graph, "dense_4_input"), 0} // the operation would look like "input_placeholder:0" in model.inputs[0] (keras)
                     {TF_GraphOperationByName(graph, "dense_2_input"), 0} // the operation would look like "input_placeholder:0" in model.inputs[0] (keras)
                     //{TF_GraphOperationByName(graph, "dense_6_input"), 0} // the operation would look like "input_placeholder:0" in model.inputs[0] (keras)
                     //{TF_GraphOperationByName(graph, "dense_4_input"), 0} // the operation would look like "input_placeholder:0" in model.inputs[0] (keras)
    
                         ),
                             // Defines an output operation - this must correspond to the name specified in python API
                                 output
                                     (
                                          //{TF_GraphOperationByName(graph, "dense_7/Sigmoid"), 0} // the operation would look like "output_value/BiasAdd:0" in model.outputs[0] (keras)
                                          //{TF_GraphOperationByName(graph, "dense_10/Sigmoid"), 0} // the operation would look like "output_value/BiasAdd:0" in model.outputs[0] (keras)
                                          //{TF_GraphOperationByName(graph, "dense_61/Sigmoid"), 0} // the operation would look like "output_value/BiasAdd:0" in model.outputs[0] (keras)
                                          //{TF_GraphOperationByName(graph, "dense_5/Sigmoid"), 0} // the operation would look like "output_value/BiasAdd:0" in model.outputs[0] (keras)
                                          {TF_GraphOperationByName(graph, "dense_3/Sigmoid"), 0} // the operation would look like "output_value/BiasAdd:0" in model.outputs[0] (keras)
                                          //{TF_GraphOperationByName(graph, "dense_7/Sigmoid"), 0} // the operation would look like "output_value/BiasAdd:0" in model.outputs[0] (keras)
                                          //{TF_GraphOperationByName(graph, "dense_5/Sigmoid"), 0} // the operation would look like "output_value/BiasAdd:0" in model.outputs[0] (keras)
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

    // Store the species composition according to the species index
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

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialized (by default 'active' is true)
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

            // Check if the species file is provided, if not set inactive
            // and NO_WRITE
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

}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void Foam::TDACChemistryModel<ReactionThermo, ThermoType>::omega
(
    const scalarField& c, // Contains all species even when mechRed is active
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
    const scalarField& c, // Contains all species even when mechRed is active
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

    // Find the matrix element and element position for the rhs
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
        // When using DAC, the ODE solver submit a reduced set of species the
        // complete set is used and only the species in the simplified mechanism
        // are updated
        this->c_ = completeC_;

        // Update the concentration of the species in the simplified mechanism
        // the other species remain the same and are used only for third-body
        // efficiencies
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

    // Constant pressure
    // dT/dt = ...
    scalar rho = 0;
    for (label i=0; i<this->c_.size(); i++)
    {
        const scalar W = this->specieThermo_[i].W();
        rho += W*this->c_[i];
    }

    scalar cp = 0;
    for (label i=0; i<this->c_.size(); i++)
    {
        // cp function returns [J/kmol/K]
        cp += this->c_[i]*this->specieThermo_[i].cp(p, T);
    }
    cp /= rho;

    // When mechanism reduction is active
    // dT is computed on the reduced set since dcdt is null
    // for species not involved in the simplified mechanism
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

        // ha function returns [J/kmol]
        const scalar hi = this->specieThermo_[si].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    dcdt[this->nSpecie_] = -dT;

    // dp/dt = ...
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

    // If the mechanism reduction is active, the computed Jacobian
    // is compact (size of the reduced set of species)
    // but according to the information of the complete set
    // (i.e. for the third-body efficiencies)

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

    // The species derivatives of the temperature term are partially computed
    // while computing dwdc, they are completed hereunder:
    scalar cpMean = 0;
    scalar dcpdTMean = 0;
    forAll(this->c_, i)
    {
        cpMean += this->c_[i]*cpi[i]; // J/(m^3 K)
        // Already multiplied by rho
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
                dTdt += hi[i]*dcdt[si]; // J/(m^3 s)
            }
        }
        else
        {
            dTdt += hi[i]*dcdt[i]; // J/(m^3 s)
        }
    }
    dTdt /= -cpMean; // K/s
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
        J(this->nSpecie_, i) += cpi[si]*dTdt; // J/(mol s)
        J(this->nSpecie_, i) /= -cpMean;    // K/s / (mol/m^3)
    }

    // ddT of dTdt
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
    // Increment counter of time-step
    timeSteps_++;

    const bool reduced = mechRed_->active();

    label nAdditionalEqn = (tabulation_->variableTimeStep() ? 1 : 0);

    basicSpecieMixture& composition = this->thermo().composition();

    // CPU time analysis
    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();
    scalar reduceMechCpuTime_ = 0;
    scalar addNewLeafCpuTime_ = 0;
    scalar growCpuTime_ = 0;
    scalar solveChemistryCpuTime_ = 0;
    scalar searchISATCpuTime_ = 0;

    this->resetTabulationResults();

    // Average number of active species
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

    // Composition vector (Yi, T, p)
    scalarField phiq(this->nEqns() + nAdditionalEqn);

    scalarField Rphiq(this->nEqns() + nAdditionalEqn);
    
    

    label num_chem_cells = 0;
    forAll(rho, celli)
    {
        num_chem_cells++;
    }
    
    const int n_always_inds = 12;
    int always_inds_arr[n_always_inds] = {1, 2, 3, 4, 5, 6, 7, 13, 15, 17, 34, 47};
    labelField always_inds(n_always_inds);
    for (label i=0; i<n_always_inds; i++)
    {  
        always_inds[i] = always_inds_arr[i];
    }
    
    
    const int n_never_inds = 1;
    int never_inds_arr[n_never_inds] = {48};
    labelField never_inds(n_never_inds);
    for (label i=0; i<n_never_inds; i++)
    {  
        never_inds[i] = never_inds_arr[i];
    }
    
    
    const int n_vari_inds = 40;
    int vari_inds_arr[n_vari_inds] = {0, 8, 9, 10, 11, 12, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52};
    labelField vari_inds(n_vari_inds);
    for (label i=0; i<n_vari_inds; i++)
    {  
        vari_inds[i] = vari_inds_arr[i];
    }
    
    
    const int n_input_inds = 12;
    int input_inds[n_input_inds-2] = {13, 5, 4, 1, 14, 3, 15, 2, 12, 9};
    
    float input_maxes[n_input_inds] = {2050.0, 1.0, 0.25157232704402227, 0.16594424796519067, 0.020616671507983172, 0.02995083692452817, 0.12241510069592972, 0.20790020790020763, 0.07873180240865373, 0.01777963214031066, 0.010632139480001737, 6.506121535281325e-05};
    float input_mins[n_input_inds] = {290.0, 1.0, 0.0, 4.1275912666522404e-32, 3.4510516978126006e-24, 9.310360588297054e-17, 5.097113825829686e-39, 0.0014552980509695756, 2.263388605568541e-39, 6.024805233091142e-22, 0.0, 0.0};
    
    float input_vals[num_chem_cells][n_input_inds];
    float output_vals[num_chem_cells][n_vari_inds];

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
        
        //if (!(tabulation_->active() && tabulation_->retrieve(phiq, Rphiq)))
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
        
        
    SCOPE_EXIT{ TF_DeleteStatus(status); }; // Auto-delete on scope exit.   
    auto options = TF_NewSessionOptions();
    SCOPE_EXIT{ TF_DeleteSessionOptions(options); }; // Auto-delete on scope exit.
    auto sess = TF_NewSession(graph, options, status);


    SCOPE_EXIT{ tf_utils::DeleteSession(sess); }; // Auto-delete on scope exit.
int num_cells = num_chem_cells;
Info<<"Number of chem cells:"<<num_chem_cells<<endl;

   //   TF_DeleteSessionOptions(options);
// Batch-wise implementation to avoid memory leak problem in TensorFlow TF_SessionRun

int batch_size = 1000000;

std::vector<float> output_data; // create a vector to hold the output data for all batches

for (int i = 0; i < num_cells; i += batch_size) {
  // create a sub-array for this batch
  int num_cells_batch = std::min(batch_size, num_cells - i);
  float input_vals_batch[num_cells_batch][n_input_inds];
  float output_vals_batch[num_cells_batch][n_vari_inds];
  for (int j = 0; j < num_cells_batch; j++) {
    for (int k = 0; k < n_input_inds; k++) {
      input_vals_batch[j][k] = input_vals[i+j][k];
    }
  }

  // create input and output tensors for this batch
  const std::vector<std::int64_t> input_dims_batch = {num_cells_batch, n_input_inds};
  const std::vector<std::int64_t> output_dims_batch = {num_cells_batch, n_vari_inds};
  auto input_tensor_batch = tf_utils::CreateTensor(TF_FLOAT,
                                                 input_dims_batch.data(), input_dims_batch.size(),
                                                 &input_vals_batch, num_cells_batch * n_input_inds * sizeof(float));
  auto output_tensor_batch = tf_utils::CreateTensor(TF_FLOAT,
                                                  output_dims_batch.data(), output_dims_batch.size(),
                                                  &output_vals_batch, num_cells_batch * n_vari_inds * sizeof(float));
  
  SCOPE_EXIT{ tf_utils::DeleteTensor(input_tensor_batch); }; // Auto-delete on scope exit.
  SCOPE_EXIT{ tf_utils::DeleteTensor(output_tensor_batch); }; // Auto-delete on scope exit.

  // call TF_SessionRun for this batch
  TF_SessionRun(sess,
                nullptr, // Run options.
                &input_ph, &input_tensor_batch, 1, // Input tensor ops, input tensor values, number of inputs.
                &output, &output_tensor_batch, 1, // Output tensor ops, output tensor values, number of outputs.
                nullptr, 0, // Target operations, number of targets.
                nullptr, // Run metadata.
                status // Output status.
                );
  
  // append the output data for this batch to the output_data array
  const auto batch_data = static_cast<float*>(TF_TensorData(output_tensor_batch));
  output_data.insert(output_data.end(), batch_data, batch_data + num_cells_batch * n_vari_inds);
}

// use the output_data variable to access the output data for all batches
const auto data = output_data.data();





// Batch-wise code ends here    






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


        // Initialise time progress
        scalar timeLeft = deltaT[celli];

        // Not sure if this is necessary
        Rphiq = Zero;

        clockTime_.timeIncrement();

        // When tabulation is active (short-circuit evaluation for retrieve)
        // It first tries to retrieve the solution of the system with the
        // information stored through the tabulation method
        if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
        {
            // Retrieved solution stored in Rphiq
            for (label i=0; i<this->nSpecie(); i++)
            {
                c[i] = rhoi*Rphiq[i]/this->specieThermo_[i].W();
            }

            searchISATCpuTime_ += clockTime_.timeIncrement();
        }
        // This position is reached when tabulation is not used OR
        // if the solution is not retrieved.
        // In the latter case, it adds the information to the tabulation
        // (it will either expand the current data or add a new stored point).
        else
        {
            // Store total time waiting to attribute to add or grow
            scalar timeTmp = clockTime_.timeIncrement();

            if (reduced)
            {
                // Reduce mechanism change the number of species (only active)
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

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                if (reduced)
                {
                    // completeC_ used in the overridden ODE methods
                    // to update only the active species
                    completeC_ = c;

                    // Solve the reduced set of ODE
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

            // If tabulation is used, we add the information computed here to
            // the stored points (either expand or add)
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

            // When operations are done and if mechanism reduction is active,
            // the number of species (which also affects nEqns) is set back
            // to the total number of species (stored in the mechRed object)
            if (reduced)
            {
                this->nSpecie_ = mechRed_->nSpecie();
            }
            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);
        }

        // Set the RR vector (used in the solver)
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
        // Every time-step, look if the tabulation should be updated
        tabulation_->update();

        // Write the performance of the tabulation
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
        // Write average number of species
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
    // Don't allow the time-step to change more than a factor of 2
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
