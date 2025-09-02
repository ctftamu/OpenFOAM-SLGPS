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

#include "GPS_ML.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryReductionMethods::GPS_ML<CompType, ThermoType>::GPS_ML
(
    const IOdictionary& dict,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    chemistryReductionMethod<CompType, ThermoType>(dict, chemistry),
     graph
    (
        tf_utils::LoadGraph("./sup_model.pb")
    ),

    // Defines an input operation - this must correspond to the name specified in python API
         input_ph
             (
                     {TF_GraphOperationByName(graph, "dense_5_input"), 0} // the operation would look like "input_placeholder:0" in model.inputs[0] (keras)
    
                         ),
                             // Defines an output operation - this must correspond to the name specified in python API
                                 output
                                     (
                                          {TF_GraphOperationByName(graph, "dense_7/Sigmoid"), 0} // the operation would look like "output_value/BiasAdd:0" in model.outputs[0] (keras)
                                                 ),
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
Foam::chemistryReductionMethods::GPS_ML<CompType, ThermoType>::~GPS_ML()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::chemistryReductionMethods::GPS_ML<CompType, ThermoType>::reduceMechanism
(
    const scalarField &c,
    const scalar T,
    const scalar p
)
{   
    
    
    
    const int n_always_inds = 13;
    int always_inds[n_always_inds] = {0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 47};
    const int n_never_inds = 5;
    int never_inds[n_never_inds] = {33, 41, 42, 43, 48};
    const int n_vari_inds = 35;
    int vari_inds[n_vari_inds] = {8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 49, 50, 51, 52};
    
    const int n_input_inds = 12;
    int input_inds[n_input_inds-2] = {13, 5, 4, 1, 14, 3, 15, 2, 12, 9};
    
    float input_maxes[n_input_inds] = {3294.3137, 100.0, 0.1282, 0.1806, 0.02904, 0.03321, 0.09929, 0.1976, 0.06597, 0.02169, 0.01106, 8.2008e-05};
    float input_mins[n_input_inds] = {1299.9992, 1.0, 6.3475e-19, 1.7356e-13, 1.5037e-12, 7.2028e-12, 1.9753e-19, 0.0005175, 2.4630e-21, 1.4075e-12, 5.3383e-18, 7.4433e-23};
    
    
    
    float input_vals[1][n_input_inds];
    float output_vals[1][n_vari_inds] = {0};
    scalar csum = 0.0;
    for (label i=0; i<53; i++) {
        csum += c[i];
    }


scalar mol_frac = 0.0;

for (label i=2; i<n_input_inds; i++) {



       mol_frac = c[input_inds[i-2]]/csum;

       input_vals[0][i] = (mol_frac - input_mins[i])/(input_maxes[i] - input_mins[i]);
}
    input_vals[0][0] = (T-input_mins[0])/(input_maxes[0] - input_mins[0]);
    
    input_vals[0][1] = ((p/101325)-input_mins[1])/(input_maxes[1] - input_mins[1]);
    
//  TF_Tensor* output_tensor = nullptr;
   
  //auto data = static_cast<float*>(TF_TensorData(output_tensor));   
    
//    auto graph = tf_utils::LoadGraph("sup_model.pb");
//    auto input_ph = TF_Output{TF_GraphOperationByName(graph, "dense_5_input"), 0};
//   auto output = TF_Output{TF_GraphOperationByName(graph, "dense_7/Sigmoid"), 0};
   

    auto status = TF_NewStatus();
 SCOPE_EXIT{ TF_DeleteStatus(status); }; // Auto-delete on scope exit.   
 auto options = TF_NewSessionOptions();
SCOPE_EXIT{ TF_DeleteSessionOptions(options); }; // Auto-delete on scope exit.
    auto sess = TF_NewSession(graph, options, status);


SCOPE_EXIT{ tf_utils::DeleteSession(sess); }; // Auto-delete on scope exit.

//   TF_DeleteSessionOptions(options);
    
   int num_cells = 1;
    const std::vector<std::int64_t> input_dims = {num_cells, n_input_inds};
    
    const std::vector<std::int64_t> output_dims = {num_cells, n_vari_inds};
    
    //TF_Tensor* output_tensor = nullptr;
    auto input_tensor = tf_utils::CreateTensor(TF_FLOAT,
                                           input_dims.data(), input_dims.size(),
                                          &input_vals, num_cells*n_input_inds*sizeof(float));
   
   auto output_tensor = tf_utils::CreateTensor(TF_FLOAT,
                                           output_dims.data(), output_dims.size(),
                                          &output_vals, num_cells*n_vari_inds*sizeof(float));

   SCOPE_EXIT{ tf_utils::DeleteTensor(input_tensor); }; // Auto-delete on scope exit.
    SCOPE_EXIT{ tf_utils::DeleteTensor(output_tensor); }; // Auto-delete on scope exit.


    TF_SessionRun(sess,
                nullptr, // Run options.
                &input_ph, &input_tensor, 1, // Input tensor ops, input tensor values, number of inputs.
                &output, &output_tensor, 1, // Output tensor ops, output tensor values, number of outputs.
                nullptr, 0, // Target operations, number of targets.
                nullptr, // Run metadata.
                status // Output status.
                );


//tf_utils::PrintOp(graph);
//tf_utils::PrintTensorInfo(graph, "dense_7/Sigmoid");



//float data[12] = TF_TensorData(output_tensor);
const auto data = static_cast<float*>(TF_TensorData(output_tensor));
  
//    Info<<"Output_tensor:"<<TF_TensorData(output_tensor)<<endl;

 //   Info<<"here is some output8"<<endl; 
  //  Info<<"Status: "<<TF_GetCode(status)<<endl;
   // Info<<"Status: "<<TF_Message(status)<<endl;

/*
   if (TF_GetCode(status) != TF_OK) {
    std::cout << "Error run session";
    }
    
    TF_CloseSession(sess, status);
    if (TF_GetCode(status) != TF_OK) {
        std::cout << "Error close session";
    }
    
    TF_DeleteSession(sess, status);
    if (TF_GetCode(status) != TF_OK) {
        std::cout << "Error delete session";
    }

    
//   auto data = static_cast<float*>(TF_TensorData(output_tensor));
//   Info<<"Size of data: "<<sizeof(data)<<endl;
//    Info<<"Size of output tensor: "<<sizeof(TF_TensorData(output_tensor))<<endl;
//    Info<<"Size of c: "<<sizeof(c)<<endl;
//    Info<<"Size of inputs: "<<sizeof(input_vals[0])<<endl;
    
    
    
    tf_utils::DeleteTensor(input_tensor);
	tf_utils::DeleteTensor(output_tensor);
    TF_DeleteStatus(status);
//}

// auto data = static_cast<float*>(TF_TensorData(output_tensor));

*/
    int speciesNumber = 0;
    Info<<"INPUT: ";
    for (label i=0; i<n_input_inds; i++)
    {
        Info<<input_vals[0][i]<<", ";
    }
    Info<<endl;
    Info<<"OUTPUT: ";
    for (label i=0; i<n_vari_inds; i++)
    {
        Info<<data[i]<<", ";
    }
    Info<<endl;
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
    
Info<<"End of ML"<<endl;   
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
    Info<<this->nSpecie_<<" species reduced to "<<this->NsSimp_<<endl;
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

Info<<"End of GPS"<<endl;
}


// ************************************************************************* //
