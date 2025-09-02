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
    //Info<<c[13]<<endl;
    //Info<<T<<endl;
    //Info<<p<<endl;
    
    
    
    
    //Info<<"here is some output1"<<endl;
    //const int n_always_inds = 13;
    //int always_inds[n_always_inds] = {0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 47};
    const int n_always_inds = 13;
    int always_inds[n_always_inds] = {0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 47};
    //const int n_never_inds = 5;
    //int never_inds[n_never_inds] = {33, 41, 42, 43, 48};
    const int n_never_inds = 8;
    int never_inds[n_never_inds] = {33, 39, 41, 42, 43, 44, 45, 48};
    //const int n_vari_inds = 35;
    //int vari_inds[n_vari_inds] = {8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 49, 50, 51, 52};
    const int n_vari_inds = 32;
    int vari_inds[n_vari_inds] = {8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 40, 46, 49, 50, 51, 52};
    
    //const int n_input_inds = 12;
    //int input_inds[n_input_inds-2] = {13, 5, 4, 1, 14, 3, 15, 2, 12, 9};
    const int n_input_inds = 12;
    int input_inds[n_input_inds-2] = {13, 5, 4, 1, 14, 3, 15, 2, 12, 9};
    
    //float input_maxes[n_input_inds] = {3294.3137, 100.0, 0.1282, 0.1806, 0.02904, 0.03321, 0.09929, 0.1976, 0.06597, 0.02169, 0.01106, 8.2008e-05};
    //float input_mins[n_input_inds] = {1299.9992, 1.0, 6.3475e-19, 1.7356e-13, 1.5037e-12, 7.2028e-12, 1.9753e-19, 0.0005175, 2.4630e-21, 1.4075e-12, 5.3383e-18, 7.4433e-23};
    float input_maxes[n_input_inds] = {3219.9131346732624, 98.88285621310023, 0.12812228934173156, 0.18619097575556737, 0.02879275852922004, 0.035414855595551846, 0.09861926058935107, 0.1976103193196484, 0.07393701039839022, 0.02310716016769376, 0.012957167247582675, 0.00011822748907861104};
    float input_mins[n_input_inds] = {1300.0604102031048, 0.10013999234648284, 2.679008231681836e-19, 6.495569289946552e-33, 2.829678300694169e-25, 1.1796952977710393e-17, 5.945326919960085e-41, 0.00011732549819055488, 1.718204883434161e-41, 2.145338540125709e-23, 6.315723958055864e-18, 2.8714414299479102e-34};
    
    //Info<<"here is some output2"<<endl;
    
    
    
    float input_vals[1][n_input_inds];
    //Info<<"here is some output2.1"<<endl;
    float output_vals[n_vari_inds];
    //Info<<"here is some output2.2"<<endl;
    scalar csum = 0.0;
    for (label i=0; i<53; i++) {
        csum += c[i];
    }
    
    for (label i=2; i<n_input_inds; i++) {
        //Info<<"Success at loop #"<<i+1<<" at start"<<endl;
        scalar mol_frac = c[input_inds[i-2]]/csum;
        //Info<<"Success at loop #"<<i+1<<" at end"<<endl;
        //Info<<input_vals<<endl;
        //if (mol_frac>input_mins[i] && mol_frac<input_maxes[i]) {
        //for (label j=0; j<sizeof(c); j++) {
        input_vals[0][i] = (mol_frac - input_mins[i])/(input_maxes[i] - input_mins[i]);
        //}
        //}
        
        //else if (mol_frac<input_mins[i]) {
        //    for (label j=0; j<sizeof(c); j++) {
        //        input_vals[0][i] = 0.0;
        //    }
        //}
        //
        //else {
        //    for (label j=0; j<sizeof(c); j++) {
        //        input_vals[0][i] = 1.0;
        //    }
        //}
    }
    //for (label j=0; j<sizeof(T); j++) {
    //        input_vals[n_input_inds-1][j] = (&T[j]-input_mins[n_input_inds-1])/(input_maxes[n_input_inds-1]-input_mins[n_input_inds-1]);
    //        input_vals[n_input_inds][j] = (&p[j]-input_mins[n_input_inds])/(input_maxes[n_input_inds]-input_mins[n_input_inds]);
    //    }
    
    //Info<<"here is some output2.3"<<endl;
    //if (T>input_mins[0] && T<input_maxes[0]) {
    input_vals[0][0] = (T-input_mins[0])/(input_maxes[0] - input_mins[0]);
    //    Info<<"T: "<<input_mins[0]<<", "<<input_vals[0][0]<<", "<<input_maxes[0]<<endl;
    //}
    
    //else if (T<input_mins[0]) {
    //    input_vals[0][0] = 0.0;
    //}
    
    //else {
    //    input_vals[0][0] = 1.0;
    //}
    //Info<<"here is some output2.4"<<endl;
    //if ((p/101325)>input_mins[1] && (p/101325)<input_maxes[1]) {
    input_vals[0][1] = ((p/101325)-input_mins[1])/(input_maxes[1] - input_mins[1]);
    //}
    
    //else if ((p/101325)<input_mins[0]) {
    //    input_vals[0][1] = 0.0;
    //}
    
    //else {
    //    input_vals[0][1] = 1.0;
    //}
    //for (label i=0; i<n_input_inds; i++) {
    //    Info<<input_vals[i]<<endl;
    //}
    
    Info<<"here is some output3"<<endl;
    
    auto graph = tf_utils::LoadGraph("model_c_500_a_0.001_l_16_dt_200.pb");
    //auto graph = tf_utils::LoadGraph("sup_model.pb");
    Info<<"here is some output3.1"<<endl;
    //Info<<"here is some output3.1"<<endl;
    //Info<<"here is some output3.1"<<endl;
    //Info<<"here is some output3.1"<<endl;
    //Info<<"here is some output3.1"<<endl;
    auto input_ph = TF_Output{TF_GraphOperationByName(graph, "dense_input_2"), 0};
    //auto input_ph = TF_Output{TF_GraphOperationByName(graph, "dense_5_input"), 0};
    Info<<"here is some output3.2"<<endl;
    //Info<<"here is some output3.2"<<endl;
    //Info<<"here is some output3.2"<<endl;
    //Info<<"here is some output3.2"<<endl;
    //Info<<"here is some output3.2"<<endl;
    auto output = TF_Output{TF_GraphOperationByName(graph, "dense_1/Sigmoid"), 0};
    //auto output = TF_Output{TF_GraphOperationByName(graph, "dense_7/Sigmoid"), 0};
    
    Info<<"here is some output4"<<endl;
    
    auto status = TF_NewStatus();
    auto options = TF_NewSessionOptions();
    auto sess = TF_NewSession(graph, options, status);
    TF_DeleteSessionOptions(options);
    
    Info<<"here is some output5"<<endl;
    
    const std::vector<std::int64_t> input_dims = {1, n_input_inds};
    
    Info<<"here is some output6"<<endl;
    
    //Info<<"Input array: "<<input_vals<<endl;
    
    TF_Tensor* output_tensor = nullptr;
    auto input_tensor = tf_utils::CreateTensor(TF_FLOAT,
                                           input_dims.data(), input_dims.size(),
                                          &input_vals, 1*(n_input_inds)*sizeof(float));
    Info<<"here is some output7"<<endl;                                      
    //Info<<"Input tensor: "<<input_tensor<<endl;
    for (label i=0; i<2; i++) {
        TF_SessionRun(sess,
                    nullptr, // Run options.
                    &input_ph, &input_tensor, 1, // Input tensor ops, input tensor values, number of inputs.
                    &output, &output_tensor, 1, // Output tensor ops, output tensor values, number of outputs.
                    nullptr, 0, // Target operations, number of targets.
                    nullptr, // Run metadata.
                    status // Output status.
                    );
    }
    Info<<"here is some output8"<<endl; 
    Info<<"Status: "<<TF_GetCode(status)<<endl;
    Info<<"Status: "<<TF_Message(status)<<endl;
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
    
    //Info<<"Output tensor: "<<output_tensor<<endl;
    Info<<"here is some output9"<<endl;
    auto data = static_cast<float*>(TF_TensorData(output_tensor));
    Info<<"Size of data: "<<sizeof(data)<<endl;
    Info<<"Size of output tensor: "<<sizeof(TF_TensorData(output_tensor))<<endl;
    Info<<"Size of c: "<<sizeof(c)<<endl;
    Info<<"Size of inputs: "<<sizeof(input_vals[0])<<endl;
    
    
    
    //Info<<"Output data: "<<data<<endl;
    Info<<"here is some output10"<<endl;
    tf_utils::DeleteTensor(input_tensor);
	tf_utils::DeleteTensor(output_tensor);
    TF_DeleteStatus(status);
    
    int speciesNumber = 0;
    Info<<"here is some output11"<<endl;
    Info<<"INPUT: "<<endl;
    for (label i=0; i<n_input_inds+4; i++)
    {
        Info<<input_vals[0][i]<<", ";
    }
    Info<<endl;
    Info<<"OUTPUT: ,";
    for (label i=0; i<n_vari_inds+4; i++)
    {
        Info<<data[i]<<", ";
    }
    Info<<endl;
    //Info<<"here is some output9"<<endl;
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
    
    //for (label i=0; i<this->nSpecie_; i++)
    //{
    //    Info<<"Species Active: "<<this->activeSpecies_[i]<<endl;
    //}
//    int num_inputs = 2;
//  int num_outputs = 4;
//   int num_cells = 101;

//       auto graph = tf_utils::LoadGraph("./alpha.pb") ;
//        auto input_ph = TF_Output{TF_GraphOperationByName(graph, "dense_13_input"), 0};
//        auto output = TF_Output{TF_GraphOperationByName(graph, "dense_24/BiasAdd"), 0};
//
//if (input_ph.oper == nullptr) {
//    std::cout << "Can't init input_ph" << std::endl;
//  }

//  if (output.oper == nullptr) {
//    std::cout << "Can't init output" << std::endl;
//  }
//Info<<"graph:"<<graph<<endl;



//auto status = TF_NewStatus();
//  auto options = TF_NewSessionOptions();
//  auto sess = TF_NewSession(graph, options, status);
//  TF_DeleteSessionOptions(options);
//if (TF_GetCode(status) != TF_OK) {
//Info<<"Graph not loaded successfully"<<endl;  

//}

//float output_vals[num_cells][4];
  
//       float input_vals[num_cells][num_inputs];
//float input_vals2[num_cells][num_inputs];
//const std::vector<std::int64_t> input_dims = {num_cells, num_inputs};
//const std::vector<std::int64_t> input_dims2 = {num_cells,num_inputs};
//  double mean_array[num_inputs+num_outputs];
// double std_array[num_inputs+num_outputs] ;
//scalar currentTime = 5;
//scalar mean1 = 50.0;
//scalar mean2 = 0.00200875;
//scalar mean3 = 0.8495495147524752;
//scalar mean4 = 0.0016548139603960396;
//scalar mean5 = -0.020853931683168316;
//scalar mean6 = -5.395049504950494e-07;

//scalar std1 = 29.154759474226502;
//scalar std2 = 0.0005566328270652627;
//scalar std3 = 0.25465459284723363;
//scalar std4 = 2.799978925236133e-05;
//scalar std5 =  4.8795400145759155e-05;
//scalar std6 = 5.020291105078218e-05;

//  for (int i=0; i<=num_cells; i++)
//  {
//input_vals2[i][0] = (i - mean1)/std1;
//Info<<"input_valsTAG:"<<input_vals2[i][0];
//input_vals2[i][1] = (currentTime - mean2)/std2;
//Info<<"Time:"<<input_vals2[i][1];
//}


  





// Set up TF C API stuff
//    TF_Tensor* output_tensor = nullptr;
    


//auto input_tensor = tf_utils::CreateTensor(TF_FLOAT,
//                                           input_dims2.data(), input_dims2.size(),
//                                          &input_vals2, num_cells*num_inputs*sizeof(float));
//                                         


   
  


//Info<<"input tensors:"<<input_tensor<<endl<<"output tensors:"<<output_tensor<<endl;
//Info<<"session:"<<sess<<endl<<"status:"<<status<<endl;
//    TF_SessionRun(sess,
//                nullptr, // Run options.
//                &input_ph, &input_tensor, 1, // Input tensor ops, input tensor values, number of inputs.
//                &output, &output_tensor, 1, // Output tensor ops, output tensor values, number of outputs.
//                nullptr, 0, // Target operations, number of targets.
//                nullptr, // Run metadata.
//                status // Output status.
//                );
//Info<<"TF-6"<<endl;

//Info<<"Output_tensor"<<output_tensor<<endl;


//if (TF_GetCode(status) != TF_OK) {
//    std::cout << "Error run session";
//  }

//  TF_CloseSession(sess, status);
//  if (TF_GetCode(status) != TF_OK) {
//    std::cout << "Error close session";
//  }

//  TF_DeleteSession(sess, status);
//  if (TF_GetCode(status) != TF_OK) {
//    std::cout << "Error delete session";
//  }


// auto data = static_cast<float*>(TF_TensorData(output_tensor));


//for (int i = 0; i < num_cells; i++)
//	{
//
//		output_vals[i][0] = data[i]*std3 + mean3; // Funnel changes back into OF - row major order	
//                output_vals[i][1] = data[num_cells + i]*std4 + mean4; 
//                 output_vals[i][2] = data[num_cells*2 + i]*std5 + mean5;
//		 output_vals[i][3] = data[num_cells*3 + i]*std6 + mean6;
//
//     }
//
/*
Info << "alpha(min/max):"<<min(output_vals[][0])<<max(output_vals[][0])<<endl;
Info << "Ux(min/max):"<<min(output_vals[][1])<<max(output_vals[][1])<<endl;
Info << "Uy(min/max):"<<min(output_vals[][2])<<max(output_vals[][2])<<endl;
Info << "Uz(min/max):"<<min(output_vals[][3])<<max(output_vals[][3])<<endl;
*/
// Normalize the output data
//
//
//	tf_utils::DeleteTensor(input_tensor);
//	tf_utils::DeleteTensor(output_tensor);
//    TF_DeleteStatus(status);
//    
//int speciesNumber = 0;
//    
//    for (label i=0; i<this->nSpecie_; i++)
//    {
//        this->activeSpecies_[i] = false;
//    }
//    for (label i=0; i<10; i++)
//    {
//        this->activeSpecies_[i] = true;
//        speciesNumber++;
//    }

    
    
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
}


// ************************************************************************* //