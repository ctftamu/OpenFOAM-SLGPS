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
    const scalar p,
    const scalarField &data,
    const labelField &always_inds,
    const labelField &never_inds,
    const labelField &vari_inds
)
{   
}

template<class CompType, class ThermoType>
void Foam::chemistryReductionMethods::GPS_ML<CompType, ThermoType>::reduceMechanism
(
    const scalarField &c,
    const scalar T,
    const scalar p
)
{   
    Info<<"Temperature: ,"<<T<<endl;
    Info<<"Pressure: ,"<<p<<endl;
    
    

    //const int n_always_inds = 13;
    //int always_inds[n_always_inds] = {0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 47};
    //for train_data_TDAC_a_0.001_less_nev
    //const int n_always_inds = 18;
    //int always_inds[n_always_inds] = {0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 47};
    //for train_data_LOW_TEMP_a_0.001
    //const int n_always_inds = 10;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 12, 13, 15, 47};
    //for train_data_LOW_TEMP_a_0.005
    //const int n_always_inds = 9;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 12, 13, 15, 47};
    //for train_data_LOW_TEMP_a_0.0001
    //const int n_always_inds = 9;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 13, 15, 47};
    //for train_data_LOW_TEMP_a_0.000001
    //const int n_always_inds = 9;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 13, 15, 47};
    //for train_data_LOW_TEMP_a_0.000000001
    //const int n_always_inds = 10;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 13, 15, 34, 47};
    //for train_data_LOW_TEMP_a_0.000000000001
    //const int n_always_inds = 12;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 7, 13, 15, 17, 34, 47};
    //for train_data_EXTRA_PROD_a_0.000000000001
    //const int n_always_inds = 12;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 7, 13, 15, 17, 34, 47};
    //for train_data_EXTRA_PROD_a_0.0000000000001
    //const int n_always_inds = 13;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 7, 13, 15, 16, 17, 34, 47};
    //for train_data_EXTRA_PROD_a_0.0000000000005
    //const int n_always_inds = 12;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 7, 13, 15, 17, 34, 47};
    //for train_data_EXTRA_PROD_a_0.000000000001_v2
    //const int n_always_inds = 12;
    //int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 7, 13, 15, 17, 34, 47};
    //for train_data_EXTRA_PROD_a_0.00000000000075
    const int n_always_inds = 12;
    int always_inds[n_always_inds] = {1, 2, 3, 4, 5, 6, 7, 13, 15, 17, 34, 47};
    
    //const int n_never_inds = 5;
    //int never_inds[n_never_inds] = {33, 41, 42, 43, 48};
    //const int n_never_inds = 8;
    //int never_inds[n_never_inds] = {33, 39, 41, 42, 43, 44, 45, 48};
    //for train_data_TDAC_a_0.001_less_nev
    //const int n_never_inds = 7;
    //int never_inds[n_never_inds] = {32, 33, 38, 41, 43, 44, 48};
    //for train_data_LOW_TEMP_a_0.001
    //const int n_never_inds = 10;
    //int never_inds[n_never_inds] = {32, 33, 37, 38, 41, 42, 43, 44, 45, 48};
    //for train_data_LOW_TEMP_a_0.005
    //const int n_never_inds = 10;
    //int never_inds[n_never_inds] = {32, 33, 37, 38, 41, 42, 43, 44, 45, 48};
    //for train_data_LOW_TEMP_a_0.0001
    //const int n_never_inds = 3;
    //int never_inds[n_never_inds] = {33, 42, 48};
    //for train_data_LOW_TEMP_a_0.000001
    //const int n_never_inds = 1;
    //int never_inds[n_never_inds] = {48};
    //for train_data_LOW_TEMP_a_0.000000001
    //const int n_never_inds = 1;
    //int never_inds[n_never_inds] = {48};
    //for train_data_LOW_TEMP_a_0.000000000001
    //const int n_never_inds = 1;
    //int never_inds[n_never_inds] = {48};
    //for train_data_EXTRA_PROD_a_0.000000000001
    //const int n_never_inds = 1;
    //int never_inds[n_never_inds] = {48};
    //for train_data_EXTRA_PROD_a_0.0000000000001
    //const int n_never_inds = 1;
    //int never_inds[n_never_inds] = {48};
    //for train_data_EXTRA_PROD_a_0.0000000000005
    //const int n_never_inds = 1;
    //int never_inds[n_never_inds] = {48};
    //for train_data_EXTRA_PROD_a_0.000000000001_v2
    //const int n_never_inds = 1;
    //int never_inds[n_never_inds] = {48};
    //for train_data_EXTRA_PROD_a_0.00000000000075
    const int n_never_inds = 1;
    int never_inds[n_never_inds] = {48};
    
    
    //const int n_vari_inds = 35;
    //int vari_inds[n_vari_inds] = {8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 49, 50, 51, 52};
    //const int n_vari_inds = 32;
    //int vari_inds[n_vari_inds] = {8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 40, 46, 49, 50, 51, 52};
    //for train_data_TDAC_a_0.001_less_nev
    //const int n_vari_inds = 28;
    //int vari_inds[n_vari_inds] = {7, 8, 9, 10, 19, 20, 21, 22, 23, 24, 27, 28, 29, 30, 31, 34, 35, 36, 37, 39, 40, 42, 45, 46, 49, 50, 51, 52};
    //for train_data_LOW_TEMP_a_0.001
    //const int n_vari_inds = 33;
    //int vari_inds[n_vari_inds] = {0, 7, 8, 9, 10, 11, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 34, 35, 36, 39, 40, 46, 49, 50, 51, 52};
    //for train_data_LOW_TEMP_a_0.005
    //const int n_vari_inds = 34;
    //int vari_inds[n_vari_inds] = {0, 6, 7, 8, 9, 10, 11, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 34, 35, 36, 39, 40, 46, 49, 50, 51, 52};
    //for train_data_LOW_TEMP_a_0.0001
    //const int n_vari_inds = 41;
    //int vari_inds[n_vari_inds] = {0, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 49, 50, 51, 52};
    //for train_data_LOW_TEMP_a_0.000001
    //const int n_vari_inds = 43;
    //int vari_inds[n_vari_inds] = {0, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52};
    //for train_data_LOW_TEMP_a_0.000000001
    //const int n_vari_inds = 42;
    //int vari_inds[n_vari_inds] = {0, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52};
    //for train_data_LOW_TEMP_a_0.00000000001
    //const int n_vari_inds = 40;
    //int vari_inds[n_vari_inds] = {0, 8, 9, 10, 11, 12, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52};
    //for train_data_EXTRA_PROD_a_0.000000000001
    //const int n_vari_inds = 40;
    //int vari_inds[n_vari_inds] = {0, 8, 9, 10, 11, 12, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52};
    //for train_data_EXTRA_PROD_a_0.0000000000001
    //const int n_vari_inds = 39;
    //int vari_inds[n_vari_inds] = {0, 8, 9, 10, 11, 12, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52};
    //for train_data_EXTRA_PROD_a_0.0000000000005
    //const int n_vari_inds = 40;
    //int vari_inds[n_vari_inds] = {0, 8, 9, 10, 11, 12, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52};
    //for train_data_EXTRA_PROD_a_0.000000000001_v2
    //const int n_vari_inds = 40;
    //int vari_inds[n_vari_inds] = {0, 8, 9, 10, 11, 12, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52};
    //for train_data_EXTRA_PROD_a_0.00000000000075
    const int n_vari_inds = 40;
    int vari_inds[n_vari_inds] = {0, 8, 9, 10, 11, 12, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52};
    
    const int n_input_inds = 12;
    int input_inds[n_input_inds-2] = {13, 5, 4, 1, 14, 3, 15, 2, 12, 9};
    
    //for sup_model.pb
    //float input_maxes[n_input_inds] = {3294.3137, 100.0, 0.1282, 0.1806, 0.02904, 0.03321, 0.09929, 0.1976, 0.06597, 0.02169, 0.01106, 8.2008e-05};
    //float input_mins[n_input_inds] = {1299.9992, 1.0, 6.3475e-19, 1.7356e-13, 1.5037e-12, 7.2028e-12, 1.9753e-19, 0.0005175, 2.4630e-21, 1.4075e-12, 5.3383e-18, 7.4433e-23};
    
    //for model_c_500_a_0.001_l_16_dt_200.pb
    //float input_maxes[n_input_inds] = {3219.9131346732624, 98.88285621310023, 0.12812228934173156, 0.18619097575556737, 0.02879275852922004, 0.035414855595551846, 0.09861926058935107, 0.1976103193196484, 0.07393701039839022, 0.02310716016769376, 0.012957167247582675, 0.00011822748907861104};
    //float input_mins[n_input_inds] = {1300.0604102031048, 0.10013999234648284, 2.679008231681836e-19, 6.495569289946552e-33, 2.829678300694169e-25, 1.1796952977710393e-17, 5.945326919960085e-41, 0.00011732549819055488, 1.718204883434161e-41, 2.145338540125709e-23, 6.315723958055864e-18, 2.8714414299479102e-34};
    
    //for train_data_TDAC_a_0.001_less_nev
    //float input_maxes[n_input_inds] = {2046.654767, 109342.576, 0.24961178507773144, 0.18258268881610723, 0.0063433596064096, 0.00688823373426578, 0.05625323171802815, 0.2072921783034408, 0.07792360967677478, 0.0030995467818457126, 0.003786400395665176, 3.6955261602451327e-06};
    //float input_mins[n_input_inds] = {291.0049001, 89405.12668, 0.0, 8.841391508717879e-58, 7.456887893341725e-61, 0.0, 2.621564557981277e-59, 0.000916465197709732, 3.9586762769402496e-58, 0.0, 0.0, 0.0};
    
    //for train_data_LOW_TEMP_a_0.001
    //float input_maxes[n_input_inds] = {2050.0, 1.0, 0.12820512820512625, 0.1549800889370777, 0.022980452669296724, 0.02953415018868222, 0.09919159769982865, 0.19762845849802277, 0.04696833435880292, 0.01739863744998275, 0.008647235313766525, 6.977407128991841e-05};
    //float input_mins[n_input_inds] = {290.0, 1.0, 4.080078113982859e-18, 2.8185517989202128e-31, 7.020553381952472e-24, 3.3907900710632095e-16, 6.957341601925623e-39, 0.00537787127960494, 3.0572978425022195e-39, 9.425046590255887e-22, 5.116936109578897e-17, 7.125253238441416e-32};
    
    //for train_data_LOW_TEMP_a_0.005
    //float input_maxes[n_input_inds] = {2050.0, 1.0, 0.2396166134185276, 0.15466071220546077, 0.02117195253762832, 0.028976087087382085, 0.11571049258720222, 0.19762845849802277, 0.04696833435880292, 0.01739863744998275, 0.010764833433680747, 6.873226144522943e-05};
    //float input_mins[n_input_inds] = {290.0, 1.0, 4.080078113982859e-18, 2.8185517989202128e-31, 7.020553381952472e-24, 3.3907900710632095e-16, 6.957341601925623e-39, 0.00537787127960494, 3.0572978425022195e-39, 9.425046590255887e-22, 5.116936109578897e-17, 7.125253238441416e-32};
    
    //for train_data_LOW_TEMP_a_0.0001
    //float input_maxes[n_input_inds] = {2050.0, 1.0, 0.25157232704402227, 0.1532121710249013, 0.020616671507983172, 0.02995083692452817, 0.11656211250668806, 0.20790020790020763, 0.047979620436186046, 0.01777963214031066, 0.010632139480001737, 6.506121535281325e-05};
    //float input_mins[n_input_inds] = {290.0, 1.0, 0.0, 4.1275912666522404e-32, 3.4510516978126006e-24, 9.310360588297054e-17, 5.097113825829686e-39, 0.0014552980509695756, 2.263388605568541e-39, 6.024805233091142e-22, 0.0, 0.0};
    
    //for train_data_LOW_TEMP_a_0.000001
    //float input_maxes[n_input_inds] = {2050.0, 1.0, 0.25157232704402227, 0.1532121710249013, 0.020616671507983172, 0.02995083692452817, 0.11656211250668806, 0.20790020790020763, 0.047979620436186046, 0.01777963214031066, 0.010632139480001737, 6.506121535281325e-05};
    //float input_mins[n_input_inds] = {290.0, 1.0, 0.0, 4.1275912666522404e-32, 3.4510516978126006e-24, 9.310360588297054e-17, 5.097113825829686e-39, 0.0014552980509695756, 2.263388605568541e-39, 6.024805233091142e-22, 0.0, 0.0};
    
    //for train_data_LOW_TEMP_a_0.000000001
    //float input_maxes[n_input_inds] = {2050.0, 1.0, 0.25157232704402227, 0.1532121710249013, 0.020616671507983172, 0.02995083692452817, 0.11656211250668806, 0.20790020790020763, 0.047979620436186046, 0.01777963214031066, 0.010632139480001737, 6.506121535281325e-05};
    //float input_mins[n_input_inds] = {290.0, 1.0, 0.0, 4.1275912666522404e-32, 3.4510516978126006e-24, 9.310360588297054e-17, 5.097113825829686e-39, 0.0014552980509695756, 2.263388605568541e-39, 6.024805233091142e-22, 0.0, 0.0};
    
    //for train_data_LOW_TEMP_a_0.000000000001
    //float input_maxes[n_input_inds] = {2050.0, 1.0, 0.25157232704402227, 0.1532121710249013, 0.020616671507983172, 0.02995083692452817, 0.11656211250668806, 0.20790020790020763, 0.047979620436186046, 0.01777963214031066, 0.010632139480001737, 6.506121535281325e-05};
    //float input_mins[n_input_inds] = {290.0, 1.0, 0.0, 4.1275912666522404e-32, 3.4510516978126006e-24, 9.310360588297054e-17, 5.097113825829686e-39, 0.0014552980509695756, 2.263388605568541e-39, 6.024805233091142e-22, 0.0, 0.0};
    
    //for train_data_EXTRA_PROD_a_0.000000000001
    //float input_maxes[n_input_inds] = {2050.0, 1.0, 0.25157232704402227, 0.1629839892138927, 0.020616671507983172, 0.02995083692452817, 0.12280170740217664, 0.20790020790020763, 0.0712743841602297, 0.01777963214031066, 0.010632139480001737, 6.506121535281325e-05};
    //float input_mins[n_input_inds] = {290.0, 1.0, 0.0, 4.1275912666522404e-32, 3.4510516978126006e-24, 9.310360588297054e-17, 5.097113825829686e-39, 0.0010507965128367855, 2.263388605568541e-39, 6.024805233091142e-22, 0.0, 0.0};
    
    //for train_data_EXTRA_PROD_a_0.0000000000001
    //float input_maxes[n_input_inds] = {2050.0, 1.0, 0.25157232704402227, 0.16594424796519067, 0.020616671507983172, 0.02995083692452817, 0.12241510069592972, 0.20790020790020763, 0.07873180240865373, 0.01777963214031066, 0.010632139480001737, 6.506121535281325e-05};
    //float input_mins[n_input_inds] = {290.0, 1.0, 0.0, 4.1275912666522404e-32, 3.4510516978126006e-24, 9.310360588297054e-17, 5.097113825829686e-39, 0.0014552980509695756, 2.263388605568541e-39, 6.024805233091142e-22, 0.0, 0.0};
    
    //for train_data_EXTRA_PROD_a_0.000000000001_v2
    //float input_maxes[n_input_inds] = {2050.0, 1.0, 0.25157232704402227, 0.16594424796519067, 0.020616671507983172, 0.02995083692452817, 0.12241510069592972, 0.20790020790020763, 0.07873180240865373, 0.01777963214031066, 0.010632139480001737, 6.506121535281325e-05};
    //float input_mins[n_input_inds] = {290.0, 1.0, 0.0, 4.1275912666522404e-32, 3.4510516978126006e-24, 9.310360588297054e-17, 5.097113825829686e-39, 0.0014552980509695756, 2.263388605568541e-39, 6.024805233091142e-22, 0.0, 0.0};
    
    //for train_data_EXTRA_PROD_a_0.00000000000075
    float input_maxes[n_input_inds] = {2050.0, 1.0, 0.25157232704402227, 0.16594424796519067, 0.020616671507983172, 0.02995083692452817, 0.12241510069592972, 0.20790020790020763, 0.07873180240865373, 0.01777963214031066, 0.010632139480001737, 6.506121535281325e-05};
    float input_mins[n_input_inds] = {290.0, 1.0, 0.0, 4.1275912666522404e-32, 3.4510516978126006e-24, 9.310360588297054e-17, 5.097113825829686e-39, 0.0014552980509695756, 2.263388605568541e-39, 6.024805233091142e-22, 0.0, 0.0};
    
    float input_vals[1][n_input_inds];
    float output_vals[1][n_vari_inds] = {0};
    scalar csum = 0.0;
    for (label i=0; i<53; i++) {
        csum += c[i];
    }

Info<<"Input Concentrations: ,";
for (label i=0; i < n_input_inds-2; i++) {
    Info<<c[input_inds[i]]<<", ";
}
Info<<endl;

Info<<"Concentrations: ,";
for (label i=0; i<53; i++) {
    Info<<c[i]<<", ";
}
Info<<endl;


Info<<"Input Mole Fractions: ";

scalar mol_frac = 0.0;

for (label i=2; i < n_input_inds; i++) {



       mol_frac = c[input_inds[i-2]]/csum;
       Info<<mol_frac<<", ";
       input_vals[0][i] = (mol_frac - input_mins[i])/(input_maxes[i] - input_mins[i]);
}
    input_vals[0][0] = (T-input_mins[0])/(input_maxes[0] - input_mins[0]);
    
    //input_vals[0][1] = ((p/101325)-input_mins[1])/(input_maxes[1] - input_mins[1]);
    input_vals[0][1] = 0.0;
    
Info<<endl; 


//for (label i=2; i<n_input_inds; i++) {



//       mol_frac = c[input_inds[i-2]]/csum;
       
//    if ((mol_frac - input_mins[i])/(input_maxes[i] - input_mins[i]) < 0) {
//       input_vals[0][i] = 0;
//    }
    
//    else if ((mol_frac - input_mins[i])/(input_maxes[i] - input_mins[i]) > 1) {
//       input_vals[0][i] = 1;
//    }
    
//    else {
//       input_vals[0][i] = (mol_frac - input_mins[i])/(input_maxes[i] - input_mins[i]);
//    }
//}



//    if ((T-input_mins[0])/(input_maxes[0] - input_mins[0]) < 0) {
//       input_vals[0][0] = 0;
//    }
    
//    else if ((T-input_mins[0])/(input_maxes[0] - input_mins[0]) > 1) {
//       input_vals[0][0] = 1;
//    }
    
//    else {
//       input_vals[0][0] = (T-input_mins[0])/(input_maxes[0] - input_mins[0]);
//    }
    
    
    
//    if (((p/101325)-input_mins[1])/(input_maxes[1] - input_mins[1]) < 0) {
//       input_vals[0][1] = 0;
//    }
    
//    else if (((p/101325)-input_mins[1])/(input_maxes[1] - input_mins[1]) > 1) {
//       input_vals[0][1] = 1;
//    }
    
//    else {
//       input_vals[0][1] = ((p/101325)-input_mins[1])/(input_maxes[1] - input_mins[1]);
//    }
    
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
    label speciesNumber = 0;
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
