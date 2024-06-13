#include "controller.h"
#include <iostream>
using namespace std;


controller::controller()
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        // Pos PID
        Kp_pos[i] = 0.0;
        Kd_pos[i] = 0.0;

        error_pos[i] = 0.0;
        error_old_pos[i] = error_pos[i];
        error_dot_pos[i] = 0.0;
        error_dot_old_pos[i] = error_dot_pos[i];
        PID_output_pos[i] = 0;

        // Vel PID
        Kp_vel[i] = 0.0;
        Kd_vel[i] = 0.0;

        error_vel[i] = 0.0;
        error_old_vel[i] = error_vel[i];
        error_dot_vel[i] = 0.0;
        error_dot_old_vel[i] = error_dot_vel[i];
        PID_output_vel[i] = 0;

        //Admittance
        deltaPos[i] = 0.0;
        deltaPos_old[i] = deltaPos[i];
        deltaPos_old2[i] = deltaPos_old[i];

        //RWDOB
        rhs_dob[i] = 0.0;
        rhs_dob_old[i] = rhs_dob[i];
        lhs_dob[i] = 0.0;
        lhs_dob_old[i] = lhs_dob[i];
        lhs_dob_LPF[i] = 0.0;
        lhs_dob_LPF_old[i] = lhs_dob_LPF[i];
        rhs_dob_LPF[i] = 0.0;
        rhs_dob_LPF_old[i] = rhs_dob_LPF[i];
        tauDist_hat[i] = 0.0;

        //RWFOB
        rhs_fob[i] = 0.0;
        rhs_fob_old[i] = rhs_fob[i];
        lhs_fob_LPF[i] = 0.0;
        lhs_fob_LPF_old[i] = lhs_fob_LPF[i];
        rhs_fob_LPF[i] = 0.0;
        rhs_fob_LPF_old[i] = rhs_fob_LPF[i];
        tauExt_hat[i] = 0.0;
        forceExt_hat[i] = 0.0;
        forceExt_hat_old[i] = forceExt_hat[i];
        forceExt_hat_old2[i] = forceExt_hat_old[i];


        

    }
    
};
controller::~controller(){}

void controller::pid_gain_pos(double kp, double kd, double cut_off)
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        Kp_pos[i] = kp;
        Kd_pos[i] = kd;
        cut_off_D_pos = cut_off;

    }

};

void controller::pid_gain_vel(double kp, double kd, double cut_off)
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        Kp_vel[i] = kp;
        Kd_vel[i] = kd;
        cut_off_D_vel = cut_off;

    }

};

void controller::ctrl_update()
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        //PID pos
        error_old_pos[i] = error_pos[i];
        error_dot_old_pos[i] = error_dot_pos[i];
        //PID vel
        error_old_vel[i] = error_vel[i];
        error_dot_old_vel[i] = error_dot_vel[i];

        //admittance
        deltaPos_old2[i] = deltaPos_old[i];
        deltaPos_old[i] = deltaPos[i];

        //DOB
        rhs_dob_old[i] = rhs_dob[i];
        lhs_dob_old[i] = lhs_dob[i];
        rhs_dob_LPF_old[i] = rhs_dob_LPF[i];
        lhs_dob_LPF_old[i] = lhs_dob_LPF[i];

        //FOB
        rhs_fob_LPF_old[i] = rhs_fob_LPF[i];
        lhs_fob_LPF_old[i] = lhs_fob_LPF[i];
        forceExt_hat_old2[i] = forceExt_hat_old[i];
        forceExt_hat_old[i] = forceExt_hat[i];

        
    }

};

Vector2d controller::PID_pos(StateModel_* state_model)
{
    for (int i = 0; i < NDOF_LEG; i++) // Error를 state 모델에 넣을 필요 있는지 생각해봐야함. error는 여기에 있어도 됨. //error들 update 해줘야함
    {
        error_pos[i] = state_model->posRW_ref[i] - state_model->posRW[i];
        error_old_pos[i] = state_model->posRW_ref_old[i] - state_model->posRW_old[i];
        
        error_dot_pos[i] = tustin_derivative(error_pos[i], error_old_pos[i], error_dot_old_pos[i], cut_off_D_pos);
        
        // DOB 끄면 PID만 사용해야하는데 state model에 넣지 않아도 되는지 생각해봐야함.
        PID_output_pos[i] = Kp_pos[i] * error_pos[i] + Kd_pos[i] * error_dot_pos[i]; // 이걸 return을 사용하면?
    }
    return PID_output_pos;

};

void controller::PID_vel(StateModel_* state_model)
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        error_vel[i] = state_model->velRW_ref[i] - state_model->velRW[i]; //error is no vector
        error_old_vel[i] = state_model->velRW_ref_old[i] - state_model->velRW_old[i];
        
        error_dot_vel[i] = tustin_derivative(error_vel[i], error_old_vel[i], error_dot_old_vel[i], cut_off_D_vel);
        
        // DOB 끄면 PID만 사용해야하는데 state model에 넣지 않아도 되는지 생각해봐야함.
        PID_output_vel[i] = Kp_vel[i] * error_vel[i] + Kd_vel[i] * error_dot_vel[i]; // 이걸 return을 사용하면?
    }
}; // negative velocity PID feedback

void controller::admittanceCtrl(StateModel_* state_model, double m , double b, double k, int flag)
{
    // 현재 omega_n, zeta, k 로 tunning 하고 있는데, 변환식을 통해 아래에 적어주면 된다
    ad_M = m;
    ad_B = b;
    ad_K = k;

    double c1 = 4 * ad_M + 2 * ad_B * Ts + ad_K * pow(Ts, 2);
    double c2 = -8 * ad_M + 2 * ad_K * pow(Ts, 2);
    double c3 = 4 * ad_M - 2 * ad_B * Ts + ad_K * pow(Ts, 2);
   
    //Delta pos는? force_hat? -> state_model에 있어야하나?
    deltaPos[0] =
        (pow(Ts, 2) * forceExt_hat[0] + 2 * pow(Ts, 2) * forceExt_hat_old[0] +
            pow(Ts, 2) * forceExt_hat_old2[0] - c2 * deltaPos_old[0] - c3 * deltaPos_old2[0]) / c1;

        if (flag == true)
            state_model->posRW_ref[0] = state_model->posRW_ref[0] + deltaPos[0];
};

void controller::DOBRW(StateModel_* state_model, double cut_off ,int flag)
{
    cut_off_dob = cut_off;
    
    //printf(" tau = %f ,%f \n", state_model->qddot_bi_tustin[0],state_model->qddot_bi_tustin[1]);
    lhs_dob = state_model->tau_bi;
    rhs_dob = state_model->Lamda_nominal_DOB * state_model->qddot_bi_tustin;
    // printf(" tau_bi = %f ,%f \n", state_model->tau_bi[0],state_model->tau_bi[1]);
    // printf(" qddot_bi_tustin = %f ,%f \n", state_model->qddot_bi_tustin[0],state_model->qddot_bi_tustin[1]);
    // printf(" Lamda_nominal_DOB(0) = %f ,%f \n", state_model->Lamda_nominal_DOB(0,0),state_model->Lamda_nominal_DOB(0,1));
    // printf(" Lamda_nominal_DOB(1) = %f ,%f \n", state_model->Lamda_nominal_DOB(1,0),state_model->Lamda_nominal_DOB(1,1));
    
    if (state_model->time < 0.00001)
    {
        // printf(" tau_bi = %f ,%f \n", state_model->tau_bi[0],state_model->tau_bi[1]);
        // printf(" qddot_bi_tustin = %f ,%f ", state_model->qddot_bi_tustin[0],state_model->qddot_bi_tustin[1]);
        // printf(" qddot_bi = %f ,%f \n", state_model->qddot_bi[0],state_model->qddot_bi[1]);
        //printf(" Lamda_nominal_DOB(0) = %f ,%f ", state_model->Lamda_nominal_DOB(0,0),state_model->Lamda_nominal_DOB(0,1));
        //printf(" Lamda_nominal_DOB(1) = %f ,%f \n", state_model->Lamda_nominal_DOB(1,0),state_model->Lamda_nominal_DOB(1,1));
        //printf("rhs_dob: %f, rhs_dob_old: %f, rhs_dob_LPF: %f, rhs_dob_LPF_old: %f, cut_off_dob: %f\n", rhs_dob[0], rhs_dob_old[0], rhs_dob_LPF[0],rhs_dob_LPF_old[0], cut_off_dob);
        //printf("lhs_dob: %f, lhs_dob_old: %f, lhs_dob_LPF: %f, lhs_dob_LPF_old: %f, cut_off_dob: %f\n", lhs_dob[0], lhs_dob_old[0], lhs_dob_LPF[0],lhs_dob_LPF_old[0], cut_off_dob);

    }

    if (flag == true)
    {
        for (int i = 0; i < NDOF_LEG; i++)
        {
            lhs_dob_LPF[i] = lowpassfilter(lhs_dob[i], lhs_dob_old[i], lhs_dob_LPF_old[i], cut_off_dob); 
            rhs_dob_LPF[i] = lowpassfilter(rhs_dob[i], rhs_dob_old[i], rhs_dob_LPF_old[i], cut_off_dob); 
            // lhs_dob_LPF[i] = Highpassfilter(lhs_dob[i], lhs_dob_old[i], , cut_off_dob); 
            // //rhs_dob_LPF[i] = Highpassfilter(rhs_dob[i], rhs_dob_old[i], 0, cut_off_dob);
            // rhs_dob_LPF[i] = (Ts * (rhs_dob[i] + rhs_dob_old[i])- (Ts - 2 * time_const)*rhs_dob_LPF_old[i] ) / (Ts + 2 * time_const); //
            tauDist_hat[i] = lhs_dob_LPF[i] - rhs_dob_LPF[i];
        }
        //printf("rhs_dob: %f, rhs_dob_old: %f, rhs_dob_LPF_old: %f, time_const: %f\n", rhs_dob[0], rhs_dob_old[0], rhs_dob_LPF_old[0], time_const);
        //printf(" rhs_dob = %f ,%f \n", rhs_dob[0],rhs_dob[1]);
        //printf(" rhs_dob_old = %f ,%f \n", rhs_dob_old[0],rhs_dob_old[1]);
        //printf(" rhs_dob_LPF = %f ,%f \n", rhs_dob_LPF[0],rhs_dob_LPF[1]);
        //printf(" rhs_dob_LPF_old = %f ,%f \n", rhs_dob_LPF_old[0],rhs_dob_LPF_old[1]);
        // printf(" cut_off_dob = %f \n", cut_off_dob);


    }
    else
    {
        for (int i = 0; i < NDOF_LEG; i++)
            tauDist_hat[i] = 0;
    }
    printf(" tauDist_hat = %f ,%f \n", tauDist_hat[0],tauDist_hat[1]);
    state_model->tau_bi = state_model->tau_bi + tauDist_hat;

}; // Rotating Workspace DOB

void controller::FOBRW(StateModel_* state_model, double cut_off)
{
    cut_off_fob = 2*pi*cut_off;
    
    // Corioli & Gravity term 만들어 놓음 필요하면 쓰면 됩니닷

    rhs_fob = state_model->Lamda_nominal_FOB * state_model->qddot_bi_tustin_old;
    rhs_fob_old = state_model->Lamda_nominal_FOB * state_model->qddot_bi_tustin;

    for (int i = 0; i < NDOF_LEG; i++)
    {
        lhs_fob_LPF[i] = lowpassfilter(state_model->tau_bi[i], state_model->tau_bi_old[i], lhs_fob_LPF_old[i], cut_off_fob);
        rhs_fob_LPF[i] = lowpassfilter(rhs_fob[i], rhs_fob_old[i], rhs_fob_LPF_old[i], cut_off_fob);

        tauExt_hat[i] = rhs_fob_LPF[i] - lhs_fob_LPF[i];

    }
    forceExt_hat = state_model->jacbRW_trans_inv * tauExt_hat;
    
    
         
} // Rotating WorkspaceForce Observer



