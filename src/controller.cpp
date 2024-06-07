#include "controller.h"



controller::controller()
{
    // sate_model에 있다가 controller class로 이동한 얘들은 전부 초기화를 여기서 해줘야하고 함수 내에서 update 해줘야함. PID에서는 없을 거 같고 DOB,FOB, Admittance에 있을 거 같음.
    // RW PD Controller
    // param_tuning->freq_cut_D = 150;

    // param_tuning->Kp_pos[0] = 700;
    // param_tuning->Kp_pos[1] = 700;

    // param_tuning->Kd_pos[0] = 80;
    // param_tuning->Kd_pos[1] = 80;

    // // RWDOB & RWFOB cutoff frequency
    // param_tuning->freq_cut_Qd = 50;
    // param_tuning->freq_cut_Qf = 50;

    // param_tuning->delta_default = 0.0001;
    // param_tuning->zeta = 1;

    // param_tuning->Ma = 0.25 * param_model->m_trunk; //m_total�� ��� think. mass ���� 
    // param_tuning->Ka = param_tuning->Ma * g / param_tuning->delta_default;
    // param_tuning->Ba = 2 * param_tuning->zeta * sqrt(param_tuning->Ma * param_tuning->Ka);

    // double lhs_dob[NDOF_LEG] = { 0 }, lhs_dob_old[NDOF_LEG] = { 0 };
    // double rhs_dob[NDOF_LEG] = { 0 }, rhs_dob_old[NDOF_LEG] = { 0 };

    //double rhs_fob[NDOF_LEG] = { 0 }, rhs_fob_old[NDOF_LEG] = { 0 };

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

void controller::ctrlupdate()
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        //PID pos
        error_dot_old_pos[i] = error_dot_pos[i];
        //PID vel
        error_dot_old_vel[i] = error_dot_vel[i];

        //admittance
        deltaPos_old2[i] = deltaPos_old[i];
        deltaPos_old[i] = deltaPos[i];

        //DOB
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
    cut_off_dob = 1/(2*pi*cut_off);
    
    lhs_dob = state_model->tau_bi;
    lhs_dob_old = state_model->tau_bi_old;

    rhs_dob = state_model->Lamda_nominal_DOB * state_model->qddot_bi_tustin;
    rhs_dob = state_model->Lamda_nominal_DOB * state_model->qddot_bi_tustin_old;
    
    if (flag == true)
    {
        for (int i = 0; i < NDOF_LEG; i++)
        {
            lhs_dob_LPF[i] = lowpassfilter(lhs_dob[i], lhs_dob_old[i], lhs_dob_LPF_old[i], cut_off_dob);
            rhs_dob_LPF[i] = lowpassfilter(rhs_dob[i], rhs_dob_old[i], rhs_dob_LPF_old[i], cut_off_dob);

            tauDist_hat[i] = -rhs_dob_LPF[i] + lhs_dob_LPF[i];
        }

    }
    else
    {
        for (int i = 0; i < NDOF_LEG; i++)
            tauDist_hat[i] = 0;
    }
    state_model->tau_bi = state_model->tau_bi + tauDist_hat;

}; // Rotating Workspace DOB

void controller::FOBRW(StateModel_* state_model, double cut_off)
{
    cut_off_fob = 1/(2*pi*cut_off);
    
    rhs_fob = state_model->Lamda_nominal_FOB * state_model->qddot_bi_tustin_old;
    rhs_fob_old = state_model->Lamda_nominal_FOB * state_model->qddot_bi_tustin;

    for (int i = 0; i < NDOF_LEG; i++)
    {
        lhs_fob_LPF[i] = lowpassfilter(state_model->tau_bi[i], state_model->tau_bi_old[i], lhs_fob_LPF_old[i], cut_off_fob);
        rhs_fob_LPF[i] = lowpassfilter(rhs_fob[i], rhs_fob_old[i], rhs_fob_LPF_old[i], cut_off_fob);

        tauExt_hat[i] = rhs_fob_LPF[i] - lhs_fob_LPF[i];
        //lowpassfilter(est_torque_ext[i], est_torque_ext_old[i], &torque_LPF[i], &torque_LPF_old[i], cutoff_freq);
    }
    forceExt_hat = state_model->jacbRW_trans_inv * tauExt_hat;
    
    
         
} // Rotating WorkspaceForce Observer


