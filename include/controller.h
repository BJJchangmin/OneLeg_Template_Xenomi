#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__

#include "globVariable.h"
#include "kinematics.h"
#include "trajectory.h"
#include "filter.h"


int loop_con = 0;
void controller_init(const mjModel* m, mjData* d, ParamModel* param_model, ParamTuning* param_tuning)
{
    // RW PD Controller
    param_tuning->freq_cut_D = 150;

    param_tuning->Kp_pos[0] = 700;
    param_tuning->Kp_pos[1] = 700;

    param_tuning->Kd_pos[0] = 80;
    param_tuning->Kd_pos[1] = 80;

    // RWDOB & RWFOB cutoff frequency
    param_tuning->freq_cut_Qd = 50;
    param_tuning->freq_cut_Qf = 50;

    param_tuning->delta_default = 0.0001;
    param_tuning->zeta = 1;

    param_tuning->Ma = 0.25 * param_model->m_trunk; //m_total�� ��� think. mass ���� 
    param_tuning->Ka = param_tuning->Ma * g / param_tuning->delta_default;
    param_tuning->Ba = 2 * param_tuning->zeta * sqrt(param_tuning->Ma * param_tuning->Ka);
    //printf("%f \n", param_tuning->Ka);

    // d->ctrl[0] = m->key_ctrl[0];
    // d->ctrl[1] = m->key_ctrl[1];
}

void admittanceCtrl(ParamModel* param_model, ParamTuning* param_tuning, StateModel* state_model, int flag)
{
    
    double c1 = 4 * param_tuning->Ma_new + 2 * param_tuning->Ba_new * Ts + param_tuning->Ka_new * pow(Ts, 2);
    double c2 = -8 * param_tuning->Ma_new + 2 * param_tuning->Ka_new * pow(Ts, 2);
    double c3 = 4 * param_tuning->Ma_new - 2 * param_tuning->Ba_new * Ts + param_tuning->Ka_new * pow(Ts, 2);
   
    
    state_model->deltaPos[0] =
        (pow(Ts, 2) * state_model->forceExt_hat[0] + 2 * pow(Ts, 2) * state_model->forceExt_hat_old[0] +
            pow(Ts, 2) * state_model->forceExt_hat_old2[0] - c2 * state_model->deltaPos_old[0] - c3 * state_model->deltaPos_old2[0]) / c1;
    //printf("delta : %f, fxR : %f, fxR_old : %f, fxR_old2 : %f \n", state_model->deltaPos[0], state_model->forceExt_hat[0], state_model->forceExt_hat_old[0], state_model->forceExt_hat_old2[0]);

        if (flag == true)
            state_model->posRW_ref[0] = state_model->posRW_ref[0] + state_model->deltaPos[0];
    

    
} // Admittance control

void posFeedbackPD(ParamTuning* param_tuning, StateModel* state_model,double t)
{
    
    for (int i = 0; i < NDOF_LEG; i++)
    {
        state_model->error_pos[i] = state_model->posRW_ref[i] - state_model->posRW[i];
        state_model->error_pos_old[i] = state_model->posRW_ref_old[i] - state_model->posRW_old[i];
        
        state_model->error_dot_pos[i] = tustin_derivative(state_model->error_pos[i], state_model->error_pos_old[i], state_model->error_dot_pos_old[i], param_tuning->freq_cut_D);
        
        state_model->ctrl_input_RW[i] = param_tuning->Kp_pos[i] * state_model->error_pos[i] + param_tuning->Kd_pos[i] * state_model->error_dot_pos[i];
    }
    
} // negative position PID feedback

void velFeedbackPD(ParamTuning* param_tuning, StateModel* state_model)
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        state_model->error_vel[i] = state_model->velRW_ref[i] - state_model->velRW[i];
        state_model->error_vel_old[i] = state_model->velRW_ref_old[i] - state_model->velRW_old[i];
        state_model->error_dot_vel_old[i] = state_model->error_dot_vel[i];

        state_model->error_dot_vel[i] = tustin_derivative(state_model->error_vel[i], state_model->error_vel_old[i], state_model->error_dot_vel_old[i], param_tuning->freq_cut_D);

        state_model->ctrl_input_RW[i] = param_tuning->Kp_vel[i] * state_model->error_vel[i] + param_tuning->Kd_vel[i] * state_model->error_dot_vel[i];
    }
} // negative velocity PID feedback

void distObserverRW(ParamModel* param_model, ParamTuning* param_tuning, StateModel* state_model, int flag)
{
    double matA[NDOF_LEG * NDOF_LEG] = { 0 };
    matA[0] = 1;
    matA[1] = -1;
    matA[2] = -1;
    matA[3] = 1;

    double matB[NDOF_LEG * NDOF_LEG] = { 0 };
    matB[0] = 1;
    matB[1] = 1;
    matB[2] = 1;
    matB[3] = 1;

    double matC[NDOF_LEG * NDOF_LEG] = { 0 };
    mju_add(matC, matA, matB, 4);
    //printf("%f %f\n", matC[0], matC[1]);
    //printf("%f %f\n\n", matC[2], matC[3]);

    double matD[NDOF_LEG * NDOF_LEG] = { 0 };
    mju_scl(matD, matA, 0.5 * param_model->JzzR_thigh, 4);

    double matE[NDOF_LEG * NDOF_LEG] = { 0 };
    mju_scl(matE, matB, 0.5 * param_model->JzzR_shank, 4);

    double matF[NDOF_LEG * NDOF_LEG] = { 0 };
    mju_add(matF, matD, matE, 4);

    double lhs_dob[NDOF_LEG] = { 0 }, lhs_dob_old[NDOF_LEG] = { 0 };
    double rhs_dob[NDOF_LEG] = { 0 }, rhs_dob_old[NDOF_LEG] = { 0 };
    //printf("%f %f \n\n", state_model->tau_bi[0], state_model->tau_bi_old[0]);

    mju_mulMatVec(lhs_dob, matC, state_model->tau_bi, 2, 2);
    mju_mulMatVec(lhs_dob_old, matC, state_model->tau_bi_old, 2, 2);

    mju_mulMatVec(rhs_dob, matF, state_model->qddot_bi_tustin, 2, 2);
    mju_mulMatVec(rhs_dob_old, matF, state_model->qddot_bi_tustin_old, 2, 2);

    if (flag == 1)
    {
        for (int i = 0; i < NDOF_LEG; i++)
        {
            state_model->lhs_dob_LPF[i] = lowpassfilter(lhs_dob[i], lhs_dob_old[i], state_model->lhs_dob_LPF_old[i], param_tuning->freq_cut_Qd);
            state_model->rhs_dob_LPF[i] = lowpassfilter(rhs_dob[i], rhs_dob_old[i], state_model->rhs_dob_LPF_old[i], param_tuning->freq_cut_Qd);
            state_model->tauDist_hat[i] = -state_model->rhs_dob_LPF[i] + state_model->lhs_dob_LPF[i];
        }
        mju_scl(state_model->tauDist_hat, state_model->tauDist_hat, 0.5, NDOF_LEG);
        //printf("%f %f \n", tauDist_hat[0], tauDist_hat[1]);
        mju_add(state_model->tau_bi, state_model->tau_bi, state_model->tauDist_hat, NDOF_LEG);
    }
    else
    {
        for (int i = 0; i < NDOF_LEG; i++)
        {
            state_model->tauDist_hat[i] = 0;
        }
    }
} // Rotating Workspace DOB

void forceObserverRW(mjData* d, ParamModel* param_model, ParamTuning* param_tuning, StateModel* state_model)
{
    // Coriolis & Gravity 
    double h[NDOF_LEG] = { 0 }, h_old[NDOF_LEG] = { 0 };

    h[0] = -param_model->m_shank * param_model->d_shank * param_model->L * sin(state_model->q[1]) * pow(state_model->qdot_bi[1], 2)
        - g * (param_model->m_thigh * param_model->d_thigh + param_model->m_shank * param_model->L) * cos(state_model->q_bi[0]);
    h_old[0] = -param_model->m_shank * param_model->d_shank * param_model->L * sin(state_model->q_old[1]) * pow(state_model->qdot_bi_old[1], 2)
        - g * (param_model->m_thigh * param_model->d_thigh + param_model->m_shank * param_model->L) * cos(state_model->q_bi_old[0]);

    h[1] = param_model->m_shank * param_model->d_shank * param_model->L * sin(state_model->q[1]) * pow(state_model->qdot_bi[0], 2)
        - g * param_model->m_shank * param_model->d_shank * cos(state_model->q_bi[1]);
    h_old[1] = param_model->m_shank * param_model->d_shank * param_model->L * sin(state_model->q_old[1]) * pow(state_model->qdot_bi_old[0], 2)
        - g * param_model->m_shank * param_model->d_shank * cos(state_model->q_bi_old[1]);
    //printf("%f %f \n", param_model->MatInertia_bi[0], param_model->MatInertia_bi[1]);
    //printf("%f %f \n\n", param_model->MatInertia_bi[2], param_model->MatInertia_bi[3]);

    //h[0] = h[0] - h[1];

    double rhs_fob[NDOF_LEG] = { 0 }, rhs_fob_old[NDOF_LEG] = { 0 };

    mju_mulMatVec(rhs_fob_old, param_model->MatInertia_bi, state_model->qddot_bi_tustin_old, 2, 2);
    mju_mulMatVec(rhs_fob, param_model->MatInertia_bi, state_model->qddot_bi_tustin, 2, 2);

    //mju_sub(rhs_fob_old, rhs_fob_old, h_old, NDOF_LEG);
    //mju_sub(rhs_fob, rhs_fob, h, NDOF_LEG);

    //mju_add(rhs_fob_old, rhs_fob_old, h_old, NDOF_LEG);
    //mju_add(rhs_fob, rhs_fob, h, NDOF_LEG);

    for (int i = 0; i < NDOF_LEG; i++)
    {
        state_model->lhs_fob_LPF[i] = lowpassfilter(state_model->tau_bi[i], state_model->tau_bi_old[i], state_model->lhs_fob_LPF_old[i], param_tuning->freq_cut_Qf);
        state_model->rhs_fob_LPF[i] = lowpassfilter(rhs_fob[i], rhs_fob_old[i], state_model->rhs_fob_LPF_old[i], param_tuning->freq_cut_Qf);

        state_model->tauExt_hat[i] = state_model->rhs_fob_LPF[i] - state_model->lhs_fob_LPF[i];
        //lowpassfilter(est_torque_ext[i], est_torque_ext_old[i], &torque_LPF[i], &torque_LPF_old[i], cutoff_freq);
    }
    //ó�� 8���� Ȯ���ϱ� ���ؼ� �ִ� �ڵ�
    while (loop_con< 8) {
        //printf("%f  %f\n", state_model->tauExt_hat[0], state_model->tauExt_hat[1]);
        //printf(" %f \n", state_model->q[1]); //q[1]�� 0���� ���°� ������ ��� 0�� ������ Ȯ���غ�����
        
        break;
    }
    
    mju_mulMatVec(state_model->forceExt_hat, state_model->jacbRW_trans_inv, state_model->tauExt_hat, 2, 2);
    
    while (loop_con < 80000) {
        state_model->forceExt_hat[0] = 0;
        state_model->forceExt_hat[1] = 0;
        break;
    }
    
      

    loop_con += 1;
    //printf("fxR_hat : %f, fyR_hat : %f \n", lhs_fob_LPF[0], lhs_fob_LPF[1]);
    //printf("fxR_hat : %f, fyR_hat : %f \n", rhs_fob[0], rhs_fob[1]); 
    //printf("fxR_hat : %f, fyR_hat : %f \n\n", rhs_fob_LPF[0], rhs_fob_LPF[1]);
    //printf("fxR_hat : %f, fyR_hat : %f \n", tauExt_hat[0], tauExt_hat[1]);
    //printf("fxR_hat : %f, fyR_hat : %f \n", forceExt_hat[0], forceExt_hat[1]);
} // Rotating WorkspaceForce Observer

#endif // !__CONTROLLER_H__
