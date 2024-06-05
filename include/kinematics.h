#ifndef __KINEMATICS_H__
#define __KINEMATICS_H__

#include "globVariable.h"
#include "filter.h"
//#include "state.h"

void state_init(const mjModel* m, mjData* d, StateModel* state_model, ParamModel* param_model)
{
    
    state_model->q[0] = d->qpos[1];
    state_model->q[1] = d->qpos[2];

    state_model->q_bi[0] = d->qpos[1];
    state_model->q_bi[1] = d->qpos[1] + d->qpos[2];

    // RW coordinates initialization
    state_model->r0 = 2 * param_model->L * cos((state_model->q_bi[1] - state_model->q_bi[0]) / 2);

    state_model->posRW[0] = 2 * param_model->L * cos((state_model->q_bi[1] - state_model->q_bi[0]) / 2);
    state_model->posRW[1] = (state_model->q_bi[0] + state_model->q_bi[1]) / 2;

    state_model->posRW_ref[0] = 2 * param_model->L * cos((state_model->q_bi[1] - state_model->q_bi[0]) / 2);
    state_model->posRW_ref[1] = (state_model->q_bi[0] + state_model->q_bi[1]) / 2;

    //Mg Trajectory
    state_model->step = 0;
    state_model->step_old = state_model->step;
    state_model->step_old2 = state_model->step_old;
    
    state_model->out_acc = 0;
    state_model->out_acc_old = state_model->out_acc;
    state_model->out_acc_old2 = state_model->out_acc_old;

    state_model->out_pos = 0;
    state_model->out_pos_old = state_model->out_pos;
    state_model->out_pos_old2 = state_model->out_pos_old;

    state_model->touch_sensor = 0;


    //printf("%f \n", state_model->r0);
    for (int i = 0; i < NDOF_LEG; i++)
    {
        // Joint coordinates [k-1] values
        state_model->q_bi_old[i] = state_model->q_bi[i];

        state_model->qdot_bi[i] = 0.;
        state_model->qdot_bi_tustin[i] = 0.;
        state_model->qdot_bi_tustin_old[i] = state_model->qdot_bi_tustin[i];

        state_model->qddot_bi[i] = 0.;
        state_model->qddot_bi_tustin[i] = 0.;
        state_model->qddot_bi_tustin_old[i] = state_model->qddot_bi_tustin[i];

        state_model->tau_bi[i] = 0.;
        state_model->tau_bi_old[i] = state_model->tau_bi[i];

        // RW coordinates [k-1] values
        state_model->posRW_old[i] = state_model->posRW[i];
        state_model->posRW_ref_old[i] = state_model->posRW_ref[i];
        state_model->posRW_ref_old2[i] = state_model->posRW_ref_old[i];

        state_model->velRW[i] = 0.;
        state_model->velRW_old[i] = state_model->velRW[i];
        state_model->velRW_ref[i] = 0.;
        state_model->velRW_ref_old[i] = state_model->velRW_ref[i];

        // Feedback controller
        state_model->error_pos[i] = 0.;
        state_model->error_pos_old[i] = state_model->error_pos[i];
        state_model->error_dot_pos[i] = 0.;
        state_model->error_dot_pos_old[i] = state_model->error_dot_pos[i];

        state_model->error_vel[i] = 0.;
        state_model->error_vel_old[i] = state_model->error_vel[i];
        state_model->error_dot_vel[i] = 0.;
        state_model->error_dot_vel_old[i] = state_model->error_dot_vel[i];

        state_model->ctrl_input_RW[i] = 0.;
        state_model->ctrl_input_RW_old[i] = state_model->ctrl_input_RW[i];

        // RWDOB
        state_model->lhs_dob_LPF[i] = 0.;
        state_model->lhs_dob_LPF_old[i] = state_model->lhs_dob_LPF[i];
        state_model->rhs_dob_LPF[i] = 0.;
        state_model->rhs_dob_LPF_old[i] = state_model->rhs_dob_LPF[i];

        state_model->tauDist_hat[i] = 0.;

        // RWFOB
        state_model->lhs_fob_LPF[i] = 0.;
        state_model->lhs_fob_LPF_old[i] = state_model->lhs_fob_LPF[i];
        state_model->rhs_fob_LPF[i] = 0.;
        state_model->rhs_fob_LPF_old[i] = state_model->rhs_fob_LPF[i];

        state_model->tauExt_hat[i] = 0.;
        state_model->forceExt_hat[i] = 0.;
        state_model->forceExt_hat_old[i] = state_model->forceExt_hat[i];
        state_model->forceExt_hat_old2[i] = state_model->forceExt_hat_old[i];

        // Admittance
        state_model->deltaPos[i] = 0.;
        state_model->deltaPos_old[i] = state_model->deltaPos[i];
        state_model->deltaPos_old2[i] = state_model->deltaPos_old[i];

        // Mg Trajectory
        state_model->tau_ff[i]=0.;
    }
    // printf("%f, %f \n", state_model->q_bi_old[0], state_model->q_bi_old[1]);
}

void state_update(StateModel* state_model)
{

    //Mg Trajectory
    state_model->step_old2 = state_model->step_old;
    state_model->step_old = state_model->step;

    state_model->out_acc_old2 = state_model->out_acc_old;
    state_model->out_acc_old = state_model->out_acc;

    state_model->out_pos_old2 = state_model->out_pos_old;
    state_model->out_pos_old = state_model->out_pos;

    //printf("%f, %f %f \n\n", state_model->forceExt_hat[0], state_model->forceExt_hat_old[0], state_model->forceExt_hat_old2[0]);
    for (int i = 0; i < NDOF_LEG; i++)
    {
        // Joint
        state_model->q_old[i] = state_model->q[i];
        state_model->q_bi_old[i] = state_model->q_bi[i]; //okay
        state_model->qdot_bi_old[i] = state_model->qdot_bi[i];
        state_model->qdot_bi_tustin_old[i] = state_model->qdot_bi_tustin[i]; //okay
        state_model->qddot_bi_tustin_old[i] = state_model->qddot_bi_tustin[i]; //okay

        // Feedback - RW Kinematics
        state_model->posRW_old[i] = state_model->posRW[i];
        state_model->velRW_old[i] = state_model->velRW[i];
        state_model->posRW_ref_old2[i] = state_model->posRW_ref_old[i];
        state_model->posRW_ref_old[i] = state_model->posRW_ref[i];
        state_model->velRW_ref_old[i] = state_model->velRW_ref[i];

        // control input
        state_model->ctrl_input_RW_old[i] = state_model->ctrl_input_RW[i];
        state_model->tau_bi_old[i] = state_model->tau_bi[i];

        // RWDOB
        state_model->rhs_dob_LPF_old[i] = state_model->rhs_dob_LPF[i];
        state_model->lhs_dob_LPF_old[i] = state_model->lhs_dob_LPF[i];

        // RWFOB
        state_model->lhs_fob_LPF_old[i] = state_model->lhs_fob_LPF[i];
        state_model->rhs_fob_LPF_old[i] = state_model->rhs_fob_LPF[i];

        // Admittance
        state_model->deltaPos_old2[i] = state_model->deltaPos_old[i];
        state_model->deltaPos_old[i] = state_model->deltaPos[i];

        state_model->forceExt_hat_old2[i] = state_model->forceExt_hat_old[i];
        state_model->forceExt_hat_old[i] = state_model->forceExt_hat[i];

        state_model->error_dot_pos_old[i] = state_model->error_dot_pos[i];
        state_model->error_dot_vel_old[i] = state_model->error_dot_vel[i];
        
    }
    //printf("%f, %f \n", state_model->lhs_fob_LPF[0], state_model->lhs_fob_LPF_old[0]);
}

void model_param_cal(const mjModel* m, mjData* d, ParamModel* param_model, StateModel* state_model)
{
    /* Trunk Parameters */
    param_model->m_hip = 2.5;
    param_model->m_trunk_front = 0.;
    param_model->m_trunk_rear = 0.;
    param_model->m_trunk = 4 * param_model->m_hip + param_model->m_trunk_front + param_model->m_trunk_rear;

    /* Leg Parameters */
    param_model->L = 0.25;
    param_model->d_thigh = 0.11017; // local position of CoM of thigh
    param_model->d_shank = 0.12997; // local position of CoM of shank
    // printf("d_thigh : %f, d_shank : %f \n", param_model->d_thigh, param_model->d_shank);

    param_model->m_thigh = 1.017; // mass of thigh link
    param_model->m_shank = 0.143; // mass of shank link
    param_model->m_leg = param_model->m_thigh + param_model->m_shank;
    param_model->m_total = param_model->m_trunk + 4 * param_model->m_leg;
    // printf("m_thigh : %f, m_shank : %f \n", param_model->m_thigh, param_model->m_shank);

    param_model->Izz_thigh = 0.0057;     // MoI of thigh w.r.t. CoM
    param_model->Izz_shank = 8.0318e-04; // MoI of shank w.r.t. CoM
    // printf("Izz_thigh : %f, Izz_shank : %f \n", param_model->Izz_thigh, param_model->Izz_shank);

    param_model->Jzz_thigh =
        param_model->Izz_thigh + param_model->m_thigh * pow(param_model->d_thigh, 2); // MoI of thigh w.r.t. HFE
    param_model->Jzz_shank =
        param_model->Izz_shank + param_model->m_shank * pow(param_model->d_shank, 2); // MoI of thigh w.r.t. KFE
    // printf("Jzz_thigh : %f, Jzz_shank : %f \n", param_model->Jzz_thigh, param_model->Jzz_shank);

    double M1 = param_model->Jzz_thigh + param_model->m_shank * pow(param_model->L, 2);
    double M2 = param_model->m_shank * param_model->d_shank * param_model->L * cos(state_model->q[1]);
    double M12 = param_model->Jzz_shank;

    param_model->MatInertia_bi[0] = M1;
    param_model->MatInertia_bi[1] = M12;
    param_model->MatInertia_bi[2] = M12;
    param_model->MatInertia_bi[3] = M2;
     //printf("%f %f \n", M1, M12);
     //printf("%f \n\n", param_model->MatInertia_bi[3]);

    param_model->JzzR_thigh = param_model->Jzz_thigh + param_model->m_shank * pow(param_model->L, 2) -
        2 * param_model->m_shank * param_model->d_shank * param_model->L * cos(state_model->q[1]);
    param_model->JzzR_couple =
        param_model->Jzz_thigh + param_model->m_shank * pow(param_model->L, 2) - param_model->Jzz_shank;
    param_model->JzzR_shank = param_model->Jzz_thigh + param_model->m_shank * pow(param_model->L, 2) +
        2 * param_model->m_shank * param_model->d_shank * param_model->L * cos(state_model->q[1]);
    // printf("JzzR_thigh : %f, JzzR_shank : %f, JzzR_couple : %f \n", param_model->JzzR_thigh, param_model->JzzR_shank,
    // param_model->JzzR_couple);

    param_model->MatInertia_RW[0] =
        param_model->JzzR_thigh / (4 * pow(param_model->L, 2) * pow(sin(state_model->q[1] / 2), 2));
    param_model->MatInertia_RW[1] = param_model->JzzR_couple / (2 * pow(param_model->L, 2) * sin(state_model->q[1]));
    param_model->MatInertia_RW[2] = param_model->JzzR_couple / (2 * pow(param_model->L, 2) * sin(state_model->q[1]));
    param_model->MatInertia_RW[3] =
        param_model->JzzR_shank / (4 * pow(param_model->L, 2) * pow(cos(state_model->q[1] / 2), 2));
} // param_model parameter

void sensor_measure(const mjModel* m, mjData* d, StateModel* state_model, ParamTuning* param_tuning)
{
    //Sensor ����� ��� ��Ī�Ǵ����� XML������ simulate�� �־ data load �Ѵ�, �׷���


    /*** (Serial) Joint position ***/
    // state_model->q[0] = d->qpos[0];  // (relative) HFE angle
    // state_model->q[1] = d->qpos[1];  // (relative) KFE angle

    state_model->q[0] = d->sensordata[6]; // (relative) HFE angle
    state_model->q[1] = d->sensordata[7]; // (relative) KFE angle
    state_model->touch_sensor = d->sensordata[8];

    /*** Biarticular Transformation ***/
    // state_model->q_bi[0] = d->qpos[0];              // (absolute) HFE angle
    // state_model->q_bi[1] = d->qpos[0] + d->qpos[1]; // (absolute) KFE angle

    state_model->q_bi[0] = d->sensordata[6];                    // (absolute) HFE angle
    state_model->q_bi[1] = d->sensordata[6] + d->sensordata[7]; // (absolute) KFE angle

    // printf("q1 : %f, q2 : %f \n", d->qpos[0], d->qpos[1]);
    // printf("qm : %f, qb : %f \n\n", q_bi[0], q_bi[1]);

    state_model->qdot_bi[0] = d->qvel[1];
    state_model->qdot_bi[1] = d->qvel[1] + d->qvel[2];
    // printf("qd1 : %f, qd2 : %f \n", d->qvel[0], d->qvel[1]);
    // printf("qdm : %f, qdb : %f \n\n", qdot_bi[0], qdot_bi[1]);

    state_model->qddot_bi[0] = d->qacc[1];
    state_model->qddot_bi[1] = d->qacc[1] + d->qacc[2];
    // printf("qdd1 : %f, qdd2 : %f \n", d->qacc[0], d->qacc[1]);S
    // printf("qddm : %f, qddb : %f \n\n", qddot_bi[0], qddot_bi[1]);

    for (int i = 0; i < NDOF_LEG; i++)
    {
        state_model->qdot_bi_tustin[i] =
            tustin_derivative(state_model->q_bi[i], state_model->q_bi_old[i], state_model->qdot_bi_tustin_old[i],
                param_tuning->freq_cut_D);
        state_model->qddot_bi_tustin[i] =
            tustin_derivative(state_model->qdot_bi_tustin[i], state_model->qdot_bi_tustin_old[i],
                state_model->qddot_bi_tustin_old[i], param_tuning->freq_cut_D);
    }
    // printf("qdot_bi[0]: %f, qdot_bi[1]: %f \n qdot_bi_tust[0]: %f, qdot_bi_tust[1]: %f \n\n",
    // state_Model_FL.qdot_bi[0],
    //        state_Model_FL.qdot_bi[1], state_Model_FL.qdot_bi_tustin[0], state_Model_FL.qdot_bi_tustin[1]);
    // printf("qddot_bi[0]: %f, qddot_bi[1]: %f \n qddot_bi_tust[0]: %f, qddot_bi_tust[1]: %f \n\n",
    //        state_Model_FL.qddot_bi[0], state_Model_FL.qddot_bi[1], state_Model_FL.qddot_bi_tustin[0],
    //        state_Model_FL.qddot_bi_tustin[1]);
}

void jacobianRW(ParamModel* param_model, StateModel* state_model)
{
    /*** Rotating Workspace ***/
    state_model->jacbRW[0] = param_model->L * sin(state_model->q[1] / 2);
    state_model->jacbRW[1] = -param_model->L * sin(state_model->q[1] / 2);
    state_model->jacbRW[2] = param_model->L * cos(state_model->q[1] / 2);
    state_model->jacbRW[3] = param_model->L * cos(state_model->q[1] / 2);
    // printf("%f %f \n", JacobianRW[0][0], JacobianRW[0][1]);
    // printf("%f %f \n\n", JacobianRW[1][0], JacobianRW[1][1]);

    state_model->jacbRW_trans[0] = param_model->L * sin(state_model->q[1] / 2);
    state_model->jacbRW_trans[1] = param_model->L * cos(state_model->q[1] / 2);
    state_model->jacbRW_trans[2] = -param_model->L * sin(state_model->q[1] / 2);
    state_model->jacbRW_trans[3] = param_model->L * cos(state_model->q[1] / 2);

    state_model->jacbRW_trans_inv[0] = 1 / (2 * param_model->L * sin(state_model->q[1] / 2));
    state_model->jacbRW_trans_inv[1] = -1 / (2 * param_model->L * sin(state_model->q[1] / 2));
    state_model->jacbRW_trans_inv[2] = 1 / (2 * param_model->L * cos(state_model->q[1] / 2));
    state_model->jacbRW_trans_inv[3] = 1 / (2 * param_model->L * cos(state_model->q[1] / 2));
}

void fwdKinematics_cal(ParamModel* param_model, StateModel* state_model)
{
    state_model->posRW[0] = 2 * param_model->L * cos((state_model->q_bi[1] - state_model->q_bi[0]) / 2); // r
    state_model->posRW[1] = (state_model->q_bi[0] + state_model->q_bi[1]) / 2;                           // qr

    mju_mulMatVec(state_model->velRW, state_model->jacbRW, state_model->qdot_bi_tustin, 2, 2);
    // printf("r : %f, qr : %f \n", posRW[0], posRW[1]);
    // printf("xRdot : %f, yRdot : %f \n\n", velRW[0], velRW[1]);
}

#endif // !__KINEMATICS_H__

