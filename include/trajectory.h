#ifndef __TRAJECTORY_H__
#define	__TRAJECTORY_H__

#include "globVariable.h"

void trajectory_init(ParamSquat *squat, ParamJump *jump, ParamTuning* param_tuning)
{
    // Squat ���� ����
    squat->T_pause = 1.5;
    squat->freq_squat = 1;
    squat->r0 = 0.3536;
    squat->rc = 0.2;

    //Jumping parameter
    param_tuning->K_thrust = 500;
    param_tuning->zeta_thrust = 0;
    param_tuning->K_land = 1000;
    param_tuning->zeta_land = 1;

    jump->r0 = 0.3536; // initail pose 
    jump->rc = 0.20;
    jump->rt = 0.40;
    
    jump->T_stand = 2;
    jump->T_crouch = 1;
    jump->T_pause = 0.5;
    jump->T_land = 1;
    jump->T_recover = 1;
    jump->qd_max = 20;
}


void trajectory_Squat(double t, ParamModel* param_model, ParamSquat* squat, StateModel* state_model)
{
    double deltaR = squat->r0 - squat->rc;

    double T_squat = 1 / squat->freq_squat;
    double T_period = T_squat + squat->T_pause;

    double t_norm = t - T_period * floor(t / T_period); // nominalized time
    double t1 = squat->T_pause;
    double t2 = t1 + T_squat;

    if (0 <= t_norm && t_norm < t1)
    {
        state_model->posRW_ref[0] = squat->r0;
        state_model->posRW_ref[1] = pi / 2;
    }
    else if (t1 <= t_norm && t_norm < t2)
    {
        state_model->posRW_ref[0] = squat->r0 - 0.5 * deltaR * (1 - cos(2 * pi * squat->freq_squat * (t_norm - t1)));
        state_model->posRW_ref[1] = pi / 2;
        state_model->velRW_ref[0] = -deltaR * (pi * squat->freq_squat) * sin(2 * pi * squat->freq_squat * (t_norm - t1));
    }
} // Squat

void trajectory_Jumping(double t, ParamModel* param_model, ParamTuning* param_tuning, ParamJump* jump, StateModel* state_model, int mode_admitt)
{
    double deltaR = jump->rt - jump->rc;

    double vel_max = (2 * param_model->L * pow(jump->qd_max, 2) * pow(sin(acos(jump->rt / (2 * param_model->L))), 2)) / g;

    double T_thrust = 3 * deltaR / vel_max;
    double T_peak = vel_max / g;    
    double T_flight = 0.7;

    double t0 = jump->T_stand;          // crouching timing
    double t1 = t0 + jump->T_crouch;    // pause timing
    double t2 = t1 + jump->T_pause;     // thrust timing
    double t3 = t2 + T_thrust;          // flight timing
    double t4 = t3 + T_flight;          // (Expected) landing timing
    double t5 = t4 + jump->T_land;      // (Expected) recovery timing
    double t6 = t5 + jump->T_recover;   // (Expected) stance timing
    double T_period = t6;               // whole jumping period

    double t_norm = t - T_period * floor(t / T_period); // nominalized time

    if (mode_admitt == 1)
    {
        param_tuning->Ma_new = param_tuning->Ma;
        param_tuning->Ka_new = param_tuning->K_thrust * param_tuning->Ka;
        param_tuning->Ba_new = 2 * param_tuning->zeta_thrust * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
    }

    if (t < t0)    // stance phase
    {
        state_model->posRW_ref[0] = jump->r0;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->Ka;
            param_tuning->Ba_new = param_tuning->Ba;
        }
    }
    else if (t0 <= t && t < t1)   // crouch phase
    {
        state_model->posRW_ref[0] = 0.5 * (jump->r0 + jump->rc) + 0.5 * (jump->r0 - jump->rc) * cos((pi / jump->T_crouch) * (t_norm - t0));
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = -0.5 * (pi / jump->T_crouch) * (jump->r0 - jump->rc) * sin((pi / jump->T_crouch) * (t_norm - t0));
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->Ka;
            param_tuning->Ba_new = param_tuning->Ba;
        }
    }
    else if (t1 <= t && t < t2)   // pause phase
    {
        state_model->posRW_ref[0] = jump->rc;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_thrust* param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_thrust * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (t2 <= t && t < t3)   // thrust phase
    {
        double a0 = jump->rc;
        double a1 = .0;
        double a2 = .0;
        double a3 = deltaR / (T_thrust * T_thrust * T_thrust);

        state_model->posRW_ref[0] = a0 + a1 * (t_norm - t2) + a2 * pow(t_norm - t2, 2) + a3 * (t_norm - t2) * pow(t_norm - t2, 2);
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = a1 + 2 * a2 * (t_norm - t2) + 3 * a3 * pow(t_norm - t2, 2);

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_thrust* param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_thrust * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (t3 <= t && t < 3.7)    // flight phase
    {
        state_model->posRW_ref[0] = jump->rt;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_land * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_land * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (3.7 <= t && t < 8)   // landing phase
    {
        state_model->posRW_ref[0] = jump->rt;//0.5 * (jump->rc + jump->rt) + 0.5 * (jump->rt - jump->rc) * cos((pi / jump->T_land) * (t_norm - t4));
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = -0.5 * (pi / jump->T_land) * (jump->rt - jump->rc) * sin((pi / jump->T_land) * (t_norm - t4));
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = 1500; //critically damper system �� �� Ka= 15000 �� �� �׳��� ���� �׸��� �׶� ���� 382
            param_tuning->Ba_new = 2 * 0.95 * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (8 <= t && t < t6)   // recovery phase
    {
        state_model->posRW_ref[0] = 0.5 * (jump->r0 + jump->rc) + 0.5 * (jump->rc - jump->r0) * cos((pi / jump->T_recover) * (t_norm - t5));
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = -0.5 * (pi / jump->T_recover) * (jump->rc - jump->r0) * sin((pi / jump->T_recover) * (t_norm - t5));
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else
    {
        state_model->posRW_ref[0] = jump->r0;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_land * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_land * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
} // Jumping

void trajectory_Jumping_Noadmittance(double t, ParamModel* param_model, ParamTuning* param_tuning, ParamJump* jump, StateModel* state_model, int mode_admitt)
{
    double deltaR = jump->rt - jump->rc;

    double vel_max = (2 * param_model->L * pow(jump->qd_max, 2) * pow(sin(acos(jump->rt / (2 * param_model->L))), 2)) / g;

    double T_thrust = 3 * deltaR / vel_max;
    double T_peak = vel_max / g;
    double T_flight = 2 * T_peak;

    double t0 = jump->T_stand;          // crouching timing
    double t1 = t0 + jump->T_crouch;    // pause timing
    double t2 = t1 + jump->T_pause;     // thrust timing
    double t3 = t2 + T_thrust;          // flight timing
    double t4 = t3 + T_flight;          // (Expected) landing timing
    double t5 = t4 + jump->T_land;      // (Expected) recovery timing
    double t6 = t5 + jump->T_recover;   // (Expected) stance timing
    double T_period = t6;               // whole jumping period

    double t_norm = t - T_period * floor(t / T_period); // nominalized time


    if (mode_admitt == 1)
    {
        param_tuning->Ma_new = param_tuning->Ma;
        param_tuning->Ka_new = param_tuning->K_thrust * param_tuning->Ka;
        param_tuning->Ba_new = 2 * param_tuning->zeta_thrust * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
    }

    if (t < t0)    // stance phase
    {
        state_model->posRW_ref[0] = jump->r0;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->Ka;
            param_tuning->Ba_new = param_tuning->Ba;
        }
    }
    else if (t0 <= t && t < t1)   // crouch phase
    {
        state_model->posRW_ref[0] = 0.5 * (jump->r0 + jump->rc) + 0.5 * (jump->r0 - jump->rc) * cos((pi / jump->T_crouch) * (t_norm - t0));
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = -0.5 * (pi / jump->T_crouch) * (jump->r0 - jump->rc) * sin((pi / jump->T_crouch) * (t_norm - t0));
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->Ka;
            param_tuning->Ba_new = param_tuning->Ba;
        }
    }
    else if (t1 <= t && t < t2)   // pause phase
    {
        state_model->posRW_ref[0] = jump->rc;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_thrust * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_thrust * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (t2 <= t && t < t3)   // thrust phase
    {
        double a0 = jump->rc;
        double a1 = .0;
        double a2 = .0;
        double a3 = deltaR / (T_thrust * T_thrust * T_thrust);

        state_model->posRW_ref[0] = a0 + a1 * (t_norm - t2) + a2 * pow(t_norm - t2, 2) + a3 * (t_norm - t2) * pow(t_norm - t2, 2);
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = a1 + 2 * a2 * (t_norm - t2) + 3 * a3 * pow(t_norm - t2, 2);

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_thrust * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_thrust * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (t3 <= t && t < 3.81115)// flight phase
    {
        state_model->posRW_ref[0] = jump->rt;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_land * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_land * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (3.81115 <= t && t < 20)   // landing phase

    {
        int t_TD = t;
        double v_0 = -3.0;
        double k = 360;
        double m = 2.5;
        double alpha = sqrt(k / m);
        double ts = 0.0001;

        //Mg Trajectory�� ���� �ʱ� setting �ϴ� �κ�
        double mg_m = 2.5;
        double mg_k = 5000; //����μ� ������ K�� 2000
        double mg_b = 2*2*sqrt(mg_k / mg_m); //12.65
        
        
        double c1 = 4 * mg_m + 2 * mg_b * ts + mg_k * pow(ts,2);
        double c2 = -8 * mg_m + 2 * mg_k * pow(ts, 2);
        double c3 = 4 * mg_m - 2 * mg_b * ts + mg_k * pow(ts, 2);
        state_model->step = -450;//-mg_m*g;

        
        //���� following
        //double FF_Trajectory = -m * alpha * exp(-alpha * (t - 3.8112)) * v_0 * (-alpha * (t - 3.8112) + 2);
        //state_model->tau_bi[0] = state_model->tau_bi[0] + state_model->jacbRW_trans[0] * FF_Trajectory;
        //state_model->tau_bi[1] = state_model->tau_bi[1] + state_model->jacbRW_trans[2] * FF_Trajectory;
        //state_model->posRW_ref[0] = 0.45 + v_0 * (t - 3.8112) * exp(-alpha * (t - 3.8112));
        //state_model->posRW_ref[1] = pi / 2.;
        
        //Mg Trajectory
        state_model->out_acc = (4 * state_model->step - 8 * state_model->step_old + 4 * state_model->step_old2 - c2 * state_model->out_acc_old - c3 * state_model->out_acc_old2) / c1;
        state_model->out_pos = (pow(ts, 2) * state_model->step +2* pow(ts, 2) * state_model->step_old + pow(ts, 2) * state_model->step_old2 - c2 * state_model->out_pos_old - c3 * state_model->out_pos_old2) / c1;
        //printf("leg: %f %f %f\n", state_model->out_acc, state_model->out_acc_old, state_model->out_acc_old2);
        //printf("leg: %f %f %f\n", c1, c2, c3);
        double Fx_Trajectory = param_model->MatInertia_RW[0] * state_model->out_acc;
        double Fy_Trajectory = param_model->MatInertia_RW[2] * state_model->out_acc;
        
        //printf("%f \n", state_model->q[1]);
        
        state_model->tau_ff[0]= state_model->jacbRW_trans[0] * Fx_Trajectory + state_model->jacbRW_trans[1] * Fy_Trajectory;
        state_model->tau_ff[1] = state_model->jacbRW_trans[2] * Fx_Trajectory + state_model->jacbRW_trans[3] * Fy_Trajectory;
        
        //printf("%f %f\n", state_model->jacbRW_trans[0], state_model->jacbRW_trans[1]);
        //printf("%f %f\n\n", state_model->jacbRW_trans[2], state_model->jacbRW_trans[3]);
        //printf("Trajectory: %f %f\n", state_model->tau_bi[0], state_model->tau_bi[1]);

        //printf("%f \n", state_model->posRW[0]);
        state_model->posRW_ref[0] = jump->rt + state_model->out_pos;
        //printf("%f \n", state_model->posRW_ref[0]);
        state_model->posRW_ref[1] = pi / 2.;

        state_model->velRW_ref[0] = -0.5 * (pi / jump->T_land) * (jump->rt - jump->rc) * sin((pi / jump->T_land) * (t_norm - t4));
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = 15000; //critically damper system �� �� Ka= 15000 �� �� �׳��� ���� �׸��� �׶� ���� 382
            param_tuning->Ba_new = 2 * 1 * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }

    else
    {
        state_model->posRW_ref[0] = jump->r0;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_land * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_land * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
} // Jumping

void trajectory_Landing_tau(double t, ParamModel* param_model, ParamTuning* param_tuning, ParamJump* jump, StateModel* state_model, int mode_admitt)
{
    double deltaR = jump->rt - jump->rc;

    double vel_max = (2 * param_model->L * pow(jump->qd_max, 2) * pow(sin(acos(jump->rt / (2 * param_model->L))), 2)) / g;

    double T_thrust = 3 * deltaR / vel_max;
    double T_peak = vel_max / g;
    double T_flight = 2 * T_peak;

    double t0 = jump->T_stand;          // crouching timing
    double t1 = t0 + jump->T_crouch;    // pause timing
    double t2 = t1 + jump->T_pause;     // thrust timing
    double t3 = t2 + T_thrust;          // flight timing
    double t4 = t3 + T_flight;          // (Expected) landing timing
    double t5 = t4 + jump->T_land;      // (Expected) recovery timing
    double t6 = t5 + jump->T_recover;   // (Expected) stance timing
    double T_period = t6;               // whole jumping period

    double t_norm = t - T_period * floor(t / T_period); // nominalized time


    if (mode_admitt == 1)
    {
        param_tuning->Ma_new = param_tuning->Ma;
        param_tuning->Ka_new = param_tuning->K_thrust * param_tuning->Ka;
        param_tuning->Ba_new = 2 * param_tuning->zeta_thrust * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
    }

    if (t < t0)    // stance phase
    {
        state_model->posRW_ref[0] = jump->r0;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->Ka;
            param_tuning->Ba_new = param_tuning->Ba;
        }
    }
    else if (t0 <= t && t < t1)   // crouch phase
    {
        state_model->posRW_ref[0] = 0.5 * (jump->r0 + jump->rc) + 0.5 * (jump->r0 - jump->rc) * cos((pi / jump->T_crouch) * (t_norm - t0));
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = -0.5 * (pi / jump->T_crouch) * (jump->r0 - jump->rc) * sin((pi / jump->T_crouch) * (t_norm - t0));
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->Ka;
            param_tuning->Ba_new = param_tuning->Ba;
        }
    }
    else if (t1 <= t && t < t2)   // pause phase
    {
        state_model->posRW_ref[0] = jump->rc;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_thrust * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_thrust * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (t2 <= t && t < t3)   // thrust phase
    {
        double a0 = jump->rc;
        double a1 = .0;
        double a2 = .0;
        double a3 = deltaR / (T_thrust * T_thrust * T_thrust);

        state_model->posRW_ref[0] = a0 + a1 * (t_norm - t2) + a2 * pow(t_norm - t2, 2) + a3 * (t_norm - t2) * pow(t_norm - t2, 2);
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = a1 + 2 * a2 * (t_norm - t2) + 3 * a3 * pow(t_norm - t2, 2);

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_thrust * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_thrust * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (t3 <= t && t < 3.81115)// flight phase
    {
        state_model->posRW_ref[0] = jump->rt;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_land * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_land * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
    else if (3.81115 <= t && t < 20)   // landing phase

    {
        int t_TD = t;
        double ts = 0.0001;

        //tau Trajectory�� ���� �ʱ� setting �ϴ� �κ�
        double k = 0.3;
        double trajectory_freq = 7;
        double f = 2 * pi * trajectory_freq;
        double tau = 1 / f;

        double mg_m = tau*tau;
        double mg_k = 1; //����μ� ������ K�� 2000
        double mg_b = 2 * tau; //12.65


        double c1 = 4 * mg_m + 2 * mg_b * ts + mg_k * pow(ts, 2);
        double c2 = -8 * mg_m + 2 * mg_k * pow(ts, 2);
        double c3 = 4 * mg_m - 2 * mg_b * ts + mg_k * pow(ts, 2);
        state_model->step = -1;//-mg_m*g;


        //tau Trajectory
        state_model->out_acc = (k * 4 * state_model->step - k * 8 * state_model->step_old + k * 4 * state_model->step_old2 - c2 * state_model->out_acc_old - c3 * state_model->out_acc_old2) / c1;
        state_model->out_pos = (k*pow(ts, 2) * state_model->step + k * 2 * pow(ts, 2) * state_model->step_old + k * pow(ts, 2) * state_model->step_old2 - c2 * state_model->out_pos_old - c3 * state_model->out_pos_old2) / c1;
        //printf("leg: %f %f %f\n", state_model->out_acc, state_model->out_acc_old, state_model->out_acc_old2);
        //printf("leg: %f %f %f\n", c1, c2, c3);
        double Fx_Trajectory = param_model->MatInertia_RW[0] * state_model->out_acc;
        double Fy_Trajectory = param_model->MatInertia_RW[2] * state_model->out_acc;

        //printf("%f \n", state_model->q[1]);

        state_model->tau_ff[0] = state_model->jacbRW_trans[0] * Fx_Trajectory + state_model->jacbRW_trans[1] * Fy_Trajectory;
        state_model->tau_ff[1] = state_model->jacbRW_trans[2] * Fx_Trajectory + state_model->jacbRW_trans[3] * Fy_Trajectory;

        //printf("%f %f\n", state_model->jacbRW_trans[0], state_model->jacbRW_trans[1]);
        //printf("%f %f\n\n", state_model->jacbRW_trans[2], state_model->jacbRW_trans[3]);
        //printf("Trajectory: %f %f\n", state_model->tau_bi[0], state_model->tau_bi[1]);

        //printf("%f \n", state_model->posRW[0]);
        state_model->posRW_ref[0] = jump->rt + state_model->out_pos;
        //printf("%f \n", state_model->posRW_ref[0]);
        state_model->posRW_ref[1] = pi / 2.;

        state_model->velRW_ref[0] = -0.5 * (pi / jump->T_land) * (jump->rt - jump->rc) * sin((pi / jump->T_land) * (t_norm - t4));
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = 15000; //critically damper system �� �� Ka= 15000 �� �� �׳��� ���� �׸��� �׶� ���� 382
            param_tuning->Ba_new = 2 * 1 * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }

    else
    {
        state_model->posRW_ref[0] = jump->r0;
        state_model->posRW_ref[1] = pi / 2.;
        state_model->velRW_ref[0] = .0;
        state_model->velRW_ref[1] = .0;

        if (mode_admitt == 1)
        {
            param_tuning->Ma_new = param_tuning->Ma;
            param_tuning->Ka_new = param_tuning->K_land * param_tuning->Ka;
            param_tuning->Ba_new = 2 * param_tuning->zeta_land * sqrt(param_tuning->Ma_new * param_tuning->Ka_new);
        }
    }
} // Jumping


void trajectory_Hold(StateModel* state_model)
{
    state_model->posRW_ref[0] = 0.3536;
    state_model->posRW_ref[1] = pi / 2;
} // Hold stance

#endif // !__TRAJECTORY_H__
