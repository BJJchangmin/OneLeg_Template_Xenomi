#ifndef __GLOBVARIABLE_H__
#define __GLOBVARIABLE_H__

#define NDOF_LEG 2		// system #(DoF)

const double Ts = 0.0001; // sampling period
const double g = 9.81;    // gravitational accel.
const double pi = 3.141592;

/* Model Variables */
typedef struct modelParam
{
    /* Trunk Parameter */
    double m_hip;         // mass of hip torso
    double m_trunk_front; // mass of front trunk
    double m_trunk_rear;  // mass of rear trunk
    double m_trunk;       // total mass of trunk
    double m_total;       // total robot mass

    /* Leg Parameter */
    double L; // leg length : thigh and shank links' length are assumed to be the same

    double m_thigh; // mass of thigh link
    double m_shank; // mass of shank link
    double m_leg;   // mass of leg

    double d_thigh; // CoM pos of thigh w.r.t HFE
    double d_shank; // CoM pos of shank w.r.t KFE

    double Izz_thigh; // MoI(z) of thigh w.r.t its CoM
    double Izz_shank; // MoI(z) of shank w.r.t its CoM

    double Jzz_thigh; // MoI(z) of thigh w.r.t HFE
    double Jzz_shank; // MoI(z) of shank w.r.t KFE

    double JzzR_thigh;
    double JzzR_shank;
    double JzzR_couple;

    double MatInertia_bi[NDOF_LEG * NDOF_LEG];
    double MatInertia_RW[NDOF_LEG * NDOF_LEG];
} ParamModel;

typedef struct modelState
{
    /* Trunk States */
    /* Joint Coordinates */
    double q[NDOF_LEG];    // Serial Coordinates
    double q_old[NDOF_LEG]; // gravity,coriolis ��� ������ ����

    double q_bi[NDOF_LEG]; // biarticular joint angle
    double q_bi_old[NDOF_LEG];

    double qdot_bi[NDOF_LEG];  // biarticular joint angular vel (sensor)
    double qdot_bi_old[NDOF_LEG];
    double qddot_bi[NDOF_LEG]; // biarticular joint angular acc (sensor)

    double qdot_bi_tustin[NDOF_LEG]; // biarticular joint angular vel (derivative)
    double qdot_bi_tustin_old[NDOF_LEG];

    double qddot_bi_tustin[NDOF_LEG]; // biarticular joint angular acc (derivative)
    double qddot_bi_tustin_old[NDOF_LEG];

    double tau_bi[NDOF_LEG]; // (Biarticular) joint torques
    double tau_bi_old[NDOF_LEG];

    /* PD Controller */
    double error_pos[NDOF_LEG]; 
    double error_pos_old[NDOF_LEG];
    double error_dot_pos[NDOF_LEG];
    double error_dot_pos_old[NDOF_LEG];

    double error_vel[NDOF_LEG];
    double error_vel_old[NDOF_LEG];
    double error_dot_vel[NDOF_LEG];
    double error_dot_vel_old[NDOF_LEG];

    //double ctrl_kd_old_pos[NDOF_LEG]; // Pd ctrl old output

    /* Rotating Workspace Coordinates */
    double r0;              // initial leg length
    double posRW[NDOF_LEG]; // RW position
    double posRW_old[NDOF_LEG];
    double posRW_ref[NDOF_LEG]; // RW position reference
    double posRW_ref_old[NDOF_LEG];
    double posRW_ref_old2[NDOF_LEG];
    // double error_pos[NDOF_LEG];
    // double error_pos_old[NDOF_LEG];

    double velRW[NDOF_LEG]; // RW velocity
    double velRW_old[NDOF_LEG];
    double velRW_ref[NDOF_LEG]; // RW velocity reference
    double velRW_ref_old[NDOF_LEG];
    // double error_vel[NDOF_LEG];
    // double error_vel_old[NDOF_LEG];

    double ctrl_input_RW[NDOF_LEG]; // control input
    double ctrl_input_RW_old[NDOF_LEG];

    /* Jacobian (Rotating Workspace) */
    double jacbRW[NDOF_LEG * NDOF_LEG];
    double jacbRW_trans[NDOF_LEG * NDOF_LEG];
    double jacbRW_trans_inv[NDOF_LEG * NDOF_LEG];

    /* Rotating Workspace Disturbance Observer */
    double tauDist_hat[NDOF_LEG];

    double lhs_dob_LPF[NDOF_LEG];
    double lhs_dob_LPF_old[NDOF_LEG * NDOF_LEG];

    double rhs_dob_LPF[NDOF_LEG];
    double rhs_dob_LPF_old[NDOF_LEG * NDOF_LEG];

    /* Rotating Workspace Force Observer */
    double tauExt_hat[NDOF_LEG]; // estimated external torque

    double forceExt_hat[NDOF_LEG]; // estimated external force (e.g. GRF)
    double forceExt_hat_old[NDOF_LEG];
    double forceExt_hat_old2[NDOF_LEG];

    double lhs_fob_LPF[NDOF_LEG];
    double lhs_fob_LPF_old[NDOF_LEG];
    double rhs_fob_LPF[NDOF_LEG];
    double rhs_fob_LPF_old[NDOF_LEG];

    // admittacne control
    double deltaPos[NDOF_LEG]; // r-direction deformation
    double deltaPos_old[NDOF_LEG];
    double deltaPos_old2[NDOF_LEG];


    // Mg Trajectory
    double step_old2;
    double step_old;
    double step;

    double out_acc_old2;
    double out_acc_old;
    double out_acc;
    
    double out_pos_old2;
    double out_pos_old;
    double out_pos;
    
    double touch_sensor;
    double tau_ff[NDOF_LEG];

    double time;

} StateModel;

/* Controller Variables */
typedef struct paramTuning
{
    double Kp_pos[NDOF_LEG];
    double Kd_pos[NDOF_LEG];
    double Ki_pos[NDOF_LEG];

    double Kp_vel[NDOF_LEG];
    double Kd_vel[NDOF_LEG];
    double Ki_vel[NDOF_LEG];

    double freq_cut_D;
    double freq_cut_Qd;
    double freq_cut_Qf;

    double delta_default;
    double zeta;
    double Ma;
    double Ba;
    double Ka;

    double zeta_thrust;
    double zeta_land;
    double Ma_new;
    double Ba_new;
    double Ka_new;

    double K_thrust; // Admittance Stiffness Gain in Thrust
    double K_land; // Admittance Stiffness Gain in Landing

} ParamTuning;

/* Trajectory Variables */
typedef struct jumpParam
{
    double r0;
    double rc; // crouching leg length
    double rt; // maximum stretched length

    double qd_max; // max. joint vel (rad/s)

    double T_stand;   // stand phase duration
    double T_crouch;  // crouching phase duration
    double T_pause;   // pause time
    double T_land;    // landing phase duration
    double T_recover; // recovery phase duration
} ParamJump;

typedef struct squatParam
{
    double r0; // initial leg length
    double rc; // crouching leg length

    double freq_squat; // squat frequency
    double T_pause;    // pause duration
} ParamSquat;

#endif // !__GLOBVARIABLE_H__