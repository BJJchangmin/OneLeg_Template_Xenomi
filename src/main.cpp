#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> //for bool
//#include<unistd.h> //for usl eep
#include <math.h>
//#include <resource.h>

#include <GLFW/glfw3.h>
#include <mujoco/mujoco.h>

#include "controller.h"
#include "dataLogging.h"
#include "animation.h"


mjvFigure figPosRW;         // RW position tracking plot
mjvFigure figFOB;           // RWFOB GRF estimation plot
mjvFigure figTrunkState;    // Trunk state plot

double simEndtime = 100;	// Simulation End Time

ParamModel  param_Model_FL;
ParamTuning param_Tuning_FL;
StateModel  state_Model_FL;

ParamJump   param_Jump;
ParamSquat  param_Squat;



/***************** Main Controller *****************/
void mycontroller(const mjModel* m, mjData* d)
{
   
    /* Controllers */
    int flag_DOB = 1;           // flag for switching ON/OFF RWDOB
    int flag_admitt = 0;        // flag for switching ON/OFF admittance control
    double time_run = d->time;
    
    
    //admittanceCtrl(&param_Model_FL, &param_Tuning_FL, &state_Model_FL, flag_admitt);    // Admittance control

    posFeedbackPD(&param_Tuning_FL, &state_Model_FL,d->time);       // RW position feedback

    mju_mulMatVec(state_Model_FL.tau_bi, state_Model_FL.jacbRW_trans, state_Model_FL.ctrl_input_RW, 2, 2);  // Jacobian transpose for joint torque implementation
    
    state_Model_FL.tau_bi[0] = state_Model_FL.tau_bi[0] + state_Model_FL.tau_ff[0];
    state_Model_FL.tau_bi[1] = state_Model_FL.tau_bi[1] + state_Model_FL.tau_ff[1];

    distObserverRW(&param_Model_FL, &param_Tuning_FL, &state_Model_FL, flag_DOB);// Rotating Workspace Disturbance Observer (RWDOB)
    forceObserverRW(d, &param_Model_FL, &param_Tuning_FL, &state_Model_FL); // Rotating Workspace Force Observer (RWFOB)
    
   // Torque input
    d->ctrl[0] = state_Model_FL.tau_bi[0] + state_Model_FL.tau_bi[1] ;
    d->ctrl[1] = state_Model_FL.tau_bi[1] ;
        
    

    if (loop_index % data_frequency == 0) {     // loop_index�� data_frequency�� ���� �������� 0�̸� �����͸� ����
        save_data(m, d, &state_Model_FL);
    }
    loop_index += 1;
}


/***************** Main Function *****************/
int main(int argc, const char** argv)
{
    // activate software
    mj_activate("mjkey.txt");

    // load and compile model
    char error[1000] = "Could not load binary model";

    // check command-line arguments
    if (argc < 2)
        m = mj_loadXML(filename, 0, error, 1000);

    else
        if (strlen(argv[1]) > 4 && !strcmp(argv[1] + strlen(argv[1]) - 4, ".mjb"))
            m = mj_loadModel(argv[1], 0);
        else
            m = mj_loadXML(argv[1], 0, error, 1000);
    if (!m)
        mju_error_s("Load model error: %s", error);


    // make data
    d = mj_makeData(m);

    // Initialize GLFW
    if (!glfwInit())
        mju_error("Could not initialize GLFW");

    // create window, make OpenGL context current, request v-sync
    GLFWwindow* window = glfwCreateWindow(1244, 700, "Demo", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // initialize visualization data structures
    mjv_defaultCamera(&cam);
    mjv_defaultOption(&opt);
    mjv_defaultScene(&scn);
    mjr_defaultContext(&con);
    mjv_makeScene(m, &scn, 2000);                // space for 2000 objects
    mjr_makeContext(m, &con, mjFONTSCALE_150);   // model-specific context

    // install GLFW mouse and keyboard callbacks
    glfwSetKeyCallback(window, keyboard);
    glfwSetCursorPosCallback(window, mouse_move);
    glfwSetMouseButtonCallback(window, mouse_button);
    glfwSetScrollCallback(window, scroll);

    double arr_view[] = {-88.95, -17.5, 1.8, 0.04, 0.000000, 0.27};
    cam.azimuth = arr_view[0];
    cam.elevation = arr_view[1];
    cam.distance = arr_view[2];
    cam.lookat[0] = arr_view[3];
    cam.lookat[1] = arr_view[4];
    cam.lookat[2] = arr_view[5];
    fid = fopen(datapath, "w");
    init_save_data();

   

    // Initialization
    mju_copy(d->qpos, m->key_qpos + 0 * m->nq, m->nq);
    model_param_cal(m, d, &param_Model_FL, &state_Model_FL);
    state_init(m, d, &state_Model_FL, &param_Model_FL);
    //trajectory_init(&param_Squat, &param_Jump, &param_Tuning_FL);
    controller_init(m, d, &param_Model_FL, &param_Tuning_FL);
    
    // custom controller
    mjcb_control = mycontroller;
    //printf(" %f  %f \n", d->ctrl[0], d->ctrl[1]);
    /***************** Simulation Loop *****************/
    // use the first while condition if you want to simulate for a period.
    int i = 0;
    while (!glfwWindowShouldClose(window))
    {
        // advance interactive simulation for 1/60 sec
        // Assuming MuJoCo can simulate faster than real-time, which it usually can,
        // this loop will finish on time for the next frame to be rendered at 60 fps.
        // Otherwise add a cpu timer and exit this loop when it is time to render.

        mjtNum simstart = d->time;
        //printf(" %f  %f \n", d->ctrl[0], d->ctrl[1]);
        state_Model_FL.time = d->time;
        while (d->time - simstart < 1.0 / 60.0)
        {
            /* Trajectory Generation */
            int cmd_motion_type = 5;
            int mode_admitt = 1;
            
            if (cmd_motion_type == 0)   // Squat
            {
                trajectory_Squat(d->time, &param_Model_FL, &param_Squat, &state_Model_FL);
            }
            else if (cmd_motion_type == 1)  // Jumping
            {
                trajectory_Jumping_Noadmittance(d->time, &param_Model_FL, &param_Tuning_FL, &param_Jump, &state_Model_FL, mode_admitt);

            }
            else if (cmd_motion_type == 2)  // Jumping
            {
                trajectory_Landing_tau(d->time, &param_Model_FL, &param_Tuning_FL, &param_Jump, &state_Model_FL, mode_admitt);
            }
            
            else
            {  
                trajectory_Hold(&state_Model_FL);  // Hold stance
                
            }

            printf("ref: %f \n", state_Model_FL.posRW_ref[1]);

            sensor_measure(m, d, &state_Model_FL, &param_Tuning_FL); // get joint sensor data & calculate biarticular angles
            model_param_cal(m, d, &param_Model_FL, &state_Model_FL); // calculate model parameters
            jacobianRW(&param_Model_FL, &state_Model_FL);            // calculate RW Jacobian
            fwdKinematics_cal(&param_Model_FL, &state_Model_FL);     // calculate RW Kinematics
           
            mj_step(m, d);
           
           
            state_update(&state_Model_FL); 

        }

        if (d->time >= simEndtime) {
            fclose(fid);
            break;
        }
        //printf("%f \n", state_Model_FL.deltaPos[0]);
        // get framebuffer viewport
        mjrRect viewport = { 0, 0, 0, 0 };
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

        


        // update scene and render
        //opt.frame = mjFRAME_WORLD;
        //cam.lookat[0] = d->qpos[0];
        //cam.lookat[1] = 0;
        //cam.lookat[2] = 0;
        mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
        mjr_render(viewport, &scn, &con);
        //printf("{%f, %f, %f, %f, %f, %f};\n",cam.azimuth,cam.elevation, cam.distance,cam.lookat[0],cam.lookat[1],cam.lookat[2]);

        // swap OpenGL buffers (blocking call due to v-sync)
        glfwSwapBuffers(window);
        
        // process pending GUI events, call GLFW callbacks
        glfwPollEvents();
    }

    // free visualization storage
    mjv_freeScene(&scn);
    mjr_freeContext(&con);

    // free MuJoCo model and data, deactivate
    mj_deleteData(d);
    mj_deleteModel(m);
    mj_deactivate();

    // terminate GLFW (crashes with Linux NVidia drivers)
#if defined(__APPLE__) || defined(_WIN32)
    glfwTerminate();
#endif

    return 1;
}

