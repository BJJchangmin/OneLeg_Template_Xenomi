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
#include "kinematics.h"
#include "trajectory.h"


mjvFigure figPosRW;         // RW position tracking plot
mjvFigure figFOB;           // RWFOB GRF estimation plot
mjvFigure figTrunkState;    // Trunk state plot

double simEndtime = 100;	// Simulation End Time
StateModel_  state_Model_FL;
controller ctrl_FL; // other class is in main loop
kinematics kin_FL;
trajectory tra_FL;


/***************** Main Controller *****************/
void mycontroller(const mjModel* m, mjData* d)
{
   
    /* Controllers */
    int flag_DOB = 1;           // flag for switching ON/OFF RWDOB
    int flag_admitt = 1;        // flag for switching ON/OFF admittance control
    double time_run = d->time;
    
    //Admittance Control
    ctrl_FL.admittanceCtrl(&state_Model_FL,5,2,5000, flag_admitt); //parameter(omega_n,zeta,k)
    
    // PID Control
    ctrl_FL.pid_gain_pos(200, 10, 150); //(kp,kd,freq)
    state_Model_FL.tau_bi = state_Model_FL.jacbRW_trans * ctrl_FL.PID_pos(&state_Model_FL); // RW position feedback
    
    // DOB control
    state_Model_FL.tau_bi = state_Model_FL.tau_bi + ctrl_FL.DOBRW(&state_Model_FL, 150, flag_DOB);
    
    // Force Observer
    ctrl_FL.FOBRW(&state_Model_FL, 100); // Rotating Workspace Force Observer (RWFOB)
    
   // Torque input Biarticular
    d->ctrl[0] = state_Model_FL.tau_bi[0] + state_Model_FL.tau_bi[1] ;
    d->ctrl[1] = state_Model_FL.tau_bi[1];
        
    

    if (loop_index % data_frequency == 0) {  
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
    mju_copy(d->qpos, m->key_qpos + 0 * m->nq, m->nq); // 얘가 젤 처음에 와야함
    //d->qpos[1] = 0.546812;
    //d->qpos[2] = 2.59478;
    
    
    
    kin_FL.model_param_cal(m, d, &state_Model_FL); // state init is before. Caution Error.
    kin_FL.state_init(m,d, &state_Model_FL); 
    
    // custom controller
    mjcb_control = mycontroller;
    

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
            int cmd_motion_type = 1;
            int mode_admitt = 1;
            
            if (cmd_motion_type == 0)   // Squat
            {
                tra_FL.Squat(d->time, &state_Model_FL);
            }
            else
            {  
                tra_FL.Hold(&state_Model_FL);  // Hold stance
                
            }
            
            kin_FL.sensor_measure(m, d, &state_Model_FL); // get joint sensor data & calculate biarticular angles
            kin_FL.model_param_cal(m, d,&state_Model_FL); // calculate model parameters
            kin_FL.jacobianRW(&state_Model_FL);            // calculate RW Jacobian
            
            if (d->time < 10)
            {
                //printf("ref: %f \n", state_Model_FL.posRW_ref[1]);
                //printf(" qddot_bi_tustin(0) = %f ,%f ", state_Model_FL.qddot_bi_tustin[0],state_Model_FL.qddot_bi_tustin[1]);
                //printf(" qddot_bi(1) = %f ,%f \n", state_Model_FL.qddot_bi[0],state_Model_FL.qddot_bi[1]);
    
            }
            
            
            kin_FL.fwdKinematics_cal(&state_Model_FL);     // calculate RW Kinematics
           
            mj_step(m, d);
           
            kin_FL.state_update(&state_Model_FL);
            ctrl_FL.ctrl_update();

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

