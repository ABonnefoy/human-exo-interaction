import numpy as np
from numpy import nan
from numpy.linalg import norm as norm
import matplotlib.pyplot as plt
import plot_utils as plut
import time
import pinocchio as pin
import tsid
import gepetto.corbaserver
import subprocess
import os
import ur5_conf as conf
from example_robot_data.robots_loader import getModelPath
from pinocchio.robot_wrapper import RobotWrapper
from mat4py import loadmat
import math



#Fichiers perso
import load #fichier load.py
import marches
import read_data_kinematics
import read_data_moments


###### DEFINITIONS DE FONCTIONS ######

def degrees2radians(degrees):
    return degrees*math.pi/180

def radians2degrees(radians):
    return radians*180/math.pi

####### FIN DE DEFINITION DE FONCTION ######


#Ne pas oublier de lancer dans un terminal autre "gepetto-gui"

print("".center(conf.LINE_WIDTH,'#'))
print(" Creating an EXOSKELETON - TWINS")
print("".center(conf.LINE_WIDTH,'#'), '\n')

# AFFICHAGE
DISPLAY=True

# PATH & URDF
path="/home/dmsm/s.otmani/Documents/Exosquelettes/"
urdf=path+"test_urdf.urdf"
path_human="/home/dmsm/s.otmani/Documents/Exosquelettes/"
urdf_human = path_human +"URDF_Human.urdf"

#Affichage des exosquelettes et de l'humain
if __name__ == '__main__':
    human_display =load.HumanNew()
    robot_display = load.Exo()
    if DISPLAY:
        human_display.robot.initViewer(loadModel=True)
        human_display.robot.display(human_display.robot.q0)
        robot_display.robot.initViewer(loadModel=True)
        robot_display.robot.display(robot_display.robot.q0)


# Creation du robot TSID
human_tsid = tsid.RobotWrapper(human_display.df_path, [str(human_display.model_path)], False)
robot_tsid = tsid.RobotWrapper(robot_display.df_path, [str(robot_display.model_path)], False)

# Creation du robot pinocchio
human_pinocchio = pin.RobotWrapper.BuildFromURDF(urdf_human,package_dirs=None,root_joint=None,verbose=False)
robot_pinocchio = pin.RobotWrapper.BuildFromURDF(urdf,package_dirs=None,root_joint=None,verbose=False)



# Configuration 
model = human_tsid.model()
p=human_pinocchio.q0
model = robot_tsid.model()
q=robot_pinocchio.q0


human_display.robot.display(p)
robot_display.robot.display(q)
time.sleep(10)
######################################################################### READ DATA (kinematics, moments, events) ##########################################################################################

studied_walk=1


## Read data
data_kinematics= read_data_kinematics.read_data("/home/dmsm/s.otmani/Documents/Données_marche/fichiers txt/marche sans attelle "+str(studied_walk)+".txt", "/home/dmsm/s.otmani/Documents/Données_marche/fichiers txt/marche sans attelle "+str(studied_walk)+".csv")
right=[data_kinematics.hip_flexion_r, data_kinematics.knee_angle_r, data_kinematics.ankle_angle_r]
left=[data_kinematics.hip_flexion_l, data_kinematics.knee_angle_l, data_kinematics.ankle_angle_l]

#Read Moments
data_moment=read_data_moments.read_data_moments("/home/dmsm/s.otmani/Documents/Données_marche/InverseDynamics/résultats sto & txt/marche sans attelle "+str(studied_walk)+".txt","/home/dmsm/s.otmani/Documents/Données_marche/InverseDynamics/résultats sto & txt/marche sans attelle  "+str(studied_walk)+".csv")
right_moments=[data_moment.hip_flexion_r_moment, data_moment.knee_angle_r_moment, data_moment.ankle_angle_r_moment]
left_moments=[data_moment.hip_flexion_l_moment, data_moment.knee_angle_l_moment, data_moment.ankle_angle_l_moment]

#Walk construction
marche=marches.Marche('marche'+str(studied_walk),data_kinematics.time,left,right,data_moment.time,left_moments,right_moments)


#Read events
path_events="/home/dmsm/s.otmani/Documents/Données_marche/dataEvents/"
events= loadmat(path_events+'marche sans attelle '+str(studied_walk)+'Events.mat')
marche.toesoff(data_kinematics.time,events["Left_Foot_Off"],events["Right_Foot_Off"])
marche.heelstrike(data_kinematics.time,events["Left_Foot_Strike"],events["Right_Foot_Strike"])
#marche.printMarche() #print l'ensemble des données de la marche
[duration,begin,end]=marche.cycle('right')



#### Exo data

studied_walk=6

## Read data
data_kinematics_exo= read_data_kinematics.read_data("/home/dmsm/s.otmani/Documents/Données_marche/fichiers txt/marche"+str(studied_walk)+".txt", "/home/dmsm/s.otmani/Documents/Données_marche/fichiers txt/marche"+str(studied_walk)+".csv")
right_exo=[data_kinematics_exo.hip_flexion_r, data_kinematics_exo.knee_angle_r, data_kinematics_exo.ankle_angle_r]
left_exo=[data_kinematics_exo.hip_flexion_l, data_kinematics_exo.knee_angle_l, data_kinematics_exo.ankle_angle_l]

#Read Moments
data_moment_exo=read_data_moments.read_data_moments("/home/dmsm/s.otmani/Documents/Données_marche/InverseDynamics/résultats sto & txt/marche "+str(studied_walk)+".txt","/home/dmsm/s.otmani/Documents/Données_marche/InverseDynamics/résultats sto & txt/marche "+str(studied_walk)+".csv")
right_moments_exo=[data_moment_exo.hip_flexion_r_moment, data_moment_exo.knee_angle_r_moment, data_moment_exo.ankle_angle_r_moment]
left_moments_exo=[data_moment_exo.hip_flexion_l_moment, data_moment_exo.knee_angle_l_moment, data_moment_exo.ankle_angle_l_moment]

#Walk construction
marche_exo=marches.Marche('marche'+str(studied_walk),data_kinematics_exo.time,left_exo,right_exo,data_moment_exo.time,left_moments_exo,right_moments_exo)


#Read events
path_events_exo="/home/dmsm/s.otmani/Documents/Données_marche/dataEvents/"
events_exo= loadmat(path_events_exo+'marche '+str(studied_walk)+'Events.mat')
marche_exo.toesoff(data_kinematics_exo.time,events_exo["Left_Foot_Off"],events_exo["Right_Foot_Off"])
marche_exo.heelstrike(data_kinematics_exo.time,events_exo["Left_Foot_Strike"],events_exo["Right_Foot_Strike"])
marche_exo.heelstrike(data_kinematics_exo.time,events_exo["Left_Foot_Strike"],events_exo["Right_Foot_Strike"])


#####


##################################################################### MAKE THE HUMAN WALK #####################################################################################################

#Joints Human:

#  1 : pelvis bending
# 31 : left hip
# 33 : left_knee
# 36 : left_ankle
# 40 : right_hip
# 42 : right_knee
# 45 : right_ankle


#Joints exo
# 0  : right_hip
# 1  : right_knee
# 2  : right_ankle
# 3  : left_hip
# 4  : right_knee
# 5  : right_ankle

########################## Decalage du COM ######################################

#Utilisation du freeflyer


########################## Test définition d'une tache de posture ###############

DISPLAY_N = 25
PLOT_JOINT_POS = 1

# Tasks weights
w_com = 1.0                       
w_posture = 1                  
w_forceRef = 1e-5                 

# Tasks gains
kp_contact = 10.0                 
kp_com = 10000     
kp_posture = 10000.0
kd_posture=300.0

# Tasks gains exo
kp_contact = 10.0                 
kp_com = 10000     
kp_posture_exo = 40020.0
kd_posture_exo=30

# Tasks priority levels
level_com = 1
level_contact = 0
level_posture = 1

# Definition of the goal exo
left_hip_exo=degrees2radians(marche_exo.left[0])
left_knee_exo=degrees2radians(marche_exo.left[1])
left_ankle_exo=degrees2radians(marche_exo.left[2]) #cheville
right_hip_exo=degrees2radians(marche_exo.right[0])
right_knee_exo=degrees2radians(marche_exo.right[1])
right_ankle_exo=degrees2radians(marche_exo.right[2])

# Definition of the goal human
# left_hip=degrees2radians(marche.left[0])
# left_knee=degrees2radians(marche.left[1])
# left_ankle=degrees2radians(marche.left[2]) #cheville
# right_hip=degrees2radians(marche.right[0])
# right_knee=degrees2radians(marche.right[1])
# right_ankle=degrees2radians(marche.right[2])

#Init human
p=np.zeros(48)
p[31]=left_hip_exo[0]
p[33]=left_knee_exo[0]
p[36]=left_ankle_exo[0]
p[40]=right_hip_exo[0]
p[42]=right_knee_exo[0]
p[45]=right_ankle_exo[0]


#Init exo
q=np.zeros(6)
q[0]=right_hip_exo[0]
q[1]=right_knee_exo[0]
q[2]=right_ankle_exo[0]
q[3]=left_hip_exo[0]
q[4]=left_knee_exo[0]
q[5]=left_ankle_exo[0]

human_display.robot.display(p)
robot_display.robot.display(q)


# Dynamics Problem initialization
t = 0.0 # time
pos_com=human_tsid.com(human_tsid.data())
pos_com_exo=robot_tsid.com(robot_tsid.data())


v = np.zeros(human_tsid.nv)
v_exo= np.zeros(robot_tsid.nv)

invdyn = tsid.InverseDynamicsFormulationAccForce("tsid", human_tsid, False)
invdyn_exo = tsid.InverseDynamicsFormulationAccForce("tsid-exo", robot_tsid, False)

invdyn.computeProblemData(t, p, v)
invdyn_exo.computeProblemData(t, q, v_exo)

data = invdyn.data()
model = human_tsid.model()

data_exo = invdyn_exo.data()
model_exo = robot_tsid.model()

#Posture Task Human
p_goal=np.zeros(48)
p_goal[31]=left_hip_exo[559]
p_goal[33]=left_knee_exo[559]
p_goal[36]=left_ankle_exo[559]
p_goal[40]=right_hip_exo[559]
p_goal[42]=right_knee_exo[559]
p_goal[45]=right_ankle_exo[559]
 
postureTask = tsid.TaskJointPosture("task-posture", human_tsid)
postureTask.setKp(kp_posture * np.ones(human_tsid.nq))
postureTask.setKd(2.0 * np.sqrt(kp_posture) * np.ones(human_tsid.nq)) 

invdyn.addMotionTask(postureTask, w_posture, level_posture, 0.0)
trajPosture = tsid.TrajectoryEuclidianConstant("traj_joint", p_goal)
samplePosture = trajPosture.computeNext() 

#Posture Task exo
q_goal=np.zeros(6)
q_goal[0]=right_hip_exo[559]
q_goal[1]=right_knee_exo[559]
q_goal[2]=right_ankle_exo[559]
q_goal[3]=left_hip_exo[559]
q_goal[4]=left_knee_exo[559]
q_goal[5]=left_ankle_exo[559]

q_goal = np.array([right_hip_exo[559],right_knee_exo[559],right_ankle_exo[559],left_hip_exo[559],left_knee_exo[559],left_ankle_exo[559]]).copy() 
postureTask_exo = tsid.TaskJointPosture("task-posture-exo", robot_tsid)
postureTask_exo.setKp(kp_posture_exo * np.ones(robot_tsid.nq))
postureTask_exo.setKd(2.0 * np.sqrt(kp_posture_exo) * np.ones(robot_tsid.nq)) 
invdyn_exo.addMotionTask(postureTask_exo, w_posture, level_posture, 0.0)
trajPosture_exo = tsid.TrajectoryEuclidianConstant("traj_joint_exo", q_goal)
samplePosture_exo = trajPosture_exo.computeNext() 


## Définition d'un couple de torsion
torsion=1000000*np.ones(48)
#torsion[31]=-10000000

# Stockage des valeurs calculées 
left_hip_computed= np.zeros(len(marche_exo.time))
left_knee_computed= np.zeros(len(marche_exo.time))
left_ankle_computed= np.zeros(len(marche_exo.time))
right_hip_computed= np.zeros(len(marche_exo.time))
right_knee_computed= np.zeros(len(marche_exo.time))
right_ankle_computed= np.zeros(len(marche_exo.time))

tau_computed= np.zeros(len(marche_exo.time))
# Solver initialization
dt=0.005
solver = tsid.SolverHQuadProgFast("qp solver")
solver.resize(invdyn.nVar, invdyn.nEq, invdyn.nIn) 

solver_exo=tsid.SolverHQuadProgFast("qp solver_exo")
solver_exo.resize(invdyn_exo.nVar, invdyn_exo.nEq, invdyn_exo.nIn) 

HQPData = invdyn.computeProblemData(t, p, v, torsion)
HQPData_exo = invdyn_exo.computeProblemData(t, q, v_exo)

sol = solver.solve(HQPData)
sol_exo = solver_exo.solve(HQPData_exo)

tau = invdyn.getActuatorForces(sol)
tau_exo = invdyn.getActuatorForces(sol_exo)

for i in range(0, len(marche_exo.time)):
    time_start = time.time()    
    
    #Human
    vector=np.zeros(48)
    vector[31]=left_hip_exo[i]
    vector[33]=left_knee_exo[i]
    vector[36]=left_ankle_exo[i]
    vector[40]=right_hip_exo[i]
    vector[42]=right_knee_exo[i]
    vector[45]=right_ankle_exo[i]

    samplePosture.value(vector)
    postureTask.setReference(samplePosture)   

    HQPData = invdyn.computeProblemData(t, p,v)
    data = invdyn.data()
    if i == 0: HQPData.print_all()

    sol = solver.solve(HQPData)
    if(sol.status!=0):
        print ("QP problem could not be solved! Error code:", sol.status)
        break
    
    tau = invdyn.getActuatorForces(sol)

    dv = invdyn.getAccelerations(sol)

    v_mean = v + 0.5*dt*dv
    v += dt*dv
    p = pin.integrate(model, p, dt*v_mean)


    #Exo
    vector_exo=np.array([right_hip_exo[i],right_knee_exo[i],right_ankle_exo[i],left_hip_exo[i],left_knee_exo[i],left_ankle_exo[i]])
    samplePosture_exo.value(vector_exo)
    postureTask_exo.setReference(samplePosture_exo)    

    HQPData_exo = invdyn_exo.computeProblemData(t, q, v_exo)
    data_exo = invdyn_exo.data()
    if i == 0: HQPData.print_all()

    sol_exo = solver_exo.solve(HQPData_exo)
    if(sol_exo.status!=0):
        print ("QP problem could not be solved! Error code:", sol_exo.status)
        break
    
    tau_exo = invdyn_exo.getActuatorForces(sol_exo)
    dv_exo = invdyn_exo.getAccelerations(sol_exo)

        
    v_mean_exo = v_exo + 0.5*dt*dv_exo
    v_exo += dt*dv_exo
    q = pin.integrate(model_exo, q, dt*v_mean_exo)
    t += dt

    ### ALlocation
    left_hip_computed[i]= q[3]
    left_knee_computed[i]= q[4]
    left_ankle_computed[i]= q[5]
    right_hip_computed[i]=  q[1]
    right_knee_computed[i]= q[2]
    right_ankle_computed[i]= q[3]
    tau_computed[i]=tau[31];
    
    if i%DISPLAY_N == 0: 
        robot_display.robot.display(q)
        human_display.robot.display(p)


    time_spent = time.time() - time_start
    if(time_spent < dt): time.sleep(dt-time_spent)
    time.sleep(0.005)

print(robot_tsid.com(data))
print(tau_computed)

if(PLOT_JOINT_POS):   

    #Côté gauche !
    plt.figure() 

    plt.subplot(221)
    plt.plot(data_kinematics_exo.time,radians2degrees(left_hip_computed))
    plt.plot(data_kinematics_exo.time,radians2degrees(left_hip_exo))
    plt.xlabel('Time (s)')
    plt.ylabel('Angles (degrees)')
    plt.title('Computed and experimental hip angular position - H.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(222)
    plt.plot(data_kinematics_exo.time,radians2degrees(left_knee_computed))
    plt.plot(data_kinematics_exo.time,radians2degrees(left_knee_exo))
    plt.xlabel('Time (s)')
    plt.ylabel('Angles (degrees)')
    plt.title('Computed and experimental knee angular position - H.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(223)
    plt.plot(data_kinematics_exo.time,radians2degrees(left_ankle_computed))
    plt.plot(data_kinematics_exo.time,radians2degrees(left_ankle_exo))
    plt.xlabel('Time (s)')
    plt.ylabel('Angles (degrees)')
    plt.title('Computed and experimental ankle angular position - H.')
    plt.legend(['Computed', 'Experimental'])


    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)
    plt.show()