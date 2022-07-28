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
from sklearn.metrics import r2_score
import seaborn as sns
import random
#from IPython import embed

#Fichiers perso
import load #fichier load.py
import marches
import read_data_kinematics
import read_data_moments
from scipy.interpolate import UnivariateSpline



###### DEFINITIONS DE FONCTIONS ######

def degrees2radians(degrees):
    return degrees*math.pi/180

def radians2degrees(radians):
    return radians*180/math.pi

def cal_torsion(S_right,S_left,D_right,D_left,pos_exo, pos_human,v_exo,v_human):
    p=pos_human
    pos_human=[p[40],p[42],p[45],p[31],p[33],p[36]]
    v=v_human
    v_human=[v[40],v[42],v[45],v[31],v[33],v[36]]
    torque_torsion_right= S_right*(pos_exo[0:2]-pos_human[0:2]) + D_right*(v_exo[0:2]-v_human[0:2])
    torque_torsion_left= S_left*(pos_exo[3:5]-pos_human[3:5]) + D_left*(v_exo[3:5]-v_human[3:5])
    torque_torsion= np.zeros(6)
    torque_torsion[0:2]=torque_torsion_right
    torque_torsion[3:5]=torque_torsion_left
    return torque_torsion      

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

studied_walk=1

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

#####
found_begin='false'
found_end='false'
for i in range(0,len(marche_exo.time)):
        if marche_exo.heelstrike_left[0]-255*(1/200) < marche_exo.time[i] and found_begin != 'true':
            begin_left_cycle_exo=i
            found_begin='true'
        elif marche_exo.heelstrike_left[1]-255*(1/200) < marche_exo.time[i] and found_end != 'true':
            end_left_cycle_exo=i
            found_end='true'


found_begin='false'
found_end='false'
for i in range(0,len(marche.time)):
        if marche.heelstrike_left[0]  < marche.time[i] and found_begin != 'true':
            begin_left_cycle_human=i
            found_begin='true'
        elif marche.heelstrike_left[1]  < marche.time[i] and found_end != 'true' and found_begin == 'true':
            end_left_cycle_human=i
            found_end='true'




##################################################################### MAKE THE HUMAN WALK #####################################################################################################

#Joints Human:

#  1 : pelvis bending
# 31 : left hip0
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
kp_posture = 5000.0
kd_posture=10.0

# Tasks gains exo
kp_contact = 10.0                 
kp_com = 10000     
kp_posture_exo = 5000.0
kd_posture_exo=30

# Tasks priority levels
level_com = 1
level_contact = 0
level_posture = 1

# Definition of the goal exo
left_hip_exo=degrees2radians(marche_exo.left[0])
left_knee_exo=degrees2radians(marche_exo.left[1])
left_ankle_exo=0*degrees2radians(marche_exo.left[2]) #cheville
right_hip_exo=degrees2radians(marche_exo.right[0])
right_knee_exo=degrees2radians(marche_exo.right[1])
right_ankle_exo=0*degrees2radians(marche_exo.right[2])

# Definition of the goal human
left_hip=degrees2radians(marche.left[0])
left_knee=degrees2radians(marche.left[1])
left_ankle=degrees2radians(marche.left[2]) #cheville
right_hip=degrees2radians(marche.right[0])
right_knee=degrees2radians(marche.right[1])
right_ankle=degrees2radians(marche.right[2])



#Rechantillonnage à la même taille que le cycle humain
#Left_hip
a=left_hip_exo[begin_left_cycle_exo:end_left_cycle_exo]
old_indices = np.arange(0,len(a))
new_length = len(left_hip[begin_left_cycle_human:end_left_cycle_human])
new_indices = np.linspace(0,len(a)-1,new_length)
spl = UnivariateSpline(old_indices,a,k=3,s=0)
left_hip_exo_cycle = spl(new_indices)
#Left_knee
a=left_knee_exo[begin_left_cycle_exo:end_left_cycle_exo]
old_indices = np.arange(0,len(a))
new_length = len(left_knee[begin_left_cycle_human:end_left_cycle_human])
new_indices = np.linspace(0,len(a)-1,new_length)
spl = UnivariateSpline(old_indices,a,k=3,s=0)
left_knee_exo_cycle = spl(new_indices)
#Right_hip
a=right_hip_exo[begin_left_cycle_exo:end_left_cycle_exo]
old_indices = np.arange(0,len(a))
new_length = len(right_hip[begin_left_cycle_human:end_left_cycle_human])
new_indices = np.linspace(0,len(a)-1,new_length)
spl = UnivariateSpline(old_indices,a,k=3,s=0)
right_hip_exo_cycle = spl(new_indices)
#Rightt_knee
a=right_knee_exo[begin_left_cycle_exo:end_left_cycle_exo]
old_indices = np.arange(0,len(a))
new_length = len(right_knee[begin_left_cycle_human:end_left_cycle_human])
new_indices = np.linspace(0,len(a)-1,new_length)
spl = UnivariateSpline(old_indices,a,k=3,s=0) 
right_knee_exo_cycle = spl(new_indices)


#Definitions de l'erreur entre la référence et le pathologique pour un même cycle.

Human_left_hip=np.zeros(len(left_hip[begin_left_cycle_human:end_left_cycle_human]))
Human_left_hip[0:len(left_hip[begin_left_cycle_human:end_left_cycle_human])]=left_hip[begin_left_cycle_human:end_left_cycle_human].copy()
Human_left_knee=np.zeros(len(left_knee[begin_left_cycle_human:end_left_cycle_human]))
Human_left_knee[0:len(left_knee[begin_left_cycle_human:end_left_cycle_human])]=left_knee[begin_left_cycle_human:end_left_cycle_human].copy()
Human_right_hip=np.zeros(len(right_hip[begin_left_cycle_human:end_left_cycle_human]))
Human_right_hip[0:len(right_hip[begin_left_cycle_human:end_left_cycle_human])]=right_hip[begin_left_cycle_human:end_left_cycle_human].copy()
Human_right_knee=np.zeros(len(right_knee[begin_left_cycle_human:end_left_cycle_human]))
Human_right_knee[0:len(right_knee[begin_left_cycle_human:end_left_cycle_human])]=right_knee[begin_left_cycle_human:end_left_cycle_human].copy()

error_ref_patho_left_hip= left_hip_exo_cycle- Human_left_hip
error_ref_patho_left_knee= left_knee_exo_cycle- Human_left_knee
error_ref_patho_right_hip= right_hip_exo_cycle- Human_right_hip
error_ref_patho_right_knee= right_knee_exo_cycle- Human_right_knee

# % de l'erreur pris en compte pour construire la marche de référence à partir de la marche patho
degr=10
degradation_knee=-degr/100
degradation_hip=-degr/100

## Définition d'un couple de torsion
S_torsion_right=0.8
D_torsion_right=0.001
# S_torsion_right=0
# D_torsion_right=0
S_torsion_left=S_torsion_right
D_torsion_left=D_torsion_right
torsion=0*np.ones(6)

#Init human
p=np.zeros(48)
p[31]=left_hip[begin_left_cycle_human]-degradation_hip*error_ref_patho_left_hip[0].copy()
p[33]=left_knee[begin_left_cycle_human]-degradation_knee*error_ref_patho_left_knee[0].copy()
p[36]=0#left_ankle_exo[begin_left_cycle_exo]
p[40]=right_hip[begin_left_cycle_human]-degradation_hip*error_ref_patho_right_hip[0].copy()
p[42]=right_knee[begin_left_cycle_human]-degradation_knee*error_ref_patho_right_knee[0].copy()
p[45]=0#right_ankle_exo[begin_left_cycle_exo]


#Init exo
q=np.zeros(6)
q[0]=right_hip[begin_left_cycle_human]-degradation_hip*error_ref_patho_right_hip[0].copy()
q[1]=right_knee[begin_left_cycle_human]-degradation_knee*error_ref_patho_right_knee[0].copy()
q[2]=0#right_ankle_exo[begin_left_cycle_exo]
q[3]=left_hip[begin_left_cycle_human]-degradation_hip*error_ref_patho_left_hip[0].copy()
q[4]=left_knee[begin_left_cycle_human]-degradation_knee*error_ref_patho_left_knee[0].copy()
q[5]=0#left_ankle_exo[begin_left_cycle_exo]

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
p_goal[31]=left_hip[559].copy()
p_goal[33]=left_knee[559].copy()
p_goal[36]=0#left_ankle_exo[559]
p_goal[40]=right_hip[559].copy()
p_goal[42]=right_knee[559].copy()
p_goal[45]=0#right_ankle_exo[559]
 
postureTask = tsid.TaskJointPosture("task-posture", human_tsid)
postureTask.setKp(kp_posture * np.ones(human_tsid.nq))
postureTask.setKd(2.0 * np.sqrt(kp_posture) * np.ones(human_tsid.nq)) 

invdyn.addMotionTask(postureTask, w_posture, level_posture, 0.0)
trajPosture = tsid.TrajectoryEuclidianConstant("traj_joint", p_goal)
samplePosture = trajPosture.computeNext() 

#Posture Task exo
q_goal=np.zeros(6)
q_goal[0]=right_hip[559].copy()
q_goal[1]=right_knee[559].copy()
q_goal[2]=0 #right_ankle_exo[559]
q_goal[3]=left_hip[559].copy()
q_goal[4]=left_knee[559].copy()
q_goal[5]=0#left_ankle_exo[559]

q_goal = np.array([right_hip[559],right_knee[559],0,left_hip[559],left_knee[559],0]).copy() 
postureTask_exo = tsid.TaskJointPosture("task-posture-exo", robot_tsid)
postureTask_exo.setKp(kp_posture_exo * np.ones(robot_tsid.nq))
postureTask_exo.setKd(2.0 * np.sqrt(kp_posture_exo) * np.ones(robot_tsid.nq)) 
invdyn_exo.addMotionTask(postureTask_exo, w_posture, level_posture, 0.0)
trajPosture_exo = tsid.TrajectoryEuclidianConstant("traj_joint_exo", q_goal)
samplePosture_exo = trajPosture_exo.computeNext() 





# Stockage des valeurs calculées 
left_hip_computed= np.zeros(len(marche.time))
left_knee_computed= np.zeros(len(marche.time))
left_ankle_computed= np.zeros(len(marche.time))
right_hip_computed= np.zeros(len(marche.time))
right_knee_computed= np.zeros(len(marche.time))
right_ankle_computed= np.zeros(len(marche.time))

left_hip_computed_human= np.zeros(len(marche.time))
left_knee_computed_human= np.zeros(len(marche.time))
left_ankle_computed_human= np.zeros(len(marche.time))
right_hip_computed_human= np.zeros(len(marche.time))
right_knee_computed_human= np.zeros(len(marche.time))
right_ankle_computed_human= np.zeros(len(marche.time))

tau_computed= np.zeros(len(marche.time))
# Solver initialization
dt=0.005
solver = tsid.SolverHQuadProgFast("qp solver")
solver.resize(invdyn.nVar, invdyn.nEq, invdyn.nIn) 

solver_exo=tsid.SolverHQuadProgFast("qp solver_exo")
solver_exo.resize(invdyn_exo.nVar, invdyn_exo.nEq, invdyn_exo.nIn) 

HQPData = invdyn.computeProblemData(t, p, v)
HQPData_exo = invdyn_exo.computeProblemData(t, q, v_exo)

sol = solver.solve(HQPData)
sol_exo = solver_exo.solve(HQPData_exo)

tau = invdyn.getActuatorForces(sol)
tau_exo = invdyn.getActuatorForces(sol_exo)


p_mes_human=p.copy()
v_mes_human=v.copy()
p_mes_exo=q.copy()
v_mes_exo=v_exo.copy()
v_mes_exo=v_exo.copy()
vector_exo_prec=v_exo.copy()
j=begin_left_cycle_exo-begin_left_cycle_human

# #Init exo
# q=np.zeros(6)
# q[0]=right_hip[begin_left_cycle_human]-degradation_hip*error_ref_patho_right_hip[0]
# q[1]=right_knee[begin_left_cycle_human]-degradation_knee*error_ref_patho_right_knee[0]
# q[2]=0#right_ankle_exo[begin_left_cycle_exo]
# q[3]=left_hip[begin_left_cycle_human]-degradation_hip*error_ref_patho_left_hip[0]
# q[4]=left_knee[begin_left_cycle_human]-degradation_knee*error_ref_patho_left_knee[0]
#q[5]=0#left_ankle_exo[begin_left_cycle_exo]

for i in range(begin_left_cycle_human, end_left_cycle_human):
    time_start = time.time()    
    
    #Human
    vector=np.zeros(48)
    vector[31]=left_hip[i].copy()
    vector[33]=left_knee[i].copy()
    vector[36]=0
    vector[40]=right_hip[i].copy()
    vector[42]=right_knee[i].copy()
    vector[45]=0

    samplePosture.value(vector)
    postureTask.setReference(samplePosture)   

    HQPData = invdyn.computeProblemData(t,p,v)
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

    #Application du couple de torsion sur l'humain
    torsion_humain=np.zeros(48)
    torsion_humain[31]=torsion[3]
    torsion_humain[33]=torsion[4]
    torsion_humain[36]=torsion[5]
    torsion_humain[40]=torsion[0]
    torsion_humain[42]=torsion[1]
    torsion_humain[45]=torsion[2]
    #print(torsion)

    dv=pin.aba(model,human_pinocchio.data,p_mes_human,v_mes_human,tau-torsion_humain)  ## Positif pour l'humain
    v_mean = v_mes_human + 0.5*dt*dv
    v_mes_human += dt*dv
    p_mes_human = pin.integrate(model, p_mes_human, dt*v_mean)

    #Exo
    new_i=i-begin_left_cycle_human

    vector_exo=np.array([right_hip[i]-(degradation_hip*error_ref_patho_right_hip[new_i]),right_knee[i]-(degradation_knee*error_ref_patho_right_knee[new_i]),0,left_hip[i]-(degradation_hip*error_ref_patho_left_hip[new_i]),left_knee[i]-(degradation_knee*error_ref_patho_left_knee[new_i]),0]).copy()
    vetor_exo_prec=vector_exo.copy()
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
    v_exo_prec=v_exo.copy()
    v_exo += dt*dv_exo
    q = pin.integrate(model_exo, q, dt*v_mean_exo)
    t += dt


    dv=pin.aba(model_exo,robot_pinocchio.data,p_mes_exo,v_mes_exo,tau_exo-torsion)  ## Negatif pour l'exo
    v_mean_exo = v_mes_exo + 0.5*dt*dv
    v_mes_exo += dt*dv
    p_mes_exo = pin.integrate(model_exo, p_mes_exo, dt*v_mean_exo)


     #### PD+ ######
    # P=3
    # D=0.1
    # tau_pd=P*(vector_exo-p_mes_exo) + D*(v_exo- v_mes_exo)
    # dv=pin.aba(model_exo,robot_pinocchio.data,p_mes_exo,v_mes_exo,tau_pd)  ## Negatif pour l'exo
    # v_mean_exo = v_mes_exo + 0.5*dt*dv
    # v_mes_exo += dt*dv
    # p_mes_exo = pin.integrate(model_exo, p_mes_exo, dt*v_mean_exo)
    # t += dt


    pin.forwardKinematics(model_exo,robot_pinocchio.data,p_mes_exo)
    pin.forwardKinematics(model,human_pinocchio.data,p_mes_human)

    ### ALlocation
    left_hip_computed[i]= p_mes_exo[3]
    left_knee_computed[i]= p_mes_exo[4]
    left_ankle_computed[i]= 0
    right_hip_computed[i]=  p_mes_exo[0]
    right_knee_computed[i]= p_mes_exo[1]
    right_ankle_computed[i]= 0
    
    left_hip_computed_human[i]= p_mes_human[31]
    left_knee_computed_human[i]= p_mes_human[33]
    left_ankle_computed_human[i]= 0
    right_hip_computed_human[i]=  p_mes_human[40]
    right_knee_computed_human[i]= p_mes_human[42]
    right_ankle_computed_human[i]= 0

    if i%DISPLAY_N == 0: 
        robot_display.robot.display(p_mes_exo)
        human_display.robot.display(p_mes_human)

    #Calcul de la torsion pour le pas suivant
        torsion= cal_torsion(S_torsion_right,S_torsion_left,D_torsion_right,D_torsion_left,p_mes_exo,p_mes_human,v_mes_exo,v_mes_human)
        # torsion[2]=0
        # torsion[5]=0
        #print(torsion)  
    
    # if i%DISPLAY_N == 0: 
    #     robot_display.robot.display(q)
    #     human_display.robot.display(p)


    time_spent = time.time() - time_start
    if(time_spent < dt): time.sleep(dt-time_spent)
    time.sleep(0.005)



