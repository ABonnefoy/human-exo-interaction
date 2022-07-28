import numpy as np
from numpy import nan
from numpy.linalg import norm as norm
import matplotlib.pyplot as plt
import time
import pinocchio as pin
import tsid
import gepetto.corbaserver
import subprocess
import os
from example_robot_data.robots_loader import getModelPath
from pinocchio.robot_wrapper import RobotWrapper
from mat4py import loadmat
import math
from IPython import embed
import sys
import plotly.graph_objects as go
from sklearn.metrics import r2_score


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

def torsion(S,D,dt,pos_ref, pos_real):
    torque_torsion= S*(pos_ref-pos_real) + D*(pos_ref-pos_real)/dt
    return torque_torsion

####### FIN DE DEFINITION DE FONCTION ######


#Ne pas oublier de lancer dans un terminal autre "gepetto-gui"

print(" Creating an EXOSKELETON - TWINS")

# AFFICHAGE
DISPLAY=True

# PATH & URD
path = os.path.dirname(os.path.realpath(__file__))
urdf = path + '/urdf_augmented.urdf'

#p[31]=left_hip
#p[33]=left_knee
#p[36]=left_ankle
#p[40]=right_hip
#p[42]=right_knee
#p[45]=right_ankle

#Affichage des exosquelettes et de l'humain
if __name__ == '__main__':
    exo_display =load.Exo()
    if DISPLAY:
        exo_display.robot.initViewer(loadModel=True)
        exo_display.robot.viewer.gui.addFloor('world/floor')


# Creation du robot TSID
robot_tsid = tsid.RobotWrapper(exo_display.df_path, [str(exo_display.model_path)], pin.JointModelFreeFlyer(), False)

# Creation du robot pinocchio
robot_pinocchio = pin.RobotWrapper.BuildFromURDF(urdf,package_dirs=None,root_joint=None,verbose=False)

# Configuration 
model = robot_tsid.model()
data = robot_tsid.data()
p=np.zeros(robot_tsid.nq)


#for i, name in enumerate(robot_tsid.model().names): print(i, name)
#for i, name in enumerate(robot_tsid.model().frames): print(i, name)
#sys.exit()

######################################################################### READ DATA (kinematics, moments, events) ##########################################################################################

studied_walk=1


## Read data
data_kinematics= read_data_kinematics.read_data(path+"/Données_marche/fichiers txt/marche"+str(studied_walk)+".txt", path + "/Données_marche/fichiers txt/marche"+str(studied_walk)+".csv")
right=[data_kinematics.hip_flexion_r, data_kinematics.knee_angle_r, data_kinematics.ankle_angle_r]
left=[data_kinematics.hip_flexion_l, data_kinematics.knee_angle_l, data_kinematics.ankle_angle_l]

#Read Moments
data_moment=read_data_moments.read_data_moments(path + "/Données_marche/InverseDynamics/résultats sto & txt/marche "+str(studied_walk)+".txt",path + "/Données_marche/InverseDynamics/résultats sto & txt/marche "+str(studied_walk)+".csv")
right_moments=[data_moment.hip_flexion_r_moment, data_moment.knee_angle_r_moment, data_moment.ankle_angle_r_moment]
left_moments=[data_moment.hip_flexion_l_moment, data_moment.knee_angle_l_moment, data_moment.ankle_angle_l_moment]

#Walk construction
marche=marches.Marche('marche'+str(studied_walk),data_kinematics.time,left,right,data_moment.time,left_moments,right_moments)


#Read events
path_events= path + "/Données_marche/dataEvents/"
events= loadmat(path_events+'marche '+str(studied_walk)+'Events.mat')
marche.toesoff(data_kinematics.time,events["Left_Foot_Off"],events["Right_Foot_Off"])
marche.heelstrike(data_kinematics.time,events["Left_Foot_Strike"],events["Right_Foot_Strike"])

#Read FirstFoot
firstFoot=marche.firstFoot()

#Events timestamps 
dt=0.005
delay = 255*dt

t_hs_left = np.array(marche.heelstrike_left) - delay
t_hs_left = np.append(t_hs_left, np.zeros(1)) # 
print('HS Left ', t_hs_left)
i_hs_left = 0

t_hs_right = np.array(marche.heelstrike_right) - delay
t_hs_right = np.append(t_hs_right, np.zeros(1))
print('HS Right ', t_hs_right)
i_hs_right = 0

t_to_left = np.array(marche.toesoff_left) - delay
t_to_left = np.append(t_to_left, np.zeros(1))
print('TO Left ', t_to_left)
i_to_left = 0

t_to_right = np.array(marche.toesoff_right) - delay
t_to_right = np.concatenate((np.array([0.18]), t_to_right, np.zeros(1)))
print('TO Right ', t_to_right)
i_to_right = 0

#sys.exit()

#Read markers
data_markers=loadmat(path + "/Données_marche/markers_data.mat")
data_markers=data_markers['output']['marker_data']['Markers']
#print(data_markers.keys())

RLM = 1e-3 * np.array(data_markers['RLM']) # /!\ l'ordre du repère est -y, x, z
LLM = 1e-3 * np.array(data_markers['LLM']) # /!\ l'ordre du repère est -y, x, z

RGT = 1e-3 * np.array(data_markers['RGT']) # /!\ l'ordre du repère est -y, x, z
LGT = 1e-3 * np.array(data_markers['LGT']) # /!\ l'ordre du repère est -y, x, z


frame_names = ['joint_right_foot_FE', 'joint_left_foot_FE']
frame_ids = [model.getFrameId('joint_right_foot_FE'), model.getFrameId('joint_left_foot_FE')]

RASIS = 1e-3 * np.array(data_markers['RASIS']) # /!\ l'ordre du repère est -y, x, z
LASIS = 1e-3 * np.array(data_markers['LASIS']) # /!\ l'ordre du repère est -y, x, z

RPSIS = 1e-3 * np.array(data_markers['RPSIS']) # /!\ l'ordre du repère est -y, x, z
LPSIS = 1e-3 * np.array(data_markers['LPSIS']) # /!\ l'ordre du repère est -y, x, z

x_com = (1/4) * (RASIS[0,1] + LASIS[0,1] + RPSIS[0,1] + LPSIS[0,1])
y_com = -(1/4) * (RASIS[0,0] + LASIS[0,0] + RPSIS[0,0] + LPSIS[0,0])
z_com = (1/4) * (RASIS[0,2] + LASIS[0,2] + RPSIS[0,2] + LPSIS[0,2])
com_target = np.array([x_com, y_com, z_com])







# Stockage des valeurs calculées 
left_hip_computed = np.zeros(len(marche.time))
left_knee_computed = np.zeros(len(marche.time))
left_ankle_computed = np.zeros(len(marche.time))
right_hip_computed = np.zeros(len(marche.time))
right_knee_computed = np.zeros(len(marche.time))
right_ankle_computed = np.zeros(len(marche.time))

left_foot_computed = np.zeros((3,len(marche.time)))
right_foot_computed = np.zeros((3,len(marche.time)))
left_foot_contact = len(marche.time) * [None]
right_foot_contact = len(marche.time) * [None]

com_computed = np.zeros((3,len(marche.time)))
com_experimental = np.zeros((3,len(marche.time)))


########################### Test définition d'une tache de posture ###############

DISPLAY_N = 10
PLOT_JOINT_POS = 1
PLOT_WALK = 1

# Tasks weights
w_com = 1           
w_posture_body = 1e-3
w_posture_legs = 1     
w_forceRef = 1e-5
w_foot = 1e-1

# Tasks gains
kp_contact = 1e3     
kp_com = 1e3 * np.array([1., 2., 2.])
kp_foot = 5e3

kp_posture_min = 3e3
kp_posture_max = 5e3
kp_posture = kp_posture_min * np.ones(robot_tsid.nq-7)
kp_posture[2]=kp_posture_max
kp_posture[3]=kp_posture_max
kp_posture[5]=kp_posture_max
kp_posture[8]=kp_posture_max
kp_posture[9]=kp_posture_max
kp_posture[11]=kp_posture_max

# Tasks priority levels
level_com = 1
level_contact = 0
level_posture = 1
level_foot = 1


#Contacts
mu = 1000                        
fMin = 5.0                          
fMax = 1000.0    
f_n = np.array([0., 0., 1.])   


# Goal trajectories
left_hip_FE=degrees2radians(marche.left[0])
left_knee_FE=degrees2radians(marche.left[1])
left_ankle_FE=degrees2radians(marche.left[2]) #cheville
right_hip_FE=degrees2radians(marche.right[0])
right_knee_FE=degrees2radians(marche.right[1])
right_ankle_FE=degrees2radians(marche.right[2])




# Dynamics Problem initialization
t = 0.0 # time



k = 0


v = np.zeros(robot_tsid.nv)
p = np.zeros(robot_tsid.nq)


p[7:]=np.array([0., 0., left_hip_FE[k],left_knee_FE[k],0.,left_ankle_FE[k],0.,0.,right_hip_FE[k],right_knee_FE[k],0.,right_ankle_FE[k]])
v = np.zeros(robot_tsid.nv)

p[3:7] = np.array([0.026, -0.061, -0.002, 0.998]) 

pin.computeAllTerms(model, data, p, v)
pin.updateFramePlacements(model, data)

com_pos = robot_tsid.com(data)
com_offset = com_target - com_pos


#p[:3]  = p[:3] + com_offset
p[:2] += com_offset[:2].copy()

'''p[0] = p[0] + (RLM[0, 1] - data.oMf[frame_ids[0]].translation[0]) # Décalage de x
p[1] = p[1] + (- RLM[0, 0] - data.oMf[frame_ids[0]].translation[1]) # Décalage de y
p[2] = p[2] + (RLM[0, 2] - data.oMf[frame_ids[0]].translation[2]) # Décalage de z'''



pin.computeAllTerms(model, data, p, v)
pin.updateFramePlacements(model, data)

com_pos = robot_tsid.com(data)
com_offset = com_target - com_pos





#print(p)
exo_display.robot.display(p)
#embed()


# Definition of the goal 

x_foot = 0.2
y_foot = 0.1
z_foot = 0.02

right_contact_point = np.ones((3,4)) * (-z_foot)
right_contact_point[0,:] = np.array([0, 0, x_foot, x_foot])
right_contact_point[1,:] = np.array([0, y_foot, 0, y_foot])

left_contact_point = np.ones((3,4)) * (-z_foot)
left_contact_point[0,:] = np.array([0, 0, x_foot, x_foot])
left_contact_point[1,:] = np.array([-y_foot, 0, -y_foot, 0])


invdyn = tsid.InverseDynamicsFormulationAccForce("tsid", robot_tsid, False)
invdyn.computeProblemData(t, p, v)
data = invdyn.data()


# Waist Task creation    
'''waistTask = tsid.TaskSE3Equality("task-waist", robot_tsid, 'base_link')
waistTask.setKp(kp_waist * np.ones(6))
waistTask.setKd(2.0 * np.sqrt(kp_waist) * np.ones(6))
waistMask = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
waistTask.setMask(waistMask)
waistTask.useLocalFrame(True)
waist_ref = robot_tsid.framePosition(data, model.getFrameId('base_link'))
trajWaist = tsid.TrajectorySE3Constant("traj_waist", waist_ref)
invdyn.addMotionTask(waistTask, w_waist, level_waist, 0.0)
sampleWaist = trajWaist.computeNext()
waistTask.setReference(sampleWaist)'''


#Posture Task
q_goal = p[7:].copy()
postureTask_legs = tsid.TaskJointPosture("task-posture_legs", robot_tsid)
postureTask_legs.setKp(kp_posture)
postureTask_legs.setKd(2.0 * np.sqrt(kp_posture)) 
postureMask_legs = np.zeros(robot_tsid.nq-7)
postureMask_legs[2]=1.
postureMask_legs[3]=1.
postureMask_legs[5]=1.
postureMask_legs[8]=1.
postureMask_legs[9]=1.
postureMask_legs[11]=1.
postureTask_legs.setMask(postureMask_legs)
invdyn.addMotionTask(postureTask_legs, w_posture_legs, level_posture, 0.0)
trajPosture_legs = tsid.TrajectoryEuclidianConstant("traj_joint_legs", q_goal)
samplePosture_legs = trajPosture_legs.computeNext()

#Posture Task
q_goal = p[7:].copy()
postureTask_body = tsid.TaskJointPosture("task-posture-body", robot_tsid)
postureTask_body.setKp(kp_posture)
postureTask_body.setKd(2.0 * np.sqrt(kp_posture)) 
postureMask_body = np.ones(robot_tsid.nq-7)
postureMask_body[2]=0.
postureMask_body[3]=0.
postureMask_body[5]=0.
postureMask_body[8]=0.
postureMask_body[9]=0.
postureMask_body[11]=0.
postureTask_body.setMask(postureMask_body)
invdyn.addMotionTask(postureTask_body, w_posture_body, level_posture, 0.0)
trajPosture_body = tsid.TrajectoryEuclidianConstant("traj_joint_body", q_goal)
samplePosture_body = trajPosture_body.computeNext()

# CoM Task creation      
comTask = tsid.TaskComEquality("task-com", robot_tsid)
comTask.setKp(kp_com) 
comTask.setKd(2 * np.sqrt(kp_com))
#comMask = np.array([1., 1., 0.])
#comTask.setMask(comMask)
invdyn.addMotionTask(comTask, w_com, level_com, 0.0)
com_target[2] = com_target[2] - com_offset[2]
trajCom = tsid.TrajectoryEuclidianConstant("traj_com", com_target)
sampleCom = trajCom.computeNext()



# Left Ankle Task creation    
leftFootTask = tsid.TaskSE3Equality("task-LeftFoot", robot_tsid, 'joint_left_foot_FE')
leftFootTask.setKp(kp_foot * np.ones(6))
leftFootTask.setKd(2.0 * np.sqrt(kp_foot) * np.ones(6))
leftFootMask = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
leftFootTask.setMask(leftFootMask)
leftFootTask.useLocalFrame(False)
left_foot_ref = data.oMf[frame_ids[1]].copy()
left_foot_ref.translation = np.array([LLM[0,1], -LLM[0,0], LLM[0,2]])
trajLeftFoot = tsid.TrajectorySE3Constant("traj-LeftFoot", left_foot_ref)
invdyn.addMotionTask(leftFootTask, w_foot, level_foot, 0.0)
sampleLeftFoot = trajLeftFoot.computeNext()
leftFootTask.setReference(sampleLeftFoot)

contacts = 2*[None]

contact_ref_RF = data.oMf[frame_ids[0]].copy()
f_n = data.oMf[frame_ids[0]].rotation[-1,:]  ### normal in local frame
contact_RF = tsid.Contact6d('joint_right_foot_FE', robot_tsid, 'joint_right_foot_FE', right_contact_point, f_n, mu, fMin, fMax)
contact_RF.setKp(kp_contact*np.ones(6))
contact_RF.setKd(2.0*np.sqrt(kp_contact)*np.ones(6))
contact_RF.setReference(contact_ref_RF)
invdyn.addRigidContact(contact_RF,w_forceRef, 1., level_contact)
contacts[0] = contact_RF

contact_ref_LF = data.oMf[frame_ids[1]].copy()
'''f_n = data.oMf[frame_ids[1]].rotation[-1,:]  ### normal in local frame
contact_LF = tsid.Contact6d('joint_left_foot_FE', robot_tsid, 'joint_left_foot_FE', left_contact_point, f_n, mu, fMin, fMax)
contact_LF.setKp(kp_contact*np.ones(6))
contact_LF.setKd(2.0*np.sqrt(kp_contact)*np.ones(6))
contact_LF.setReference(contact_ref_LF)
invdyn.addRigidContact(contact_LF,w_forceRef, 1., level_contact)
contacts[1] = contact_LF'''


# Solver initialization
solver = tsid.SolverHQuadProgFast("qp solver")
solver.resize(invdyn.nVar, invdyn.nEq, invdyn.nIn) 

HQPData = invdyn.computeProblemData(t, p, v)
HQPData.print_all()

sol = solver.solve(HQPData)
if(sol.status!=0):
        print ("QP problem could not be solved! Error code:", sol.status, 'WRONG INIT')


exo_display.robot.display(p)

i = k

### ALlocation
left_hip_computed[i] = p[7+2]
left_knee_computed[i] = p[7+3]
left_ankle_computed[i] = p[7+5]
right_hip_computed[i] =  p[7+8]
right_knee_computed[i] = p[7+9]
right_ankle_computed[i] = p[7+11]

left_foot_computed[:, i] = contact_ref_LF.translation.copy()
right_foot_computed[:, i] = contact_ref_RF.translation.copy()
left_foot_contact[i] = int(invdyn.checkContact('joint_left_foot_FE', sol))
right_foot_contact[i] = int(invdyn.checkContact('joint_right_foot_FE', sol))

com_experimental[:,i] = com_target.copy()
com_computed[:,i] = com_pos[:].copy()

#embed()


while (i<len(marche.time)-10) and (i_hs_left<3):
    time_start = time.time()    

    #break

    x_com = (1/4) * (RASIS[i,1] + LASIS[i,1] + RPSIS[i,1] + LPSIS[i,1])
    y_com = -(1/4) * (RASIS[i,0] + LASIS[i,0] + RPSIS[i,0] + LPSIS[i,0])
    z_com = (1/4) * (RASIS[i,2] + LASIS[i,2] + RPSIS[i,2] + LPSIS[i,2]) - com_offset[2]
    com_target = np.array([x_com, y_com, z_com])
    sampleCom.value(com_target)
    comTask.setReference(sampleCom)
    
    q_goal=np.array([0., 0., left_hip_FE[i],left_knee_FE[i],0.,left_ankle_FE[i],0.,0.,right_hip_FE[i],right_knee_FE[i],0.,right_ankle_FE[i]])
    samplePosture_legs.value(q_goal)
    postureTask_legs.setReference(samplePosture_legs)

    samplePosture_body.value(q_goal)
    postureTask_body.setReference(samplePosture_body)


    '''sampleWaist = trajWaist.computeNext()
    waistTask.setReference(sampleWaist)'''

    

    if (not invdyn.checkContact('joint_right_foot_FE', sol)):
        #right_foot_ref = sampleRightFoot.value().copy()
        #right_foot_ref[:3] = np.array([RLM[i,1], -RLM[i,0], RLM[i,2]])
        #sampleRightFoot.value = right_foot_ref

        right_foot_ref = data.oMf[frame_ids[0]].copy()
        right_foot_ref.translation = np.array([RLM[i,1], -RLM[i,0], RLM[i,2]])
        trajRightFoot = tsid.TrajectorySE3Constant("traj-RightFoot", right_foot_ref)
        sampleRightFoot = trajRightFoot.computeNext()

        rightFootTask.setReference(sampleRightFoot)


    if (not invdyn.checkContact('joint_left_foot_FE', sol)):
        #left_foot_ref = sampleLeftFoot.value().copy()
        #left_foot_ref[:3] = np.array([LLM[i,1], -LLM[i,0], LLM[i,2]])
        #sampleLeftFoot.value = left_foot_ref

        left_foot_ref = data.oMf[frame_ids[1]].copy()
        left_foot_ref.translation = np.array([LLM[i,1], -LLM[i,0], LLM[i,2]])
        trajLeftFoot = tsid.TrajectorySE3Constant("traj-LeftFoot", left_foot_ref)
        sampleLeftFoot = trajLeftFoot.computeNext()

        leftFootTask.setReference(sampleLeftFoot)




    HQPData = invdyn.computeProblemData(t, p, v)
    data = invdyn.data()

    sol = solver.solve(HQPData)
    if(sol.status!=0):
        print ("QP problem could not be solved! Error code:", sol.status)
        break
    
    tau = invdyn.getActuatorForces(sol)
    dv = invdyn.getAccelerations(sol)

        
    v_mean = v + 0.5*dt*dv
    v += dt*dv
    p = pin.integrate(model, p, dt*v_mean)

    robot_tsid.computeAllTerms(data, p, v)

    if i%DISPLAY_N == 0: exo_display.robot.display(p)
    
    com_pos = robot_tsid.com(data)

    contact_ref_RF = data.oMf[frame_ids[0]]
    f_n_right = data.oMf[frame_ids[0]].rotation[-1,:]

    contact_ref_LF = data.oMf[frame_ids[1]]
    f_n_left = data.oMf[frame_ids[1]].rotation[-1,:]

    if(t<t_to_right[i_to_right]<=t+dt):
        print(('Toes Off Right'))
        invdyn.removeRigidContact('joint_right_foot_FE', 0.0)
        i_to_right += 1

        rightFootTask = tsid.TaskSE3Equality("task-RightFoot", robot_tsid, 'joint_right_foot_FE')
        rightFootTask.setKp(kp_foot * np.ones(6))
        rightFootTask.setKd(2.0 * np.sqrt(kp_foot) * np.ones(6))
        rightFootMask = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        rightFootTask.setMask(rightFootMask)
        rightFootTask.useLocalFrame(False)
        right_foot_ref = data.oMf[frame_ids[0]].copy()
        right_foot_ref.translation = np.array([RLM[0,1], -RLM[0,0], RLM[0,2]])
        trajRightFoot = tsid.TrajectorySE3Constant("traj_RightFoot", right_foot_ref)
        invdyn.addMotionTask(rightFootTask, w_foot, level_foot, 0.0)
        sampleRightFoot = trajRightFoot.computeNext()
        rightFootTask.setReference(sampleRightFoot)

    if(t<t_to_left[i_to_left]<=t+dt):
        print(('Toes Off Left'))
        invdyn.removeRigidContact('joint_left_foot_FE', 0.0)
        i_to_left += 1

        leftFootTask = tsid.TaskSE3Equality("task-LeftFoot", robot_tsid, 'joint_left_foot_FE')
        leftFootTask.setKp(kp_foot * np.ones(6))
        leftFootTask.setKd(2.0 * np.sqrt(kp_foot) * np.ones(6))
        leftFootMask = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        leftFootTask.setMask(leftFootMask)
        leftFootTask.useLocalFrame(False)
        left_foot_ref = data.oMf[frame_ids[1]].copy()
        left_foot_ref.translation = np.array([LLM[0,1], -LLM[0,0], LLM[0,2]])
        trajLeftFoot = tsid.TrajectorySE3Constant("traj-LeftFoot", left_foot_ref)
        invdyn.addMotionTask(leftFootTask, w_foot, level_foot, 0.0)
        sampleLeftFoot = trajLeftFoot.computeNext()
        leftFootTask.setReference(sampleLeftFoot)



    if (t<t_hs_right[i_hs_right]<=t+dt):
        print(('Heel Strike Right'))
        contact_RF = tsid.Contact6d('joint_right_foot_FE', robot_tsid, 'joint_right_foot_FE', right_contact_point, f_n_right, mu, fMin, fMax)
        contact_RF.setKp(kp_contact*np.ones(6))
        contact_RF.setKd(2.0*np.sqrt(kp_contact)*np.ones(6))
        contact_RF.setReference(contact_ref_RF)
        invdyn.addRigidContact(contact_RF,w_forceRef, 1., level_contact)
        contacts[0] = contact_RF
        i_hs_right += 1
        invdyn.removeTask("task-RightFoot", 0.0)

    if (t<t_hs_left[i_hs_left]<=t+dt):
        print(('Heel Strike Left'))
        contact_LF = tsid.Contact6d('joint_left_foot_FE', robot_tsid, 'joint_left_foot_FE', left_contact_point, f_n_left, mu, fMin, fMax)
        contact_LF.setKp(kp_contact*np.ones(6))
        contact_LF.setKd(2.0*np.sqrt(kp_contact)*np.ones(6))
        contact_LF.setReference(contact_ref_LF)
        invdyn.addRigidContact(contact_LF,w_forceRef, 1., level_contact)
        contacts[1] = contact_LF
        i_hs_left += 1
        invdyn.removeTask("task-LeftFoot", 0.0)

    '''

    contact_ref_RF = data.oMf[frame_ids[0]]
    #f_n_right = data.oMf[frame_ids[0]].rotation[-1,:]
    f_n_right = np.array([0., 0., 1.])

    contact_ref_LF = data.oMf[frame_ids[1]]
    #f_n_left = data.oMf[frame_ids[1]].rotation[-1,:]
    f_n_left = np.array([0., 0., 1.])

    contact_len = abs(contact_ref_LF.translation[0]-contact_ref_RF.translation[0])


    if(invdyn.checkContact('joint_right_foot_FE', sol)) and (invdyn.checkContact('joint_left_foot_FE', sol)) and (com_pos[0]> (contact_ref_RF.translation[0]+leg_change_ratio*contact_len)) and ((contact_ref_LF.translation[0]>contact_ref_RF.translation[0])):
        print(('Toes Off Right'))
        invdyn.removeRigidContact('joint_right_foot_FE', 0.0)

        rightFootTask = tsid.TaskSE3Equality("task-RLM", robot_tsid, 'joint_right_foot_FE')
        rightFootTask.setKp(kp_foot * np.array([1., 1., 10., 1., 1., 1.]))
        rightFootTask.setKd(2.0 * np.sqrt(kp_foot) * np.array([1., 1., 10., 1., 1., 1.]))
        rightFootMask = np.array([1.0, 0.0, 1.0, 0.0, 0.0, 0.0])
        rightFootTask.setMask(rightFootMask)
        rightFootTask.useLocalFrame(False)
        #right_foot_ref = robot_tsid.framePosition(data, model.getFrameId('joint_right_foot_FE'))
        right_foot_ref = data.oMf[frame_ids[0]]
        right_foot_ref.translation = np.array([RLM[i,1], -RLM[i,0], RLM[i,2]])
        trajRightFoot = tsid.TrajectorySE3Constant("traj_RLM", right_foot_ref)
        invdyn.addMotionTask(rightFootTask, w_foot, level_foot, 0.0)
        sampleRightFoot = trajRightFoot.computeNext()
        rightFootTask.setReference(sampleRightFoot)

    if (invdyn.checkContact('joint_left_foot_FE', sol)) and (invdyn.checkContact('joint_right_foot_FE', sol)) and (com_pos[0]> (contact_ref_LF.translation[0]+leg_change_ratio*contact_len)) and ((contact_ref_LF.translation[0]<contact_ref_RF.translation[0])):  
        print(('Toes Off Left'))
        invdyn.removeRigidContact('joint_left_foot_FE', 0.0)

        leftFootTask = tsid.TaskSE3Equality("task-LLM", robot_tsid, 'joint_left_foot_FE')
        leftFootTask.setKp(kp_foot * np.array([1., 1., 10., 1., 1., 1.]))
        leftFootTask.setKd(2.0 * np.sqrt(kp_foot) * np.array([1., 1., 10., 1., 1., 1.]))
        leftFootMask = np.array([1.0, 0.0, 1.0, 0.0, 0.0, 0.0])
        leftFootTask.setMask(leftFootMask)
        leftFootTask.useLocalFrame(False)
        #left_foot_ref = robot_tsid.framePosition(data, model.getFrameId('joint_left_foot_FE'))
        left_foot_ref = data.oMf[frame_ids[1]]
        left_foot_ref.translation = np.array([LLM[i,1], -LLM[i,0], LLM[i,2]])
        trajLeftFoot = tsid.TrajectorySE3Constant("traj-LLM", left_foot_ref)
        invdyn.addMotionTask(leftFootTask, w_foot, level_foot, 0.0)
        sampleLeftFoot = trajLeftFoot.computeNext()
        leftFootTask.setReference(sampleLeftFoot)


    if (not invdyn.checkContact('joint_right_foot_FE', sol)) and (contact_ref_RF.translation[2]<=z_foot+security_height) and (contact_ref_LF.translation[0]<contact_ref_RF.translation[0]):
        print(('Heel Strike Right'))
        contact_RF = tsid.Contact6d('joint_right_foot_FE', robot_tsid, 'joint_right_foot_FE', right_contact_point, f_n_right, mu, fMin, fMax)
        contact_RF.setKp(kp_contact*np.ones(6))
        contact_RF.setKd(2.0*np.sqrt(kp_contact)*np.ones(6))
        contact_RF.setReference(contact_ref_RF)
        invdyn.addRigidContact(contact_RF,w_forceRef, 1., level_contact)
        contacts[0] = contact_RF

        invdyn.removeTask("task-RLM", 0.0)

    if (not invdyn.checkContact('joint_left_foot_FE', sol)) and (contact_ref_LF.translation[2]<=z_foot+security_height) and (contact_ref_LF.translation[0]>contact_ref_RF.translation[0]):
        print(('Heel Strike Left'))
        contact_LF = tsid.Contact6d('joint_left_foot_FE', robot_tsid, 'joint_left_foot_FE', left_contact_point, f_n_left, mu, fMin, fMax)
        contact_LF.setKp(kp_contact*np.ones(6))
        contact_LF.setKd(2.0*np.sqrt(kp_contact)*np.ones(6))
        contact_LF.setReference(contact_ref_LF)
        invdyn.addRigidContact(contact_LF,w_forceRef, 1., level_contact)
        contacts[1] = contact_LF

        invdyn.removeTask("task-LLM", 0.0)'''

        #embed()
    

    '''if(invdyn.checkContact('joint_right_foot_FE', sol)):
        #print('Right contact OK')
        contacts[0].setContactNormal(f_n_right)
        #contacts[0].setReference(right_foot_ref)
    if(invdyn.checkContact('joint_left_foot_FE', sol)):
        #print('Left contact OK')
        contacts[1].setContactNormal(f_n_left)
        #contacts[1].setReference(left_foot_ref)'''





    t += dt
    i += 1

    ### ALlocation
    left_hip_computed[i] = p[7+2]
    left_knee_computed[i] = p[7+3]
    left_ankle_computed[i] = p[7+5]
    right_hip_computed[i] =  p[7+8]
    right_knee_computed[i] = p[7+9]
    right_ankle_computed[i] = p[7+11]

    left_foot_computed[:, i] = contact_ref_LF.translation.copy()
    right_foot_computed[:, i] = contact_ref_RF.translation.copy()
    left_foot_contact[i] = int(invdyn.checkContact('joint_left_foot_FE', sol))
    right_foot_contact[i] = int(invdyn.checkContact('joint_right_foot_FE', sol))

    com_experimental[:,i] = com_target.copy()
    com_computed[:,i] = com_pos.copy()




    time_spent = time.time() - time_start
    if(time_spent < dt): time.sleep(dt-time_spent)

print("Nombre d'iterations: ", i, " sur ", len((marche.time)))

if(i==0): 
    i=1   
    trace_mode = 'markers'
else: 
    trace_mode = 'lines'


############################################ PLOT #######################################################################################"
if(PLOT_JOINT_POS):   

    #Côté gauche !
    plt.figure(0) 

    plt.subplot(221)
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(left_hip_computed[k:i-10]))
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(left_hip_FE[k:i-10]))
    plt.xlabel('Time (s)')
    plt.ylabel('Angles (degrees)')
    plt.title('Computed and experimental hip angular position - H.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(222)
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(left_knee_computed[k:i-10]))
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(left_knee_FE[k:i-10]))
    plt.xlabel('Time (s)')
    plt.ylabel('Angles (degrees)')
    plt.title('Computed and experimental knee angular position - H.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(223)
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(left_ankle_computed[k:i-10]))
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(left_ankle_FE[k:i-10]))
    plt.xlabel('Time (s)')
    plt.ylabel('Angles (degrees)')
    plt.title('Computed and experimental ankle angular position - H.')
    plt.legend(['Computed', 'Experimental'])


    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)


    #Côté droit !
    plt.figure(1) 

    plt.subplot(221)
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(right_hip_computed[k:i-10]))
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(right_hip_FE[k:i-10]))
    plt.xlabel('Time (s)')
    plt.ylabel('Angles (degrees)')
    plt.title('Computed and experimental hip angular position - H.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(222)
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(right_knee_computed[k:i-10]))
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(right_knee_FE[k:i-10]))
    plt.xlabel('Time (s)')
    plt.ylabel('Angles (degrees)')
    plt.title('Computed and experimental knee angular position - H.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(223)
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(right_ankle_computed[k:i-10]))
    plt.plot(data_kinematics.time[k:i-10],radians2degrees(right_ankle_FE[k:i-10]))
    plt.xlabel('Time (s)')
    plt.ylabel('Angles (degrees)')
    plt.title('Computed and experimental ankle angular position - H.')
    plt.legend(['Computed', 'Experimental'])


    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)


    
    #Pied gauche
    '''plt.figure(2) 

    plt.subplot(221)
    plt.plot(data_kinematics.time[k:i-10],left_foot_computed[0, :i-10])
    plt.plot(data_kinematics.time[k:i-10],LLM[k:i-10, 1])
    plt.xlabel('Time (s)')
    plt.ylabel('Position [mm]')
    plt.title('Computed and experimental left foot position, x-axis.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(222)
    plt.plot(data_kinematics.time[k:i-10],left_foot_computed[1, :i-10])
    plt.plot(data_kinematics.time[k:i-10],-LLM[k:i-10, 0])
    plt.xlabel('Time (s)')
    plt.ylabel('Position [mm]')
    plt.title('Computed and experimental left foot position, y-axis.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(223)
    plt.plot(data_kinematics.time[k:i-10],left_foot_computed[2, :i-10])
    plt.plot(data_kinematics.time[k:i-10],LLM[k:i-10, 2])
    plt.xlabel('Time (s)')
    plt.ylabel('Position [mm]')
    plt.title('Computed and experimental left foot position, z-axis.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(224)
    plt.plot(data_kinematics.time[k:i-10], left_foot_contact[k:i-10])
    plt.xlabel('Time (s)')
    plt.ylabel('Bool')
    plt.title('Left contact.')'''


    #Pied droit
    '''plt.figure(3) 

    plt.subplot(221)
    plt.plot(data_kinematics.time[k:i-10],right_foot_computed[0, :i-10])
    plt.plot(data_kinematics.time[k:i-10],RLM[k:i-10, 1])
    plt.xlabel('Time (s)')
    plt.ylabel('Position [mm]')
    plt.title('Computed and experimental right foot position, x-axis.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(222)
    plt.plot(data_kinematics.time[k:i-10],right_foot_computed[1, :i-10])
    plt.plot(data_kinematics.time[k:i-10],-RLM[k:i-10, 0])
    plt.xlabel('Time (s)')
    plt.ylabel('Position [mm]')
    plt.title('Computed and experimental right foot position, y-axis.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(223)
    plt.plot(data_kinematics.time[k:i-10],right_foot_computed[2, :i-10])
    plt.plot(data_kinematics.time[k:i-10],RLM[k:i-10, 2])
    plt.xlabel('Time (s)')
    plt.ylabel('Position [mm]')
    plt.title('Computed and experimental right foot position, z-axis.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(224)
    plt.plot(data_kinematics.time[k:i-10], right_foot_contact[k:i-10])
    plt.xlabel('Time (s)')
    plt.ylabel('Bool')
    plt.title('right contact.')


    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)'''
    


    #CoM
    plt.figure(4) 

    plt.subplot(221)
    plt.plot(data_kinematics.time[k:i-10],com_computed[0, k:i-10])
    plt.plot(data_kinematics.time[k:i-10],com_experimental[0, k:i-10])
    plt.xlabel('Time (s)')
    plt.ylabel('Position [mm]')
    plt.title('Computed and experimental CoM position, x-axis.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(222)
    plt.plot(data_kinematics.time[k:i-10],com_computed[1, k:i-10])
    plt.plot(data_kinematics.time[k:i-10],com_experimental[1, k:i-10])
    plt.xlabel('Time (s)')
    plt.ylabel('Position [mm]')
    plt.title('Computed and experimental CoM position, y-axis.')
    plt.legend(['Computed', 'Experimental'])

    plt.subplot(223)
    plt.plot(data_kinematics.time[k:i-10],com_computed[2, k:i-10])
    plt.plot(data_kinematics.time[k:i-10],com_experimental[2, k:i-10])
    plt.xlabel('Time (s)')
    plt.ylabel('Position [mm]')
    plt.title('Computed and experimental CoM position, z-axis.')
    plt.legend(['Computed', 'Experimental'])




    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)
    
    
    
    
    
    
    
    
    plt.show()



if(PLOT_WALK):
    trace_right_foot_exp = go.Scatter3d(x=RLM[k:i,1], y=-RLM[k:i,0], z=RLM[k:i,2],
                                    mode=trace_mode, name= 'Rigt Foot - experimental')
    trace_right_foot_comp = go.Scatter3d(x=right_foot_computed[0,k:i], y=right_foot_computed[1,k:i], z=right_foot_computed[2,k:i],
                                    mode=trace_mode, name= 'Rigt Foot - computed')

    trace_left_foot_exp = go.Scatter3d(x=LLM[k:i,1], y=-LLM[k:i,0], z=LLM[k:i,2],
                                    mode=trace_mode, name= 'Left Foot - experimental')
    trace_left_foot_comp = go.Scatter3d(x=left_foot_computed[0,k:i], y=left_foot_computed[1,k:i], z=left_foot_computed[2,k:i],
                                    mode=trace_mode, name= 'Left Foot - computed')

    trace_com_exp = go.Scatter3d(x=com_experimental[0,k:i], y=com_experimental[1,k:i], z=com_experimental[2,k:i],
                                    mode=trace_mode, name= 'CoM - experimental')
    trace_com_comp = go.Scatter3d(x=com_computed[0,k:i], y=com_computed[1,k:i], z=com_computed[2,k:i],
                                    mode=trace_mode, name= 'CoM - computed')
                                    

    x_start = [RLM[0,1], LLM[0,1], com_experimental[0,0]]
    y_start = [-RLM[0,0], -LLM[0,0], com_experimental[1,0]]
    z_start = [RLM[0,2], LLM[0,2], com_experimental[2,0]]

    '''x_hs = []
    y_hs = []
    z_hs = []
    for iter in (hs_right_list):
        x_hs.append(RLM[iter, 1])
        y_hs.append(-RLM[iter, 0])
        z_hs.append(RLM[iter, 2])
    for iter in (hs_left_list):
        x_hs.append(LLM[iter, 1])
        y_hs.append(-LLM[iter, 0])
        z_hs.append(LLM[iter, 2])

    x_to = []
    y_to = []
    z_to = []
    for iter in (to_right_list):
        x_to.append(RLM[iter, 1])
        y_to.append(-RLM[iter, 0])
        z_to.append(RLM[iter, 2])
    for iter in (to_left_list):
        x_to.append(LLM[iter, 1])
        y_to.append(-LLM[iter, 0])
        z_to.append(LLM[iter, 2])'''




    trace_start = go.Scatter3d(x=x_start, y=y_start, z=z_start,
                                    mode='markers', name= 'Start')

    '''trace_heel_strike = go.Scatter3d(x=x_hs, y=y_hs, z=z_hs,
                                    mode='markers', name= 'Heel Strike')

    trace_toes_off = go.Scatter3d(x=x_to, y=y_to, z=z_to,
                                    mode='markers', name= 'Toes Off')'''


    name = 'eye = (x:0., y:2.5, z:0.)'
    camera = dict(
        up=dict(x=0, y=0., z=1.),
        eye=dict(x=0, y=-3, z=0.3),
        center=dict(x=0, y=0., z=-0.5)
    )


    fig = go.Figure(data=[trace_right_foot_exp, trace_right_foot_comp, trace_left_foot_exp, trace_left_foot_comp, trace_com_exp, trace_com_comp, trace_start])
    fig.update_layout(scene_aspectmode='manual',
                    scene_aspectratio=dict(x=2, y=1, z=1), scene_camera=camera, legend=dict(font=dict(size=20), yanchor="top",
        y=0.8,
        xanchor="right",
        x=0.9)) #, xaxis = dict(tickfont = dict(size=20)), yaxis = dict(tickfont = dict(size=20)))
    fig.update_yaxes(tickfont=dict(size=50))
    fig.update_traces(marker=dict(size=3))
    fig.show()

print('R2 left hip = ', r2_score(left_hip_FE[k:i], left_hip_computed[k:i]))
print('R2 right hip = ', r2_score(right_hip_FE[k:i], right_hip_computed[k:i]))
print('R2 left knee = ', r2_score(left_knee_FE[k:i], left_knee_computed[k:i]))
print('R2 right knee = ', r2_score(right_knee_FE[k:i], right_knee_computed[k:i]))
print('R2 left ankle = ', r2_score(left_ankle_FE[k:i], left_ankle_computed[k:i]))
print('R2 right ankle = ', r2_score(right_ankle_FE[k:i], right_ankle_computed[k:i]))

print('R2 right foot = ', r2_score(np.array([RLM[k:i,1], -RLM[k:i,0], RLM[k:i,2]]), right_foot_computed[:,k:i]))
print('R2 left foot = ', r2_score(np.array([LLM[k:i,1], -LLM[k:i,0], LLM[k:i,2]]), left_foot_computed[:,k:i]))

print('R2 com = ', r2_score(com_experimental[:,k:i], com_computed[:,k:i]))