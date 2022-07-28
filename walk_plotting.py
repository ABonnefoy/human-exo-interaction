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



######################################################################### READ DATA (kinematics, moments, events) ##########################################################################################

studied_walk=1

# PATH & URD
path = os.path.dirname(os.path.realpath(__file__))
urdf = path + '/urdf_augmented.urdf'

## Read data
data_kinematics= read_data_kinematics.read_data(path+"/Données_marche/fichiers txt/marche sans attelle "+str(studied_walk)+".txt", path + "/Données_marche/fichiers txt/marche sans attelle "+str(studied_walk)+".csv")
right=[data_kinematics.hip_flexion_r, data_kinematics.knee_angle_r, data_kinematics.ankle_angle_r]
left=[data_kinematics.hip_flexion_l, data_kinematics.knee_angle_l, data_kinematics.ankle_angle_l]

#Read Moments
data_moment=read_data_moments.read_data_moments(path + "/Données_marche/InverseDynamics/résultats sto & txt/marche sans attelle "+str(studied_walk)+".txt",path + "/Données_marche/InverseDynamics/résultats sto & txt/marche sans attelle "+str(studied_walk)+".csv")
right_moments=[data_moment.hip_flexion_r_moment, data_moment.knee_angle_r_moment, data_moment.ankle_angle_r_moment]
left_moments=[data_moment.hip_flexion_l_moment, data_moment.knee_angle_l_moment, data_moment.ankle_angle_l_moment]

#Walk construction
marche=marches.Marche('marche'+str(studied_walk),data_kinematics.time,left,right,data_moment.time,left_moments,right_moments)


#Read events
path_events= path + "/Données_marche/dataEvents/"
events= loadmat(path_events+'marche sans attelle '+str(studied_walk)+'Events.mat')
marche.toesoff(data_kinematics.time,events["Left_Foot_Off"],events["Right_Foot_Off"])
marche.heelstrike(data_kinematics.time,events["Left_Foot_Strike"],events["Right_Foot_Strike"])
#marche.printMarche() #print l'ensemble des données de la marche
#print('Heel strike', marche.heelstrike_right)
#Read FirstFoot
firstFoot=marche.firstFoot()

#delay = 255 * 0.005
delay = 0

t_hs_left = np.array(marche.heelstrike_left) - delay
t_hs_left = np.append(t_hs_left, np.zeros(1))
#print(t_hs_left)
i_hs_left = 0

t_hs_right = np.array(marche.heelstrike_right) - delay
t_hs_right = np.append(t_hs_right, np.zeros(1))
#print(t_hs_right)
i_hs_right = 0

t_to_left = np.array(marche.toesoff_left) - delay
t_to_left = np.append(t_to_left, np.zeros(1))
#print(t_to_left)
i_to_left = 0

t_to_right = np.array(marche.toesoff_right) - delay
t_to_right = np.concatenate((np.array([0.18]), t_to_right, np.zeros(1)))
#print(t_to_right)
i_to_right = 0

to_left_list = []
to_right_list = []
hs_left_list = []
hs_right_list = []
for i, t_event in enumerate(data_kinematics.time):
    print(t_event)
    if (t_hs_left[i_hs_left]<=t_event<t_hs_left[i_hs_left]+0.005):
        hs_left_list.append(i)
        i_hs_left += 1        
    if (t_hs_right[i_hs_right]<=t_event<t_hs_right[i_hs_right]+0.005):
        hs_right_list.append(i) 
        i_hs_right += 1
    if (t_to_left[i_to_left]<=t_event<t_to_left[i_to_left]+0.005):
        to_left_list.append(i)
        i_to_left += 1
    if (t_to_right[i_to_right]<=t_event<t_to_right[i_to_right]+0.005):
        to_right_list.append(i)
        i_to_right += 1





#Read markers
data_markers=loadmat(path + "/Données_marche/markers_data_patho.mat")
data_markers=data_markers['output']['marker_data']['Markers']
#print(data_markers.keys())

RLM = 1e-3 * np.array(data_markers['RLM']) # /!\ l'ordre du repère est -y, x, z
LLM = 1e-3 * np.array(data_markers['LLM']) # /!\ l'ordre du repère est -y, x, z

RGT = 1e-3 * np.array(data_markers['RGT']) # /!\ l'ordre du repère est -y, x, z
LGT = 1e-3 * np.array(data_markers['LGT']) # /!\ l'ordre du repère est -y, x, z


RASIS = 1e-3 * np.array(data_markers['RASIS']) # /!\ l'ordre du repère est -y, x, z
LASIS = 1e-3 * np.array(data_markers['LASIS']) # /!\ l'ordre du repère est -y, x, z

RPSIS = 1e-3 * np.array(data_markers['RPSIS']) # /!\ l'ordre du repère est -y, x, z
LPSIS = 1e-3 * np.array(data_markers['LPSIS']) # /!\ l'ordre du repère est -y, x, z


x_com = (1/4) * (RASIS[:,1] + LASIS[:,1] + RPSIS[:,1] + LPSIS[:,1])
y_com = -(1/4) * (RASIS[:,0] + LASIS[:,0] + RPSIS[:,0] + LPSIS[:,0])
z_com = (1/4) * (RASIS[:,2] + LASIS[:,2] + RPSIS[:,2] + LPSIS[:,2])



trace_right_foot = go.Scatter3d(x=RLM[:,1], y=-RLM[:,0], z=RLM[:,2],
                                   mode='lines', name= 'Right Foot')
trace_left_foot = go.Scatter3d(x=LLM[:,1], y=-LLM[:,0], z=LLM[:,2],
                                   mode='lines', name= 'Left Foot')
trace_com = go.Scatter3d(x=x_com, y=y_com, z=z_com,
                                   mode='lines', name= 'CoM')

x_start = [RLM[0,1], LLM[0,1], x_com[0]]
y_start = [-RLM[0,0], -LLM[0,0], y_com[0]]
z_start = [RLM[0,2], LLM[0,2], z_com[0]]

#x_start = [RLM[0,1], LLM[0,1]]
#y_start = [-RLM[0,0], -LLM[0,0]]
#z_start = [RLM[0,2], LLM[0,2]]

x_hs = []
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
    z_to.append(LLM[iter, 2])




trace_start = go.Scatter3d(x=x_start, y=y_start, z=z_start,
                                   mode='markers', name= 'Start')

trace_heel_strike = go.Scatter3d(x=x_hs, y=y_hs, z=z_hs,
                                   mode='markers', name= 'Heel Strike')

trace_toes_off = go.Scatter3d(x=x_to, y=y_to, z=z_to,
                                   mode='markers', name= 'Toes Off')

name = 'eye = (x:0., y:2.5, z:0.)'
camera = dict(
    up=dict(x=0, y=0., z=1.),
    eye=dict(x=0, y=-3, z=0.3),
    center=dict(x=0, y=0., z=-0.5)
)




#fig = go.Figure(data=[trace_right_foot, trace_left_foot, trace_start, trace_heel_strike, trace_toes_off])
fig = go.Figure(data=[trace_right_foot, trace_left_foot, trace_com, trace_start, trace_heel_strike, trace_toes_off])
fig.update_layout(scene_aspectmode='manual',
                  scene_aspectratio=dict(x=2, y=1, z=1), scene_camera=camera, legend=dict(font=dict(size=20), yanchor="top",
    y=0.8,
    xanchor="right",
    x=0.75)) #, xaxis = dict(tickfont = dict(size=20)), yaxis = dict(tickfont = dict(size=20)))
fig.update_yaxes(tickfont=dict(size=50))
fig.update_traces(marker=dict(size=3))
fig.show()


