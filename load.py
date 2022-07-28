from pathlib import Path
from example_robot_data.robots_loader import RobotLoader
from os.path import join
from pinocchio.robot_wrapper import RobotWrapper



# Declaration des différents exosquelettes comme étant des robots pinocchio 
class Exo(RobotLoader):
    path = ''
    urdf_filename = 'test_urdf.urdf'
    urdf_subpath = ''
    model_path = '/home/dmsm/s.otmani/Documents/Exosquelettes/' #Path(__file__).parent
    free_flyer = False

class Exo2(RobotLoader):
    path = 'exo2'
    urdf_filename = 'exo2.urdf'
    urdf_subpath = 'urdf'
    model_path = '/home/dmsm/s.otmani/Documents/Exosquelettes/' #Path(__file__).parent
    free_flyer = True


class Exo_Attach(RobotLoader):
    path = 'Exo_attach.SLDASM'
    urdf_filename = 'Exo_attach.SLDASM.urdf'
    urdf_subpath = 'urdf'
    model_path = '/home/dmsm/s.otmani/Documents/Exosquelettes/' #Path(__file__).parent
    free_flyer = True

class Exo_Attach2(RobotLoader):
    path = 'Exo_attach2.SLDASM'
    urdf_filename = 'Exo_attach2.SLDASM.urdf'
    urdf_subpath = 'urdf'
    model_path = '/home/dmsm/s.otmani/Documents/Exosquelettes/' #Path(__file__).parent
    free_flyer = True



#Declaration de l'humain!
class Human(RobotLoader):
    path = 'human-gazebo-master'
    urdf_filename = 'humanSubject01_48dof (copie).urdf'
    urdf_subpath = 'humanSubject01'
    model_path = '/home/dmsm/s.otmani/Documents/Exosquelettes/' #Path(__file__).parent
    free_flyer = False

#Declaration de l'humain!
class HumanNew(RobotLoader):
    path = ''
    urdf_filename = 'URDF_Human.urdf'
    urdf_subpath = ''
    model_path = '/home/dmsm/s.otmani/Documents/Exosquelettes/' #Path(__file__).parent
    free_flyer = False

#Problème : comment avoir accès aux différentes articulations (names) --> fichier srdf?
#Problème : URDF de l'exo bug