# Read data kinematics #

import txt2csv
import pandas as pd

########## FICHIERS DE MARCHES DE SIRINE ######################

def read_data(path_txt,path_csv):   
    txt= open(path_txt,"r")
    csv= open(path_csv,"w")
    toCSV=txt2csv.txt2csv(path_txt,path_csv)
    fichier_csv_bis = open(path_csv, "r")
    data = pd.read_csv(fichier_csv_bis, skiprows=10) # on ne regarde pas le header des data kinematics
    return data
