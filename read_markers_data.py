import loadmat

########## FICHIERS DE MARCHES DE SIRINE ######################

def read_markers_data(path_markers,file):   
    markers_data= loadmat(path_events+file)
    return markers_data