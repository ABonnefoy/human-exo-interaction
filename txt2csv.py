#Fonction qui transvase le contenu d'un fichier txt à un fichier .csv que l'on a cree au préalable
import csv

def txt2csv(path_txt, path_csv):
	in_txt = csv.reader(open(path_txt, "r"), delimiter = '\t')
	out_csv = csv.writer(open(path_csv, 'w'))
	out_csv.writerows(in_txt)

