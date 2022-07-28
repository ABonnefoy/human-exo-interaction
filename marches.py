# Classe des marches

import numpy as np

class Marche:
	def __init__(self,name,time,left,right,time_moments,left_moments,right_moments):
		self.name=name
		self.time=time
		self.left=[left[0],left[1],left[2]]
		self.right=[right[0],right[1],right[2]]
		self.time_moment=time_moments
		self.moment_left= [left_moments[0],left_moments[1],left_moments[2]]
		self.moment_right= [right_moments[0],right_moments[1],right_moments[2]]
		self.toesoff_left=[0,0,0]
		self.toesoff_right=[0,0,0]
		self.heelstrike_left=[0,0,0]
		self.heelstrike_right=[0,0,0]

	def toesoff(self,time,array_toesoff_left,array_toesoff_right):
		self.toesoff_left = array_toesoff_left
		self.toesoff_right = array_toesoff_right

	def heelstrike(self,time,array_heelstrike_left,array_heelstrike_right):
		self.heelstrike_left= array_heelstrike_left
		self.heelstrike_right= array_heelstrike_right

	def numberOfCycles(self):
		return len(self.toesoff_right)

	def firstFoot(self):
		if self.heelstrike_right[0] < self.heelstrike_left[0]:
			return 'right'
		else: 
			return 'left'

	def printMarche(self):
		print("Name:%s \n Time:%s \n Right Kinematics:%s \n Left Kinematics: %s \n Right Moment : %s \nT Left Moment: %s \n ToesOff_left : %s \n ToesOff_right: %s \n HeelStrike_left: %s \n HeelStrike_right: %s" % (self.name, self.time, self.right,self.left,self.moment_right, self.moment_left, self.toesoff_left,self.toesoff_right,self.heelstrike_left,self.heelstrike_right) ) 


	def cycle(self,side):
		if (side=='right'):
			for i in range(len(self.time)):
				#print('hello')
				if round(self.heelstrike_right[0],3) == round(self.time[i],3):
					begin=i
				elif round(self.heelstrike_right[1],3) == round(self.time[i],3):
					end=i
			duration=round(self.heelstrike_right[1],3)-round(self.heelstrike_right[0],3)  
		else:
			for i in range(len(self.time)):
					if round(self.heelstrike_left[0],3) == round(self.time[i],3):
						begin=i
					elif round(self.heelstrike_left[1],3) == round(self.time[i],3):
						end=i
			duration=round(self.heelstrike_left[1],3)-round(self.heelstrike_left[0],3) 
		return [duration,begin,end]



