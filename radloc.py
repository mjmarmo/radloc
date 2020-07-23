from scipy import linalg
from PIL import Image
import numpy as np
import sys
import math
import json



# Simple Kriging (Euclidean Covariance)
class SKI(object):
	def dimension_filter(self,r_filename, r_filetype=''):
		data = np.genfromtxt(fname=r_filename, delimiter=r_filetype)
		n = len(data)
		z = 79.89
		altitude = []
		if data.shape == (n,4):
			return data
		if data.shape == (n,3):
			for i in range(n):
				altitude = np.r_[altitude, [z]]
			new_data = np.column_stack((data[:,0:2],altitude,data[:,2]))
			return new_data



	def zero_filter(self,raw_data):

		raw = raw_data
		fx = []
		fy = []
		fz = []
		fr = []
		for i in range(len(raw)):
			if raw[i,0] and raw[i,1] != 0:
						fx.append(raw[i,0])
						fy.append(raw[i,1])
						fz.append(raw[i,2])
						fr.append(raw[i,3])															
		return np.column_stack((fx,fy,fz,fr))



	def reduction_filter(self,raw_data):
		filtered_data = np.vstack(set(map(tuple, raw_data)))
		return filtered_data


	def iRange(self,min, max, r):
		while min < max:
			yield min
			min += r
	
	def __init__(self, in_fname, delimiter=''):
		self.raw_data = self.dimension_filter(r_filename=in_fname, r_filetype=delimiter)
		self.Data = self.reduction_filter(self.zero_filter(self.raw_data)) 
		self.N = len(self.Data)
		self.Exposure = self.Data[:,3]
		self.GPS = self.Data[:,0:3]
		self.GPS_xmin = self.GPS[:,0].min()
		self.GPS_xmax = self.GPS[:,0].max()
		self.GPS_ymin = self.GPS[:,1].min()
		self.GPS_ymax = self.GPS[:,1].max()
		self.GPS_zmin = self.GPS[:,2].min()
		self.GPS_zmax = self.GPS[:,2].max()	
		# Hyperparameters 
		self.center = 1
		self.epsilon = 1
		self.rbf_power = -0.5
		self.gamma = 0.01
		self.metric = 0.5
		# filtered data file
		self.filtered_fname = in_fname + 'filtered'


	def CM_solve(self):
		def cov( i,j):
			covariance = pow( (self.GPS[i,0] - self.GPS[j,0])**2 + (self.GPS[i,1] - self.GPS[j,1])**2, self.metric)
			return covariance
		CoMatrix = np.matrix([cov(i,j) for i in range(self.N) for j in range(self.N)]).reshape(self.N,self.N)	
		COM = np.matrix(CoMatrix)


		return np.linalg.inv(COM)




	def Interpolate(self,A, x, y):
		def covector(x,y,i):
			norm = pow((self.GPS[i,0] - x)**2 + (self.GPS[i,1] - y)**2, self.metric)
			return norm

		
		CoVector = np.matrix([covector(x,y,i) for i in range(self.N)]).T	
		
		Weight = np.dot(A, COV)
		
		
		prediction = np.dot(Weight, self.Exposure)
		
		return float(prediction)
			
	def vrange(self,vector, inumber):
		if vector=='Lat':
			return np.linspace(self.GPS_xmin, self.GPS_xmax, inumber)				
		if vector=='Lon':
			return np.linspace(self.GPS_ymin, self.GPS_ymax, inumber)
		


	def vInterpolate(self,X,Y):
		A = self.CM_solve()
		voki = np.vectorize(self.Interpolate)
		return voki(A,X,Y,pred_variance=pred_variance)

	def iSKI(self,resolution):
		self.oki_inv = self.CoMatrix_inverse()
		for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
			for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
				print '{} {} {}'.format(x, y, self.OrdinaryKriging(x,y,pred_variance=pred_variance))

	def gen_data(self, resolution):
		

		if ftype == 'txt':
			self.oki_inv = self.CoMatrix_inverse()
			for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
				for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
					print '{} {} {}'.format(x, y, self.OrdinaryKriging(x,y,pred_variance=pred_variance))
		
		if ftype == 'json':		
			with open(out_fname + ".metadata.json", "w") as outf:	
				for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
					for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
						outf.write(json.dumps({"rad": self.OrdinaryKriging(x,y,pred_variance=variance), "lat" : x, "lon" : y, }))			
			

	
	def gen_fdata(self, f_fname):	


		np.savetxt(f_fname , self.GPS[:,0:2])

	


	def genjson(self, out_fname, resolution):
		self.oki_inv = self.CoMatrix_inverse()
		
		with open(out_fname + "metadata.json", "w") as outf:	
			for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
				for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
					outf.write(json.dumps({"rad": self.OrdinaryKriging(x,y,pred_variance=variance), "lat" : x, "lon" : y, }))			

	def gen_json(self, out_fname, resolution):
		self.oki_inv = self.CoMatrix_inverse()
		
		with open(out_fname + ".metadata.json", "w") as outf:	
			for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
				for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
					outf.write(json.dumps({"rad": self.OrdinaryKriging(x,y,pred_variance=variance), "lat" : x, "lon" : y, }))			

	




# Ordinary Kriging (Euclidean Covariance)
class OKI(object):
	def dimension_filter(self,r_filename, r_filetype=''):
		data = np.genfromtxt(fname=r_filename, delimiter=r_filetype)
		n = len(data)
		z = 79.89
		altitude = []
		if data.shape == (n,4):
			return data
		if data.shape == (n,3):
			for i in range(n):
				altitude = np.r_[altitude, [z]]
			new_data = np.column_stack((data[:,0:2],altitude,data[:,2]))
			return new_data



	def zero_filter(self,raw_data):

		raw = raw_data
		fx = []
		fy = []
		fz = []
		fr = []
		for i in range(len(raw)):
			if raw[i,0] and raw[i,1] != 0:
						fx.append(raw[i,0])
						fy.append(raw[i,1])
						fz.append(raw[i,2])
						fr.append(raw[i,3])															
		return np.column_stack((fx,fy,fz,fr))



	def reduction_filter(self,raw_data):
		filtered_data = np.vstack(set(map(tuple, raw_data)))
		return filtered_data


	def iRange(self,min, max, r):
		while min < max:
			yield min
			min += r
	
	def __init__(self, in_fname, delimiter=''):
		self.raw_data = self.dimension_filter(r_filename=in_fname, r_filetype=delimiter)
		self.Data = self.reduction_filter(self.zero_filter(self.raw_data)) 
		self.N = len(self.Data)
		self.Exposure = self.Data[:,3]
		self.GPS = self.Data[:,0:3]
		self.GPS_xmin = self.GPS[:,0].min()
		self.GPS_xmax = self.GPS[:,0].max()
		self.GPS_ymin = self.GPS[:,1].min()
		self.GPS_ymax = self.GPS[:,1].max()
		self.GPS_zmin = self.GPS[:,2].min()
		self.GPS_zmax = self.GPS[:,2].max()	
		# Hyperparameters 
		self.center = 1
		self.epsilon = 1
		self.rbf_power = -0.5
		self.gamma = 0.01
		self.metric = 0.5
		# filtered data file
		self.filtered_fname = in_fname + 'filtered'


	def CM_solve(self):
		def cov( i,j):
			covariance = pow( (self.GPS[i,0] - self.GPS[j,0])**2 + (self.GPS[i,1] - self.GPS[j,1])**2, self.metric)
			return covariance
		CoMatrix = np.matrix([cov(i,j) for i in range(self.N) for j in range(self.N)]).reshape(self.N,self.N)	
		CM_side = np.ones((self.N,1))
		CM_corner = np.zeros((1,1))
		CM_base = np.ones((1,self.N))
		CM_body = np.vstack((CoMatrix, CM_base))
		CM_left = np.vstack((CM_side,CM_corner))
		COM = np.column_stack((CM_body,CM_left))


		return np.linalg.inv(COM)




	def Interpolate(self,A, x, y, pred_variance=False):
		def covector(x,y,i):
			norm = pow((self.GPS[i,0] - x)**2 + (self.GPS[i,1] - y)**2, self.metric)
			return norm

		
		CoVector = np.matrix([covector(x,y,i) for i in range(self.N)]).T	
		COV = np.vstack((CoVector,1))
		Weight = np.dot(A, COV)
		weight_1 = np.delete(Weight, (-1), axis=0)
		lagrange_multiplier = Weight[-1]
		prediction = np.dot(weight_1.T, self.Exposure)
		variance = np.dot(Weight.T, COV)
		if pred_variance==True:
			return float(variance)
		if pred_variance==False:	
			return np.abs((float(prediction)))
			
			
	def vrange(self,vector, inumber):
		if vector=='Lat':
			return np.linspace(self.GPS_xmin, self.GPS_xmax, inumber)				
		if vector=='Lon':
			return np.linspace(self.GPS_ymin, self.GPS_ymax, inumber)
		


	def vInterpolate(self,X,Y,pred_variance=False):
		A = self.CM_solve()
		voki = np.vectorize(self.Interpolate)
		return voki(A,X,Y,pred_variance=pred_variance)

	def iOKI(self,resolution,pred_variance=False):
		self.oki_inv = self.CoMatrix_inverse()
		for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
			for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
				print '{} {} {}'.format(x, y, self.OrdinaryKriging(x,y,pred_variance=pred_variance))

	def gen_data(self, resolution, pred_variance=False, ftype='txt'):
		

		if ftype == 'txt':
			self.oki_inv = self.CoMatrix_inverse()
			for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
				for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
					print '{} {} {}'.format(x, y, self.OrdinaryKriging(x,y,pred_variance=pred_variance))
		
		if ftype == 'json':		
			with open(out_fname + ".metadata.json", "w") as outf:	
				for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
					for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
						outf.write(json.dumps({"rad": self.OrdinaryKriging(x,y,pred_variance=variance), "lat" : x, "lon" : y, }))			
			

	
	def gen_fdata(self, f_fname):	


		np.savetxt(f_fname , self.GPS[:,0:2])

	


	def genjson(self, out_fname, resolution, variance=False):
		self.oki_inv = self.CoMatrix_inverse()
		
		with open(out_fname + "metadata.json", "w") as outf:	
			for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
				for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
					outf.write(json.dumps({"rad": self.OrdinaryKriging(x,y,pred_variance=variance), "lat" : x, "lon" : y, }))			

	def gen_json(self, out_fname, resolution, variance=False):
		self.oki_inv = self.CoMatrix_inverse()
		
		with open(out_fname + ".metadata.json", "w") as outf:	
			for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
				for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
					outf.write(json.dumps({"rad": self.OrdinaryKriging(x,y,pred_variance=variance), "lat" : x, "lon" : y, }))			

	




# Radial Basis Function Neural Network (IMQ Method)
class RBF(object):
	def dimension_filter(self,r_filename, r_filetype=''):
		data = np.genfromtxt(fname=r_filename, delimiter=r_filetype)
		n = len(data)
		z = 79.89
		altitude = []
		if data.shape == (n,4):
			return data
		if data.shape == (n,3):
			for i in range(n):
				altitude = np.r_[altitude, [z]]
			new_data = np.column_stack((data[:,0:2],altitude,data[:,2]))
			return new_data



	def zero_filter(self,raw_data):

		raw = raw_data
		fx = []
		fy = []
		fz = []
		fr = []
		for i in range(len(raw)):
			if raw[i,0] and raw[i,1] != 0:
						fx.append(raw[i,0])
						fy.append(raw[i,1])
						fz.append(raw[i,2])
						fr.append(raw[i,3])															
		return np.column_stack((fx,fy,fz,fr))
	


	def reduction_filter(self,raw_data):
		filtered_data = np.vstack(set(map(tuple, raw_data)))
		return filtered_data
	def iRange(self,min, max, r):
		while min < max:
			yield min
			min += r
	def __init__(self, filename, delimiter=''):
		self.raw_data = self.dimension_filter(r_filename=filename, r_filetype=delimiter)
		self.Data = self.reduction_filter(self.zero_filter(self.raw_data)) 
		self.N = len(self.Data)
		self.Exposure = self.Data[:,3]
		self.GPS = self.Data[:,0:3]
		self.GPS_xmin = self.GPS[:,0].min()
		self.GPS_xmax = self.GPS[:,0].max()
		self.GPS_ymin = self.GPS[:,1].min()
		self.GPS_ymax = self.GPS[:,1].max()
		self.GPS_zmin = self.GPS[:,2].min()
		self.GPS_zmax = self.GPS[:,2].max()	
		# Hyperparameters 
		self.center = 1
		self.epsilon = 1000
		self.rbf_power = -2
		self.gamma = 10
		self.length = 1
		
		
	def basis(self,i,j):
		norm = (self.GPS[i,0]-self.GPS[j,0])**2 + (self.GPS[i,1]-self.GPS[j,1])**2
		mq = self.center**2 + (self.epsilon**2)*norm
		return pow(mq, self.rbf_power)
		
	def bvect(self,x,y,i):
		norm = (self.GPS[i,0]-x)**2 + (self.GPS[i,1]-y)**2
		mq = self.center**2 + (self.epsilon**2)*norm
		return pow(mq, self.rbf_power)
		
		
		
	def CM_solve(self):
		gmq = [self.basis(i,j) for i in range(self.N) for j in range(self.N)]
		gmqmatrix = np.matrix(gmq).reshape(self.N,self.N)
		gmqi = np.linalg.inv(gmqmatrix)
		#lower = np.tril(gmqmatrix)
		#weight = linalg.solve(lower,self.Exposure, check_finite=False)
		weight = np.dot(gmqi, self.Exposure)
		return weight
		
	def Interpolate(self, x, y):
		
		gmqVector = np.matrix([self.bvect(x,y,i) for i in range(self.N)])
	
		
		
		interp = float(np.dot(np.matrix(A), gmqVector.T))
		
		return interp*self.gamma
		
	def vrange(self,vector, inumber):
		if vector=='Lat':
			return np.linspace(self.GPS_xmin, self.GPS_xmax, inumber)				
		if vector=='Lon':
			return np.linspace(self.GPS_ymin, self.GPS_ymax, inumber)
			
	def vInterpolate(self, X, Y):
		A = self.CM_solve()
		voki = np.vectorize(self.Interpolate)
		return voki(A,X,Y)

	def iRBF(self,resolution):
		for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
			for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
				print '{} {} {}'.format(x, y, self.rbf_interpolation(x,y))
	
		
	def GID(self,outfile, resolution):	
		import sys
		iData = outfile
		i_data = open(iData, 'w')
		sys.stdout = i_data
	
		self.irbf(resolution=resolution)
		
		i_data.close()
		sys.stdout = sys.__stdout__
	
	




# Inverse Distance Weighting (Shepard's Method, p=2)
class IDW(object):
	def dimension_filter(self,r_filename, r_filetype=''):
		data = np.genfromtxt(fname=r_filename, delimiter=r_filetype)
		n = len(data)
		z = 79.89
		altitude = []
		if data.shape == (n,4):
			return data
		if data.shape == (n,3):
			for i in range(n):
				altitude = np.r_[altitude, [z]]
			new_data = np.column_stack((data[:,0:2],altitude,data[:,2]))
			return new_data



	def zero_filter(self,raw_data):

		raw = raw_data
		fx = []
		fy = []
		fz = []
		fr = []
		for i in range(len(raw)):
			if raw[i,0] and raw[i,1] != 0:
						fx.append(raw[i,0])
						fy.append(raw[i,1])
						fz.append(raw[i,2])
						fr.append(raw[i,3])															
		return np.column_stack((fx,fy,fz,fr))
	


	def reduction_filter(self,raw_data):
		filtered_data = np.vstack(set(map(tuple, raw_data)))
		return filtered_data


	def iRange(self,min, max, r):
		while min < max:
			yield min
			min += r
	def __init__(self, filename, delimiter=''):
		self.raw_data = self.dimension_filter(r_filename=filename, r_filetype=delimiter)
		self.Data = self.reduction_filter(self.zero_filter(self.raw_data)) 
		self.N = len(self.Data)
		self.Exposure = self.Data[:,3]
		self.GPS = self.Data[:,0:3]
		self.GPS_xmin = self.GPS[:,0].min()
		self.GPS_xmax = self.GPS[:,0].max()
		self.GPS_ymin = self.GPS[:,1].min()
		self.GPS_ymax = self.GPS[:,1].max()
		self.GPS_zmin = self.GPS[:,2].min()
		self.GPS_zmax = self.GPS[:,2].max()	
		# Hyperparameters 
		self.center = 1
		self.epsilon = 1
		self.rbf_power = -0.5
		self.gamma = 0.01
		self.metric = 0.5	
		self.IDW_power = 2

	def CM_solve(self):
		return True
	def Interpolate(self,A, x, y):
		if A==True:
			def distance_vector(x, y, i):
				distvector = (pow(((x - self.GPS[i,0])**2 + (y -  self.GPS[i,1])**2), self.metric))**self.IDW_power
				if distvector == 0:
			 		return 1
				else: 
			 		return distvector
	
			DistVector = [distance_vector(x,y,i) for i in range(self.N)]
			Weights = (1 / np.matrix(DistVector)).T
			return float(np.dot(self.Exposure, Weights)/Weights.sum())




	def Interpolate_3D(self,A, x, y, z):
		if A==True:
			def distance_vector(x, y, z, i):
				distvector = (pow(((x - self.GPS[i,0])**2 + (y -  self.GPS[i,1])**2 + (z -  self.GPS[i,2])**2), self.metric))**self.IDW_power
				if distvector == 0:
			 		return 1
				else: 
			 		return distvector
		
			DistVector = [distance_vector(x,y,z,i) for i in range(self.N)]
			Weights = (1 / np.matrix(DistVector)).T
			return float(np.dot(self.Exposure, Weights)/Weights.sum())
 
  

	def vInterpolate(self, X, Y):
		A = self.CM_solve()
		voki = np.vectorize(self.Interpolate)
		return voki(A,X,Y)

	def iIDW(self,resolution):

		for x in self.iRange(self.GPS_xmin, self.GPS_xmax, resolution):
			for y in self.iRange(self.GPS_ymin, self.GPS_ymax, resolution):									
				for z in self.iRange(self.GPS_zmin, self.GPS_zmax, 1):
					print '{} {} {} {}'.format(x, y, z, self.Interpolate_3D(A,x,y,z))

	def vrange(self,vector, inumber):
		if vector=='Lat':
			return np.linspace(self.GPS_xmin, self.GPS_xmax, inumber)				
		if vector=='Lon':
			return np.linspace(self.GPS_ymin, self.GPS_ymax, inumber)
			
	
	
	def GID(self,outfile, resolution):	
		import sys
		iData = outfile
		i_data = open(iData, 'w')
		sys.stdout = i_data
	
		self.iIDW(resolution=resolution)
		
		i_data.close()
		sys.stdout = sys.__stdout__











# Uses interpolation algorithm for image reconstruction
class GSImage(object):
	def __init__(self, fname, iMethod, MAX_X, MAX_Y, image_fname, MIN_CR=None, MAX_CR=None, delimiter=''):
		self.Server = 'http://localhost:8000/'
		self.fname = fname
		self.image_fname = image_fname
		if iMethod == 'InverseDistanceWeighting':
			MODE = IDW
		if iMethod == 'OrdinaryKriging':
			MODE = OKI
		if iMethod == 'SimpleKriging':
			MODE = SKI
		if iMethod == 'RadialBasisFunction':
			MODE = RBF
		self.MODE = MODE(self.fname, delimiter=delimiter)
		self.MIN_LAT = self.MODE.Data[:,0].min()
		self.MAX_LAT = self.MODE.Data[:,0].max()
		self.MIN_LON = self.MODE.Data[:,1].min()
		self.MAX_LON = self.MODE.Data[:,1].max()
		self.MAX_X = MAX_X
		self.MAX_Y = MAX_Y
		self.MIN_CR = MIN_CR
		self.MAX_CR = MAX_CR		
		self.time = self.MODE.Data[:,3]



	def centroid(self):
		delta_lat = np.abs(self.MAX_LAT) - np.abs(self.MIN_LAT)
		delta_lon = np.abs(self.MAX_LON) - np.abs(self.MIN_LON)	
		Clat = self.MAX_LAT - np.abs(delta_lat)/2
		Clon = self.MIN_LON + np.abs(delta_lon)/2
		return Clat, Clon
	
	
		
	def params(self, bounds):
		if bounds == 'SW':
			return self.MIN_LAT, self.MIN_LON
		if bounds == 'NE':
			return self.MAX_LAT, self.MAX_LON
		
	def GPS_to_pixel(self, lat, lon):
	    adj_lat = float(lat) - self.MIN_LAT
	    adj_lon = float(lon) - self.MIN_LON

	    delta_lat = self.MAX_LAT - self.MIN_LAT
	    delta_lon = self.MAX_LON - self.MIN_LON

	    # x is lon, y is lat
	    # 0,0 is MIN_LON, MAX_LAT

	    lon_frac = adj_lon/delta_lon
	    lat_frac = adj_lat/delta_lat

	    x = int(lon_frac*self.MAX_X)
	    y = int((1-lat_frac)*self.MAX_Y)

	    return x,y

	def pixel_to_GPS(self, x, y):
	    delta_lat = self.MAX_LAT - self.MIN_LAT
	    delta_lon = self.MAX_LON - self.MIN_LON

	    # x is lon, y is lat
	    # 0,0 is MIN_LON, MAX_LAT

	    x_frac = float(x)/self.MAX_X
	    y_frac = float(y)/self.MAX_Y

	    lon = self.MIN_LON + x_frac*delta_lon
	    lat = self.MAX_LAT - y_frac*delta_lat


	    calc_x, calc_y = self.GPS_to_pixel(lat, lon)

	    if abs(calc_x-x) > 1 or abs(calc_y-y) > 1:
	        print( "Mismatch: %s, %s => %s %s" % (
	            x,y, calc_x, calc_y))

	    return lat, lon



	


	def greyscale(self, ival):
	    grey = int(256*float(ival)/3000)
	    return grey, grey, grey

	def color(self, ival, buckets):
		#if 0 <= ival <= 40:
		#	return (0,255,0)
		#if 40 < ival <= 80:
		#	return (255,255,0)
		#if 80<= ival < 110:
		#	return (255,125,0)
		#if ival >= 120:
		#	return (255,0,0)
			
			
				
		colors = [(255, 0, 0),
		          (255, 91, 0),
		          (255, 127, 0),
		          (255, 171, 0),
		          (255, 208, 0),
		          (255, 240, 0),
		          (255, 255, 0),
		          (218, 255, 0),
		          (176, 255, 0),
		          (128, 255, 0),
		          (0, 255, 0),
		          (0, 255, 255),
		          (0, 240, 255),
		          (0, 213, 255),
		          (0, 171, 255),
		          (0, 127, 255),
		          (0, 86, 255),
		          (0, 0, 255),
		          ]

		for price, color in zip(buckets, colors):
		    if ival > price:
		        return color
		return colors[-1]


	


	def get_image(self):
		filtered_data = self.MODE.Data[:,0:2]
		rad = {}
		iMatrix = self.MODE.CM_solve()
		for x in range(self.MAX_X):
			for y in range(self.MAX_Y):
				lat, lon = self.pixel_to_GPS(x,y)
				lat = np.asarray(lat)
				lon = np.asarray(lon)
				if self.MODE != None:
	  				irad = self.MODE.Interpolate(iMatrix, lat, lon)				
				rad[x,y] = irad
		 
		
		
		all_irad_areas = [x for x in sorted(rad.values()) if x is not None]
		total_irad_area = len(all_irad_areas)	
		buckets = []
		divisions = 17.0
		stride = total_irad_area/ (divisions + 1)
		next_i = int(stride)
		error_i = stride - next_i
		#________________________________________________		
		for i, val in enumerate(all_irad_areas):
			buckets.append(val)
			delta_i = stride + error_i
			next_i += int(delta_i)
			error_i = delta_i - int(delta_i)
		buckets.reverse()

		# color region by radiation level
		I = Image.new('RGBA', (self.MAX_X, self.MAX_Y))
		IM = I.load()
		for x in range(self.MAX_X):
			for y in range(self.MAX_Y):
				IM[x,y] = self.color(rad[x,y], buckets)
		#________________________________________________
		# add filtered points of robot path
		#for lat, lon in filtered_data:
		#	x, y = self.GPS_to_pixel(lat,lon)
		#	if 0 <= x < self.MAX_X and 0 <= y < self.MAX_Y:
		#		IM[x,y] = (0,0,0)
	
		out_fname = self.image_fname
		I.save(out_fname, "PNG")
		

		

# Overlays reconstructed image over geospatial region
class GSMap(object):
   
		def __init__(self, fname, Centroid, SWbounds, NEbounds):
			self.centerLat, self.centerLon = Centroid
			self.swLat, self.swLon = SWbounds
			self.neLat, self.neLon = NEbounds
			self.fname = fname
			
		def __str__(self):
				centerLat = self.centerLat
				centerLon = self.centerLon
				swLat = self.swLat
				swLon =  self.swLon
				neLat = self.neLat
				neLon = self.neLon
				image_path = 'http://localhost:8000/' + self.fname + '.heatmap.png'
				markersCode = "\n".join(
				    [ """function USGSOverlay(bounds, image, map) {
									  this.bounds_ = bounds;
									  this.image_ = image;
									  this.map_ = map;
									  this.div_ = null;
									  this.setMap(map);
									}
									USGSOverlay.prototype.onAdd = function() {
									  var div = document.createElement('div');
									  div.style.borderStyle = 'none';
									  div.style.borderWidth = '0px';
									  div.style.position = 'absolute';
									  var img = document.createElement('img');
									  img.src = this.image_;
									  img.style.width = '100%';
									  img.style.height = '100%';
									  img.style.position = 'absolute';
										img.style.opacity = 0.3;
										img.style.filter = 'alpha (opacity=30)'
									  div.appendChild(img);
									  this.div_ = div;
									  var panes = this.getPanes();
									  panes.overlayLayer.appendChild(div);
									};
									USGSOverlay.prototype.draw = function() {
									  var overlayProjection = this.getProjection();
									  var sw = overlayProjection.fromLatLngToDivPixel(this.bounds_.getSouthWest());
									  var ne = overlayProjection.fromLatLngToDivPixel(this.bounds_.getNorthEast());
									  var div = this.div_;
									  div.style.left = sw.x + 'px';
									  div.style.top = ne.y + 'px';
									  div.style.width = (ne.x - sw.x) + 'px';
									  div.style.height = (sw.y - ne.y) + 'px';
									};
									USGSOverlay.prototype.onRemove = function() {
									  this.div_.parentNode.removeChild(this.div_);
									  this.div_ = null;
									};"""

				    ])
				return """
				    <script src="https://maps.googleapis.com/maps/api/js?v=3.exp&sensor=false"></script>
				    <div id="map-canvas" style="height: 100%; width: 100%"></div>
				    <script type="text/javascript">
								var overlay;
								USGSOverlay.prototype = new google.maps.OverlayView();
								
				        function show_map() {{
				           var map = new google.maps.Map(document.getElementById("map-canvas"), {{
				                zoom: 18,
				                center: new google.maps.LatLng({centerLat}, {centerLon}),
												mapTypeId: google.maps.MapTypeId.SATELLITE
				            }});
									  var swBound = new google.maps.LatLng({swLat}, {swLon});
									  var neBound = new google.maps.LatLng({neLat}, {neLon});
									  var bounds = new google.maps.LatLngBounds(swBound, neBound);
									  var srcImage = '{image_path}';
									  overlay = new USGSOverlay(bounds, srcImage, map);
										
				        }}
								
				    		{markersCode}
								
								google.maps.event.addDomListener(window, 'load', show_map);
				    
						</script>
				""".format(centerLat=centerLat, centerLon=centerLon, swLat=swLat, swLon=swLon, neLat=neLat, neLon=neLon, image_path=image_path,
				           markersCode=markersCode)