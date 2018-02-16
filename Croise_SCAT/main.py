#from utils import *
from croise_scatt import *
from scipy.io.wavfile import read
import sys
import time
import glob
import cPickle
import csv
from pylab import *
from sklearn.decomposition import PCA



window_size = 2**14
J1,Q1= 5 ,8
J2,Q2 = 4,1
family_names = ['Morlet','Paul','Gammatone']
family_params = [6,2,[6.,.5]]




#f=sort(glob.glob('../Scattering-MFCWT/wav/'+names_wav[i]+'.wav'))
#_,signal = read(data_files)


def transform(x,window_size,family_names,family_params,J1,Q1,J2,Q2):
	c = PCA(1)
	#GENERATE FILTER BANK
	filter_bank1,filters_1_fft = get_filter_banks(window_size,J1,Q1,family_names,family_params)
	filter_bank2,filters_2_fft = get_filter_banks(window_size,J2,Q2,family_names,family_params)
	#INIT LISTS
	L1,L2,S1,S2,V,P =[],[],[],[],[],[]
	start_indices = xrange(0,len(x)-window_size,window_size/2)
	print "NUMBER OF WINDOWS",len(start_indices)
	for w in start_indices:
		print w
		x_window = x[w:w+window_size].astype('float32')*hanning(window_size)
		x_window /= norm(x_window)
		s1,s2,l1,l2 =  scattering_3d(x_window,family_names,filter_bank1,filter_bank2)	
		L1 = asarray(l1)
		L2 = asarray(l2)
		S1.append(s1)
                S2.append(s2)
		VV = []
		PP = []
		for kk in xrange(L1.shape[0]):#accross families
			print kk
			v_local = []
			p_local = []
			for ii in xrange(J2*Q2):
				p_local.append(c.fit_transform(L2[-1][kk,ii*J1*Q1:(ii+1)*J1*Q1,:].T).mean())
				v_local.append(c.explained_variance_[0])
			VV.append(asarray(v_local))
			PP.append(asarray(p_local))
		V.append(VV)
		P.append(PP)
	S1=asarray(S1)
	S2=asarray(S2)
	L1=asarray(L1)
	L2=asarray(L2)
	V=asarray(V)
	P=asarray(P)
	print shape(S1),shape(V)
	return S1,S2,L1,L2,V,P

data = transform(randn(2**17),window_size,family_names,family_params,J1,Q1,J2,Q2)

