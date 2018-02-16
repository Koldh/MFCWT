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
J1,Q1= 8 ,16
J2,Q2 = 8,6
family_names = ['Morlet','Paul','Gammatone']
family_params = [15,8,[6.,.5]]



f=sort(glob.glob('./*.wav'))
_,signal = read(f[0])


def transform(x,window_size,family_names,family_params,J1,Q1,J2,Q2):
	c = PCA(1)
	#GENERATE FILTER BANK
	filter_bank1,filters_1_fft = get_filter_banks(len(x),J1,Q1,family_names,family_params)
	filter_bank2,filters_2_fft = get_filter_banks(len(x),J2,Q2,family_names,family_params)
	#INIT LISTS
	L1,L2,S1,S2,V,P =[],[],[],[],[],[]
	start_indices = xrange(0,len(x)-window_size,window_size/2)
	print "NUMBER OF WINDOWS",len(start_indices)
#	s1,s2,l1,l2 =  scattering_3d(x,family_names,filter_bank1,filter_bank2)
	l1 = asarray([abs(ifft(filters_1*fft(x).reshape((1,-1)))) for filters_1 in filters_1_fft])
	l2 = asarray([concatenate([abs(ifft(fft(ll1,axis=1)*w.reshape((1,-1)))) for w in filters_2],axis=0) for ll1,filters_2 in zip(l1,filters_2_fft)])
	print shape(l1),shape(l2)
	for w in start_indices:
		print w
		S1.append(asarray([v.mean(1) for v in l1[:,:,w:w+window_size]]))
                S2.append(asarray([v.mean(1) for v in l2[:,:,w:w+window_size]]))
		VV = []
		PP = []
		for kk in xrange(l1.shape[0]):#accross families
			print kk
			v_local = []
			p_local = []
			for ii in xrange(J2*Q2):
				print ii
				p_local.append(abs(c.fit_transform(l2[kk,ii*J1*Q1:(ii+1)*J1*Q1,w:w+window_size].T)).mean())
				v_local.append(abs(c.explained_variance_[0]))
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
	print shape(S1),shape(S2),shape(L1),shape(L2),shape(V),shape(P)
	return S1,S2,l1,l2,V,P


sig = signal[:2**16]
data = transform(sig,2**11,[family_names[1]],[family_params[1]],J1,Q1,J2,Q2)
subplot(411)
plot(sig)
xlim([0,len(sig)])
subplot(412)
imshow(data[-1][:,0,:].T,aspect='auto')
subplot(413)
imshow(data[-2][:,0,:].T,aspect='auto')
subplot(414)
imshow(data[2][0],aspect='auto')
show()
