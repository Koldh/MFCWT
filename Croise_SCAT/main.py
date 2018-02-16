#from utils import *
from croise_scatt import *
from scipy.io.wavfile import read
import sys
import time
import glob
import cPickle
import csv
from pylab import *

csv_file = open('ff1010bird_metadata.csv','rb')
label_csv = csv.reader(csv_file)
data_list = vstack(list(label_csv))
names_wav = data_list[1:,0]
labels    = data_list[1:,1]
names_wav = names_wav[int(sys.argv[1]):int(sys.argv[2])]
labels    = labels[int(sys.argv[1]):int(sys.argv[2])]



window_size = 2**16
J1,Q1= 5 ,8
J2,Q2 = 4,1
family_names = ['Morlet','Paul','Gammatone']
family_params = [6,2,[6.,.5]]


print ' GENERATE FILTER BANK'
filter_bank1,filters_1_fft = get_filter_banks(window_size,J1,Q1,family_names,family_params)
filter_bank2,filters_2_fft = get_filter_banks(window_size,J2,Q2,family_names,family_params)
print 'GET CROISE SCATTERING COEFFS'
for i in xrange((len(names_wav))):
	print 'FILE NB: ' + str(i)+ ' OUT OF: ' +str(len(names_wav))
	L2,L4,S1,S2,V1,V2,RISK_1,RISK_2,FAMILY_1,FAMILY_2 =[], [],[],[], [],[],[],[],[],[]
	data_files = sort(glob.glob('../Scattering-MFCWT/wav/'+names_wav[i]+'.wav'))
        Fs,x = read(data_files)
	cpt = 0
	for w in xrange(0,len(x)-window_size,window_size-20000):
		t = time.time()	
		x_window = x[w:w+window_size].astype('float32')
		x_window /= norm(x_window)
		S1_w,S2_w,L2_w,L4_w =  scattering_3d(x_window,family_names,filter_bank1,filter_bank2)
		print time.time()-t
		S1.append(S1_w)
                S2.append(S2_w)
		#L2.append(L2_w)
		#L4.append(L4_w)
		print 'WINDOW:  '+str(cpt) + '   DATA FILE   '+names_wav[i] + '  LABEL  '+ labels[i]
		cpt +=1
	features = [S1,S2,int(labels[i])]
	f = open('./SCATTERING-FEATURES-NEW/scat_features_'+names_wav[i]+'_label_'+labels[i]+'.pkl','wb')
	cPickle.dump(features,f)
	print './scattering_features_'+names_wav[i]
	f.close()








