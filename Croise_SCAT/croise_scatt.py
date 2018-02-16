from pylab import *
import matplotlib 
import time
from scipy.interpolate import interp1d
from scipy.io.wavfile import read



matplotlib.rc('xtick', labelsize=35) 
matplotlib.rc('ytick', labelsize=30) 

class sparse_filter:
	def __init__(self,coeffs,support,mu):
		self.mu      = mu
		self.coeffs  = coeffs
		self.support = support

class sparse_vector:
	def __init__(self,coeffs,support):
		self.coeffs  = coeffs
		self.support = support 


def create_filters(N,J,Q,eps=0.001,wavelet_name='Morlet',wavelet_params=60):
	delta_w = 2*pi/N
        c       = sqrt(-2.0*log(eps))
	if wavelet_name == 'Morlet':
	        domain                  = linspace(-pi,pi,N)
	        j                       = arange(J*Q)
	        sigma                   = wavelet_params
	        lambda_                 = (.5+sigma)/(pi)*2**(j/float(Q)).reshape(-1,1)
	        filters                 = pi**(-0.25)*exp(-(domain*lambda_-sigma)**2/2)
	        filters[:,:N/2]         = 0
		filters_sparse  	= [None]*J*Q
  	        filters                 = fftshift(filters,axes=1)
		mus                     = zeros_like(lambda_)
		sigmas                  = zeros_like(lambda_)
  	        for i in xrange(shape(filters)[0]):
			#mus[i]		  = (sigma/(lambda_[i]))/(delta_w)
			#sigmas[i]	  = sigma*(1./(2*lambda_[i]**2))/delta_w
			#support    	  = [int(floor((mus[i]-sigmas[i]))),int(min(ceil((mus[i]+sigmas[i])),N/2))]
                        #filters[i,:]     /= norm(filters[i,:])
			plot(filters[i,:])
			#axvline(mus[i],color='r')
			#axvline(int(floor((mus[i]-sigmas[i]))),color='k')
                        #axvline(int(min(ceil((mus[i]+sigmas[i])),N/2)),color='k')	
			#coeffs 		  = filters[i,support[0]:support[1]]
	                #filters_sparse[i] = sparse_filter(coeffs,support,mus[i])
		show()
	elif wavelet_name == 'Gammatone':
        	alpha = wavelet_params[1]
        	m     = wavelet_params[0]
        	B = alpha
        	xi         = (pi-alpha)
        	r          = sqrt(1./2)
        	domain     = linspace(-pi,pi,N)
        	j          = arange(J*Q)
        	lambda_    = 2**(j/float(Q)).reshape(-1,1)
        	sigma      = (r**(2./m)*(1-r**(2./m))*m**2*(xi)**2)/2*(sqrt(1+B**2/((1-r**(2./m))**2*m**2*(xi)**2))-1)
        	filters   = (1j*domain*lambda_*(prod(arange(1,m))))/(sigma+1j*(domain*lambda_-xi))**m
        	filters[:,:N/2]         = 0
        	filters                 = fftshift(filters,axes=1)
                filters_sparse          = [None]*J*Q                
		mus                     = zeros_like(lambda_)
                sigmas                  = zeros_like(lambda_)
        	for i in xrange(shape(filters)[0]):
                        #mus[i]            = (1./lambda_[i])*xi/(delta_w)
                        #sigmas[i]         = sigma*(B)/(delta_w*lambda_[i])
                        #support           = [int(floor((mus[i]-sigmas[i]))),int(min(ceil((mus[i]+sigmas[i])),N/2))]
                        #filters[i,:] /= norm(filters[i,:])
                        plot(filters[i,:])
                        #coeffs            = filters[i,support[0]:support[1]]
                        #filters_sparse[i] = sparse_filter(coeffs,support,mus[i])
		show()

	elif wavelet_name ==  'Paul':
		m                 = wavelet_params
     	  	domain            = linspace(-pi,pi,N)
     	   	j                 = arange(J*Q)
      		sigma             = .8
		lambda_           = (((2*m+1)+sqrt(2*m+1))/(2*pi))*(2**(j/float(Q))).reshape(-1,1)
		filters           = (2**m)/(m*prod(sqrt(arange(1,2*m))))*((domain*lambda_)**m)*(exp(-domain*lambda_/sigma))
		filters[:,:N/2+1] 	= 0
		filters 		= fftshift(filters,axes=1)
                filters_sparse          = [None]*J*Q
                mus                     = zeros_like(lambda_)
                sigmas                  = zeros_like(lambda_)
                for i in xrange(shape(filters)[0]):
                        #mus[i]            = (2*m+1)/(2*lambda_[i]*delta_w)
                	#sigmas[i]         = sqrt(2*m+1)/(2*lambda_[i]*delta_w)
                        #support           = [int(floor((mus[i]-sigmas[i]))),int(min(ceil((mus[i]+sigmas[i])),N/2))]
                        #filters[i,:] /= norm(filters[i,:])
			plot(filters[i,:])
                        #coeffs            = filters[i,support[0]:support[1]]
                        #plot(filters[i,:])
                        #axvline(mus[i],color='r')
                        #axvline(int(floor((mus[i]-sigmas[i]))),color='k')
                        #axvline(int(min(ceil((mus[i]+sigmas[i])),N/2)),color='k')      
                        #filters_sparse[i] = sparse_filter(coeffs,support,mus[i])
		show()
        return filters#,filters_sparse




def get_filter_banks(N,J,Q,family_names,family_params):
	filter_bank1     = []
	filters_1_fft    = []
        for i in xrange(len(family_names)):
                print 'GEN FILTER BANK: ', family_names[i]
                filters_1 = create_filters(N,J,Q,0.001,family_names[i],family_params[i])
                #filter_bank1.append(bank1)
		filters_1_fft.append(filters_1)
	return filter_bank1,filters_1_fft


def get_scatt(sparse_input,N):
        N_IN = len(sparse_input)
        output = [None]*N_IN
        scalo = abs(create_representation(sparse_input,N))
        for i in xrange(N_IN):
                output[i]      = sparse_vector(fft(scalo[i])[:N/2],[0,N/2])
	return output,scalo.sum(1),scalo

def scattering_3d(x,family_names,filter_bank1,filter_bank2):
	LL2           = []
	L2 = []
	LL4           = []
	L4 = []
        first        = [sparse_vector(fft(x)[:len(x)/2],[0,len(x)/2])]
	S1,V1        = [],[]
	S2,V2        = [],[]
	## LAYER 1 FOR ALL FAMILIES
	for i in xrange(len(family_names)):
	        L1          = apply_filter_bank(first,filter_bank1[i])
		l2,s1,l2full    = get_scatt(L1,len(x))
		LL2.append(l2)
		L2.append(l2full)
		S1.append(s1)
	for i in xrange(len(family_names)):
		L3 = apply_filter_bank(LL2[i],filter_bank2[i])
		_,s2,L1 = get_scatt(L3,len(x))
		L4.append(L1)
	        S2.append(s2)
	return asarray(S1),asarray(S2),asarray(L2),asarray(L4)
	
	

def support_intersection(supp1,supp2):
	if(len(supp1)==0 or len(supp2)==0):
		return []
	elif(supp1[0]>=supp2[0]):
		if(supp1[0]>supp2[-1]):
			return []
		elif(supp1[-1]>=supp2[-1]):
			return [supp1[0],supp2[-1]]
		else:
			return [supp1[0],supp1[-1]]
	else:
		if(supp1[-1]<supp2[0]):
			return []
		elif(supp1[-1]>=supp2[-1]):
			return [supp2[0],supp2[-1]]
		else:
			return [supp2[0],supp1[-1]]



def apply_filter_bank(sparse_input,sparse_filter_bank):
	N_IN          = len(sparse_input)
	N_FILTER_BANK = len(sparse_filter_bank)
	output        = [None]*N_IN*N_FILTER_BANK
	for i in xrange(N_FILTER_BANK):
		for j in xrange(N_IN):
			if(sparse_input[j] != None):
				inter_support = support_intersection(sparse_input[j].support,sparse_filter_bank[i].support)
				if(len(inter_support)==0):
					output[i*N_IN+j]=sparse_vector([],[])
				else:
					supp   = inter_support[1]-inter_support[0]
                                        common_j = [inter_support[0]-sparse_input[j].support[0],(inter_support[0]-sparse_input[j].support[0]+supp)]
					common_i = [(inter_support[0]-sparse_filter_bank[i].support[0]),(inter_support[0]-sparse_filter_bank[i].support[0]+supp)]
					output[i*N_IN+j]=sparse_vector((sparse_input[j].coeffs[common_j[0]:common_j[1]]/(norm(sparse_input[j].coeffs[common_j[0]:common_j[1]])+0.0000000001))*sparse_filter_bank[i].coeffs[common_i[0]:common_i[1]],inter_support)
	return output






def create_phi(N,sigma,eps=0.001):
	w       = linspace(-pi,pi,N)
	phi     = exp(-w**2/(2*sigma**2))
	return fftshift(phi)



def extract_mean_scattering(sparse_input):
	N_IN     = len(sparse_input)
	features = zeros(N_IN)
	for i in xrange(N_IN):
		features[i]=sparse_input[i].coeffs[0]
	return features

def extract_dispersion_scattering(sparse_input):
        N_IN     = len(sparse_input)
        features = zeros(N_IN)
        for i in xrange(N_IN):
                features[i]=norm(sparse_input[i].coeffs[1:])
	return features


def create_representation(sparse_vector,N_in):
	N_IN        = len(sparse_vector)
	to_show     = zeros((N_IN,N_in),dtype='complex64')
        for i in xrange(N_IN):
#                print sparse_vector[i].support
		if(sparse_vector[i].support[0]==0):
                        to_show[i,:sparse_vector[i].support[1]]    = sparse_vector[i].coeffs
                        to_show[i,-sparse_vector[i].support[1]+1:] = conj(flipud(sparse_vector[i].coeffs[1:]))
		else:
	                to_show[i,sparse_vector[i].support[0]:sparse_vector[i].support[1]]    = sparse_vector[i].coeffs
#	                to_show[i,-sparse_vector[i].support[1]+1:-sparse_vector[i].support[0]+1] = conj(flipud(sparse_vector[i].coeffs))
        to_show=ifft(to_show)#[:,::2]
	return to_show



def extract_features(L):
	N_IN = len(L)
	S    = list()
	V    = list()
	for i in xrange(N_IN):
		if(L[i] != None):
			S.append(abs(L[i].coeffs[0]))
			V.append(2*norm(L[i].coeffs[1:]))
	return array(S),array(V)





#
