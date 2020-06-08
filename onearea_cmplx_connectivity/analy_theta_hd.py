#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 15:50:32 2020

@author: shni2598
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 15:45:12 2020

@author: shni2598
"""


import load_data_dict
import brian2.numpy_ as np
from brian2.only import *
import post_analysis as psa
import connection as cn
import pickle
import matplotlib as mpl
mpl.use('Agg')
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft
#%%
lattice_ext = cn.coordination.makelattice(63, 62, [0, 0])

#%%
'''
ie_ratio = 3.375*np.arange(0.9, 1.21, 0.03)
tau_s_di = [3., 6.] # ms
new_delta_gk = np.arange(0, 8.1, 1) # nS
'''
#%%
#i = 10; j = 1; k = 8
sys_argv = int(sys.argv[1])
loop_num = sys_argv
'''
loop_num = -1
for i in range(len(ie_ratio)):
    for j in range(len(tau_s_di)):
        for k in range(len(new_delta_gk)):
            loop_num += 1
            if loop_num == sys_argv:
                print(i,j,k)
                break
        else: continue
        break
    else:continue
    break
'''
#i = 2; j = 0; k = 2
#data_num = i*2*9+j*9+k
#fileadr = '/import/headnode1/shni2598/brian2/data/onearea_chg_adapt_simp_2/'
with open('data%s'%(loop_num), 'rb') as file:
    spke1 = load_data_dict.data_onegroup(pickle.load(file))

spke1.t = spke1.t*0.1*ms
#%%
chg_adapt_loca = [0, 0]
chg_adapt_range = spke1.param['chg_adapt_range']
chg_adapt_neuron = cn.findnearbyneuron.findnearbyneuron(lattice_ext, chg_adapt_loca, chg_adapt_range, 62)

#%%
spkrate1, centre_ind1, jump_size1, jump_dist1 = psa.get_rate_centre_jumpdist(spke1, starttime=0*ms, endtime=8000*ms, binforrate=10*ms, sample_interval=1*ms, slide_interval = 1, jump_interval = 15, dt=0.1*ms, show_trajectory=False)
#%%
'''
plt.figure()
stat = plt.hist2d(centre_ind1[:1000,1], centre_ind1[:1000,0], bins=np.linspace(0,62,10))
#plt.matshow(stat[0])
plt.colorbar()
#%%
plt.matshow(stat[0]/len(centre_ind1[:1000]), extent=[0,62,62,0])
plt.title('residence time of pattern at\n different locations in the network (%)')
ax = plt.gca()
ax.xaxis.set_ticks_position('bottom')
#ax.tick_params(labeltop=False)
plt.colorbar()
#%%
fig, ax = plt.subplots(1,1)
psa.plot_traj(ax, centre_ind1[1000:])
ax.set_xlim([0,62])
ax.set_ylim([0,62])
ax.set_title('trajectory of pattern after adaptation change')
'''
#%%
'''
fftcoef = fft(centre_ind1[1000:,1])
magcoef = abs(fftcoef)
fig, [ax1, ax2] = plt.subplots(1,2)
ax1.plot(np.arange(len(magcoef))[1:300]/len(magcoef)*1000, magcoef[1:300])
ax2.plot(np.arange(len(magcoef))[1:int(len(magcoef)/2)]/len(magcoef)*1000, magcoef[1:int(len(magcoef)/2)])
ax2.set_xscale('log')
ax2.set_yscale('log')
ax1.set_xlabel('Hz'); ax2.set_xlabel('Hz')
fig.suptitle('fourier transform of the x coordinate of pattern trajactory')
#spkrate1 = psa.overlap_centreandspike(centre_ind1, spkrate1, show_trajectory = False)
#anititle='change adaptation in region near centre to %.1fnS at 1000 ms\ntau_s_di: %.1fms; ie_ratio: %.3f'%(new_delta_gk[k], tau_s_di[j], ie_ratio[i])
#ani1 = psa.show_pattern(spkrate1, spkrate2=0, area_num = 1, frames = 5500, bottom_up=1, top_down=1, stimu_onset = -1, start_time = 0, anititle=anititle)
'''
#%%

rate_chgadapt_neu = spkrate1.reshape(3969,-1)[chg_adapt_neuron]
poprate_rate_chgadapt_neu = rate_chgadapt_neu.sum(0)/len(chg_adapt_neuron)/0.01
#%%
'''
import scipy.io as sio
sio.savemat('adpt_rate%s.mat'%loop_num, {'data':{'rate':poprate_rate_chgadapt_neu, 'param':{'new_delta_gk':new_delta_gk[k],'tau_s_di':tau_s_di[j],'ie_ratio':ie_ratio[i]}}})
'''
#%%
import matlab.engine
'''
#%%
eng = matlab.engine.start_matlab() # '-nodisplay'
#%%
jumpsize_mat = jump_size1.reshape(-1)
#plt.figure()
#plt.hist(jumpsize_mat)
jumpsize_mat = matlab.double(list(jumpsize_mat))
jumpsize_mat.reshape([2*len(jump_size1),1])
#%%
eng.fitdist(jumpsize_mat,'stable',nargout=0)
eng.workspace['result_fit'] = eng.fitdist(jumpsize_mat,'stable')#,nargout=0)

eng.eval("result_fit=fitdist(jumpsize_mat,'stable')", nargout=0)
paramci = eng.eval('result_fit.paramci')
fitparam = eng.eval('[result_fit.alpha result_fit.beta result_fit.gam result_fit.delta]')
fitparam = np.concatenate((np.array(fitparam),np.array(paramci)),0).T
import pickle
with open('fitparam%s'%loop_num, 'wb') as file:
    pickle.dump(fitparam, file)
'''
#%%
eng = matlab.engine.start_matlab('-nodisplay') # '-nodisplay'
#%%
mt_rate = matlab.double(list(poprate_rate_chgadapt_neu))
#cwt_coef =  eng.cwt(mt_rate, 1000.)
#%%
savename = (spke1.param['chg_adapt_range'],spke1.param['new_delta_gk'], \
 spke1.param['delta_gk'], spke1.param['tau_s_di'], \
 spke1.param['tau_s_de'],spke1.param['ie_ratio'], loop_num)
#%%
eng.figure(nargout=0)
eng.eval("set(gcf, 'Position',[10,10,900,350]);", nargout=0)
eng.cwt(mt_rate, 1000., nargout=0)
eng.caxis(matlab.double([0,35]),nargout=0)
matcmd = "saveas(gcf, 'wvt_{:.0f}_{:.0f}_{:.0f}_{:.1f}_{:.1f}_{:.3f}_{:d}.png')".format(*savename)
eng.eval(matcmd, nargout=0)
#%%
eng.quit()
#%%
plt.figure(figsize=[12,6])
plt.plot(np.arange(len(poprate_rate_chgadapt_neu)),poprate_rate_chgadapt_neu)
plt.plot([1000,1000],[0,140],'r--', label='time of adaptation change')
plt.xlabel('time(ms)')
plt.ylabel('rate(Hz)')
title = '''adpt_range:%d,new_delta_gk:%.0f,delta_gk:%.0f
tau_s_di:%.1f,tau_s_de:%.1f,ie_ratio:%.3f'''%savename[:-1]
plt.title(title)#('firing rate of neurons with decreased adaptation\nnew_delta_gk:%.1f;delta_gk:%.1f\ntau_s_di: %.1fms; ie_ratio: %.3f'%(spke1.param['new_delta_gk'], spke1.param['delta_gk'], spke1.param['tau_s_di'], spke1.param['ie_ratio']))
plt.legend()
plt.savefig('rate_chg_adapt_%.0f_%.0f_%.0f_%.1f_%.1f_%.3f_%d.png'%savename)
#%%
magf = abs(fft(poprate_rate_chgadapt_neu[1000:]-poprate_rate_chgadapt_neu[1000:].mean()))
#%%
fig, [ax1, ax2] = plt.subplots(1,2,figsize=(12,6))
fig.suptitle(title)#('Power spectrum of rate of neurons after adaptation change\nnew_delta_gk: %.1f; ndelta_gk: %.1f; tau_s_di: %.1fms; ie_ratio: %.3f'%(spke1.param['new_delta_gk'], spke1.param['delta_gk'], spke1.param['tau_s_di'], spke1.param['ie_ratio']))
#plt.figure(figsize=[12,6])
ax1.plot((np.arange(len(magf))/(len(magf))*1000)[:250], magf[:250])
ax1.set_xlabel('freqeuncy (Hz)')
ax1.set_title('Fourier Transform')
ax2.loglog((np.arange(len(magf))/(len(magf))*1000)[1:int(len(magf)/2)], magf[1:int(len(magf)/2)])
ax2.set_ylim((magf[1:].min()/10, magf[1:].max()*10))
ax2.set_title('Fourier Transform (log-log)')
ax2.set_xlabel('freqeuncy (log(Hz))')
#%%
fig.savefig('fftrate_chg_adapt_%.0f_%.0f_%.0f_%.1f_%.1f_%.3f_%d.png'%savename)
#%%
#plt.figure()
#plt.loglog((np.arange(len(magf))/(len(magf))*1000)[:2400], magf[:2400],'*')
#%%
'''
spkrate1, centre_ind1, jump_size1, jump_dist1 = psa.get_rate_centre_jumpdist(spke1, starttime=0*ms, endtime=8000*ms, binforrate=10*ms, sample_interval=1*ms, slide_interval = 1, jump_interval = 15, dt=0.1*ms, show_trajectory=False)
spkrate1 = psa.overlap_centreandspike(centre_ind1, spkrate1, show_trajectory = False)
anititle='change adaptation in region near centre to %.1fnS at 1000 ms\ndelta_gk: %.1fnS; tau_s_di: %.1fms; ie_ratio: %.3f'%(spke1.param['new_delta_gk'], spke1.param['delta_gk'], spke1.param['tau_s_di'], spke1.param['ie_ratio'])
ani1 = psa.show_pattern(spkrate1, spkrate2=0, area_num = 1, frames = 5500, bottom_up=1, top_down=1, stimu_onset = -1, start_time = 0, anititle=anititle)
#ani1.save('chg_adapt_%d_%.1f_%.1f_%.1f_%.1f_%.3f_%d.png'%savename)
'''



