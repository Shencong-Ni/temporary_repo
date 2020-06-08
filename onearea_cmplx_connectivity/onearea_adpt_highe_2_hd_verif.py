#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:33:44 2020

@author: shni2598
"""

'''
change adaptation of neurons in a small area to make this small area produce theta oscillation
high excitatory coupling strength and high inhibitory coupling strength(E-I, E-E; I-I, I-E)
'''
import brian2.numpy_ as np
from brian2.only import *
from connection import pre_process
import connection as cn
import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.pyplot as plt
#import matplotlib.animation as animation
import time
import sys
#import scipy.io as sio
import pickle

#%%
dir_cache = './cache/cache%s' %sys.argv[1]
prefs.codegen.runtime.cython.cache_dir = dir_cache
prefs.codegen.max_cache_dir_size = 20
#%%
'''
scale_ie_ratio = np.arange(0.8,1.25,0.05)
scale_ei = np.arange(0.8,1.25,0.05)
scale_ii = np.arange(0.8,1.25,0.05)
#%%
sys_argv = int(sys.argv[1])
loop_num = -1
for i in range(len(scale_ie_ratio)):
    for j in range(len(scale_ei)):
        for k in range(len(scale_ii)):
            loop_num += 1
            if loop_num == sys_argv:
                print(i,j,k)
                break
        else:continue
        break
    else:continue
    break
'''
#%%
sys_argv = int(sys.argv[1])
loop_num = -1
for chg_adapt_range in [6]*10: #np.arange(6,9):
    for ie_ratio in [3.20625]: #3.375*np.arange(0.9, 1.21, 0.05):
        for tau_s_de in [4.5]: #np.arange(4., 6.1, 0.5):
            for tau_s_di in [4.]: #np.arange(3., 6.1, 0.5): # ms
                for new_delta_gk in [2]: #np.arange(0, 3.1, 1): # nS
                    for delta_gk in [9]: #np.arange(8, 14.1, 1):
                        loop_num += 1
                        if loop_num == sys_argv:
                            print(loop_num)
                            break
                    else: continue
                    break
                else:continue
                break
            else:continue
            break
        else:continue
        break
    else:continue
    break
        
#%%
ijwd = pre_process.get_ijwd()
scale_ee = 1.4; scale_ei = 1.4
ijwd.ie_ratio = ie_ratio# 3.375 #3.375 #* scale_ie_ratio[i]
ijwd.mean_J_ee = 4*10**-3 * scale_ee#* 1.4 # usiemens
ijwd.sigma_J_ee = 1.9*10**-3 * scale_ee#* 1.4# usiemens

ijwd.generate_ijw()
ijwd.generate_d_rand()

ijwd.w_ei = 5*10**(-3) * scale_ei * 1.05#usiemens
ijwd.w_ii = 25*10**(-3) * 1.4  #usiemens
#%%
with open('ijwd%s'%loop_num, 'wb') as file:
    pickle.dump(ijwd, file)
#%%
chg_adapt_loca = [0, 0]
#chg_adapt_range = 5
chg_adapt_neuron = cn.findnearbyneuron.findnearbyneuron(ijwd.lattice_ext, chg_adapt_loca, chg_adapt_range, ijwd.width)

#%%
start_scope()

neuronmodel_e = '''
dv/dt = (1/C)*(-g_l*(v - v_l) + (-g_k)*(v - v_k) + (-g_I)*(v - v_rev_I) +(-g_E-g_E_inter)*(v - v_rev_E) + I_extnl) : volt (unless refractory)
dg_k/dt = -g_k/tau_k :siemens
delta_gk : siemens 
g_I : siemens
g_E : siemens
g_E_inter : siemens
I_extnl : amp
'''
neuronmodel_i = '''
dv/dt = (1/C)*(-g_l*(v - v_l) + (-g_I)*(v - v_rev_I) +(-g_E-g_E_inter)*(v - v_rev_E) + I_extnl) : volt (unless refractory)
g_I : siemens
g_E : siemens
g_E_inter : siemens
I_extnl : amp
'''

synapse_e = '''
w: siemens
g_E_post = w*s : siemens (summed)
ds/dt = -s/tau_s_de + rect_puls*(1 - s) : 1 (clock-driven)
rect_puls = (1/tau_s_re)*effect : Hz
effect : integer
'''
synapse_i = '''
w: siemens
g_I_post = w*s : siemens (summed)
ds/dt = -s/tau_s_di + rect_puls*(1 - s) : 1 (clock-driven)
rect_puls = (1/tau_s_ri)*effect : Hz
effect : integer
'''


synapse_e_v1v2 = '''
w: siemens
g_E_inter_post = w*s : siemens (summed)
ds/dt = -s/tau_s_de + rect_puls*(1 - s) : 1 (clock-driven)
rect_puls = (1/tau_s_re)*effect : Hz
effect : integer
'''

#synapse_ext = '''
#w: siemens
#g_ext_post = w*s : siemens (summed)
#ds/dt = -s/tau_s_de + rect_puls : 1 (clock-driven)
#rect_puls = (1/tau_s_re)*effect : Hz
#effect : 1
#
#''' group_e_v1, group_i_v1, syn_ee_v1, syn_ei_v1, syn_ie_v1, syn_ii_v1, vexi, spike_e
#%%
group_e_v1 =NeuronGroup(ijwd.Ne, model=neuronmodel_e,
                     threshold='v>v_threshold', method='euler',
                     reset='''v = v_reset
                              g_k += delta_gk''', refractory='(t-lastspike)<t_ref')

group_i_v1 =NeuronGroup(ijwd.Ni, model=neuronmodel_i,
                     threshold='v>v_threshold', method='euler',
                     reset='v = v_reset', refractory='(t-lastspike)<t_ref')

syn_ee_v1 = Synapses(group_e_v1, group_e_v1, model=synapse_e, method='euler',
                  on_pre={'up':'effect += 1', 'down': 'effect -= 1'})
syn_ei_v1 = Synapses(group_e_v1, group_i_v1, model=synapse_e, method='euler',
                  on_pre={'up':'effect += 1', 'down': 'effect -= 1'})
syn_ie_v1= Synapses(group_i_v1, group_e_v1, model=synapse_i, method='euler',
                  on_pre={'up':'effect += 1', 'down': 'effect -= 1'})
syn_ii_v1 = Synapses(group_i_v1, group_i_v1, model=synapse_i, method='euler',
                  on_pre={'up':'effect += 1', 'down': 'effect -= 1'})
#%%
syn_ee_v1.connect(i=ijwd.i_ee, j=ijwd.j_ee)
syn_ei_v1.connect(i=ijwd.i_ei, j=ijwd.j_ei)
syn_ie_v1.connect(i=ijwd.i_ie, j=ijwd.j_ie)
syn_ii_v1.connect(i=ijwd.i_ii, j=ijwd.j_ii)

#%%
syn_ee_v1.w = ijwd.w_ee*usiemens
syn_ei_v1.w = ijwd.w_ei*usiemens #5*nS
syn_ii_v1.w = ijwd.w_ii*usiemens #25*nS
syn_ie_v1.w = ijwd.w_ie*usiemens
#w_ext = 2*nS
#syn_pois_e.w = w_ext
#syn_pois_i.w = w_ext
#%%
def set_delay(syn, delay_up):
    #n = len(syn)
    syn.up.delay = delay_up*ms
    syn.down.delay = (delay_up + 1)*ms
    
    return syn   
#%%
#ijwd.generate_d_dist()
#ijwd.generate_d_rand()
syn_ee_v1 = set_delay(syn_ee_v1, ijwd.d_ee)
syn_ie_v1 = set_delay(syn_ie_v1, ijwd.d_ie)
syn_ei_v1 = set_delay(syn_ei_v1, ijwd.d_ei)
syn_ii_v1 = set_delay(syn_ii_v1, ijwd.d_ii)
#syn_pois_e = set_delay(syn_pois_e)
#syn_pois_i = set_delay(syn_pois_i)

#syn_ee_v1.up.delay = 3*ms; syn_ee_v1.down.delay = 4*ms; 
#syn_ie_v1.up.delay = 3*ms; syn_ie_v1.down.delay = 4*ms; 
#syn_ei_v1.up.delay = 3*ms; syn_ei_v1.down.delay = 4*ms; 
#syn_ii_v1.up.delay = 3*ms; syn_ii_v1.down.delay = 4*ms; 

syn_ee_v1.effect = 0
syn_ie_v1.effect = 0
syn_ei_v1.effect = 0
syn_ii_v1.effect = 0
#syn_pois_e.effect = 0
#syn_pois_i.effect = 0
#group_e_v1.v = np.random.random(Ne)*35*mV-85*mV
#group_i_v1.v = np.random.random(Ni)*35*mV-85*mV
group_e_v1.v = np.random.random(ijwd.Ne)*35*mV-85*mV
group_i_v1.v = np.random.random(ijwd.Ni)*35*mV-85*mV
group_e_v1.delta_gk = delta_gk*nS
#group_e_v1.v = np.random.random(ijwd.Ne)*10*mV-60*mV
#group_i_v1.v = np.random.random(ijwd.Ni)*10*mV-60*mV
group_e_v1.I_extnl = 0.51*nA
group_i_v1.I_extnl = 0.6*nA

#%%
#v_v1 = StateMonitor(group_e_v1, ('v'), record=True,dt=0.5*ms,name='v_v1')
spk_e1 = SpikeMonitor(group_e_v1, record = True)

#%%
net = Network(collect())
#net.add(spk_e1)
net.store('state1')
#%%
print('ie_w: %fnsiemens' %(syn_ie_v1.w[0]/nsiemens))
#Ne = 63*63; Ni = 1000;
C = 0.25*nF # capacitance
g_l = 16.7*nS # leak capacitance
v_l = -70*mV # leak voltage
v_threshold = -50*mV
v_reset = -60*mV
v_rev_I = -80*mV
v_rev_E = 0*mV
v_k = -85*mV
tau_k = 80*ms# 80*ms
#delta_gk = 10*nS #10*nS
t_ref = 4*ms # refractory period

'''
tau_s_de_list = np.arange(2,8,0.5)*ms
tau_s_di_list = np.arange(2,10,0.5)*ms
sys_argv = int(sys.argv[1])
loop_num = 0
for i in range(len(tau_s_de_list)):
    for j in range(len(tau_s_di_list)):
        loop_num += 1
        if loop_num == sys_argv:
            print(i,j)
            break
    else:
        continue
    break
'''        
      
tau_s_de = tau_s_de*ms # 5*ms
tau_s_di = tau_s_di*ms # 3*ms
tau_s_re = 1*ms
tau_s_ri = 1*ms
tic = time.perf_counter()
#seed(10)
simu_time1 = 1000*ms
simu_time2 = 7000*ms

net.run(simu_time1, report = 'text', profile=False) #,namespace={'tau_k': 80*ms}
group_e_v1.delta_gk[chg_adapt_neuron] = new_delta_gk*nS
net.run(simu_time2, report = 'text', profile=False) #,namespace={'tau_k': 80*ms}

print('%ds'%(time.perf_counter() - tic))

#%%
'''
sio.savemat('index%s.mat'%sys_argv, {'ind':spk_e1.i[:]})
spk_tstep = np.round(spk_e1.t/(0.1*ms)).astype(int)
sio.savemat('time%s.mat'%sys_argv, {'time':spk_tstep})
'''
#%
spk_tstep = np.round(spk_e1.t/(0.1*ms)).astype(int)
param = {'delta_gk':delta_gk,'new_delta_gk':new_delta_gk,'tau_s_di':tau_s_di/ms,\
         'tau_s_de':tau_s_de/ms, 'ie_ratio':ie_ratio, 'chg_adapt_range':chg_adapt_range}
data = {'param':param, 'i':spk_e1.i[:], 't':spk_tstep}

with open('data%s'%loop_num, 'wb') as file:
    pickle.dump(data, file)

#%%
'''
import post_analysis as psa
#%%
spkrate1, centre_ind1, jump_size1, jump_dist1 = psa.get_rate_centre_jumpdist(spk_e1, starttime=0*ms, endtime=8000*ms, binforrate=10*ms, sample_interval=1*ms, slide_interval = 1, jump_interval = 15, dt=0.1*ms, show_trajectory=False)
spkrate1 = psa.overlap_centreandspike(centre_ind1, spkrate1, show_trajectory = False)
anititle='change adaptation in region near centre to %.1fnS at 1000 ms\ntau_s_di: %.1fms; ie_ratio: %.3f'%(new_delta_gk, tau_s_di/ms, ie_ratio)
ani1 = psa.show_pattern(spkrate1, spkrate2=0, area_num = 1, frames = 5500, bottom_up=1, top_down=1, stimu_onset = -1, start_time = 0, anititle=anititle)
#%%
ani1.save('chg_adapt_1_%.1f_%.1f_%.3f.mp4'%(new_delta_gk, tau_s_di/ms, ie_ratio))
'''
#%%
'''
vrec1 = v_v1.v[:] #extract membrane potential data
vrec1 = vrec1.reshape(63,63,vrec1.shape[1])
vrec1 = vrec1/mV

#%%
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
fig, ax1= plt.subplots(1,1)
value1=ax1.matshow(vrec1[:,:,0])

def updatev(iii):
    value1.set_array(vrec1[:,:,iii])
    return value1

cbaxes = fig.add_axes([0.1, 0.06, 0.8, 0.03]) 
cb=fig.colorbar(value1, cax = cbaxes, orientation='horizontal') 
value1.set_clim(vmin=-85, vmax=-50)

ani=animation.FuncAnimation(fig, updatev, frames=int(vrec1.shape[2]), interval=1)   

#%%
import post_analysis
#%%
analy1 = post_analysis.analysis()

#analy1.get_spike_rate_animation

#anly_tmp = post_analysis.analysis()
spk = analy1.get_spike_rate_animation(spk_e1, start_time=0*ms, ani_time=400*ms, n_neuron = 3969, window = 10*ms, dt = 0.1*ms)
centre = analy1.get_centre_mass(spk)
spk = analy1.overlap_centreandspike(centre, spk)
#%%
fig, ax1= plt.subplots(1,1)


value1=ax1.matshow(spk[:,:,0], cmap=plt.cm.get_cmap('viridis', 5))
ax1.axis('off')
#cmap=plt.cm.get_cmap('binary', 3)

def updatev(iii):
    value1.set_array(spk[:,:,iii])
    return value1
#
cbaxes = fig.add_axes([0.1, 0.06, 0.8, 0.03]) 
cb=fig.colorbar(value1, cax = cbaxes, orientation='horizontal',ticks=[0,1,2,3,4]) 
value1.set_clim(vmin=0, vmax=4)

#fig.suptitle('ie_ratio:%s\ncentre of pattern and number of spikes per 10ms\n ' %ie_ratio)
ani=animation.FuncAnimation(fig, updatev, frames=spk.shape[2], interval=0.1)  
#ani.show()

'''
#%%
#net.restore('state1')

#%%

#%%


