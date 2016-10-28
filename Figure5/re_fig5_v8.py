# -*- coding: utf-8 -*-
'''
v7:
	delete subplot3.
v8:
	display the numbers of dtag-parasites over 15nt 
'''


import struct
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import numpy 
plt.style.use('seaborn-paper')

stepnum = 10000000
recint  = 1000 # interval for record funtion
maxlen  = 100 # rna length
records_num= int(stepnum/recint+1)

plotint=100
dpi=600
path1='pmv_influence/contrast/'
path4='pmv_influence/pmv-/step_600_0.001/'

######################### figure layout set #########
fig=plt.figure(figsize=[6400/dpi,5000/dpi])
fig.subplots_adjust(hspace=0.3, wspace=0.3)
dotnum   = int(stepnum/(recint*plotint)+1)

gs1 = gridspec.GridSpec(2, 1) #(col,row)
gs1.update( bottom=0.54,top=0.92,hspace=0.0) #gs1.update(bottom=0.48,top=0.8)
ax1 = plt.subplot(gs1[0, 0])
ax2 = plt.subplot(gs1[1, 0])

gs2 = gridspec.GridSpec(2, 1)
gs2.update( bottom=0.07,top=0.45,hspace=0.0)
ax4 = plt.subplot(gs2[0, 0])
ax5 = plt.subplot(gs2[1, 0])

sci_formater=matplotlib.ticker.ScalarFormatter(useMathText=True) # scientific notation 
sci_formater.set_powerlimits((0,0))

# for subplot2,5
bar_width= int(recint*plotint*0.9)
lenbox=[6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]#, 16, 17, 18, 19, 20]
colorb=['#D80FB1','#D40C35','#D82933','#DC4632','#E06331','#E58030','#E99D2F','#EDBA2E','#F1D72D','#F6F42C','#FFFFFF']#a
my_cmap = matplotlib.colors.ListedColormap(colorb)
CS3 = plt.contourf([[0,0],[0,0]], lenbox+[lenbox[-1]+1], cmap=my_cmap)

#### data import for subplot 1 ###
def dat_import(path):
	''' import poltag_num,para_tag_num,contl_num from dat file '''
	pol_array  =numpy.zeros(records_num)
	tag_array  =numpy.zeros(records_num)
	contl_array=numpy.zeros(records_num)
	with open(path+'pol-a-3.dat','rb')as f1: # poltag
		with open(path+'pol-a-4.dat','rb')as f2: # poltagcom
			with open(path+'pol-a-5.dat','rb')as f3: # contrl1
				with open(path+'pol-a-11.dat','rb')as f4: # contrl2
					with open(path+'pol-a-6.dat','rb')as f5: # tag-parasite
						for i in range(records_num):
							pol_array[i]  +=struct.unpack('f',f1.read(4))[0]
							pol_array[i]  +=struct.unpack('f',f2.read(4))[0]
							contl_array[i]+=struct.unpack('f',f3.read(4))[0]
							contl_array[i]+=struct.unpack('f',f4.read(4))[0]
							tag_array[i]  +=struct.unpack('f',f5.read(4))[0]
	return pol_array,tag_array,contl_array

def dat_import2(path):
	''' import dtag-parasite length distribution info from dat file '''
	dtag_array =numpy.zeros( (records_num,maxlen) )
	with open(path+'pol-a-26.dat','rb') as f1 :
	    with open(path+'pol-a-27.dat','rb') as f2 :
	        for stepi in range(records_num):
	            for leni in range(maxlen):
	                f1.read(4) # skip single taged parasites
	                f2.read(4)
	                dtag_array[stepi,leni]=struct.unpack('i',f1.read(4))[0] + struct.unpack('i',f2.read(4))[0]
	return dtag_array

def ax_plot(ax,x,pol,tag,contrl,
			xlabel='Time (Monte Carlo steps)',
			ylabel='Number of molecules'     ):
	ax.set_xlabel(xlabel,size=11)
	ax.set_ylabel(ylabel,size=11)
	ax.axis([0,stepnum,0,800])
	ax.yaxis.set_label_coords(-0.05, -0.025)
	ax.xaxis.set_major_formatter(sci_formater)
	ax.xaxis.set_visible(False)
	# mew:markeredgewidth ,mec:markeredgecolor , mfc:markerfacecolor
	ax.plot(x, pol,    marker='o',linestyle='None',mfc ='w' ,mec ='k',mew=1, markersize=3, label='Polymerase',  )
	ax.plot(x, tag,    marker='^',linestyle='None',mfc ='w' ,mec ='k',mew=1, markersize=3, label='Tag parasite' )
	ax.plot(x, contrl, marker='s',linestyle='None',mfc ='w' ,mec ='k',mew=1, markersize=3, label='Control'      )

def ax_plot2(ax2,x,dtag_array):
	bottom=numpy.zeros(dotnum)
	for i in range(len(lenbox)) :
	    if( i==len(lenbox)-1 ) : # plot bar of 16nt~ 
	        ax2.bar( x, dtag_array[::plotint,lenbox[i]-1:].sum(axis=1), width=bar_width,facecolor = colorb[i],edgecolor = 'k',align='center', bottom=bottom)
	        break
	    ax2.bar( x2, dtag_array[::plotint,lenbox[i]-1], width=bar_width,facecolor = colorb[i],edgecolor = 'k', bottom=bottom,align='center')
	    bottom+=dtag_array[::plotint,lenbox[i]-1]
	ax2.set_xlabel('Time (Monte Carlo steps)',size=11)
	ax2.xaxis.set_major_formatter(sci_formater)
	ax2.axis([0,stepnum,0,400])#dtag_array.max()])
	ax2.yaxis.get_major_ticks()[-1].label1.set_visible(False) # hide the last tick of y-axis
	ax2.tick_params(axis='x', which='major',direction='out')

def add_colorbar(main_ax,cbar_ax):
	''' add color at cbar_ax '''
	fcbar=main_ax.colorbar(CS3,cax=cbar_ax ,orientation='horizontal',drawedges=True)
	fcbar.set_ticks(list(i+0.5 for i in lenbox)+[lenbox[-1]+1])
	fcbar.set_ticklabels(lenbox[:-1]+['16~','       (nt)'], update_ticks=True)
	fcbar.ax.tick_params(axis='x', length=0)

###### subplot 1 #####
x=numpy.arange(0,stepnum+1,recint)
pol1,tag1,contrl1=dat_import(path1)
ax_plot(ax1,x[::plotint],pol1[::plotint],tag1[::plotint],contrl1[::plotint])
ax1.text(-0.06, 1.1,'A',size=11, ha='center', va='center', transform=ax1.transAxes) # A marker

#### subplot 2 ###
x2=numpy.array(list(i*recint*plotint for i in range(dotnum)))
dtag_array = dat_import2(path1)
ax_plot2(ax2,x2,dtag_array)

### colorbar for subplot 2
cbar_ax2 = fig.add_axes([0.15, 0.7, 0.25, 0.01]) # [left, bottom, width, height]
add_colorbar(fig,cbar_ax2)

###### subplot 4 #####
x4=numpy.arange(0,stepnum+1,recint)
pol4,tag4,contrl4=dat_import(path4)

ax_plot(ax4,x4[::plotint],pol4[::plotint],tag4[::plotint],contrl4[::plotint])
ax4.text(-0.06, 1.1,'B',size=11, ha='center', va='center', transform=ax4.transAxes) # B marker

#### subplot 5 ###
x5=numpy.array(list(i*recint*plotint for i in range(dotnum)))
dtag_array5 =dat_import2(path4)

ax_plot2(ax5,x5,dtag_array5)

### colorbar
cbar_ax5 = fig.add_axes([0.15, 0.23, 0.25, 0.01])
add_colorbar(fig,cbar_ax5)

plt.savefig('fig5_{}_v8.png'.format(recint*plotint),dpi=dpi)
#plt.show()
