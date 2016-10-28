# -*- coding: utf-8 -*-
import struct
import matplotlib.pyplot as plt
import matplotlib
import numpy 
plt.style.use('seaborn-paper')

stepnum = 10000000
stepnum2= 4000000
recint  = 1000 # interval for record funtion
records_num = int(stepnum /recint+1)
records_num2= int(stepnum2/recint+1)

plotint=50
dpi=600
path1='withtag/Release/'
path2="withouttag 3+9+3/Release/"
path3="without_tag_0+9+0_stat/"
path4="withtag_prb_changed@300W/Release/"

sci_formater=matplotlib.ticker.ScalarFormatter(useMathText=True) # scientific notation 
sci_formater.set_powerlimits((0,0))

def dat_import(path,records):
	pol_array  =numpy.zeros(records)
	tag_array  =numpy.zeros(records)
	contl_array=numpy.zeros(records)
	with open(path+'pol-a-3.dat','rb')as f1: # poltag
		with open(path+'pol-a-4.dat','rb')as f2: # poltagcom
			with open(path+'pol-a-5.dat','rb')as f3: # contrl1
				with open(path+'pol-a-11.dat','rb')as f4: # contrl2
					with open(path+'pol-a-6.dat','rb')as f5: # tag-parasite
						for i in range(records):
							pol_array[i]  +=struct.unpack('f',f1.read(4))[0]
							pol_array[i]  +=struct.unpack('f',f2.read(4))[0]
							contl_array[i]+=struct.unpack('f',f3.read(4))[0]
							contl_array[i]+=struct.unpack('f',f4.read(4))[0]
							tag_array[i]  +=struct.unpack('f',f5.read(4))[0]
	return pol_array,tag_array,contl_array

def dat_import_notag(path,records):
	pol_array  =numpy.zeros(records)
	tag_array  =numpy.zeros(records)
	contl_array=numpy.zeros(records)
	#－2 pol_num
	#-23 notag_par_num
	with open(path+'pol-a-3.dat','rb')as f1: # poltag
		with open(path+'pol-a-4.dat','rb')as f2: # poltagcom
			with open(path+'pol-a-5.dat','rb')as f3: # contrl1
				with open(path+'pol-a-11.dat','rb')as f4: # contrl2
					with open(path+'pol-a-6.dat','rb')as f5: # tag-parasite
						with open(path+'pol-a-2.dat','rb')as f6: # pol
							with open(path+'pol-a-28.dat','rb')as f7: # polcom
								with open(path+'pol-a-23.dat','rb')as f8: # notag-parasite
									for i in range(records):
										pol_array[i]  +=struct.unpack('f',f1.read(4))[0]
										pol_array[i]  +=struct.unpack('f',f2.read(4))[0]
										contl_array[i]+=struct.unpack('f',f3.read(4))[0]
										contl_array[i]+=struct.unpack('f',f4.read(4))[0]
										tag_array[i]  +=struct.unpack('f',f5.read(4))[0]
										pol_array[i]  +=struct.unpack('f',f6.read(4))[0]
										pol_array[i]  +=struct.unpack('f',f7.read(4))[0]
										tag_array[i]  +=struct.unpack('f',f8.read(4))[0]
	return pol_array,tag_array,contl_array

def ax_plot(ax,x,pol,tag,contrl,
			xlabel='Time (Monte Carlo steps)',
			ylabel='Number of molecules'     ):
	ax.set_xlabel(xlabel,size=11)
	ax.set_ylabel(ylabel,size=11)
	ax.xaxis.set_major_formatter(sci_formater)
	# mew:markeredgewidth ,mec:markeredgecolor , mfc:markerfacecolor
	ax.plot(x, pol,    marker='o',linestyle='None',mfc ='w' ,mec ='k',mew=1, markersize=3, label='Polymerase',  )
	ax.plot(x, tag,    marker='^',linestyle='None',mfc ='w' ,mec ='k',mew=1, markersize=3, label='Tag parasite' )
	ax.plot(x, contrl, marker='s',linestyle='None',mfc ='w' ,mec ='k',mew=1, markersize=3, label='Control'      )

fig=plt.figure(figsize=[6400/dpi,5000/dpi])
fig.subplots_adjust(hspace=0.50, wspace=0.3)
x=numpy.arange(0,stepnum+1,recint)
pol1,tag1,contrl1=dat_import(path1,records_num )
pol2,tag2,contrl2=dat_import(path2,records_num2)
pol3,tag3,contrl3=dat_import_notag(path3,records_num2)
pol4,tag4,contrl4=dat_import(path4,records_num )

ax1 =fig.add_subplot(3,1,1)
ax_plot(ax1,x[::plotint],pol1[::plotint],tag1[::plotint],contrl1[::plotint])
ax1.axis([0,6000000,0,400])
ax1.text(-0.05, 1.1,'A',size=11, ha='center', va='center', transform=ax1.transAxes) # A marker

x2=numpy.arange(0,stepnum2+1,recint)
ax2 =fig.add_subplot(3,2,3)
ax_plot(ax2,x2[::plotint],pol2[::plotint],tag2[::plotint],contrl2[::plotint])
ax2.axis([0,2500000,0,400])
ax2.text(-0.117, 1.1,'B',size=11, ha='center', va='center', transform=ax2.transAxes) # B marker

ax3 =fig.add_subplot(3,2,4)
ax_plot(ax3,x2[::plotint],pol3[::plotint],tag3[::plotint],contrl3[::plotint])
ax3.axis([0,2500000,0,400])
ax3.text(-0.117, 1.1,'C',size=11, ha='center', va='center', transform=ax3.transAxes) # C marker

ax4 =fig.add_subplot(3,1,3)
ax_plot(ax4,x[::plotint],pol4[::plotint],tag4[::plotint],contrl4[::plotint])
ax4.axis([0,6000000,0,400])
ax4.text(-0.05, 1.1,'D',size=11, ha='center', va='center', transform=ax4.transAxes) # D marker

plt.savefig('fig2_v3_{}.png'.format(recint*plotint),dpi=dpi)
plt.show()
