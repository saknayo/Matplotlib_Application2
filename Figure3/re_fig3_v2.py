import struct
import matplotlib.pyplot as plt
import matplotlib
import numpy
plt.style.use('seaborn-paper')

#################
dpi = 600
fig=plt.figure(figsize=[6400/dpi,5000/dpi])
fig.subplots_adjust(hspace=0.3, wspace=0.3)
sci_formater=matplotlib.ticker.ScalarFormatter(useMathText=True) # scientific notation 
sci_formater.set_powerlimits((0,0))
bar_color={ 'dtag_fc' :'black' , 'dtag_ec' :'black',
			'stag_fc' :(160/255,160/255,160/255)  , 'stag_ec' :'black',
			'notag_fc':'white' , 'notag_ec':'black'}
######### subplot 1 #########

####### data for subplot 1 ########
step  = 500000
PLMC  = [2E-6,50E-6,100E-6,150E-6,200E-6,250E-6,300E-6,350E-6,400E-6,450E-6,500E-6,550E-6,600E-6,650E-6,700E-6,750E-6,800E-6,850E-6,900E-6,950E-6,1000E-6]
NOTAG = numpy.zeros(len(PLMC))
STAG  = numpy.zeros(len(PLMC))
DTAG  = numpy.zeros(len(PLMC))
dirs  = ['2','50','100','150','200','250','300','350','400','450','500','550','600','650','700','750','800','850','900','950','1000']
for i in range(len(dirs)): 
	with open('non-ribozyme derived parasites/'+ dirs[i] + '/Release/pol-a-6.dat','rb')as f1: # tag-parasite (including dtaged-parasite)
		with open('non-ribozyme derived parasites/'+ dirs[i] + '/Release/pol-a-14.dat','rb')as f2: # dtag-parasite
		    with open('non-ribozyme derived parasites/'+ dirs[i] + '/Release/pol-a-23.dat','rb')as f3: # notag-rna & len > POLCOVERSEQ
			    f1.seek(int(step/1000-1)*4)			
			    f2.seek(int(step/1000-1)*4)
			    f3.seek(int(step/1000-1)*4)
			    STAG[i] =struct.unpack('f',f1.read(4))[0]
			    DTAG[i] =struct.unpack('f',f2.read(4))[0]
			    NOTAG[i]=struct.unpack('f',f3.read(4))[0]
STAG=STAG-DTAG # exclude dtags
#print(STAG,DTAG,NOTAG)

ax1=fig.add_subplot(2,2,1)
bar_width1=max(PLMC)/len(PLMC)*0.7
ax1.xaxis.set_major_formatter(sci_formater)
ax1.set_xlabel('PRL (Probability of random ligation)',size=11)
ax1.set_ylabel('Number of molecules',size=11)
ax1.text(-0.15, 1.1,'A',size=11, ha='center', va='center', transform=ax1.transAxes) # A marker
ax1.axis([ 0, 1200e-6, 0, 8000])
#ax1.xaxis.set_visible(False)
#ax1.yaxis.set_label_coords(-0.12, 0.1)
# ax1.xaxis.tick_top()
ax1.tick_params(axis='x', which='major',direction='out')

#ax1.plot(PLMC,NOTAG,'ko')
#ax1.plot(PLMC,STAG,'yo')
#ax1.plot(PLMC,DTAG,'ro')
ax1.bar( PLMC,  DTAG, width=bar_width1, align='center', facecolor =bar_color['dtag_fc']  ,edgecolor =bar_color['dtag_ec']  )
ax1.bar( PLMC,  STAG, width=bar_width1, align='center', facecolor =bar_color['stag_fc']  ,edgecolor =bar_color['stag_ec'] , bottom=DTAG )
ax1.bar( PLMC, NOTAG, width=bar_width1, align='center', facecolor =bar_color['notag_fc'] ,edgecolor =bar_color['notag_ec'], bottom=DTAG+STAG )

## subplot of ax1
ax3=fig.add_subplot(2,2,3)
ax3.xaxis.set_major_formatter(sci_formater)
ax3.axis([ 0, 1200e-6, 0, 80])
#ax3.yaxis.get_major_ticks()[-1].label1.set_visible(False) # hide the last tick of y-axis
ax3.set_xlabel('PRL (Probability of random ligation)',size=11)
ax3.set_ylabel('Number of molecules',size=11)
ax3.tick_params(axis='x', which='major',direction='out')
#ax3.text(-0.1, 1.1,'A',size=11, ha='center', va='center', transform=ax1.transAxes) # A marker
#ax1.plot(PLMC,NOTAG,'ko')
#ax1.plot(PLMC,STAG,'yo')
#ax1.plot(PLMC,DTAG,'ro')
ax3.bar( PLMC,  DTAG, width=bar_width1, align='center', facecolor =bar_color['dtag_fc']  ,edgecolor =bar_color['dtag_ec'] )
ax3.bar( PLMC,  STAG, width=bar_width1, align='center', facecolor =bar_color['stag_fc']  ,edgecolor =bar_color['stag_ec'] , bottom=DTAG )
ax3.bar( PLMC, NOTAG, width=bar_width1, align='center', facecolor =bar_color['notag_fc'] ,edgecolor =bar_color['notag_ec'], bottom=DTAG+STAG )
bbox3=ax3.get_position()
ax3.set_position((bbox3.x0,(bbox3.y0+bbox3.y1)/2,bbox3.x1-bbox3.x0,(bbox3.y1-bbox3.y0)/2)) # (left,bottom,width,height)

######## data for subplot 2 ######
stepnum = 40000
recint  = 1000 # interval for record funtion
records_num= int(stepnum/recint+1)

path = 'ribozyme-derived parasites/Release/'
tag_sb2  =numpy.zeros( records_num)#, dtype=numpy.int8 )
dtag_sb2 =numpy.zeros( records_num)#, dtype=numpy.int8 )
notag_sb2=numpy.zeros( records_num)#, dtype=numpy.int8 )
with open(path+'pol-a-6.dat','rb')as f1 : # tag-parasite (including dtaged-parasite)
    with open(path+'pol-a-14.dat','rb')as f2 : # dtag-parasite
        with open(path+'pol-a-23.dat','rb')as f3 : # notag-rna & len > POLCOVERSEQ
            for i in range(records_num):
                tag_sb2[i]  = struct.unpack('f',f1.read(4))[0] 
                dtag_sb2[i] = struct.unpack('f',f2.read(4))[0] 
                notag_sb2[i]= struct.unpack('f',f3.read(4))[0] 
                #print(tag_sb2[i],dtag_sb2[i],notag_sb2[i])
            tag_sb2=tag_sb2-dtag_sb2 # tag_sb2 only contain single taged parasite

#### subplot 2 ###
ax2=fig.add_subplot(2,2,2)
plotint=1
bar_width=recint*plotint*0.7
dotnum   = int(stepnum/(recint*plotint)+1)

ax2.set_xlabel('Time (Monte Carlo steps)',size=11)
ax2.set_ylabel('Number of Molecules',size=11)
ax2.text(-0.15, 1.1,'B',size=11, ha='center', va='center', transform=ax2.transAxes) # B marker
ax2.xaxis.set_major_formatter(sci_formater)
ax2.axis([10000,35000,0,8000])#(dtag_sb2[::plotint]+tag_sb2[::plotint]+notag_sb2[::plotint]).max()])
ax2.tick_params(axis='x', which='major',direction='out')
#ax2.xaxis.set_visible(False)
#ax2.yaxis.set_label_coords(-0.12, 0.1)

x=numpy.array(list(i*recint*plotint for i in range(dotnum)))
ax2.bar( x,  dtag_sb2[::plotint], width=bar_width, align='center', facecolor =bar_color['dtag_fc']  ,edgecolor =bar_color['dtag_ec'] )
ax2.bar( x,   tag_sb2[::plotint], width=bar_width, align='center', facecolor =bar_color['stag_fc']  ,edgecolor =bar_color['stag_ec'],
		 bottom=dtag_sb2[::plotint] )
ax2.bar( x, notag_sb2[::plotint], width=bar_width, align='center', facecolor =bar_color['notag_fc'] ,edgecolor =bar_color['notag_ec'], 
		 bottom=dtag_sb2[::plotint]+tag_sb2[::plotint] )

# subplot for ax2
ax4=fig.add_subplot(2,2,4)
ax4.set_xlabel('Time (Monte Carlo steps)',size=11)
ax4.set_ylabel('Number of Molecules',size=11)
#ax4.text(-0.1, 1.1,'B',size=11, ha='center', va='center', transform=ax2.transAxes) # B marker
ax4.xaxis.set_major_formatter(sci_formater)
ax4.axis([10000,35000,0,80])#(dtag_sb2[::plotint]+tag_sb2[::plotint]+notag_sb2[::plotint]).max()])
ax4.tick_params(axis='x', which='major',direction='out')

ax4.bar( x,  dtag_sb2[::plotint], width=bar_width, align='center', facecolor =bar_color['dtag_fc']  ,edgecolor =bar_color['dtag_ec'] )
ax4.bar( x,   tag_sb2[::plotint], width=bar_width, align='center', facecolor =bar_color['stag_fc']  ,edgecolor =bar_color['stag_ec'],
		 bottom=dtag_sb2[::plotint] )
ax4.bar( x, notag_sb2[::plotint], width=bar_width, align='center', facecolor =bar_color['notag_fc'] ,edgecolor =bar_color['notag_ec'], 
		 bottom=dtag_sb2[::plotint]+tag_sb2[::plotint] )
bbox4=ax4.get_position()
ax4.set_position((bbox4.x0,(bbox4.y0+bbox4.y1)/2,bbox4.x1-bbox4.x0,(bbox4.y1-bbox4.y0)/2)) # (left,bottom,width,height)

#####################
plt.savefig('fig3_{}.png'.format(recint*plotint),dpi=dpi) 
plt.show()