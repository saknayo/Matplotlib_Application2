# -*- coding: utf-8 -*-
'''
v3:
    display the numbers of dtag-parasites over 15nt 
'''
import struct
import matplotlib.pyplot as plt
import matplotlib
import numpy 
import numpy as np
plt.style.use('seaborn-paper')

def colorbar_index(ncolors, mappable,cbar_ax):
    #cmap = cmap_discretize(cmap, ncolors)
    #mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable,cax=cbar_ax,orientation='horizontal')
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    colorbar.set_ticklabels(range(ncolors))
    return colorbar

path='withtag/Release/'
stepnum = 6000000
recint  = 1000 # interval for record funtion
maxlen  = 100 # rna length
records_num= int(stepnum/recint+1)

dpi=600
plotint  = 50 # for display
dotnum   = int(stepnum/(recint*plotint)+1)

fig=plt.figure(figsize=[6400/dpi,5000/dpi])
fig.subplots_adjust(hspace=0.45, wspace=0.3)
ax1=fig.add_subplot(2,1,1)
ax2=fig.add_subplot(2,1,2)
sci_formater=matplotlib.ticker.ScalarFormatter(useMathText=True) # scientific notation 
sci_formater.set_powerlimits((0,0))

bar_width= int(recint*plotint*0.9)
bar_color={ 'dtag_fc' :'black' , 'dtag_ec' :'black',
            'stag_fc' :(160/255,160/255,160/255)  , 'stag_ec' :'black',
            'notag_fc':'white' , 'notag_ec':'black'} # for subplot 1
# for subplot 2
lenbox=[6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]#, 16, 17, 18, 19, 20]
#colorb=['#FF0000','#FF1493','#FF4040','#FF6A6A','#FF7F24','#FF83FA','#FFA500','#FFB6C1',
#        '#FFC125','#FFDAB9','#FFE4C4','#FFEC8B','#FFF5EE','#FFFAF0','#FFFFF0']
#colorb=list((1,(255-130)*(i/len(lenbox))/255,(255-130)*(i/len(lenbox))/255) for i in range(len(lenbox)))
#colorb=list((214/255,(39+(214-39)*(i/len(lenbox)))/255,39/255) for i in range(len(lenbox)))
#colorb=['#d300a9']+list((214/255,(39+(214-39)*(i/len(lenbox)))/255,39/255) for i in range(1,len(lenbox)))
#colorb=['#ff2acf','#b760e6',        '#981b1e','#e31c3d','#e59393','#f9dede',         '#f9a560',        '#f2d96e','#d6ad5e','#a1884e']
#colorb=['#d300a9']+list((1,(255*(i/len(lenbox)))/255,(255*(i/len(lenbox)))/255) for i in range(1,len(lenbox)))
#colorb=list( ((i/len(lenbox)),(i/len(lenbox)),(i/len(lenbox))) for i in range(len(lenbox)) )
#colorb=['#d300a9']+list( (  (255+0*(i/len(lenbox)))/255, 
#                            (0+255*(i/len(lenbox)))/255,
#                            (0+0*(i/len(lenbox)))/255
#                         )  for i in range(len(lenbox)-1))
#colorb=['#d300a9','#db0000','#c50000','#a90000','#900000', '#cc5803','#e2711d','#ff9505','#ffb627','#f6e228']
#
#colorb=['#712387','#900000','#a90000','#c50000','#db0000', '#cc5803','#e2711d','#ff9505','#ffb627','#f6e228']
#
#colorb=['#cd19d8','#b10821','#c90926','#e20a2a','#f41133', '#cc5803','#e2711d','#ff9505','#ffb627','#f6e228']

colorb=['#D80FB1','#D40C35','#D82933','#DC4632','#E06331','#E58030','#E99D2F','#EDBA2E','#F1D72D','#F6F42C','#FFFFFF']#a
#colorb=['#AA0CC5','#C50C65','#CC164D','#D22137','#D9362D','#DF6039','#E68845','#ECAE53','#F3D161','#F9F270']#b
#colorb=['#E217C8','#E00B46','#E21730','#E52A23','#E75230','#EA783D','#ED9B4B','#EFBB58','#F2D766','#F5F074']#c
my_cmap = matplotlib.colors.ListedColormap(colorb)
CS3 = plt.contourf([[0,0],[0,0]], lenbox+[lenbox[-1]+1], cmap=my_cmap)

#### data for subplot 1 ###
tag_sb1  =numpy.zeros( records_num)#, dtype=numpy.int8 )
dtag_sb1 =numpy.zeros( records_num)#, dtype=numpy.int8 )
notag_sb1=numpy.zeros( records_num)#, dtype=numpy.int8 )
with open(path+'pol-a-6.dat','rb')as f1 : # tag-parasite
    with open(path+'pol-a-14.dat','rb')as f2 : # dtag-parasite
        with open(path+'pol-a-23.dat','rb')as f3 : # notag-rna
            for i in range(records_num):
                tag_sb1[i]  = struct.unpack('f',f1.read(4))[0] 
                dtag_sb1[i] = struct.unpack('f',f2.read(4))[0] 
                notag_sb1[i]= struct.unpack('f',f3.read(4))[0] 
                #print(tag_sb1[i],dtag_sb1[i],notag_sb1[i])
            tag_sb1=tag_sb1-dtag_sb1 # tag_sb1 only contain single taged parasite


#### data for subplot 2 ###
dtag_array =numpy.zeros( (records_num,maxlen) )
with open(path+'pol-a-26.dat','rb') as f1 : #single chain parasites
    with open(path+'pol-a-27.dat','rb') as f2 : #double chain parasites
        for stepi in range(records_num):
            for leni in range(maxlen):
                f1.read(4) # skip single taged parasites
                f2.read(4)
                dtag_array[stepi,leni]=struct.unpack('i',f1.read(4))[0] + struct.unpack('i',f2.read(4))[0]

#### subplot 1: ration of parasites ###

ax1.set_xlabel('Time (Monte Carlo steps)',size=11)
ax1.set_ylabel('Number of Molecules',size=11)
ax1.text(-0.05, 1.06,'A',size=11, ha='center', va='center', transform=ax1.transAxes) # A marker
ax1.xaxis.set_major_formatter(sci_formater)
ax1.axis([0,stepnum,0,800])#(dtag_sb1[::plotint]+tag_sb1[::plotint]+notag_sb1[::plotint]).max()])
ax1.tick_params(axis='x', which='major',direction='out')

x=numpy.array(list(i*recint*plotint for i in range(dotnum)))
ax1.bar( x,  dtag_sb1[::plotint], width=bar_width, facecolor =bar_color['dtag_fc'] ,edgecolor =bar_color['dtag_ec'],align='center', bottom=numpy.zeros(dotnum))
ax1.bar( x,   tag_sb1[::plotint], width=bar_width, facecolor =bar_color['stag_fc'] ,edgecolor =bar_color['stag_ec'],align='center', bottom=dtag_sb1[::plotint] )
ax1.bar( x, notag_sb1[::plotint], width=bar_width, facecolor =bar_color['notag_fc'] ,edgecolor =bar_color['notag_ec'],
        align='center', bottom=dtag_sb1[::plotint]+tag_sb1[::plotint] )
bbox1=ax1.get_position()
ax1.set_position((bbox1.x0,bbox1.y0,bbox1.x1-bbox1.x0,3*(bbox1.y1-bbox1.y0)/4)) # (left,bottom,width,height)

#### subplot 2: length distribution of dtags ###


x=numpy.array(list(i*recint*plotint for i in range(dotnum)))
bottom=numpy.zeros(dotnum)

for i in range(len(lenbox)) :
    if( i==len(lenbox)-1 ) :
        ax2.bar( x, dtag_array[::plotint,lenbox[i]-1:].sum(axis=1), width=bar_width,facecolor = colorb[i],edgecolor = 'k',align='center', bottom=bottom)
        break
    ax2.bar( x, dtag_array[::plotint,lenbox[i]-1], width=bar_width,facecolor = colorb[i],edgecolor = 'k',align='center', bottom=bottom)
    bottom+=dtag_array[::plotint,lenbox[i]-1]

ax2.set_xlabel('Time (Monte Carlo steps)',size=11)
ax2.set_ylabel('Number of Molecules',size=11)
ax2.text(-0.05, 1.06,'B',size=11, ha='center', va='center', transform=ax2.transAxes) # B marker
ax2.xaxis.set_major_formatter(sci_formater)
ax2.axis([0,stepnum,0,100])#dtag_array.max()])
ax2.tick_params(axis='x', which='major',direction='out')
bbox2=ax2.get_position()
ax2.set_position((bbox2.x0,bbox2.y0+(bbox2.y0+bbox2.y1)/4,bbox2.x1-bbox2.x0,3*(bbox2.y1-bbox2.y0)/4)) # (left,bottom,width,height)

### colorbar
#fig.subplots_adjust(bottom=0.8,right=0.1)
#cbar_ax = fig.add_axes([0.2, 0.03, 0.6, 0.01]) # [left, bottom, width, height]
cbar_ax = fig.add_axes([0.15, 0.45, 0.25, 0.01]) # [left, bottom, width, height]
fcbar=fig.colorbar(CS3,cax=cbar_ax ,orientation='horizontal',drawedges=True)
fcbar.set_ticks(list(i+0.5 for i in lenbox)+[lenbox[-1]+1])
fcbar.set_ticklabels(lenbox[:-1]+['16~','       (nt)'], update_ticks=True)
fcbar.ax.tick_params(axis='x', length=0)
#fcbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar

########################
plt.savefig('fig4_{}_v3.png'.format(recint*plotint),dpi=dpi)



