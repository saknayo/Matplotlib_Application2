from PIL import Image
from PIL import ImageDraw 
from PIL import ImageFont
import re,os

posp_=re.compile('pos\d=[(](?P<l>\d+),(?P<u>\d+),(?P<r>\d+),(?P<d>\d+)[)]')
with open('config.ini','r',encoding='utf-8')as cf:
	s=posp_.search(cf.readline()).group('l','u','r','d')
	pos1=tuple(int(i) for i in s)
	dir1=re.search("dir\d=(?P<dir>\w+)",cf.readline()).group('dir')
	s=posp_.search(cf.readline()).group('l','u','r','d')
	pos2=tuple(int(i) for i in s)
	dir2=re.search("dir\d=(?P<dir>\w+)",cf.readline()).group('dir')

imgbox1={}
imglist1=[]
stepbox=[]
for i in os.listdir(dir1) :
	if i.find('.jpg')!=-1 or i.find('.png')!=-1 :
		pass
	else:
		continue
	stepbox.append(re.findall('(\d+)[.]jpg',i)[0])
	imp=dir1+'/'+i
	with Image.open(imp) as im :
		im=im.crop(pos1)
		imgbox1[i]=im
		imglist1.append(im)

imgbox2={}
imglist2=[]
for i in os.listdir(dir2) :
	if i.find('.jpg')!=-1 or i.find('.png')!=-1 :
		pass
	else:
		continue
	imp=dir2+'/'+i
	with Image.open(imp) as im :
		im=im.crop(pos2)
		imgbox2[i]=im
		imglist2.append(im)
print(sorted(imgbox1))
print(stepbox)
newlength=(pos1[2]-pos1[0])
newwideth=(pos1[3]-pos1[1])
newim= Image.new("RGB", (newlength*len(imglist1)+400,newwideth*2+400), "white")
imdraw=ImageDraw.Draw(newim)
#font = ImageFont.truetype("Times New Roman.ttf", 64)
for i in range(len(imglist1)):
	#imdraw.text((200+(i)*newlength, 120),stepbox[i],(0,0,0),font=font)
	newim.paste( imglist1[i], (200+i*newlength,200) )
	if i<len(imglist2):
		newim.paste( imglist2[i], (200+i*newlength,200+newwideth) )

newim.save('newv1.png')