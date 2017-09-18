# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import random

particleNumber=256
boxsize=8.0
eps=1.0
sigma=1.0
promien=0.2
deltat=0.0001
temp=2.5

delta=2


class czastka:
	def __init__(self,promien,pos,vel):
		self.promien = promien
		self.r=pos # polożenie
		self.v=vel # prędkosc

'''def sila(r):
	f_v=np.zeros(2)
	dlugosc_r=np.linalg.norm(r)
	fv=np.zeros(2)
	fv=-r/dlugosc_r
	fv=fv*G*M*m/(dlugosc_r**2)
	return fv'''


nx = int(np.sqrt(particleNumber))
ny = int(np.sqrt(particleNumber))
particles = []
for i in range(nx):
	for j in range(ny):
		polozenie = np.array([i*delta, j*delta])
		predkosc=np.array([(random.random()-1./2),(random.random() -1./2)])
		particles.append(czastka(promien,polozenie,predkosc) )

sumv=0.0
sumv2=0.0
for p in particles:
	sumv=sumv+p.v
sumv=sumv/particleNumber # prędkość środka masy
for p in particles:
	p.v=(p.v-sumv) # teraz środek masy spoczywa
for p in particles:
	sumv2=sumv2+np.dot(p.v,p.v)/2.0
sumv2=sumv2/particleNumber # średnia energia kinetyczna
fs=np.sqrt(temp/sumv2) # czynnik skalujący, temp - żądana temperatura
for p in particles:
	p.v=p.v*fs # skalujemy


t=0
tmax = 30
dt = 1

while t<tmax:
	#rysowanie
	import matplotlib.pyplot as plt
	from matplotlib.patches import Circle
	if (t%1==0): # co 100-na klatka
		plt.clf() # wyczyść obrazek
		F = plt.gcf() # zdefiniuj nowy
		for i in range(particleNumber): # pętla po cząstkach
			p = particles[i]
			a = plt.gca() # ‘get current axes’ (to add smth to them)
			cir = Circle((p.r[0],p.r[1]), radius=p.promien) # kółko tam gdzie jest cząstka
			a.add_patch(cir) # dodaj to kółko do rysunku
			plt.plot() # narysuj
		plt.xlim((0,boxsize)) # obszar do narysowania
		plt.ylim((0,boxsize))
		F.set_size_inches((6,6)) # rozmiar rysunku
		nStr=str(t) #nagraj na dysk – numer pliku z 5 cyframi, na początku zera, np 00324.png
		nStr=nStr.rjust(5,'0')
		plt.title("Symulacja gazu Lennarda-Jonesa, krok "+nStr)
		plt.savefig('img'+nStr+'.png')
	#oblicz siły
	#aktualizuj prędkości i położenia
	for p in particles:
		p.r= p.r+p.v*0.1*dt
		
		if p.r[0] > boxsize:
			p.r[0] -= boxsize
		if p.r[0] < 0:
			p.r[0] += boxsize
		if p.r[1] >boxsize:
			p.r[1] -= boxsize
		if p.r[1] < 0:
			p.r[1] += boxsize


	t=t+dt
