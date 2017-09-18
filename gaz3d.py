# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import random
from matplotlib.patches import Circle

# WARTOSCI POCZATKOWE
particleNumber=4
boxsize=16.0
eps=1.0
sigma=1.0
promien=0.5
deltat=0.001
temp=0.1
b2=boxsize/2
b=boxsize
delta=2

tab_t = []

# KLASA CZASTKA
class czastka:
	def __init__(self,promien,pos,vel):
		self.promien = promien	#promien
		self.r=pos # polożenie
		self.v=vel # prędkosc

# FUNKCJA OBLICZAJACA SILE
def sila(i,particles):
	fv = np.zeros(2)
	for p in particles:
		r = particles[i].r-p.r
		if r[0] > b2: # b2 – połowa pudełka b2=b/2
			r[0] =r[0] -b # przesuwamy współrzędną x wektora r_vect
		elif r[0] <-b2:
			r[0] =r[0]+ b # b – bok pudełka
		if r[1] > b2: # to samo dla y
			r[1] =r[1] - b
		elif r[1] <-b2:
			r[1] =r[1] +b
		dlugosc_r=np.linalg.norm(r)
		if dlugosc_r > 0:
			f_temp = -r/dlugosc_r
			if dlugosc_r < 2.5:
					fv = fv + f_temp*(-48.*((1/dlugosc_r)**14-0.5*(1/dlugosc_r)**8))
	return fv

# PROCEDURA RYSUJACA
def rysuj():
	if (en%100==0): # co 100-na klatka
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
		nStr=str(en) #nagraj na dysk – numer pliku z 5 cyframi, na początku zera, np 00324.png
		nStr=nStr.rjust(5,'0')
		plt.title("Symulacja gazu Lennarda-Jonesa, krok "+nStr)
		plt.savefig('img'+nStr+'.png')


# Tworzenie czastek
nx = int(np.sqrt(particleNumber))
ny = int(np.sqrt(particleNumber))
particles = []
for i in range(nx):
	for j in range(ny):
		polozenie = np.array([i*delta+1, j*delta+1])
		predkosc = np.array([(random.random()-1./2),(random.random() -1./2)])
		particles.append(czastka(promien,polozenie,predkosc) )
# Unieruchomienie srodka masy i ustawienie temp
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

#warunki pocz
t=0
en = 0
tmax = 10

#Euler na start
vm05 = []
vp05 = []
for i in range(particleNumber):
	vm05.append([0,0])
	vp05.append([0,0])
for i in range(particleNumber):
	f = sila(i,particles)
	vm05[i] = particles[i].v-0.5*f*deltat
######

# Glowna petla programu
while t<tmax:
	# rysowanie
	rysuj()
	# Leapfrog
	for i in range(particleNumber):
		r = particles[i].r
		p = particles[i].v	# p=v
		f = sila(i,particles)
		vp05[i]=vm05[i]+f*deltat
		p=0.5*(vm05[i]+vp05[i])
		r=r+vp05[i]*deltat
		vm05[i]=vp05[i]
	 	particles[i].r = r
		particles[i].v = p
	#update położenia
	for p in particles:
		if p.r[0] > boxsize or p.r[0] < boxsize:
			p.r[0] = p.r[0]%boxsize
		if p.r[1] > boxsize or p.r[1] < boxsize:
			p.r[1] = p.r[1]%boxsize
	#temperatura
	sumv2 =0.0
	for p in particles:
		sumv2=sumv2+np.dot(p.v,p.v)/2.0
	sumv2=sumv2/particleNumber # średnia energia kinetyczna
	T=sumv2		# warto podzielic przez stala bolzmana	
	print T
	tab_t.append(T)
	#next step
	t=t+deltat
	en=en+1
#plt.show()
#plt.ylim("T")
#plt.plot(tab_t)
#plt.show()
