import math
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

f1 = open("dane/data_spf1.txt", "r")
f2 = open("dane/data_spf2.txt", "r")
gx1 = open("dane/sp2_F_GX3_1.jpg", 'r')
gx2 = open("dane/sp2_F_GX3_2.jpg", 'r')

f1_data = [[],[],[],[],[],[],[]]
f2_data = [[],[],[],[],[],[],[]]
gx1_data = [[],[],[],[],[],[],[]]
gx2_data = [[],[],[],[],[],[],[]]

def GenAB(tab):
    #przyjmuje U jako tab[0] i a jako tab[1]
    A = tab[1][0]/(tab[0][0]-tab[0][1])+tab[1][1]/(tab[0][1]-tab[0][0])
    B = tab[1][0]*tab[0][1]/(tab[0][1]-tab[0][0])+tab[1][1]*tab[0][0]/(tab[0][0]-tab[0][1])
    return [A, B]
    
def FindIndex(tab, row, val):
    indx = 0
    for i in tab[row]:
        if(i<val):
            indx+=1
        else:
            break
    return int(indx)

def GenAverage(tab, start, end, row):
    st = FindIndex(tab, row[0], start)
    ed = FindIndex(tab, row[0], end)
    return np.average(tab[row[1]][st:ed])
    
def Integral(tab, start, end, row):
    st = FindIndex(tab, row[0], start)
    ed = FindIndex(tab, row[0], end)
    sum = 0
    for i in range(st, ed):
        sum+=0.5*(tab[row[1]][i]+tab[row[1]][i+1])*(tab[row[0]][i+1]-tab[row[0]][i])
    return sum    

for line in f1:
	ll = line.split("\t")
	for i in range(len(ll)):
		f1_data[i].append(float(ll[i]))

for line in f2:
	ll = line.split("\t")
	for i in range(len(ll)):
		f2_data[i].append(float(ll[i]))

for line in gx1:
	ll = line.split("\t")
	for i in range(len(ll)):
		gx1_data[i].append(float(ll[i]))
		
for line in gx2:
	ll = line.split("\t")
	for i in range(len(ll)):
		gx2_data[i].append(float(ll[i]))
		
f1.close()
f2.close()
gx1.close()
gx2.close()

t_scale_wait = 2.37#czas oczekiwania na pomiar skali

#skalowanie przyspieszenia
f1_file_scale = open("output/f1_skale.txt", 'w')
f1_file_scale.write("oś\tA\tB\tU_avg_min\tU_avg_max\n")
f1_a_x_scale = [[GenAverage(f1_data, 0, 2.37, [0,1]), GenAverage(f1_data, 4.5, 7.21, [0,1])], [-9.81, 9.81]]
ab = GenAB(f1_a_x_scale)
f1_file_scale.write("x\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}\n".format(ab[0], ab[1], GenAverage(f1_data, 0, 2.37, [0,1]), GenAverage(f1_data, 4.5, 7.21, [0,1])))
f1_a_y_scale = [[GenAverage(f1_data, 8.67, 10.40, [0,2]), GenAverage(f1_data, 12.84, 16.43, [0,2])], [-9.81, 9.81]]
ab = GenAB(f1_a_y_scale)
f1_file_scale.write("y\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}\n".format(ab[0], ab[1], GenAverage(f1_data, 8.67, 10.40, [0,2]), GenAverage(f1_data, 12.84, 16.43, [0,2])))
f1_a_z_scale = [[GenAverage(f1_data, 18.08, 19.81, [0,3]), GenAverage(f1_data, 23.24, 27.0, [0,3])], [-9.81, 9.81]]
ab = GenAB(f1_a_z_scale)
f1_file_scale.write("z\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}\n".format(ab[0], ab[1], GenAverage(f1_data, 18.08, 19.81, [0,3]), GenAverage(f1_data, 23.24, 27.0, [0,3])))
f1_file_scale.close()

#skalowanie giroskopów
f2_file_scale = open("output/f2_skala.txt", 'w')
f2_times = [[[0.0, 11.35], [0.0, 11.35], [14.78, 31.46]], [[39.17, 45.91], [19.0, 28.27], [1.88, 9.42]]]
f2_Uoi_w = []#[GenAverage(f2_data, f2_times[0][0], 11.35, [0,4]), GenAverage(f2_data, 0.0, 11.35, [0,5]), GenAverage(f2_data, 14.78, 31.46, [0,6])]
for i in range(len(f2_times[0])):
    f2_Uoi_w.append(GenAverage(f2_data, f2_times[0][i][0], f2_times[0][i][1], [0, i+4]))
f2_file_scale.write("U_ox: {0:.4f}\nU_oy: {1:.4f}\nU_oz: {2:.4f}\n".format(f2_Uoi_w[0], f2_Uoi_w[1], f2_Uoi_w[2]))
f2_rescaled_vals = []
for i in range(4):
    f2_rescaled_vals.append([])
    for o in range(len(f2_data[4])):
        if(i>0):
            f2_rescaled_vals[i].append(f2_data[i+3][o]-f2_Uoi_w[i-1])
        else:
            f2_rescaled_vals[i].append(f2_data[i][o])
        #print(f2_rescaled_vals[i][o])
f2_ki = []
for i in range(1, 4):
    print("{0}, {1}\n".format(f2_times[1][i-1][0], f2_times[1][i-1][1]))
    intg = Integral(f2_rescaled_vals, f2_times[1][i-1][0], f2_times[1][i-1][1], [0, i])
    print(intg)
    f2_ki.append(90.0/intg)
    f2_file_scale.write("k{0}: {1:.4f}\n".format(i, f2_ki[i-1]))
f2_file_scale.close()

#przyspieszenie
plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(f1_data[0], f1_data[1], label="x")
plt.plot(f1_data[0], f1_data[2], label="y")
plt.plot(f1_data[0], f1_data[3], label="z")
#plt.plot([f1_data[0][indx], f1_data[0][indx]], [-300, 300], label = "granica pomiaru do skalowania")
plt.grid(True)
plt.title("Wykres odczytanego napięcia na IMU-ZAIOL-01 w funkcji czasu (składowa \"$a$\").")
plt.ylabel("U [V]")
plt.xlabel("t [s]")
plt.ylim(1.85, 3)
plt.legend(loc = "lower right")
plt.savefig("output/IMU-ZAIOL-01_1_a_U_t.png")

#obrót
plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(f2_data[0], f2_data[4], label="oś x")
plt.plot(f2_data[0], f2_data[5], label="oś y")
plt.plot(f2_data[0], f2_data[6], label="oś z")
plt.grid(True)
plt.title("Wykres odczytanego napięcia na IMU-ZAIOL-01 w funkcji czasu (składowa \"$\omega$\").")
plt.ylabel(r"w $\left [ \frac{m}{s} \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/IMU-ZAIOL-01_2_U_w_t.png")
plt.show()