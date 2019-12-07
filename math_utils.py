from math import pow, sqrt, exp

#import matplotlib.pyplot as plt

def average(l):
	if len(l) == 0:
		return 0
	total = 0
	for i in l:
		total += i
	return total / len(l)

def variance(l):
	a = average(l)
	return average([pow(x-a, 2) for x in l])

def std_dev(l):
	return sqrt(variance(l))

def net_change(l):
	return max(l) - min(l)

def sortIndex(l, i):
	return sorted(l, key=lambda x: x[i])

# ---- Statistical Analysis ---- #

# Corrected sample standard deviation
def sample_std_dev(l):
	avg = average(l)
	return sqrt(sum([(x-avg)**2 for x in l]) / (len(l)-1))

# Standard Error of the Mean
def std_error(l, name=""):
	# Calculate normalized correlation function
	C = []
	avg = average([x*x for x in l])
	if avg == 0:
		print("Average is 0")
		return std_dev(l) / sqrt(len(l))

	for dt in range(len(l)/2):
		correlate = [l[i] * l[i+dt] for i in range(len(l)-dt)]
		C.append([dt, average(correlate) / avg])
	# Search for (dt,C) pair closest to C=1/e
	tau = C[0][0]
	c_tau = C[0][1]
	for x in C:
		dt = x[0]
		c = x[1]
		if abs(c - exp(-1)) < abs(c_tau - exp(-1)):
			tau = dt
			c_tau = c
	# Compute std error
	print("Correlation time: " + str(tau))
	N = len(l) / tau
	print("\tN = " + str(N))

	#plt.plot([x[0] for x in C], [x[1] for x in C])
	#plt.xlabel("dt")
	#plt.ylabel("C(dt)")
	#plt.title("Correlation function (tau=" + str(tau) + ")")
	#plt.ylim(0, 1)
	#plt.hlines(exp(-1), C[0][0], C[-1][0])
	#plt.axvline(x=tau)
	#plt.savefig("crl_" + str(name) + ".png")
	#plt.clf()

	return std_dev(l) / sqrt(N)

# Relative Standard Error
# 	Expressed as decimal 0 <= rse <= 1
def rel_std_error(l):
	return std_error(l) / avg(l)
