from math import sqrt, acos, cos, sin, pi

def make_vec(point1, point2):
	return [point2[i] - point1[i] for i in range(len(point1))]

def norm(v):
	total = 0
	for x in v:
		total += x*x
	return sqrt(total)

def normalize(v):
	n = norm(v)
	return [x / n for x in v]

def dot(v1, v2):
	total = 0
	for i in range(len(v1)):
		total += v1[i] * v2[i]
	return total

def angle(v1, v2):
	try:
		d = dot(v1,v2)
		if d == 0:
			return 0
		return acos(d / (norm(v1) * norm(v2)))
	except ValueError:
		print("v1: " + str(v1))
		print("v2: " + str(v2))
		print("dot: " + str(dot(v1,v2)))
		print("n(v1): " + str(norm(v1)))
		print("n(v2): " + str(norm(v2)))
		print("arg: " + str(dot(v1,v2) / (norm(v1)*norm(v2))))
		exit()


def cross(v1, v2):
	out = [0,0,0]
	out[0] = v1[1] * v2[2] - v1[2] * v2[1]
	out[1] = v1[2] * v2[0] - v1[0] * v2[2]
	out[2] = v1[0] * v2[1] - v1[1] * v2[0]
	return out

def dihedral(v1, v2, v3):
	n1 = cross(v1, v2)
	n2 = cross(v2, v3)
	return angle(n1, n2)

# Multiply matrix times a vector
def matrix_mult(m, v):
	out = [0 for x in m]
	for i in range(len(m)):
		for j in range(len(m[i])):
			out[i] += m[i][j] * v[j]
	return out

# Multiply matrix times a matrix
def mat_mult(m1, m2):
	if (len(m1[0]) != len(m2)):
		print("ERROR: incorrect dimensions ({a}x{b}) * ({c}x{d})".format(a=len(m1), b=len(m1[0]), c=len(m2), d=len(m2[0])))
		return

	out = [[0 for x in range(len(m2[0]))] for y in range(len(m1))]
	for i in range(len(m1)):
		for j in range(len(m2[0])):
			for k in range(len(m2)):
				out[i][j] += m1[i][k] * m2[k][j]
	return out

def rotate(v, axis, angle):
	return matrix_mult(rotation_matrix(axis, angle), v)

def rotation_matrix(axis, angle):
	c = cos(angle)
	s = sin(angle)
	n = norm(axis)
	u0 = axis[0] / n
	u1 = axis[1] / n
	u2 = axis[2] / n
	rotmat = [[c+u0*u0*(1-c), u0*u1*(1-c)-u2*s, u0*u2*(1-c)+u1*s],
		  [u0*u1*(1-c)+u2*s, c+u1*u1*(1-c), u1*u2*(1-c)-u0*s],
		  [u0*u2*(1-c)-u1*s, u1*u2*(1-c)+u0*s, c+u2*u2*(1-c)]]
	return rotmat

def midpoint(v1, v2):
	return [(v1[i] + v2[i]) / 2 for i in range(len(v1))]

def vec_sum(vecs):
	total = [0 for i in vecs[0]]
	for v in vecs:
		for i in range(len(vecs[0])):
			total[i] += v[i]
	return total

def vec_add(v1, v2):
	return [v1[i] + v2[i] for i in range(len(v1))]

def vec_sub(v1, v2):
	return [v1[i] - v2[i] for i in range(len(v1))]

def scalar(a, v):
	return [a * x for x in v]

def vecToMatrix(v):
	return [[x] for x in v]

def transpose(matrix):
	return [[matrix[i][j] for i in range(len(matrix))] for j in range(len(matrix[0]))]
		

def mat_inverse(matrix, tol=0.0001):
	# LUPDecompose & LUPInverse from wikipedia LU_decomposition, adapted for python
	N = len(matrix)
	A = [[matrix[i][j] for j in range(N)] for i in range(N)]
	P = [i for i in range(N+1)]

	for i in range(N):
		maxA = 0
		imax = 1

		for k in range(i, N):
			absA = abs(A[k][i])
			if (abs > maxA) :
				maxA = absA
				imax = k

		if (maxA < tol):
			print("matrix is degenerate, cannot invert")
			printMatrix(matrix)
			return

		if (imax != i):
			j = P[i]
			P[i] = P[imax]
			P[imax] = j

			temp = A[i]
			A[i] = A[imax]
			A[imax] = temp

			P[N] += 1

		for j in range(i+1, N):
			A[j][i] /= A[i][i]

			for k in range(i+1, N):
				A[j][k] -= A[j][i] * A[i][k]
	
	IA = [[0 for j in range(N)] for i in range(N)]
	for j in range(N):
		for i in range(N):
			if (P[i] == j): IA[i][j] = 1
			else:		IA[i][j] = 0
			
			for k in range(i):
				IA[i][j] -= A[i][k] * IA[k][j]

		for i in reversed(range(N)):
			for k in range(i+1, N):
				IA[i][j] -= A[i][k] * IA[k][j]
			IA[i][j] = IA[i][j] / A[i][i]
	
	return IA

def resid(f, betas, data):
	r = [x[1] - f(x[0], betas) for x in data]
	return dot(r,r)

def leastSquaresReg(f, init_betas, data, partials, iterations=100):
	betas = [b for b in init_betas]
	for s in range(iterations):
		J = [[partial(x[0], betas) for partial in partials] for x in data]
		r = [x[1] - f(x[0], betas) for x in data]
		m = mat_inverse(mat_mult(transpose(J), J))
		if m is None:
			return betas
		m = mat_mult(m, transpose(J))
		d = matrix_mult(m, r)

		alpha = norm(d)	# maximum alpha value
		d = normalize(d)
		c = 0.5		# controls tolerance
		tau = 0.5	# controls rate of search
		gradS = vec_sum([[-2*r[i]*J[i][j] for j in range(len(partials))] for i in range(len(r))])
		t = -c * dot(d, gradS)
		S_old = resid(f, betas, data)
		S_new = resid(f, vec_add(betas, scalar(alpha, d)), data)
		while (S_old < S_new + alpha*t):
			alpha *= tau
			S_new = resid(f, vec_add(betas, scalar(alpha, d)), data)
		
		betas = vec_add(betas, scalar(alpha, d))
	return betas

def backtrack_line_search(f, x, p, gradf, alpha_max, c, tau):
	alpha = alpha_max
	t = -c * dot(p, gradf)
	f_old = f(x)
	f_new = f(vec_add(x, scalar(alpha, p)))
	print("\t\tf_old: " + str(f_old))
	print("\t\tf_new + at: " + str(f_new + alpha*t))
	while (f_old < f_new + alpha * t):
		alpha *= tau
		f_new = f(vec_add(x, scalar(alpha, p)))
		print("\t\tf_new + at: " + str(f_new + alpha*t))
	return alpha


def printMatrix(mat):
	for line in mat:
		out = ""
		for x in line:
			out += str(x) + " "
		print(out)
