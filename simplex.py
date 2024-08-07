from math import inf
from sortedcontainers import SortedSet
import numpy as np

# assume we are in SEF
# max c^T x
# s.t. Ax = b, x >= 0

z = 0

# is infeasible
c = [1]
A = [[0]]
b = [1]

# has optimal solution
# c = [4, 3, 0, 0, 0, 0]
# A = [[1, 0, 1, 0, 0, 0],
# 	 [0, 1, 0, 1, 0, 0],
# 	 [1, 1, 0, 0, 1, 0],
# 	 [1, 1, 0, 0, 0, -1]]
# b = [2, 3, 4, 1]

# is unbounded
# c = [10, -2, -8, 4, -2, -1]
# A = [[4, 1, -2, -1, -2, 0],
# 	 [1, -1, 2, -4, 0, 2]]
# b = [-10, 4]

n = len(c)
m = len(b)

assert len(c) == n
assert len(A) == m
assert len(b) == m
for i in range(m):
	assert len(A[i]) == n

def get_cf(A, b, c, B):
	
	A_B = []

	for j in range(m):
		row = []
		for i in range(n + m):
			if i + 1 in B:
				row.append(A[j][i])
		A_B.append(row)
	
	A_B = np.array(A_B)

	A_B_inverse = np.linalg.inv(A_B)
	A_B_inverse_transpose = A_B_inverse.T
	c_B = []

	for i in range(len(c)):
		if i + 1 in B:
			c_B.append(c[i])
	c_B = np.array(c_B).T

	y = A_B_inverse_transpose @ c_B

	return [c - y @ A, y.T @ b, A_B_inverse @ A, A_B_inverse @ b, y]

def phase_one():
	slack_A = []
	flipped = set()
	for j in range(m):
		negative = b[j] < 0
		flipped.add(j)
		row = []
		for i in range(n):
			if negative:
				row.append(-A[j][i])
			else:
				row.append(A[j][i])
		
		# we must add m basic vars
		for _ in range(m):
			if _ == j:
				if negative:
					row.append(-1)
				else:
					row.append(1)
			else:
				row.append(0)

		slack_A.append(row)

	# initial feasible basis
	B = SortedSet(list(range(n + 1, n + m + 1)))
	slack_c = []

	for _ in range(n):
		slack_c.append(0)
	
	for _ in range(m):
		slack_c.append(-1)

	slack_A = np.array(slack_A)
	slack_c = np.array(slack_c)
	slack_b = np.array(b.copy())
	y = [0] * n

	[slack_c_temp, slack_z_temp, slack_A_temp, slack_b_temp, y] = get_cf(slack_A, slack_b, slack_c, B)
	# print('initial: ')
	# print(slack_z_temp)
	# print(slack_c_temp)
	# print(slack_A_temp)
	# print(slack_b_temp)
	# print()


	it = 1
	while True:
		# print(f'iteration {it}')
		# print(f'Basis is {B}')
		# print(slack_z_temp)
		# print(slack_c_temp)
		# print(slack_A_temp)
		# print(slack_b_temp)
		# print()

		k = -1
		for i, x in enumerate(slack_c_temp):
			if x > 0:
				k = i
				break

		if k == -1:
			break
		
		# print(f'x_{k + 1} should enter the basis')

		r = -1
		min_ratio = inf


		# aux LP should always have an optimal solution, no need to check for unbounded
		for i in range(m):
			if slack_A_temp[i][k] <= 0:
				continue

			ratio = slack_b_temp[i] / slack_A_temp[i][k]

			if ratio < min_ratio:
				min_ratio = slack_b_temp[i] / slack_A_temp[i][k]
				r = i
		r = B[r] - 1
		# print(f'x_{r + 1} should leave the basis')
		# print()
		
		B.remove(r + 1)
		B.add(k + 1)
		it += 1

		[slack_c_temp, slack_z_temp, slack_A_temp, slack_b_temp, y] = get_cf(slack_A, slack_b, slack_c, B)
	
	if slack_z_temp == 0:
		print('the aux problem has an optimal basis!')
		print(f'we can start phase 2 with basis {B}')
		return (B, True)
	else:
		print('LP IS INFEASIBLE')
		print(f'Certificate: y={y}')
		return (None, False)
		

def phase_two(B):
	y = [0] * n
	[c_temp, z_temp, A_temp, b_temp, y] = get_cf(A, b, c, B)

	it = 1
	

	while True:
		k = -1 

		for i, x in enumerate(c_temp):
			if x > 0:
				k = i
				break

		if k == -1:
			break
		
		# print(f'x_{k + 1} should enter the basis')

		r = -1
		min_ratio = inf

		found_positive = False
		for i in range(m):
			if A_temp[i][k] <= 0:
				continue

			found_positive = True

			ratio = b_temp[i] / A_temp[i][k]

			if ratio < min_ratio:
				min_ratio = b_temp[i] / A_temp[i][k]
				r = i

		if not found_positive:
			x_bar = [0] * n
			idx = 0
			for i in B:
				x_bar[i - 1] = b_temp[idx]
				idx += 1
			
			d = [0] * n
			d[k] = 1.0
			
			idx = 0
			for i in range(m):
				d[B[idx] - 1] = -A_temp[i][k]
				idx += 1

			print('LP IS UNBOUNDED')
			print(f'Certificate: x_bar={x_bar}, d={d}')
			return
		
		r = B[r] - 1
		# print(f'x_{r + 1} should leave the basis')
		# print()
		
		B.remove(r + 1)
		B.add(k + 1)
		it += 1

		[c_temp, z_temp, A_temp, b_temp, y] = get_cf(A, b, c, B)
	
	x_bar = [0] * (n + m)
	idx = 0
	for i in B:
		x_bar[i - 1] = b_temp[idx]
		idx += 1

	print('LP IS OPTIMAL')
	print(f'Certificate: y={y}')
	print(f'optimal solution: {x_bar}')
	print(f'optimal value: {z_temp}')

	

if __name__ == "__main__":
	feasible_basis, feasible = phase_one()
	if feasible:
		phase_two(feasible_basis)