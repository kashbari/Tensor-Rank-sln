def qsort(L):
	if L  == []:
		return []
	else:
		p = L[0]
		l = qsort([x for x in L[1:] if x < p])
		g = qsort([x for x in L[1:] if x >= p])
		return l + [p] + g

def perm_parity(L):
	exp = -1
	for i in range(0,len(L)):
		for j in range(i+1, len(L)):
			if L[i]<L[j]:
				exp += 1
			else:
				exp = exp
	p = (-1)**exp
	return p

def perm(a,b):
	if qsort(a) == qsort(b):
		c = perm_parity(a)*perm_parity(b)
	else:
		c = 0
	return c

a = [1,2,3]
b = [1,2,2]

print(perm(a,b))
