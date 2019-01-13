class Esq:
	def __init__(self, a1, a2, b1, b2):
		self.a1 = a1;
		self.a2 = a2;
		self.b1 = b1;
		self.b2 = b2;

def gfunc(ag):
	return lambda x1,x2: ag*(x1+x2)

ficflex = [
	lambda x1, x2: (1 - x1)*(1 - x2),
	lambda x1, x2: x1*(1 - x2),
	lambda x1, x2: x1*x2,
	lambda x1, x2: (1 - x1)*x2
]

def fi(i,e,x1,x2):
	return ficflex[i-1]((x1-e.b1)/e.a1,(x2-e.b2)/e.a2)

def df_jdx_i(e,j,i):
	dfi = [
		[-1/e.a1,1/e.a1,1/e.a1,-1/e.a1],
		[-1/e.a2,-1/e.a2,1/e.a2,1/e.a2]
	]
	return dfi[i-1][j-1]
	
def b_ij(e,i,j):
	res = 0
	for ix in range(1,3):
		res += df_jdx_i(e,i,ix)*df_jdx_i(e,j,ix)
	return e.a1*e.a2*res

def onNeumann(p1,p2):
	return p1 != (0,0) && p2 != (0,0)

def l_i(e,g,i):
	points = [(0,0),(1,0),(1,1),(0,1)]
	res = 0;
	for j in range(0,len(points)):
		p1 = points[j] - (e.b1, e.b2)
		p2 = points[(j+1)%len(points)] - (e.b1, e.b2)
		if(onNeumann(p1,p2)):
			ps = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2)
			res += g(p[0], p[1])*fi(e,i,p[0],p[1])*(p1.x==p2.x ? e.a2 : e.a1)
	return res	

E = [Esq(1,1,-1,0),Esq(1,1,0,0),Esq(1,1,0,-1)]
ag = float(input())
g = gfunc(ag)

fiToE = []

