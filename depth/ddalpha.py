
import numpy as np
from ctypes import *
import math
import scipy.spatial as scsp
import sklearn.covariance as sk
import scipy.special as scspecial
import sys
import os


if sys.platform=='linux' or sys.platform=='darwin'::
	for i in sys.path:
		ddalpha_exact=glob.glob(i+'/depth/'+'ddalpha.so')
		ddalpha_approx=glob.glob(i+'/depth/'+'ddalpha_wrapper.so')
	

	libr=CDLL(ddalpha_exact[0])
	libRom=CDLL(ddalpha_approx[0])
if sys.platform=='nt':
	site_packages = next(p for p in sys.path if 'site-packages' in p)
	print(site_packages)
	os.add_dll_directory(site_packages+"\depth")
	libr=CDLL(r""+site_packages+"\depth\ddalpha.dll")
	libRom=CDLL(r""+site_packages+"\depth\ddalpha_wrapper.dll")





def MCD_fun(data,alpha,NeedLoc=False):
    cov = sk.MinCovDet(support_fraction=alpha).fit(data)
    if NeedLoc:return([cov.covariance_,cov.location_])
    else:return(cov.covariance_)

    
def longtoint(k):
  limit = 2000000000
  k1 = int(k/limit)
  k2 = int(k - k1*limit)
  return np.array([k1,k2])



def halfspace(x, data,numDirections=1000,exact=True,method="recursive",
						solver = "neldermead",
						NRandom = 100,
						option = 1,
						n_refinements = 10,
						sphcap_shrink = 0.5,
						alpha_Dirichlet = 1.25,
						cooling_factor = 0.95,
						cap_size = 1,
						start = "mean",
						space = "sphere",
						line_solver = "goldensection",
						bound_gc = True):

	if exact:
		if (method =="recursive" or method==1):
			method=1
		elif (method =="plane" or method==2):
			method=2
		elif (method =="line" or method==3):
			method=3
		else:
			print("Wrong argument, method=str(recursive) or str(plane) or str(line)")
			print("recursive by default")
			method=3
		
		
		points_list=data.flatten()
		objects_list=x.flatten()
		points=(c_double*len(points_list))(*points_list)
		objects=(c_double*len(objects_list))(*objects_list)
		k=numDirections

		points=pointer(points)

		objects=pointer(objects)
		numPoints=pointer(c_int(len(data)))
		numObjects=pointer(c_int(len(x)))
		dimension=pointer(c_int(len(data[0])))
		algNo=pointer((c_int(method)))
		depths=pointer((c_double*len(x))(*np.zeros(len(x))))
	
		libr.HDepthEx(points,objects, numPoints,numObjects,dimension,algNo,depths)
	
		res=np.zeros(len(x))
		for i in range(len(x)):
			res[i]=depths[0][i]
		return res
	else:	
		notion="halfspace"
		res=depth_approximation(x,data,notion,solver ,NRandom ,option ,n_refinements ,
			sphcap_shrink ,alpha_Dirichlet ,cooling_factor,cap_size ,start ,space ,line_solver ,bound_gc )
	
	
		return res
	
	
def zonoid(x, data,seed=0,exact=True,solver = "neldermead",
						NRandom = 100,
						option = 1,
						n_refinements = 10,
						sphcap_shrink = 0.5,
						alpha_Dirichlet = 1.25,
						cooling_factor = 0.95,
						cap_size = 1,
						start = "mean",
						space = "sphere",
						line_solver = "goldensection",
						bound_gc = True):
	if exact:
		points_list=data.flatten()
		objects_list=x.flatten()
		points=(c_double*len(points_list))(*points_list)
		objects=(c_double*len(objects_list))(*objects_list)

		points=pointer(points)
		objects=pointer(objects)
		numPoints=pointer(c_int(len(data)))
		numObjects=pointer(c_int(len(x)))
		dimension=pointer(c_int(len(data[0])))
		seed=pointer((c_int(seed)))
		depths=pointer((c_double*len(x))(*np.zeros(len(x))))

		libr.ZDepth(points,objects, numPoints,numObjects,dimension,seed,depths)

		res=np.zeros(len(x))
		for i in range(len(x)):
			res[i]=depths[0][i]
		return res
	else:
		notion="zonoid"
		option=1
		res=depth_approximation(x,data,notion,solver ,NRandom ,option ,n_refinements ,
			sphcap_shrink ,alpha_Dirichlet ,cooling_factor,cap_size ,start ,space ,line_solver ,bound_gc )
	
	
		return res
		
		
		
def mahalanobis(x, data,exact="True",mah_estimate="moment",mah_parMcd = 0.75, 
						method="recursive",
						solver = "neldermead",
						NRandom = 100,
						option = 1,
						n_refinements = 10,
						sphcap_shrink = 0.5,
						alpha_Dirichlet = 1.25,
						cooling_factor = 0.95,
						cap_size = 1,
						start = "mean",
						space = "sphere",
						line_solver = "goldensection",
						bound_gc = True):
						
						
	if exact:
	

		points_list=data.flatten()
		objects_list=x.flatten()
		
		points=(c_double*len(points_list))(*points_list)
		objects=(c_double*len(objects_list))(*objects_list)

		points=pointer(points)
		objects=pointer(objects)
		numPoints=pointer(c_int(len(data)))
		numObjects=pointer(c_int(len(x)))
		dimension=pointer(c_int(len(data[0])))
		PY_MatMCD=MCD_fun(data,mah_parMcd)
		PY_MatMCD=PY_MatMCD.flatten()
		mat_MCD=pointer((c_double*len(PY_MatMCD))(*PY_MatMCD))

		
		
		depths=pointer((c_double*len(x))(*np.zeros(len(x))))

		libr.MahalanobisDepth(points,objects,numPoints,numObjects,dimension,mat_MCD,depths)

		res=np.zeros(len(x))
		for i in range(len(x)):
			res[i]=depths[0][i]
		return res
	else:
	
		notion="mahalanobis"
		option=1
		res=depth_approximation(x,data,notion,solver ,NRandom ,option ,n_refinements ,
			sphcap_shrink ,alpha_Dirichlet ,cooling_factor,cap_size ,start ,space ,line_solver ,bound_gc )
	
	
		return res
	
	
	
def simplical(x, data,exact=1,k=0.05,seed=0):
	points_list=data.flatten()
	objects_list=x.flatten()
	points=(c_double*len(points_list))(*points_list)
	objects=(c_double*len(objects_list))(*objects_list)
	points=pointer(points)
	objects=pointer(objects)
	
	numPoints=pointer(c_int(len(data)))
	numObjects=pointer(c_int(len(x)))
	dimension=pointer(c_int(len(data[0])))
	seed=pointer((c_int(seed)))
	exact=pointer((c_int(exact)))#exact=1
	if k<=0:
		print("k must be positive")
		print("k=1")
		k=scspecial.comb(len(data),len(data[0]),exact=True)*k
		k=pointer((c_int*2)(*longtoint(k)))
	elif k<=1:
		k=scspecial.comb(len(data),len(data[0]),exact=True)*k
		k=pointer((c_int*2)(*longtoint(k)))
	else:
		k=pointer((c_int*2)(*longtoint(k)))
		
	
	depths=pointer((c_double*len(x))(*np.zeros(len(x))))


	libr.SimplicialDepth(points,objects, numPoints,numObjects,dimension,seed,exact,k,depths)


	res=np.zeros(len(x))
	for i in range(len(x)):
		res[i]=depths[0][i]
	return res
	

	
	
def potential(x, data, pretransform = "1Mom", kernel="EDKernel" ,mah_parMcd=0.75):

	if(kernel=="GKernel" or kernel==2):
		kernel=2
	elif(kernel=="EKernel" or kernel==3):
		kernel=3
	elif(kernel=="TriangleKernel" or kernel ==4):
		kernel=4
	else:
		kernel = 1
	
            
            
            
            
	if (pretransform == "1Mom" or pretransform == "NMom"):
		[mu,B_inv,cov]=Maha_moment(data)
	elif (pretransform == "1MCD" or pretransform == "NMCD"):
		[mu,B_inv,cov]=Maha_mcd(data, mah_parMcd)
	data=Maha_transform(data,mu,B_inv)
	x =Maha_transform(x,mu,B_inv)

	points_list=data.flatten()
	objects_list=x.flatten()
	points=(c_double*len(points_list))(*points_list)
	objects=(c_double*len(objects_list))(*objects_list)
	points=pointer(points)
	points2=pointer(objects)
	
	numPoints=pointer(c_int(len(data)))
	numpoints2=pointer(c_int(len(x)))
	dimension=pointer(c_int(len(data[0])))
	
	KernelType=pointer(c_int(kernel))
	ignoreself=pointer(c_int(0))
	classes=pointer((c_int(1)))
	kernel_bandwidth=pointer(c_double(math.pow(len(data),-2/(len(data[0])+4))))
	depth=pointer((c_double*len(x))(*np.zeros(len(x))))

	libr.PotentialDepthsCount(points,numPoints,dimension,classes,numPoints,points2,numpoints2,KernelType,kernel_bandwidth,ignoreself,depth)
	res=np.zeros(len(x))
	for i in range(len(x)):
		res[i]=depth[0][i]
	return res
			
			

	

def Maha_moment (x):
	x=np.transpose(x)
	mu =np.mean(x,axis=1)
	cov=np.cov(x)
	w,v=np.linalg.eig(cov)
	B_inv=np.linalg.inv(np.matmul(v,np.diag(np.sqrt(w))))
	return ([mu,B_inv,cov])


def Maha_mcd(x, alpha =0.5):
	[cov,mu] = MCD_fun(x,alpha,1)
	w,v=np.linalg.eig(cov)
	B_inv=np.linalg.inv(np.matmul(v,np.diag(np.sqrt(w))))
	return ([mu,B_inv,cov])


def Maha_transform (x, mu, B_inv): 
	return(np.transpose(np.matmul(B_inv,np.transpose(x-mu))))


    






def count_convexes(objects,points,cardinalities, seed = 0):

	tmp_x=points.flatten()
	tmp_x=pointer((c_double*len(tmp_x))(*tmp_x))
	dimension=pointer(c_int(len(points[0])))
	numClasses=pointer(c_int(1))
	tmp_objects=objects.flatten()
	tmp_objects=pointer((c_double*len(tmp_objects))(*tmp_objects))
	PY_numObjects=len(objects)
	numObjects=pointer(c_int(PY_numObjects))
	tmp_cardinalities=pointer(c_int(cardinalities))
	seed=pointer(c_int(seed))
	length=PY_numObjects*1
	init_zeros=np.zeros(length,dtype=int)
	isInConv=pointer((c_int*length)(*init_zeros))

	
	libr.IsInConvexes(tmp_x,dimension,tmp_cardinalities,numClasses,tmp_objects,numObjects,seed,isInConv)
	res=np.zeros(length)
	for i in range(length):
		res[i]=isInConv[0][i]
	res.reshape(PY_numObjects,1)
	return res
  	
  			
def is_in_convex(x, data, cardinalities, seed = 0):
	res=count_convexes(x, data, cardinalities, seed)
	return res 
	
def qhpeeling(x, data):
	points_list=data.flatten()
	objects_list=x.flatten()
	nrow_data=len(data)
	depths=np.zeros(len(x))
	tmpData=data
	for i in range(nrow_data):
		if (len(tmpData)<(len(data[0])*(len(data[0])+1)+0.5)):
			break
		tmp=is_in_convex(x,tmpData,len(tmpData))
		depths+=tmp
		tmp_conv=scsp.ConvexHull(tmpData)
		tmpData=np.delete(tmpData,np.unique(np.array(tmp_conv.simplices)),0)
	depths=depths/nrow_data
	return depths

def simplicalVolume(x,data,exact=0,k=0.05,mah_estimate="moment", mah_parMCD=0.75,seed=0):
	points_list=data.flatten()
	objects_list=x.flatten()
	if (mah_estimate == "none"):
		useCov = 0
		covEst =np.eye(len(data[0])).flatten()
	elif (mah_estimate == "moment"):
		useCov = 1
		covEst=np.cov(np.transpose(data))
    
	elif (mah_estimate == "MCD") :
		useCov = 2
		covEst = MCD_fun(data, mah_parMCD)
	else:
		print("Wrong argument \"mah.estimate\", should be one of \"moment\", \"MCD\", \"none\"")
		print("moment is use")
		useCov = 1
		covEst=np.cov(data)
        
	points=(c_double*len(points_list))(*points_list)
	objects=(c_double*len(objects_list))(*objects_list)

	points=pointer(points)
	objects=pointer(objects)
	numPoints=pointer(c_int(len(data)))
	numObjects=pointer(c_int(len(x)))
	dimension=pointer(c_int(len(data[0])))
	
	seed=pointer((c_int(seed)))
	exact=pointer((c_int(exact)))
	if k<=0:
		print("k must be positive")
		print("k=1")
		k=scspecial.comb(len(data),len(data[0]),exact=True)*k
		k=pointer((c_int*2)(*longtoint(k)))
	elif k<=1:
		k=scspecial.comb(len(data),len(data[0]),exact=True)*k
		k1=k
		k=pointer((c_int*2)(*longtoint(k)))
	else:
		k=pointer((c_int*2)(*longtoint(k)))
		
	
	
	useCov=pointer(c_int(useCov))
	covEst=covEst.flatten()
	covEst=pointer((c_double*len(covEst))(*covEst))
        
	depths=pointer((c_double*len(x))(*np.zeros(len(x))))

	libr.OjaDepth(points,objects,numPoints,numObjects,dimension,seed, exact, k, useCov, covEst, depths)

	res=np.zeros(len(x))
	for i in range(len(x)):
		res[i]=depths[0][i]
	return res



def L2(x, data,mah_estimate='moment',mah_parMcd=0.75):
	points_list=data.flatten()
	objects_list=x.flatten()
	
	if mah_estimate=='none':
		sigma=np.eye(len(data[0]))
	else:
		if mah_estimate=='moment':
			cov=np.cov(np.transpose(data))
		elif mah_estimate=='MCD':
			cov=MCD_fun(data, mah_parMcd)
		else :
			print("Wrong argument \"mah.estimate\", should be one of \"moment\", \"MCD\", \"none\"")
			print("moment is used")
			cov=np.cov(np.transpose(data))
			
		if np.sum(np.isnan(cov))==0:
			sigma=np.linalg.inv(cov)
		else:
			print("Covariance estimate not found, no affine-invariance-adjustment")
			sigma=np.eye(len(data))
	
	depths=(-1)*np.ones(len(x))
	for i in range(len(x)):
		tmp1=(x[i]-data)
		tmp2=np.matmul(tmp1,sigma)
		tmp3=np.sum(tmp2 * tmp1,axis=1)
		depths[i]=1/(1 + np.mean(np.sqrt(tmp3)))
	return depths


def BetaSkeleton(x, data, beta = 2, distance = "Lp", Lp_p = 2, mah_estimate = "moment", mah_parMcd = 0.75):
	points_list=data.flatten()
	objects_list=x.flatten()
	if (distance == "Mahalanobis"):
		code = 5
		if (mah_estimate == "none"):
			sigma = np.eye(len(data[0]))
		else:
			if(mah_estimate == "moment"):
				cov = np.cov(np.transpose(data))
			elif (mah_estimate == "MCD"):
				cov = MCD_fun(data, mah_parMcd)
			else:
				print("Wrong argument \"mah_estimate\", should be one of \"moment\", \"MCD\", \"none\"")
			
			if (np.sum(np.isnan(cov)) == 0):
				sigma = np.linalg.inv(cov)
			else:
				sigma = np.eye(len(data[0]))
				print("Covariance estimate not found, no affine-invariance-adjustment")
	else:
		sigma = np.zeros(1)
		if (distance== "Lp"):
			code=4
			if (Lp_p == 1):
				code=1
			if (Lp_p == 2):
				code = 2
			if (Lp_p==math.inf and Lp_p > 0):
				code = 3
		else:stop("Argument \"distance\" should be either \"Lp\" or \"Mahalanobis\"")

	points=pointer((c_double*len(points_list))(*points_list))
	objects=pointer((c_double*len(objects_list))(*objects_list))
	numPoints=pointer(c_int(len(data)))
	numObjects=pointer(c_int(len(x)))
	dimension=pointer(c_int(len(data[0])))
	beta=[beta]
	
	beta=pointer((c_double*1)(*beta))
	code=pointer(c_int(code))
	Lp_p=[Lp_p]
	Lp_p=pointer((c_double*1)(*Lp_p))
	sigma=pointer((c_double*len(sigma.flatten()))(*sigma.flatten()))
	depth=pointer((c_double*len(x))(*np.zeros(len(x))))

	
	
	libr.BetaSkeletonDepth(points, objects, numPoints, numObjects, dimension, beta, code, Lp_p, sigma, depth)
    	
	res=np.zeros(len(x))
	for i in range(len(x)):
		res[i]=depth[0][i]
	return res



def spatial(x, data,mah_estimate='moment',mah_parMcd=0.75):
        depths_tab=[]

        if mah_estimate=='none':
                print('none')
                lambda1=np.eye(len(data))
        elif mah_estimate=='moment':
                print('moment')
                cov=np.cov(np.transpose(data))
        elif mah_estimate=='MCD':
                print('mcd')
        if np.sum(np.isnan(cov))==0:
                w,v=np.linalg.eig(cov)
                lambda1=np.linalg.inv(np.matmul(v,np.diag(np.sqrt(w))))#invàconfirmer
        else:
                lambda1=np.eye(len(data))

        depths=np.repeat(-1,len(x),axis=0)
        for i in range(len(x)):
                interm=[]
                tmp1_ter=np.transpose(x[i]-data)
                tmp1=np.transpose(np.matmul(lambda1,tmp1_ter))
                tmp1_bis=np.sum(tmp1,axis=1)
                for elements in tmp1_bis:
                        if elements==0:
                                interm.append(False)
                        if elements!=0:
                                interm.append(True)
                
                interm=np.array(interm)
                tmp1=tmp1[interm]
                tmp2=1/np.sqrt(np.sum(np.power(tmp1,2),axis=1))
                tmp3=np.zeros([len(tmp1),len(tmp1[0])])
                tmp1=np.transpose(tmp1)
                for jj in range(len(tmp1)):
                        tmp3[:,jj]=tmp2*(tmp1[:][jj])
                tmp4=np.sum(tmp3,axis=0)/len(data)
                tmp5=np.power((tmp4),2)
                tmp6=np.sum(tmp5)
                depths_tab.append(1-np.sqrt(tmp6))
        return depths_tab





def projection(x, data,sym, 
						solver = "neldermead",
						NRandom = 100,
						option = 1,
						n_refinements = 10,
						sphcap_shrink = 0.5,
						alpha_Dirichlet = 1.25,
						cooling_factor = 0.95,
						cap_size = 1,
						start = "mean",
						space = "sphere",
						line_solver = "goldensection",
						bound_gc = True):

	if((sym=="projection") or (sym=="aprojection")):
		option=1
		res=depth_approximation(x,data,sym,solver ,NRandom ,option ,n_refinements ,
			sphcap_shrink ,alpha_Dirichlet ,cooling_factor,cap_size ,start ,space ,line_solver ,bound_gc )
		return res
	
	else:
		return "Wrong argument: sym=\"projection\" or sym=\"aprojection\" " 





def get_ext_filename(): # Return suffix from shared library based on OS 
    from distutils.sysconfig import get_config_var
    ext_suffix = get_config_var('EXT_SUFFIX')
    return ext_suffix

						
						
						
def depth_approximation(z,
						X,
						notion = "halfspace",
						solver = "neldermead",
						NRandom = 100,
						option = 1,
						n_refinements = 10,
						sphcap_shrink = 0.5,
						alpha_Dirichlet = 1.25,
						cooling_factor = 0.95,
						cap_size = 1,
						start = "mean",
						space = "sphere",
						line_solver = "goldensection",
						bound_gc = True):
	"""
	depth_approximation(z,
						X,
						notion = "halfspace",
						solver = "neldermead",
						NRandom = 100,
						option = 1,
						n_refinements = 10,
						sphcap_shrink = 0.5,
						alpha_Dirichlet = 1.25,
						cooling_factor = 0.95,
						cap_size = 1,
						start = "mean",
						space = "sphere",
						line_solver = "goldensection",
						bound_gc = True)

    Compute data depth approximation based on the weak projection property.

    Parameters
    ----------
    z : array_like
        Points whose depth is to be calculated, each row contains a d-variate point.
        Should have the same dimension as `X`.
    X : array_like
        Data where each row contains a d-variate point, w.r.t. which the depth is to be calculated.
    notion : {'halfspace', 'mahalanobis', 'zonoid', 'projection', 'aprojection'}, optional
        Which depth will be computed.
    solver : {'simplegrid', 'refinedgrid', 'simplerandom', 'refinedrandom', 'coordinatedescent', 'randomsimplices', 'neldermead', 'simulatedannealing'}, optional
        The type of solver used to approximate the depth.
    NRandom : int, optional
    	The total number of iterations to compute the depth. Some solvers are converging
    	faster so they are run several time to achieve `NRandom` iterations.
    option : int, optional
        If `option` = 1, only approximated depths are returned.
        If `option` = 2, depths calculated at every iteration are also returned.
        If `option` = 3, random directions used to project depths are also returned
        with indices of converging for the solver selected.
    n_refinements : int, optional
        For `solver` = `refinedrandom` or `refinedgrid`, set the maximum of iteration for 
        computing the depth of one point.
    sphcap_shrink : float, optional
        For `solver` = `refinedrandom` or `refinedgrid`, it's the shrinking of the spherical cap.
    alpha_Dirichlet : float, optional
        For `solver` = `randomsimplices`. it's the parameter of the Dirichlet distribution.
    cooling_factor : float, optional
        For `solver` = `randomsimplices`, it's the cooling factor.
    cap_size : float, optional
        For `solver` = `simulatedannealing` or `neldermead`, it's the size of the spherical cap.
    start : {'mean', 'random'}, optional
        For `solver` = `simulatedannealing` or `neldermead`, it's the method used to compute the first depth.
    space : {'sphere', 'euclidean'}, optional
        For `solver` = `coordinatedescent` or `neldermead`, it's the type of spacecin which
        the solver is running.
    line_solver : {'uniform', 'goldensection'}, optional
        For `solver` = `coordinatedescent`, it's the line searh strategy used by this solver.
    bound_gc : bool, optional
        For `solver` = `neldermead`, it's ``True`` if the search is limited to the closed hemisphere.

    Returns
    -------
    type
        Depending on `option`, additional results are returned.
    depths : array_like
        The approximate depth for the points `z`.
    depths_iter : array_like
        Every iteration of depths for each point in `z`.
        The minimum value is returned in `depths`.
    directions : array_like
        Every direction used to approximate depths
    ind_convergence : array_like
    	The `solver` can be run several times and `ind_convergence` returns
    	the indices for which he converged.

    Raises
    ------
    ValueError
        Because you didn't choose a valid parameter for the function.

    References
    ----------
    Cite the relevant literature, e.g. [1]_.  You may also cite these
    references in the notes section above.

    .. [1]
	    O. Rainer Dyckerhoff, Pavlo Mozharovskyi, Stanislav Nagy.
	    Approximate computation of projection depths.
	    Computational Statistics and Data Analysis, Elsevier, In press, 157, pp.107166.
	    <https://doi.org/10.1016/j.csda.2020.107166>

    Examples
    --------
    >>> import numpy as np
    >>> np.seed(1)
    >>> n = 100
    >>> d = 3
    >>> mean = np.zeros(d)
    >>> cov = np.eye(d)
    >>> X = np.random.multivariate_normal(mean, cov, n)
    >>> z = np.random.multivariate_normal(mean, cov, 20)
    >>> depths, depths_iter, directions, ind_convergence = depth.depth_approximation(z, X,
    ... notion = "halfspace", solver = "neldermead", NRandom = 100, option = 3, cap_size = 1,
    ... start = "mean", space = "sphere", bound_gc = True)
    >>> print(depths)
    [0.13 0.   0.   0.   0.   0.02 0.1  0.01 0.2  0.07
    0.04 0.07 0.   0.16 0.11 0.4  0.03 0.05 0.01 0.01]

	"""

	depth_indice = check_depth(notion)
	check_space(space)
	solver_indice = check_solver(solver, space)
	start_indice = check_start(start)
	line_solver_indice = check_line_solver(line_solver)
	check_bound(bound_gc)

	try:
		n, d = X.shape
	except ValueError:
		n = X.shape[0]
		d = 1
	
	n_z = z.shape[0]
	depths = np.empty(n_z, dtype=np.double)
	depths_iter = np.empty((n_z, NRandom), dtype=np.double)
	directions = np.full(((n_z, NRandom, d)), -1, dtype=np.double)
	directions_card = np.full(((n_z, NRandom)), -1, dtype=np.int32)	# Keep the size of every split of directions which converged
																	# Initialization with at most NRandom split of directions
#return_flag = 

	libRom.depth_approximation(
		c_void_p(z.ctypes.data),
		c_void_p(X.ctypes.data),
		c_int(depth_indice),
		c_int(solver_indice),
		c_int(NRandom),
		c_int(option),
		c_int(n_refinements),
		c_double(sphcap_shrink),
		c_double(alpha_Dirichlet),
		c_double(cooling_factor),
		c_double(cap_size),
		c_int(start_indice),
		c_int(line_solver_indice),
		c_int(bound_gc),
		c_int(n),
		c_int(d),
		c_int(n_z),
		c_void_p(depths.ctypes.data),
		c_void_p(depths_iter.ctypes.data),
		c_void_p(directions.ctypes.data),
		c_void_p(directions_card.ctypes.data)
		)
	
	# Resize and clear array of every directions used
	directions = directions.tolist()
	for i in range(n_z):
		for j in range(NRandom):
			if(directions[i][j].count(-1) != 0):
				directions[i] = directions[i][:j] # Clear -1 values
				break
	
	# Fill indices for every start of convergence
	ind_convergence = []
	for i in range(n_z):
		if(solver == "refinedgrid" or solver == "refinedrandom"): # Return indices of refinements step
			ind_convergence = np.arange(0, NRandom, NRandom//n_refinements)[:ceil(len(directions[0])/(NRandom/n_refinements))].tolist()
		else:
			ind_bin = directions_card[i, ~(directions_card[i] == -1)] # Clear every -1 value 
			ind_bin_cumsum = np.cumsum(ind_bin)
			ind_convergence.append((ind_bin_cumsum - ind_bin).tolist())


	if(option == 3):return depths, depths_iter, directions, ind_convergence
	elif(option == 2):return depths, depths_iter
	else:return depths



def check_depth(depth):
	all_depths = ["mahalanobis", "halfspace", "zonoid", "projection", "aprojection"]
	if (depth not in all_depths):
		raise ValueError("Depths approximation is available only for depths in %s, got %s."%(all_depths, depth))
	else:
		return all_depths.index(depth)

def check_solver(solver, space):
	all_solvers = ["simplegrid", "refinedgrid", "simplerandom", "refinedrandom",
				"coordinatedescent", "randomsimplices", "neldermead", "simulatedannealing"]
	if solver not in all_solvers:
		raise ValueError("Depths approximation supports only solvers in %s, got %s."%(all_solvers, solver))
	else:
		if(solver == "coordinatedescent" and space == "sphere"):return 8 # Indice of the solver in ProjectionDepths
		elif(solver == "coordinatedescent" and space == "euclidean"):return all_solvers.index("coordinatedescent")
		elif(solver == "neldermead" and space == "sphere"):return 9 # Indice of the solver in ProjectionDepths
		elif(solver == "neldermead" and space == "euclidean"):return all_solvers.index("neldermead")
		else:
			return all_solvers.index(solver)

def check_start(start):
	all_start = ["mean", "random"]
	if (start not in all_start):
		raise ValueError("Only start available are in %s, got %s."%(all_start, start))
	else:
		return all_start.index(start)

def check_space(space):
	all_space = ["sphere", "euclidean"]
	if (space not in all_space):
		raise ValueError("Only space available are in %s, got %s."%(all_space, space))

def check_line_solver(line_solver):
	all_line_solver = ["uniform", "goldensection"]
	if (line_solver not in all_line_solver):
		raise ValueError("Only line_solver available are in %s, got %s."%(all_line_solver, line_solver))
	else:
		return all_line_solver.index(line_solver)

def check_bound(bound):
	all_bound = [True, False]
	if (bound not in all_bound):
		raise ValueError("Only bound option available are in %r, got %r."%(all_bound, bound))

if __name__ == "__main__":
    print('main')
    
    
    









halfspace.__doc__= """

Description
	Calculates the exact or random Tukey (=halfspace, location) depth (Tukey, 1975) of points w.r.t. a
	multivariate data set.

Usage
	depth.halfspace(x, data, exact, method, num.directions = 1000, seed = 0)

Arguments
x 			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.

data 			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.

exact			The type of the used method. The default is exact=F, which leads to approx-
			imate computation of the Tukey depth. For exact=F, method="Sunif.1D"
			is used by default. If exact=T, the Tukey depth is computed exactly, with
			method="recursive" by default.

method			For exact=F, if method="Sunif.1D" (by default), the Tukey depth is computed
			approximately by being minimized over univariate projections (see Details be-
			low).
			For exact=T, the Tukey depth is calculated as the minimum over all combina-
			tions of k points from data (see Details below). In this case parameter method
			specifies k, with possible values 1 for method="recursive" (by default), d − 2
			for method="plane", d − 1 for method="line".
			The name of the method may be given as well as just parameter exact, in which
			case the default method will be used.

num.directions 		Number of random directions to be generated (for method="Sunif.1D"). The
			algorithmic complexity is linear in the number of observations in data, given
			the number of directions.

seed			The random seed. The default value seed=0 makes no changes (for method="Sunif.1D").


"""
	
	
	
	
zonoid.__doc__= """

Description
	Calculates the zonoid depth of points w.r.t. a multivariate data set.

Usage
	depth.zonoid(x, data, seed = 0)

Arguments
x 			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.

data 			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.

seed 			the random seed. The default value seed=0 makes no changes.

"""


mahalanobis.__doc__= """

Description
	Calculates the Mahalanobis depth of points w.r.t. a multivariate data set.

Usage
	depth.Mahalanobis(x, data, mah.estimate = "moment", mah.parMcd = 0.75)

Arguments
x 			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.


data 			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.

mah.estimate 		is a character string specifying which estimates to use when calculating the Ma-
			halanobis depth; can be "moment" or "MCD", determining whether traditional
			moment or Minimum Covariance Determinant (MCD) (see covMcd) estimates
			for mean and covariance are used. By default "moment" is used.

mah.parMcd		is the value of the argument alpha for the function covMcd; is used when
			mah.estimate = "MCD".

"""




simplical.__doc__ = """

Description
	Calculates the simplicial depth of points w.r.t. a multivariate data set.

Usage
	depth.simplicial(x, data, exact = F, k = 0.05, seed = 0)

Arguments
x 			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.

data 			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.

exact 			exact=F (by default) implies the approximative algorithm, considering k sim-
			plices, exact=T implies the exact algorithm.

k 			Number (k > 1) or portion (if 0 < k < 1) of simplices that are considered if
			exact=F. If k > 1, then the algorithmic complexity is polynomial in d but is
			independent of the number of observations in data, given k. If 0 < k < 1,
			then the algorithmic complexity is exponential in the number of observations in
			data, but the calculation precision stays approximately the same.

seed 			the random seed. The default value seed=0 makes no changes.

"""



potential.__doc__="""

Description
	Calculate the potential of the points w.r.t. a multivariate data set. The potential is the kernel-
	estimated density multiplied by the prior probability of a class. Different from the data depths, a
	density estimate measures at a given point how much mass is located around it.

Usage
	depth.potential (x, data, pretransform = "1Mom",
	kernel = "GKernel", kernel.bandwidth = NULL, mah.parMcd = 0.75)

Arguments
x 			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.

data 			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.

pretransform 		The method of data scaling.
			NULL to use the original data,
			1Mom or NMom for scaling using data moments,
			1MCD or NMCD for scaling using robust data moments (Minimum Covariance De-
			terminant (MCD) ).

kernel			"EDKernel" for the kernel of type 1/(1+kernel.bandwidth*EuclidianDistance2(x,
			y)),
			"GKernel" [default and recommended] for the simple Gaussian kernel,
			"EKernel" exponential kernel: exp(-kernel.bandwidth*EuclidianDistance(x, y)),
			"VarGKernel" variable Gaussian kernel, where kernel.bandwidth is propor-
			tional to the depth.zonoid of a point.

kernel.bandwidth	the single bandwidth parameter of the kernel. If NULL - the Scott’s rule of thumb
			is used.

mah.parMcd		is the value of the argument alpha for the function covMcd; is used when
			pretransform = "*MCD".


"""
	
	
	
	

is_in_convex.__doc__= """

Description
	Checks the belonging to at least one of class convex hulls of the training sample.
Usage

	is.in.convex(x, data, cardinalities, seed = 0)

Arguments
x 			Matrix of objects (numerical vector as one object) whose belonging to convex
			hulls is to be checked; each row contains a d-variate point. Should have the
			same dimension as data.

data 			Matrix containing training sample where each row is a d-dimensional object,
			and objects of each class are kept together so that the matrix can be thought of
			as containing blocks of objects, representing classes.

cardinalities 		Numerical vector of cardinalities of each class in data, each entry corresponds
			to one class.

seed 			the random seed. The default value seed=0 makes no changes.

"""



qhpeeling.__doc__= """

Description
	Calculates the convex hull peeling depth of points w.r.t. a multivariate data set.

Usage
	depth.qhpeeling(x, data)

Arguments
x 			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.

data 			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.



"""


simplicalVolume.__doc__="""

Description
	Calculates the simpicial volume depth of points w.r.t. a multivariate data set.

Usage
	depth.simplicialVolume(x, data, exact = F, k = 0.05, mah.estimate = "moment",
	mah.parMcd = 0.75, seed = 0)

Arguments
x 			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.

data 			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.

exact			exact=F (by default) implies the approximative algorithm, considering k sim-
			plices, exact=T implies the exact algorithm.

k			Number (k > 1) or portion (if 0 < k < 1) of simplices that are considered if
			exact=F. If k > 1, then the algorithmic complexity is polynomial in d but is
			independent of the number of observations in data, given k. If 0 < k < 1,
			then the algorithmic complexity is exponential in the number of observations in
			data, but the calculation precision stays approximately the same.

mah.estimate 		A character string specifying affine-invariance adjustment; can be "none", "moment"
			or "MCD", determining whether no affine-invariance adjustemt or moment or
			Minimum Covariance Determinant (MCD) (see covMcd) estimates of the co-
			variance are used. By default "moment" is used.

mah.parMcd 		The value of the argument alpha for the function covMcd; is used when, mah.estimate = "MCD".


seed 			The random seed. The default value seed=0 makes no changes.

"""


L2.__doc__=""" 

Description
	Calculates the L2-depth of points w.r.t. a multivariate data set.

Usage
	depth.L2(x, data, mah.estimate = "moment", mah.parMcd = 0.75)

Arguments

x 			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.

data 			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.

mah.estimate 		is a character string specifying which estimates to use when calculating sample
			covariance matrix; can be "none", "moment" or "MCD", determining whether
			traditional moment or Minimum Covariance Determinant (MCD) (see covMcd)
			estimates for mean and covariance are used. By default "moment" is used. With
			"none" the non-affine invariant version of the L2-depth is calculated	

mah.parMcd		is the value of the argument alpha for the function covMcd; is used when
			mah.estimate = "MCD".

"""



BetaSkeleton.__doc__= """ 

Description
	Calculates the beta-skeleton depth of points w.r.t. a multivariate data set.

Usage
	depth.betaSkeleton(x, data, beta = 2, distance = "Lp", Lp.p = 2,
	mah.estimate = "moment", mah.parMcd = 0.75)

Arguments
x			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.

data 			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.

beta 			The paremeter defining the positionning of the balls’ centers, see Yang and
			Modarres (2017) for details. By default (together with other arguments) equals
			2, which corresponds to the lens depth, see Liu and Modarres (2011).

distance		A character string defining the distance to be used for determining inclusion
			of a point into the lens (influence region), see Yang and Modarres (2017) for
			details. Possibilities are "Lp" for the Lp-metric (default) or "Mahalanobis" for
			the Mahalanobis distance adjustment.

Lp.p 			A non-negative number defining the distance’s power equal 2 by default (Eu-
			clidean distance); is used only when distance = "Lp".

mah.estimate 		A character string specifying which estimates to use when calculating sample
			covariance matrix; can be "none", "moment" or "MCD", determining whether
			traditional moment or Minimum Covariance Determinant (MCD) (see covMcd)
			estimates for mean and covariance are used. By default "moment" is used. Is
			used only when distance = "Mahalanobis".

mah.parMcd 		The value of the argument alpha for the function covMcd; is used when distance
			= "Mahalanobis" and mah.estimate = "MCD".

"""



spatial.__doc__=""" 

Description
	Calculates the spatial depth of points w.r.t. a multivariate data set.

Usage
	depth.spatial(x, data, mah.estimate = "moment", mah.parMcd = 0.75)

Arguments
x 			Matrix of objects (numerical vector as one object) whose depth is to be calcu-
			lated; each row contains a d-variate point. Should have the same dimension as
			data.

data			Matrix of data where each row contains a d-variate point, w.r.t. which the depth
			is to be calculated.

mah.estimate 		is a character string specifying which estimates to use when calculating sample
			covariance matrix; can be "none", "moment" or "MCD", determining whether
			traditional moment or Minimum Covariance Determinant (MCD) (see covMcd)
			estimates for mean and covariance are used. By default "moment" is used. With
			"none" the non-affine invariant version of Spatial depth is calculated

mah.parMcd 		is the value of the argument alpha for the function covMcd; is used when
			mah.estimate = "MCD".


"""



