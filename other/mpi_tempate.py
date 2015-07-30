from mpi4py import MPI

N_MODELS_TOTAL = 100 #number of models to run
nproc = MPI.COMM_WORLD.Get_size()   	# number of processes
my_rank = MPI.COMM_WORLD.Get_rank()   	# The number/rank of this process

my_node = MPI.Get_processor_name()    	# Node where this MPI process runs

# total number of models to run	
n_models = N_MODELS_TOTAL / nproc		# number of models for each thread
remainder = N_MODELS_TOTAL - ( n_models * nproc )	# the remainder. e.g. your number of models may 
# little trick to spread remainder out among threads
# if say you had 19 total models, and 4 threads
# then n_models = 4, and you have 3 remainder
# this little loop would distribute these three 
if remainder < my_rank + 1:
	my_extra = 0
	extra_below = remainder
else:
	my_extra = 1
	extra_below = my_rank

# where to start and end your loops for each thread
my_nmin = (my_rank * n_models) + extra_below
my_nmax = my_nmin + n_models + my_extra

# total number you actually do
ndo = my_nmax - my_nmin
#print 'My_Rank:', my_rank
#print 'my_node:', my_node
#print 'Time', TI.time()*(my_rank+1*np.pi)

for i in range( my_nmin, my_nmax):
    print i
	#do your loops in here	

#uncomment this at the end of your script
cur.close()
# always call this when finishing up
MPI.Finalize()