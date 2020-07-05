# Python 3 implementation of counting 
# number of possible hexagonal walks 

depth = 16
ways = [[[0 for i in range(17)] 
			for i in range(17)] 
			for i in range(17)] 

def preprocess(list, steps): 
	
	# We initialize our origin with 1 
	ways[0][8][8] = 1

	# For each N = 1 to 14, we traverse in 
	# all possible direction. Using this 3D 
	# array we calculate the number of ways 
	# at each step and the total ways for a 
	# given step shall be found at ways[step 
	# number][8][8] because all the steps after 
	# that will be used to trace back to the 
	# original point index 0:0 according to the image. 
	for N in range(1, 16, 1): 
		for i in range(1, depth, 1): 
			for j in range(1, depth, 1): 
				ways[N][i][j] = (ways[N - 1][i][j + 1] +
								ways[N - 1][i][j - 1] +
								ways[N - 1][i + 1][j] +
								ways[N - 1][i - 1][j] +
								ways[N - 1][i + 1][j - 1] +
								ways[N - 1][i - 1][j + 1]) 

		# This array stores the number of ways 
		# possible for a given step 
		list[N] = ways[N][8][8] 

	print("Number of walks possible is/are", 
								list[steps]) 
print(ways)

# Driver Code 
if __name__ == '__main__': 
	list = [0 for i in range(16)] 
	steps = 4
	
	# Preprocessing all possible ways 
	preprocess(list, steps) 
	
# This code is contributed by 
# Surendra_Gangwar 
